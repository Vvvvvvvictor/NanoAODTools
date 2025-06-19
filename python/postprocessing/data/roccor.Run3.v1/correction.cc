#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif
#include "rapidjson/document.h"
#pragma GCC diagnostic pop
#include "rapidjson/filereadstream.h"
#include "rapidjson/error/en.h"
#include <optional>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <cstdlib> // std::abort
#include "correction.h"
#if __has_include(<zlib.h>)
#include <zlib.h>
#include "gzfilereadstream.h"
#define WITH_ZLIB 1
#endif

using namespace correction;

//! helper function for parsing flow behavior from string
_FlowBehavior parse_flowbehavior(const rapidjson::Value& flowbehavior) {
  if ( flowbehavior == "clamp" ) {
    return _FlowBehavior::clamp;
  }
  else if ( flowbehavior == "error" ) {
    return _FlowBehavior::error;
  }
  else if ( flowbehavior == "wrap" ) {
    return _FlowBehavior::wrap;
  }
  else {
    return _FlowBehavior::value;
  }
}

class correction::JSONObject {
  public:
    JSONObject(rapidjson::Value::ConstObject&& json) : json_(json) { }
    // necessary to force use of const Value::GetObject() method
    // must check if json is an object in calling code!
    JSONObject(const rapidjson::Document& json) : json_(json.GetObject()) { }

    template<typename T>
    T getRequired(const char * key) const {
      const auto it = json_.FindMember(key);
      if ( it != json_.MemberEnd() ) {
        if ( it->value.template Is<T>() ) {
          return it->value.template Get<T>();
        } else {
          throw std::runtime_error(
              "Encountered invalid type for required attribute '"
              + std::string(key) + "'");
        }
      }
      throw std::runtime_error(
          "Object missing required attribute '"
          + std::string(key) + "'");
    }

    const rapidjson::Value& getRequiredValue(const char * key) const {
      const auto it = json_.FindMember(key);
      if ( it != json_.MemberEnd() ) {
        return it->value;
      }
      throw std::runtime_error(
          "Object missing required attribute '"
          + std::string(key) + "'");
    }

    template<typename T>
    std::optional<T> getOptional(const char * key) const {
      const auto it = json_.FindMember(key);
      if ( it != json_.MemberEnd() ) {
        if ( it->value.template Is<T>() ) {
          return it->value.template Get<T>();
        } else if ( it->value.IsNull() ) {
          return std::nullopt;
        } else {
          throw std::runtime_error(
              "Encountered invalid type for optional attribute '"
              + std::string(key) + "'");
        }
      }
      return std::nullopt;
    }

    // escape hatch
    inline auto FindMember(const char * key) const { return json_.FindMember(key); }
    inline auto MemberEnd() const { return json_.MemberEnd(); }

  private:
    rapidjson::Value::ConstObject json_;
};

template<>
std::string_view JSONObject::getRequired<std::string_view>(const char * key) const {
  const auto it = json_.FindMember(key);
  if ( it != json_.MemberEnd() ) {
    if ( it->value.IsString() ) {
      return std::string_view(it->value.GetString(), it->value.GetStringLength());
    } else {
      throw std::runtime_error(
          "Encountered invalid type for required attribute '"
          + std::string(key) + "'");
    }
  }
  throw std::runtime_error(
      "Object missing required attribute '"
      + std::string(key) + "'");
}

namespace {
  Content resolve_content(const rapidjson::Value& json, const Correction& context) {
    if ( json.IsNumber() ) { return json.GetDouble(); }
    else if ( json.IsObject() && json.HasMember("nodetype") ) {
      auto obj = JSONObject(json.GetObject());
      auto type = obj.getRequired<std::string_view>("nodetype");
      if ( type == "binning" ) { return Binning(obj, context); }
      else if ( type == "multibinning" ) { return MultiBinning(obj, context); }
      else if ( type == "category" ) { return Category(obj, context); }
      else if ( type == "formula" ) { return Formula(obj, context); }
      else if ( type == "formularef" ) { return FormulaRef(obj, context); }
      else { throw std::runtime_error("Unrecognized Content object nodetype"); }
    }
    throw std::runtime_error("Invalid Content node type");
  }

  struct node_evaluate {
    double operator() (double node) { return node; }

    template <class Node>
    double operator() (const Node &node) {
      return node.evaluate(values);
    }

    const std::vector<Variable::Type>& values;
  };

  std::size_t find_bin_idx(Variable::Type value_variant,
                           const std::variant<_UniformBins, _NonUniformBins> &bins_,
                           const _FlowBehavior &flow,
                           std::size_t variableIdx,
                           const char *name)
  {
    double value = std::visit([](auto&& arg) -> double {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, int>) return static_cast<double>(arg);
      else if constexpr (std::is_same_v<T, double>) return arg;
      else throw std::logic_error("I should not have ever seen a string");
      }, value_variant);
    if ( auto *bins = std::get_if<_UniformBins>(&bins_) ) { // uniform binning
      if (value < bins->low || value >= bins->high) {
        switch (flow) {
          case _FlowBehavior::value:
            return bins->n; // the default value is stored at the end of the content array, after the last bin
          case _FlowBehavior::clamp:
            return value < bins->low ? 0 : bins->n - 1; // assuming we always have at least 1 bin
          case _FlowBehavior::wrap:
            break;
          case _FlowBehavior::error:
            const std::string belowOrAbove = value < bins->low ? "below" : "above";
            auto msg = "Index " + belowOrAbove + " bounds in " + name + " for input argument " + std::to_string(variableIdx) + " value: " + std::to_string(value);
            throw std::runtime_error(std::move(msg));
        }
      }

      double norm_value = ((value - bins->low) / (bins->high - bins->low));
      if (flow == _FlowBehavior::wrap) {
        norm_value -= std::floor(norm_value);
      }
      std::size_t binIdx = bins->n * norm_value;
      return binIdx;
    }

    // otherwise we have non-uniform binning
    using namespace std::string_literals;
    const auto bins = std::get<_NonUniformBins>(bins_);
    if ( flow == _FlowBehavior::wrap ) {
      double low = bins[0];
      double high = bins[bins.size() - 1];
      double norm_value = (value - low) / (high - low);
      norm_value -= std::floor(norm_value);
      value = low + norm_value * (high - low);
    }

    auto it = std::upper_bound(std::begin(bins), std::end(bins), value);
    if ( it == std::begin(bins) ) { // underflow
      if ( flow == _FlowBehavior::value ) {
        return bins.size() - 1; // the default value is stored at the end of the content array, after the last bin
      }
      else if ( flow == _FlowBehavior::error ) {
        throw std::runtime_error("Index below bounds in "s + name + " for input argument " + std::to_string(variableIdx) + " value: " + std::to_string(value));
      }
      else if ( flow == _FlowBehavior::wrap ) {
        throw std::logic_error("I should not have ever seen an underflow");
      }
      else { // clamp
        it++;
      }
    }
    else if ( it == std::end(bins) ) { // overflow
      if ( flow == _FlowBehavior::value ) {
        return bins.size() - 1;
      }
      else if ( flow == _FlowBehavior::error ) {
        throw std::runtime_error("Index above bounds in "s + name + " for input argument " + std::to_string(variableIdx) + " value: " + std::to_string(value));
      }
      else if ( flow == _FlowBehavior::wrap ) {
        throw std::logic_error("I should not have ever seen an overflow");
      }
      else { // clamp
        it--;
      }
    }

    // -1 because upper_bound returns the edge _after_ the bin we are interested in
    const std::size_t binIdx = std::distance(std::begin(bins), it) - 1;
    return binIdx;
  }

  size_t find_input_index(const std::string_view name, const std::vector<Variable> &inputs) {
    size_t idx = 0;
    for (const auto& var : inputs) {
      if ( name == var.name() ) return idx;
      idx++;
    }
    throw std::runtime_error("Error: could not find variable " + std::string(name) + " in inputs");
  }

  double parse_edge(const rapidjson::Value& edge) {
    if ( edge.IsDouble() ) {
      return edge.GetDouble();
    } else if ( edge.IsString() ) {
      std::string_view str = edge.GetString();
      if ((str == "inf") || (str == "+inf")) return std::numeric_limits<double>::infinity();
      else if (str == "-inf") return -std::numeric_limits<double>::infinity();
    }
    throw std::runtime_error("Invalid edge type");
  }

  std::vector<double> parse_bin_edges(const rapidjson::Value::ConstArray& edges) {
    std::vector<double> result;
    result.reserve(edges.Size());
    for (const auto& edge : edges) {
      double val = parse_edge(edge);
      if ( result.size() > 0 && result.back() >= val ) {
        throw std::runtime_error("binning edges are not monotone increasing");
      }
      result.push_back(val);
    }
    return result;
  }
} // end of anonymous namespace

Variable::Variable(const JSONObject& json) :
  name_(json.getRequired<const char *>("name")),
  description_(json.getOptional<const char*>("description").value_or(""))
{
  auto type = json.getRequired<std::string_view>("type");
  if (type == "string") { type_ = VarType::string; }
  else if (type == "int") { type_ = VarType::integer; }
  else if (type == "real") { type_ = VarType::real; }
  else { throw std::runtime_error("Unrecognized variable type"); }
}

std::string Variable::typeStr() const {
  if ( type_ == VarType::string ) { return "string"; }
  else if ( type_ == VarType::integer ) { return "int"; }
  else if ( type_ == VarType::real ) { return "real"; }
  return "";
}

void Variable::validate(const Type& t) const {
  if ( std::holds_alternative<std::string>(t) ) {
    if ( type_ != VarType::string ) {
      throw std::runtime_error("Input " + name() + " has wrong type: got string expected " + typeStr());
    }
  }
  else if ( std::holds_alternative<int>(t) ) {
    if ( type_ != VarType::integer ) {
      throw std::runtime_error("Input " + name() + " has wrong type: got int expected " + typeStr());
    }
  }
  else if ( std::holds_alternative<double>(t) ) {
    if ( type_ != VarType::real ) {
      throw std::runtime_error("Input " + name() + " has wrong type: got real-valued expected " + typeStr());
    }
  }
}

// Variable Variable::from_string(const char * data) {
//   rapidjson::Document json;
//   rapidjson::ParseResult ok = json.Parse(data);
//   if (!ok) {
//     throw std::runtime_error(
//         std::string("JSON parse error: ") + rapidjson::GetParseError_En(ok.Code())
//         + " at offset " + std::to_string(ok.Offset())
//         );
//   }
//   if ( ! json.IsObject() ) { throw std::runtime_error("Expected Variable object"); }
//   return Variable(json);
// }

Formula::Formula(const JSONObject& json, const Correction& context, bool generic)
  : Formula(json, context.inputs(), generic) {}

Formula::Formula(const JSONObject& json, const std::vector<Variable>& inputs, bool generic) :
  expression_(json.getRequired<const char *>("expression")),
  generic_(generic)
{
  auto parser_type = json.getRequired<std::string_view>("parser");
  if (parser_type == "TFormula") { type_ = FormulaAst::ParserType::TFormula; }
  else if (parser_type == "numexpr") {
    type_ = FormulaAst::ParserType::numexpr;
    throw std::runtime_error("numexpr formula parser is not yet supported");
  }
  else { throw std::runtime_error("Unrecognized formula parser type"); }

  std::vector<size_t> variableIdx;
  for (const auto& item : json.getRequired<rapidjson::Value::ConstArray>("variables")) {
    auto idx = find_input_index(item.GetString(), inputs);
    if ( inputs[idx].type() != Variable::VarType::real ) {
      throw std::runtime_error("Formulas only accept real-valued inputs, got type "
          + inputs[idx].typeStr() + " for variable " + inputs[idx].name());
    }
    variableIdx.push_back(idx);
  }

  std::vector<double> params;
  if ( auto items = json.getOptional<rapidjson::Value::ConstArray>("parameters") ) {
    for (const auto& item : *items) {
      params.push_back(item.GetDouble());
    }
  }

  ast_ = std::make_unique<FormulaAst>(FormulaAst::parse(type_, expression_, params, variableIdx, !generic));
}

// Formula::Ref Formula::from_string(const char * data, std::vector<Variable>& inputs) {
//   rapidjson::Document json;
//   rapidjson::ParseResult ok = json.Parse(data);
//   if (!ok) {
//     throw std::runtime_error(
//         std::string("JSON parse error: ") + rapidjson::GetParseError_En(ok.Code())
//         + " at offset " + std::to_string(ok.Offset())
//         );
//   }
//   if ( ! json.IsObject() ) { throw std::runtime_error("Expected Formula object"); }
//   return std::make_shared<Formula>(json, inputs);
// }

double Formula::evaluate(const std::vector<Variable::Type>& values) const {
  if ( generic_ ) {
    throw std::runtime_error("Generic formulas must be evaluated with parameters");
  }
  return ast_->evaluate(values, {});
}

double Formula::evaluate(const std::vector<Variable::Type>& values, const std::vector<double>& params) const {
  return ast_->evaluate(values, params);
}

FormulaRef::FormulaRef(const JSONObject& json, const Correction& context) {
  formula_ = context.formula_ref(json.getRequired<int>("index"));
  for (const auto& item : json.getRequired<rapidjson::Value::ConstArray>("parameters")) {
    parameters_.push_back(item.GetDouble());
  }
}

double FormulaRef::evaluate(const std::vector<Variable::Type>& values) const {
  return formula_->evaluate(values, parameters_);
}

Binning::Binning(const JSONObject& json, const Correction& context)
{
  const auto& content = json.getRequired<rapidjson::Value::ConstArray>("content");

  // set bins_
  const auto &edgesObj = json.getRequiredValue("edges");
  if ( edgesObj.IsArray() ) { // non-uniform binning
    std::vector<double> edges = parse_bin_edges(edgesObj.GetArray());
    if ( edges.size() != content.Size() + 1 ) {
      throw std::runtime_error("Inconsistency in Binning: number of content nodes does not match binning");
    }
    _NonUniformBins bins{std::move(edges)};
    bins_ = std::move(bins);
  } else if ( edgesObj.IsObject() ) { // UniformBinning
    const JSONObject uniformBins{edgesObj.GetObject()};
    const auto n = uniformBins.getRequired<uint32_t>("n");
    if ( n == 0 ) {
      throw std::runtime_error("Error when processing Binning with UniformBinning: number of bins is zero");
    }
    if ( n != content.Size() ) {
      throw std::runtime_error("Inconsistency in Binning: number of content nodes does not match binning");
    }
    const auto low = uniformBins.getRequired<double>("low");
    const auto high = uniformBins.getRequired<double>("high");
    bins_ = _UniformBins{n, low, high};
  } else {
    throw std::runtime_error ("Error when processing Binning: edges are neither an array nor a UniformBinning object");
  }

  variableIdx_ = find_input_index(json.getRequired<std::string_view>("input"), context.inputs());
  if ( context.inputs().at(variableIdx_).type() == Variable::VarType::string ) {
    throw std::runtime_error("Binning cannot use string inputs as binning variables");
  }
  Content default_value{0.};
  const auto& flowbehavior = json.getRequiredValue("flow");
  flow_ = parse_flowbehavior(flowbehavior);
  if (flow_ == _FlowBehavior::value) {
      default_value = resolve_content(flowbehavior, context);
  }

  // set bin contents
  for (size_t i=0; i < content.Size(); ++i)
    contents_.push_back(resolve_content(content[i], context));
  contents_.push_back(std::move(default_value));
}

double Binning::evaluate(const std::vector<Variable::Type>& values) const
{
  std::size_t binIdx = find_bin_idx(values[variableIdx_], bins_, flow_, variableIdx_, "Binning");
  const Content& child = contents_[binIdx];
  return std::visit(node_evaluate{values}, child);
}

MultiBinning::MultiBinning(const JSONObject& json, const Correction& context)
{
  const auto& inputs = json.getRequired<rapidjson::Value::ConstArray>("inputs");

  const auto& edges = json.getRequired<rapidjson::Value::ConstArray>("edges");
  axes_.reserve(edges.Size());
  size_t idx {0};
  for (const auto& dimension : edges) {
    const auto& input = inputs[idx];
    if ( dimension.IsArray() ) { // non-uniform binning
      std::vector<double> dim_edges = parse_bin_edges(dimension.GetArray());
      if ( ! input.IsString() ) { throw std::runtime_error("invalid multibinning input type"); }
      size_t variableIdx = find_input_index(input.GetString(), context.inputs());
      if ( context.inputs().at(variableIdx).type() == Variable::VarType::string ) {
        throw std::runtime_error("MultiBinning cannot use string inputs as binning variables");
      }
      axes_.push_back({variableIdx, 0, _NonUniformBins(std::move(dim_edges))});
    } else if ( dimension.IsObject() ) { // UniformBinning
      const JSONObject uniformBins{dimension.GetObject()};
      const auto n = uniformBins.getRequired<uint32_t>("n");
      if ( n == 0 ) {
        auto msg = "Error when processing MultiBinning: number of bins for dimension " + std::to_string(idx) + " is zero";
        throw std::runtime_error(std::move(msg));
      }
      const auto low = uniformBins.getRequired<double>("low");
      const auto high = uniformBins.getRequired<double>("high");
      size_t variableIdx = find_input_index(input.GetString(), context.inputs());
      if ( context.inputs().at(variableIdx).type() == Variable::VarType::string ) {
        throw std::runtime_error("MultiBinning cannot use string inputs as binning variables");
      }
      axes_.push_back({variableIdx, 0, _UniformBins{n, low, high}});
    } else {
      auto msg = "Error when processing MultiBinning: edges for dimension " + std::to_string(idx) + " are neither an array nor a UniformBinning object";
      throw std::runtime_error (std::move(msg));
    }
    ++idx;
  }

  const auto& content = json.getRequired<rapidjson::Value::ConstArray>("content");
  size_t stride {1};
  --idx; // now corresponds to the last dimension
  for (auto it=axes_.rbegin(); it != axes_.rend(); ++it) {
    it->stride = stride;
    stride *= nbins(idx);
    --idx;
  }
  content_.reserve(content.Size() + 1); // + 1 for default value
  for (const auto& item : content) {
    content_.push_back(resolve_content(item, context));
  }
  if ( content_.size() != stride ) {
    throw std::runtime_error("Inconsistency in MultiBinning: number of content nodes does not match binning");
  }

  const auto& flowbehavior = json.getRequiredValue("flow");
  flow_ = parse_flowbehavior(flowbehavior);
  if (flow_ == _FlowBehavior::value) {
      content_.push_back(resolve_content(flowbehavior, context));
  }
}

double MultiBinning::evaluate(const std::vector<Variable::Type>& values) const
{
  size_t idx {0};
  size_t localidx {0};
  size_t dim {0};

  for (const auto& [variableIdx, stride, edgesVariant] : axes_) {
    localidx = find_bin_idx(values[variableIdx], edgesVariant, flow_, variableIdx, "MultiBinning");
    if ( localidx == nbins(dim) ) // find_bin_idx is indicating we need to return the default value
      return std::visit(node_evaluate{values}, content_.back());
    idx += localidx * stride;
    ++dim;
  }

  const Content& child = content_.at(idx);
  return std::visit(node_evaluate{values}, child);
}

size_t MultiBinning::nbins(size_t dimension) const
{
  if ( const auto *bins = std::get_if<_UniformBins>(&axes_[dimension].bins) )
    return bins->n; // using uniform bins

  // otherwise we must have non-uniform bins
  const auto &bins = std::get<_NonUniformBins>(axes_[dimension].bins);
  return bins.size() - 1;
}

Category::Category(const JSONObject& json, const Correction& context)
{
  variableIdx_ = find_input_index(json.getRequired<std::string_view>("input"), context.inputs());
  const auto& variable = context.inputs()[variableIdx_];
  if ( variable.type() == Variable::VarType::string ) {
    map_ = StrMap();
  } // (default-constructed as IntMap)
  for (const auto& kv_pair : json.getRequired<rapidjson::Value::ConstArray>("content"))
  {
    if ( ! (kv_pair.IsObject() && kv_pair.HasMember("key") && kv_pair.HasMember("value")) ) {
      throw std::runtime_error("Expected CategoryItem object");
    }
    if ( kv_pair["key"].IsString() ) {
      if ( variable.type() != Variable::VarType::string ) {
        throw std::runtime_error("Category got a key of type string, but its input is type " + variable.typeStr());
      }
      std::get<StrMap>(map_).try_emplace(kv_pair["key"].GetString(), resolve_content(kv_pair["value"], context));
    }
    else if ( kv_pair["key"].IsInt() ) {
      if ( variable.type() != Variable::VarType::integer ) {
        throw std::runtime_error("Category got a key of type int, but its input is type " + variable.typeStr());
      }
      std::get<IntMap>(map_).try_emplace(kv_pair["key"].GetInt(), resolve_content(kv_pair["value"], context));
    }
    else {
      throw std::runtime_error("Invalid key type in Category");
    }
  }

  const auto def = json.FindMember("default");
  if ( def != json.MemberEnd() && ! def->value.IsNull() ) {
    default_ = std::make_unique<Content>(resolve_content(def->value, context));
  }
}

double Category::evaluate(const std::vector<Variable::Type>& values) const {
  const Content* child = nullptr;
  if ( auto pval = std::get_if<std::string>(&values[variableIdx_]) ) {
    try {
      child = &std::get<StrMap>(map_).at(*pval);
    } catch (std::out_of_range& ex) {
      if ( default_ ) {
        child = default_.get();
      }
      else {
        throw std::out_of_range("Index not available in Category for input argument " + std::to_string(variableIdx_) + " val: " + *pval);
      }
    }
  }
  else if ( auto pval = std::get_if<int>(&values[variableIdx_]) ) {
    try {
      child = &std::get<IntMap>(map_).at(*pval);
    } catch (std::out_of_range& ex) {
      if ( default_ ) {
        child = default_.get();
      }
      else {
        throw std::out_of_range("Index not available in Category for input argument " + std::to_string(variableIdx_) + " val: " + std::to_string(*pval));
      }
    }
  } else {
    throw std::runtime_error("Invalid variable type");
  }

  return std::visit(node_evaluate{values}, *child);
}

Correction::Correction(const JSONObject& json) :
  name_(json.getRequired<const char *>("name")),
  description_(json.getOptional<const char*>("description").value_or("")),
  version_(json.getRequired<int>("version")),
  output_(json.getRequired<rapidjson::Value::ConstObject>("output"))
{
  if ( output_.type() != Variable::VarType::real ) { throw std::runtime_error("Outputs can only be real-valued"); }
  for (const auto& item : json.getRequired<rapidjson::Value::ConstArray>("inputs")) {
    if ( ! item.IsObject() ) { throw std::runtime_error("invalid input item type"); }
    inputs_.emplace_back(item.GetObject());
  }
  if ( const auto& items = json.getOptional<rapidjson::Value::ConstArray>("generic_formulas") ) {
    for (const auto& item : *items) {
      if ( ! item.IsObject() ) { throw std::runtime_error("invalid generic_formulas item type"); }
      formula_refs_.push_back(std::make_shared<Formula>(item.GetObject(), *this, true));
    }
  }

  data_ = resolve_content(json.getRequiredValue("data"), *this);
  initialized_ = true;
}

double Correction::evaluate(const std::vector<Variable::Type>& values) const {
  if ( ! initialized_ ) {
    throw std::logic_error("Not initialized");
  }
  if ( values.size() > inputs_.size() ) {
    throw std::runtime_error("Too many inputs");
  }
  else if ( values.size() < inputs_.size() ) {
    throw std::runtime_error("Insufficient inputs");
  }
  for (size_t i=0; i < inputs_.size(); ++i) {
    inputs_[i].validate(values[i]);
  }
  return std::visit(node_evaluate{values}, data_);
}

std::unique_ptr<CorrectionSet> CorrectionSet::from_file(const std::string& fn) {
  rapidjson::Document json;
  FILE* fp = fopen(fn.c_str(), "rb");
  if ( fp == nullptr ) {
    throw std::runtime_error("Failed to open file: " + fn);
  }
  constexpr unsigned char magicref[2] = {0x1f, 0x8b};
  unsigned char magic[2];
  if (fread(magic, sizeof *magic, 2, fp) != 2) {
    throw std::runtime_error("Failed to read file magic: " + fn);
  }
  rewind(fp);
  char readBuffer[65536];
  rapidjson::ParseResult ok;
  if (memcmp(magic, magicref, sizeof(magic)) == 0) {
    fclose(fp);
#ifdef WITH_ZLIB
    gzFile_s* fpz = gzopen(fn.c_str(), "r");
    rapidjson::GzFileReadStream is(fpz, readBuffer, sizeof(readBuffer));
    ok = json.ParseStream<rapidjson::kParseNanAndInfFlag>(is);
    gzclose(fpz);
#else
    throw std::runtime_error("Gzip-compressed JSON files are only supported if ZLIB is found when the package is built");
#endif
  } else {
    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    ok = json.ParseStream<rapidjson::kParseNanAndInfFlag>(is);
    fclose(fp);
  }
  if (!ok) {
    throw std::runtime_error(
        std::string("JSON parse error: ") + rapidjson::GetParseError_En(ok.Code())
        + " at offset " + std::to_string(ok.Offset())
        );
  }
  if ( ! json.IsObject() ) { throw std::runtime_error("Expected CorrectionSet object"); }
  return std::make_unique<CorrectionSet>(json);
}

// std::unique_ptr<CorrectionSet> CorrectionSet::from_string(const char * data) {
//   rapidjson::Document json;
//   rapidjson::ParseResult ok = json.Parse<rapidjson::kParseNanAndInfFlag>(data);
//   if (!ok) {
//     throw std::runtime_error(
//         std::string("JSON parse error: ") + rapidjson::GetParseError_En(ok.Code())
//         + " at offset " + std::to_string(ok.Offset())
//         );
//   }
//   if ( ! json.IsObject() ) { throw std::runtime_error("Expected CorrectionSet object"); }
//   return std::make_unique<CorrectionSet>(json);
// }

CorrectionSet::CorrectionSet(const JSONObject& json) {
  schema_version_ = json.getRequired<int>("schema_version");
  if ( schema_version_ > evaluator_version ) {
    throw std::runtime_error("Evaluator is designed for schema v" + std::to_string(evaluator_version) + " and is not forward-compatible");
  }
  else if ( schema_version_ < evaluator_version ) {
    throw std::runtime_error("Evaluator is designed for schema v" + std::to_string(evaluator_version) + " and is not backward-compatible");
  }
  description_ = json.getOptional<const char*>("description").value_or("");
  for (const auto& item : json.getRequired<rapidjson::Value::ConstArray>("corrections")) {
    if ( ! item.IsObject() ) { throw std::runtime_error("Expected Correction object"); }
    auto corr = std::make_shared<Correction>(item.GetObject());
    if ( corrections_.find(corr->name()) != corrections_.end() ) {
      throw std::runtime_error("Duplicate Correction name: " + corr->name());
    }
    corrections_[corr->name()] = corr;
  }
}

bool CorrectionSet::validate() {
  // TODO: validate with https://rapidjson.org/md_doc_schema.html
  return true;
}
