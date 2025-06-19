#include <string>
#include <memory>
#include <cmath>
#include "correction.h"
#include "correction.cc"

// Global CorrectionSet pointer used by MuonScaRe.cc
static std::shared_ptr<correction::CorrectionSet> cset_ptr;
static correction::CorrectionSet* cset;

#include "MuonScaRe.cc"  // include the core C++ implementation of scale and resolution

// Wrapper class to provide a RoccoR-like interface for Run3 muon corrections
class MuonScaReWrapper {
public:
    // Constructor loads the JSON correction parameters
    MuonScaReWrapper(const std::string &jsonFile) {
        // load CorrectionSet and set global pointer
        cset_ptr = correction::CorrectionSet::from_file(jsonFile);
        cout << "Loaded json file: " << jsonFile << endl;
        cset = cset_ptr.get();
    }

    // Data scale correction
    double kScaleDT(int charge, double pt, double eta, double phi) const {
        return pt_scale(true, pt, eta, phi, charge);
    }
    // Uncertainty on data scale correction
    double kScaleDTerror(int charge, double pt, double eta, double phi) const {
        double nom = kScaleDT(charge, pt, eta, phi);
        double up  = pt_scale_var(pt, eta, phi, charge, std::string("up"));
        return fabs(up - nom);
    }

    // MC resolution smearing
    double kSmearMC(int charge, double pt, double eta, double phi, int nTrkLayers, double rndm) const {
        double scaled = pt_scale(false, pt, eta, phi, charge);
        return pt_resol(scaled, eta, nTrkLayers);
    }
    // Uncertainty on MC smearing
    double kSmearMCerror(int charge, double pt, double eta, double phi, int nTrkLayers, double rndm) const {
        double scaled = pt_scale(false, pt, eta, phi, charge);
        double smeared = pt_resol(scaled, eta, nTrkLayers);
        double up = pt_resol_var(scaled, smeared, eta, std::string("up"));
        return fabs(up - smeared);
    }
};
