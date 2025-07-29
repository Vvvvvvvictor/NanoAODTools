#include <string>
#include <memory>
#include <cmath>
#include <tuple>
#include "correction.h"
#include "correction.cc"
// get correction code from https://github.com/cms-nanoAOD/correctionlib

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

    // MC resolution smearing - returns the corrected pt
    double kScaleMC(int charge, double pt, double eta, double phi) const {
        return pt_scale(false, pt, eta, phi, charge);
    }

    double kSmearMC(double scaled, double eta, int nTrkLayers) const {
        return pt_resol(scaled, eta, nTrkLayers);
    }
    
    // MC scale uncertainty
    double kSmearMCScaleErr(int charge, double smeared, double eta, double phi) const {
        double scale_up = pt_scale_var(smeared, eta, phi, charge, std::string("up"));
        return fabs(scale_up - smeared);
    }
    
    // MC smear uncertainty
    double kSmearMCSmearErr(int charge, double scaled, double smeared, double eta) const {
        double smear_up = pt_resol_var(scaled, smeared, eta, std::string("up"));
        return fabs(smear_up - smeared);
    }
};
