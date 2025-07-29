from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import ROOT
import os
import random
ROOT.PyConfig.IgnoreCommandLineOptions = True


def mk_safe(fct, *args):
    try:
        result = fct(*args)
        return result
    except Exception as e:
        if any('Error in function boost::math::erf_inv' in str(arg)
               for arg in e.args):
            print(
                'WARNING: catching exception and returning -1. Exception arguments: %s'
                % e.args)
            return -1.
        else:
            print('ERROR in mk_safe: %s' % e.args)
            raise e


class muonScaleResProducer(Module):
    def __init__(self, rc_dir, rc_corrections, dataYear):
        self.rc_dir = rc_dir  # Store rc_dir as instance variable
        p_postproc = '%s/src/PhysicsTools/NanoAODTools/python/postprocessing' % os.environ[
            'CMSSW_BASE']
        p_roccor = p_postproc + '/data/' + self.rc_dir
        # Use RoccoR for Run2, MuonScaReWrapper for Run3
        if self.rc_dir.startswith('roccor.Run3'):
            # Load MuonScaReWrapper
            p_helper = '%s/MuonScaReWrapper.cc' % p_roccor
            print('Loading Run3 C++ helper from ' + p_helper)
            ROOT.gROOT.ProcessLine('.L ' + p_helper)
            # Pass JSON correction file
            self._roccor = ROOT.MuonScaReWrapper(p_roccor + '/' + rc_corrections)
        else:
            # Original Run2 RoccoR interface
            if "/RoccoR_cc.so" not in ROOT.gSystem.GetLibraries():
                p_helper = '%s/RoccoR.cc' % p_roccor
                print('Loading C++ helper from ' + p_helper)
                ROOT.gROOT.ProcessLine('.L ' + p_helper)
            self._roccor = ROOT.RoccoR(p_roccor + '/' + rc_corrections)

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        # self.out.branch("Muon_pt", "F", lenVar="nMuon")
        self.out.branch("Muon_corrected_pt", "F", lenVar="nMuon")
        # For Run3, add separate scale and smear uncertainty branches
        if hasattr(self, 'rc_dir') and self.rc_dir.startswith('roccor.Run3'):
            self.out.branch("Muon_scaleUp_pt", "F", lenVar="nMuon")
            self.out.branch("Muon_scaleDown_pt", "F", lenVar="nMuon")
            self.out.branch("Muon_smearUp_pt", "F", lenVar="nMuon")
            self.out.branch("Muon_smearDown_pt", "F", lenVar="nMuon")
        else:
            self.out.branch("Muon_correctedUp_pt", "F", lenVar="nMuon")
            self.out.branch("Muon_correctedDown_pt", "F", lenVar="nMuon")
        self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        muons = Collection(event, "Muon")
        if self.is_mc:
            genparticles = Collection(event, "GenPart")
        roccor = self._roccor
        # Run3 branch: call MuonScaReWrapper interface
        if self.rc_dir.startswith('roccor.Run3'):
            # pt = []
            pt_corr = []
            pt_scale_err = []
            pt_smear_err = []
            for mu in muons:
                if self.is_mc:
                    # MC: smearing - use separate function calls instead of tuple
                    u1 = random.uniform(0.0, 1.0)
                    scale = mk_safe(roccor.kScaleMC, mu.charge, mu.bsConstrainedPt, mu.eta, mu.phi)
                    smear = mk_safe(roccor.kSmearMC, scale, mu.eta, mu.nTrackerLayers)

                    scale_err = mk_safe(roccor.kSmearMCScaleErr, mu.charge, smear, mu.eta, mu.phi)
                    smear_err = mk_safe(roccor.kSmearMCSmearErr, mu.charge, scale, smear, mu.eta)

                    corr = smear

                else:
                    # Data: scaling
                    corr = mk_safe(roccor.kScaleDT, mu.charge, mu.pt, mu.eta, mu.phi)
                    scale_err = mk_safe(roccor.kScaleDTerror, mu.charge, mu.pt, mu.eta, mu.phi)
                    smear_err = 0.0  # No smearing for data
                
                # pt.append(mu.pt)
                pt_corr.append(corr)
                pt_scale_err.append(scale_err)
                pt_smear_err.append(smear_err)

            # self.out.fillBranch("Muon_pt", pt)
            self.out.fillBranch("Muon_corrected_pt", pt_corr)
            pt_scale_up = list(
                max(pt_corr[imu] + pt_scale_err[imu], 0.0)
                for imu, mu in enumerate(muons))
            pt_scale_down = list(
                max(pt_corr[imu] - pt_scale_err[imu], 0.0)
                for imu, mu in enumerate(muons))
            pt_smear_up = list(
                max(pt_corr[imu] + pt_smear_err[imu], 0.0)
                for imu, mu in enumerate(muons))
            pt_smear_down = list(
                max(pt_corr[imu] - pt_smear_err[imu], 0.0)
                for imu, mu in enumerate(muons))
            self.out.fillBranch("Muon_scaleUp_pt", pt_scale_up)
            self.out.fillBranch("Muon_scaleDown_pt", pt_scale_down)
            self.out.fillBranch("Muon_smearUp_pt", pt_smear_up)
            self.out.fillBranch("Muon_smearDown_pt", pt_smear_down)
        else:
            # Original Run2 processing
            genparticles = Collection(event, "GenPart")
            pt_corr = []
            pt_err = []
            for mu in muons:
                genIdx = mu.genPartIdx
                if genIdx >= 0 and genIdx < len(genparticles):
                    genMu = genparticles[genIdx]
                    pt_corr.append(mu.pt *
                                   mk_safe(roccor.kSpreadMC, mu.charge, mu.pt,
                                           mu.eta, mu.phi, genMu.pt))
                    pt_err.append(mu.pt *
                                  mk_safe(roccor.kSpreadMCerror, mu.charge,
                                          mu.pt, mu.eta, mu.phi, genMu.pt))
                else:
                    u1 = random.uniform(0.0, 1.0)
                    pt_corr.append(
                        mu.pt * mk_safe(roccor.kSmearMC, mu.charge, mu.pt,
                                        mu.eta, mu.phi, mu.nTrackerLayers, u1))
                    pt_err.append(
                        mu.pt * mk_safe(roccor.kSmearMCerror, mu.charge, mu.pt,
                                        mu.eta, mu.phi, mu.nTrackerLayers, u1))

            self.out.fillBranch("Muon_corrected_pt", pt_corr)
            pt_corr_up = list(
                max(pt_corr[imu] + pt_err[imu], 0.0)
                for imu, mu in enumerate(muons))
            pt_corr_down = list(
                max(pt_corr[imu] - pt_err[imu], 0.0)
                for imu, mu in enumerate(muons))
            self.out.fillBranch("Muon_correctedUp_pt", pt_corr_up)
            self.out.fillBranch("Muon_correctedDown_pt", pt_corr_down)
        return True


muonScaleRes2016_UL16PreVFP = lambda: muonScaleResProducer('roccor.Run2.v5',
                                                'RoccoR2016aUL.txt', 2016)
muonScaleRes2016_UL16PostVFP = lambda: muonScaleResProducer('roccor.Run2.v5',
                                                'RoccoR2016bUL.txt', 2016)
muonScaleRes2017 = lambda: muonScaleResProducer('roccor.Run2.v5',
                                                'RoccoR2017UL.txt', 2017)
muonScaleRes2018 = lambda: muonScaleResProducer('roccor.Run2.v5',
                                                'RoccoR2018UL.txt', 2018)
# Add Run3 muon correction lambdas
muonScaleResRun3_2022 = lambda: muonScaleResProducer('roccor.Run3.v1', '2022_Summer22.json', 2022)
muonScaleResRun3_2022EE = lambda: muonScaleResProducer('roccor.Run3.v1', '2022_Summer22EE.json', 2022)
muonScaleResRun3_2023 = lambda: muonScaleResProducer('roccor.Run3.v1', '2023_Summer23.json', 2023)
muonScaleResRun3_2023BPix = lambda: muonScaleResProducer('roccor.Run3.v1', '2023_Summer23BPix.json', 2023)
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon
# run3 
#https://muon-wiki.docs.cern.ch/guidelines/corrections/#__tabbed_6_2
# muon scale and resolution for Run3
# https://gitlab.cern.ch/cms-muonPOG/muonscarekit
