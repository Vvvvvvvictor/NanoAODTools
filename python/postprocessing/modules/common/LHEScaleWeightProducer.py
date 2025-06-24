import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import numpy
class LHEScaleWeightProducer(Module):
    def __init__(self,year="2016"):
        self.year=year
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("LHEScaleWeight_Zero", "F")
        self.out.branch("LHEScaleWeight_One", "F")
        self.out.branch("LHEScaleWeight_Two", "F")
        self.out.branch("LHEScaleWeight_Three", "F")
        self.out.branch("LHEScaleWeight_Four", "F")
        self.out.branch("LHEScaleWeight_Five", "F")
        self.out.branch("LHEScaleWeight_Six", "F")
        self.out.branch("LHEScaleWeight_Seven", "F")
        self.out.branch("LHEScaleWeight_Eight", "F")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        # lheScaleWeights = Collection(event, "LHEScaleWeight")
        lheScaleWeights = event["LHEScaleWeight"]
        if len(lheScaleWeights)==8:
            self.out.fillBranch("LHEScaleWeight_Zero", lheScaleWeights[0])
            self.out.fillBranch("LHEScaleWeight_One", lheScaleWeights[1])
            self.out.fillBranch("LHEScaleWeight_Two", lheScaleWeights[2])
            self.out.fillBranch("LHEScaleWeight_Three", lheScaleWeights[3])
            self.out.fillBranch("LHEScaleWeight_Four", lheScaleWeights[4])
            self.out.fillBranch("LHEScaleWeight_Five", lheScaleWeights[5])
            self.out.fillBranch("LHEScaleWeight_Six", lheScaleWeights[6])
            self.out.fillBranch("LHEScaleWeight_Seven", lheScaleWeights[7])
            self.out.fillBranch("LHEScaleWeight_Eight", 0.)

        if len(lheScaleWeights)==9:
            self.out.fillBranch("LHEScaleWeight_Zero", lheScaleWeights[0])
            self.out.fillBranch("LHEScaleWeight_One", lheScaleWeights[1])
            self.out.fillBranch("LHEScaleWeight_Two", lheScaleWeights[2])
            self.out.fillBranch("LHEScaleWeight_Three", lheScaleWeights[3])
            self.out.fillBranch("LHEScaleWeight_Four", lheScaleWeights[4])
            self.out.fillBranch("LHEScaleWeight_Five", lheScaleWeights[5])
            self.out.fillBranch("LHEScaleWeight_Six", lheScaleWeights[6])
            self.out.fillBranch("LHEScaleWeight_Seven", lheScaleWeights[7])
            self.out.fillBranch("LHEScaleWeight_Eight", lheScaleWeights[8])
        return True
LHEScaleSF_UL16PreVFP  = lambda : LHEScaleWeightProducer(year="2016")
LHEScaleSF_UL16PostVFP = lambda : LHEScaleWeightProducer(year="2016post")
LHEScaleSF_UL17         = lambda : LHEScaleWeightProducer(year="UL17")
LHEScaleSF_UL18         = lambda : LHEScaleWeightProducer(year="UL18")
LHEScaleSF_UL22PreEE = lambda : LHEScaleWeightProducer(year="2022")
LHEScaleSF_UL22PostEE = lambda : LHEScaleWeightProducer(year="2022EE")
LHEScaleSF_UL23PreBPIX = lambda : LHEScaleWeightProducer(year="2023")
LHEScaleSF_UL23PostBPIX = lambda : LHEScaleWeightProducer(year="2023BPix")