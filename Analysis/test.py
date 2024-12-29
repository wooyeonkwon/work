import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("Analysis")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
directory_path = '/data1/users/dndus0107/AnalysisResults/Muon0/crab_MuonSkimming_Run2024C0/241109_142235/0000/'
file_list = os.listdir(directory_path)
file_paths = [f'file://{directory_path}{filename}' for filename in file_list]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*file_paths)
)



process.Analysis = cms.EDAnalyzer('Analysis',
    muons = cms.InputTag("muons"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    vertices = cms.InputTag("offlinePrimaryVertices"),  
    hltPath = cms.string("HLT_IsoMu24_v")
)

process.p = cms.Path(process.Analysis)
