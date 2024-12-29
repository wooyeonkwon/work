import FWCore.ParameterSet.Config as cms
import os
import gc
gc.collect()
process = cms.Process("Analysis")
# Configure the MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000 # Report every 10000 events
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
directory_path = '/data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/'
file_list = os.listdir(directory_path)
#file_paths = [f'file://{directory_path}{filename}' for filename in file_list]
file_paths = ['file:///data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_1-1.root',
              'file:///data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_1-2.root',
              'file:///data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_1-3.root',
              'file:///data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_1-4.root',
              'file:///data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_1-5.root'
              ]
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*file_paths)
)

# Enable multithreading 
process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(0),
    numberOfStreams = cms.untracked.uint32(2)
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("/data1/users/dndus0107/AnalysisResults/processed_data/data_C.root")
)

process.Analysis = cms.EDAnalyzer('Analysis',
    muons = cms.InputTag("muons"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    vertices = cms.InputTag("offlinePrimaryVertices"),  
    hltPath = cms.string("HLT_IsoMu24_v")
)
process.p = cms.Path(process.Analysis)
process.options.wantSummary = True
