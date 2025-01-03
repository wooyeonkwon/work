import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("AnalysisMC")
# Configure the MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # Report every 10000 events
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
directory_path = '/data1/users/dndus0107/AnalysisResults/Muon/crab_Run2022D-27Jun2023-v2/241231_153948/1000/'
file_list = os.listdir(directory_path)
file_paths = [f'file://{directory_path}{filename}' for filename in file_list]
#file_paths = ['file:///data1/users/dndus0107/public/skimmed_data_2022E.root']

#file_paths = ['file:///data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_1-1.root',
#              'file:///data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_1-2.root',
#              'file:///data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_1-3.root',
#              'file:///data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_1-4.root',
#              'file:///data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_1-5.root'
#              ]
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*file_paths)
)

# Enable multithreading 
process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("/data1/users/dndus0107/AnalysisResults/processed_data/Analysis_Data_22D.root")
)

process.AnalysisMC = cms.EDAnalyzer('AnalysisMC',
    muons = cms.InputTag("muons")
)
process.p = cms.Path(process.AnalysisMC)  
process.options.wantSummary = True
