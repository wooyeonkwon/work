import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Analysis")

# Configure the MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")

# Reduce the verbosity of the MessageLogger
process.MessageLogger.cerr.FwkReport.reportEvery = 10000  # Report every 10000 events

# Optional: Suppress all other messages except for errors
#process.MessageLogger.cerr.threshold = 'ERROR'
#process.MessageLogger.cerr.default.limit = 10

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

directory_path = '/eos/cms/store/relval/CMSSW_14_0_7/Muon0/AOD/140X_dataRun3_Prompt_v2_2024C_HCALCalibChecks_RelVal_2024C-v2/2580000/'
file_list = os.listdir(directory_path)
file_paths = [f'file://{directory_path}{filename}' for filename in file_list]
#file_names_string = ',\n'.join([f'"{file_path}"' for file_path in file_paths])

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*file_paths)
)

# Enable multithreading
process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),  # Number of threads
    numberOfStreams = cms.untracked.uint32(0)   # Number of streams
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("data3.root")
)

process.Analysis = cms.EDAnalyzer('Analysis',
    muons = cms.InputTag("muons")
)

process.p = cms.Path(process.Analysis)
