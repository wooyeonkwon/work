import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

# Configure the MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000  # Report every 10000 events

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

# Enable multithreading
process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("data_C.root")
)

process.Analysis = cms.EDAnalyzer('Analysis',
    muons = cms.InputTag("muons")
)

process.p = cms.Path(process.Analysis)
