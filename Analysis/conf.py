import FWCore.ParameterSet.Config as cms
import os

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
    fileName = cms.string("data.root")
)

process.Analysis = cms.EDAnalyzer('Analysis',
    muons = cms.InputTag("muons"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    vertices = cms.InputTag("offlinePrimaryVertices"),  
    hltPath = cms.string("HLT_IsoMu24_v")
)

process.p = cms.Path(process.Analysis)
