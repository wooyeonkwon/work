import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("AnalysisTnp")

# Configure the MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000  # Report every 10000 events

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


file_paths = ['file:///pnfs/knu.ac.kr/data/cms/store/user/wkwon/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_99.root'
              ]
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*file_paths)
)

# Enable multithreading 
process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0)
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string("data.root")
)

process.AnalysisTnp = cms.EDAnalyzer('AnalysisTnp',
    muons=cms.InputTag("muons"),
    triggerResults=cms.InputTag("TriggerResults", "", "HLT"),
    vertices=cms.InputTag("offlinePrimaryVertices"),
    hltPath=cms.string("HLT_IsoMu24_v")
)

process.p = cms.Path(process.AnalysisTnp)


process.options.wantSummary = True
