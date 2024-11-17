import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("skim")

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

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:mc_skimmed_data.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    ),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep recoMuons_muons__*',
        'keep recoVertexs_offlinePrimaryVertices__*',
        'keep edmTriggerResults_TriggerResults__HLT',
        'keep edmEventAuxiliary_*_*_*',
        #'keep recoTracks_generalTracks_*_*',
        #'keep recoTracks_globalMuons_*_*',
        #'keep recoTracks_standAloneMuons_*_*',
        'keep triggerTriggerEvent_hltTriggerSummaryAOD__*'
        # comment out, if real data
        #,'keep recoGenParticles_genParticles__*'
    )
)

process.p = cms.Path()
process.e = cms.EndPath(process.out)

#process.MessageLogger = cms.Service("MessageLogger",
#    destinations = cms.untracked.vstring('cout'),
#    cout = cms.untracked.PSet(
#        threshold = cms.untracked.string('INFO')
#    )
#)

process.options.wantSummary = True
