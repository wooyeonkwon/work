import FWCore.ParameterSet.Config as cms

process = cms.Process("skim")

# Configure the MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000  # Report every 10000 events

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
file_paths = ['file:///data1/users/dndus0107/public/bf1be0eb-6045-492f-a31f-03972d22b8af.root'
              ]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*file_paths)
)

# Enable multithreading 
process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(1)
)

process.skim = cms.EDFilter('skim',
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    hltPath = cms.string("HLT_IsoMu24_v")
)

process.p = cms.Path(process.skim)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:skimmed_data_22EE.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    ),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep recoMuons_muons__*',
        'keep edmEventAuxiliary_*_*_*',
        # comment out, if real data
        #,'keep recoGenParticles_genParticles__*'
    )
)

process.e = cms.EndPath(process.out)

process.options.wantSummary = True
