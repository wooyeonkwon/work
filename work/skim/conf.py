import FWCore.ParameterSet.Config as cms

process = cms.Process("skim")

# Configure the MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

# MC truth matching for muons
process.muonMCMatch = cms.EDProducer("MCMatcher",
    src         = cms.InputTag("muons"),
    matched     = cms.InputTag("genParticles"),
    mcPdgId     = cms.vint32(13),
    checkCharge = cms.bool(True),
    maxDeltaR   = cms.double(0.3),
    maxDPtRel   = cms.double(0.5),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(True)
)


# Skim filter
process.skim = cms.EDFilter('skim',
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    hltPath = cms.string("HLT_IsoMu24_v"),
    matchedMuons = cms.InputTag("muonMCMatch")
)

process.p = cms.Path(
    process.muonMCMatch +  # MC truth matching
    process.skim  # Trigger filter
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:skimmed_mc.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    ),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep recoMuons_muons__*',
        'keep recoMuons_muonMCMatch__*',  # Store matched muon information
        'keep edmEventAuxiliary_*_*_*',
        'keep recoGenParticles_genParticles__*',
        'keep genWeights_genWeight__*'
    )
)

process.e = cms.EndPath(process.out)

process.options.wantSummary = True
