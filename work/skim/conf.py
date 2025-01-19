import FWCore.ParameterSet.Config as cms

process = cms.Process("skim")

# Configure the MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring()
)

# MC truth matching for muons
#process.muonMCMatch = cms.EDProducer("MCMatcher",
#    src         = cms.InputTag("muons"),       # reco muons
#    matched     = cms.InputTag("genParticles"),  # gen particles
#    mcPdgId     = cms.vint32(13),             # Muon PDG ID
#    checkCharge = cms.bool(True),             # Check charge match
#    maxDeltaR   = cms.double(0.3),            # ΔR threshold
#    maxDPtRel   = cms.double(0.5),            # Relative pT difference
#    resolveAmbiguities = cms.bool(True),      # One-to-one matching
#    resolveByMatchQuality = cms.bool(True),   # Match by ΔR
#    mcStatus    = cms.vint32(1)               # Match only stable particles
#)

# Skim filter (trigger filtering only)
process.skim = cms.EDFilter('skim',
    triggerResults=cms.InputTag("TriggerResults", "", "HLT"),
    hltPath=cms.string("HLT_IsoMu24_v")  # No matchedMuons here
)

# Path configuration
process.p = cms.Path(
#    process.muonMCMatch +  # MC truth matching
    process.skim           # Trigger filter
)

# Output module to store matched information
process.out = cms.OutputModule("PoolOutputModule",
    fileName=cms.untracked.string('file:skimmed_data.root'),
    SelectEvents=cms.untracked.PSet(
        SelectEvents=cms.vstring('p')
    ),
    outputCommands=cms.untracked.vstring(
        'drop *',
        'keep recoMuons_muons__*',         # Store all reco muons
        'keep recoMuons_muonMCMatch__*',   # Store MC matching results
        'keep edmEventAuxiliary_*_*_*',
        # Keep gen-level particles
#        'keep recoGenParticles_genParticles__*',  
#        'keep genWeights_genWeight__*',
#        'keep recoGenParticlesedmAssociation_muonMCMatch__skim*',
#        'keep GenEventInfoProduct_generator__*'
    )
)

process.e = cms.EndPath(process.out)

process.options.wantSummary = True