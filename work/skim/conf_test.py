import FWCore.ParameterSet.Config as cms

process = cms.Process("skim")

# Configure the MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000  # Report every 10000 events

#directory_path = '/data1/users/dndus0107/AnalysisResults/Muon/crab_Run2022D-27Jun2023-v2/241231_153948/0000/'
#file_list = os.listdir(directory_path)
#file_paths = [f'file://{directory_path}{filename}' for filename in file_list]
file_paths = ["file:///data1/users/dndus0107/raw_data_samples/00656c16-fb80-415d-b50b-15fbe6824f22.root"]
#file_paths = ['file:///data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_1-1.root',
#              'file:///data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_1-2.root',
#              'file:///data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_1-3.root',
#              'file:///data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_1-4.root',
#              'file:///data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/skimmed_data_1-5.root'
#              ]

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(*file_paths)
)

# MC truth matching for muons
process.muonMCMatch = cms.EDProducer("MCMatcher",
    src         = cms.InputTag("muons"),       # reco muons
    matched     = cms.InputTag("genParticles"),  # gen particles
    mcPdgId     = cms.vint32(13),             # Muon PDG ID
    checkCharge = cms.bool(True),             # Check charge match
    maxDeltaR   = cms.double(0.3),            # ΔR threshold
    maxDPtRel   = cms.double(0.5),            # Relative pT difference
    resolveAmbiguities = cms.bool(True),      # One-to-one matching
    resolveByMatchQuality = cms.bool(True),   # Match by ΔR
    mcStatus    = cms.vint32(1)               # Match only stable particles
)

# Skim filter (trigger filtering only)
process.skim = cms.EDFilter('skim',
    triggerResults=cms.InputTag("TriggerResults", "", "HLT"),
    hltPath=cms.string("HLT_IsoMu24_v")  # No matchedMuons here
)

# Path configuration
process.p = cms.Path(
    process.muonMCMatch +  # MC truth matching
    process.skim           # Trigger filter
)

# Output module to store matched information
process.out = cms.OutputModule("PoolOutputModule",
    fileName=cms.untracked.string('file:skimmed_mc.root'),
    SelectEvents=cms.untracked.PSet(
        SelectEvents=cms.vstring('p')
    ),
    outputCommands=cms.untracked.vstring(
        'drop *',
        'keep recoMuons_muons__*',         # Store all reco muons
        'keep recoMuons_muonMCMatch__*',   # Store MC matching results
        'keep edmEventAuxiliary_*_*_*',
        'keep recoGenParticles_genParticles__*',  # Keep gen-level particles
        'keep genWeights_genWeight__*'
    )
)

process.e = cms.EndPath(process.out)

process.options.wantSummary = True