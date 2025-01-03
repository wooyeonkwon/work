import FWCore.ParameterSet.Config as cms

process = cms.Process("skim")

# Configure the MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000  # Report every 10000 events

directory_path = '/data1/users/dndus0107/AnalysisResults/Muon/crab_Run2022D-27Jun2023-v2/241231_153948/0000/'
file_list = os.listdir(directory_path)
#file_paths = [f'file://{directory_path}{filename}' for filename in file_list]
file_paths = ["file:///data1/users/dndus0107/public/00656c16-fb80-415d-b50b-15fbe6824f22.root"]
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
    numberOfThreads = cms.untracked.uint32(0),
    numberOfStreams = cms.untracked.uint32(0)
)

process.skim = cms.EDFilter('skim',
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    hltPath = cms.string("HLT_IsoMu24_v")
)

process.p = cms.Path(process.skim)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test_skimmed_data.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    ),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep recoMuons_muons__*',
        'keep edmEventAuxiliary_*_*_*'
        # comment out, if real data
        ,'keep recoGenParticles_genParticles__*'
        ,'keep GenEventInfoProduct_generator__*'
    )
)

process.e = cms.EndPath(process.out)

process.options.wantSummary = True
