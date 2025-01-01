import FWCore.ParameterSet.Config as cms

process = cms.Process("skim")

# Configure the MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000  # Report every 10000 events

directory_path = '/data1/users/dndus0107/AnalysisResults/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v4/241117_152227/0000/'
file_list = os.listdir(directory_path)
file_paths = [f'file://{directory_path}{filename}' for filename in file_list]
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
    fileName = cms.untracked.string('file:/data1/users/dndus0107/AnalysisResults/skimmed_data/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/skimmed_mc22EE.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    ),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep recoMuons_muons__*',
        'keep edmEventAuxiliary_*_*_*',
        # comment out, if real data
        'keep recoGenParticles_genParticles__*'
    )
)

process.e = cms.EndPath(process.out)

process.options.wantSummary = True
