import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Analysis")

# Configure the MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000  # Report every 10000 events

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#file_paths =  ['file:///eos/cms/store/mc/Run3Summer22DRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v1/70000/00775a3b-7386-44ac-8ade-e3aafbdc0261.root',
#               'file:///eos/cms/store/mc/Run3Summer22DRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v1/70000/00656c16-fb80-415d-b50b-15fbe6824f22.root',
#               'file:///eos/cms/store/mc/Run3Summer22DRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v1/70000/0066b911-997b-4e7a-b322-442102f0a30e.root']
file_paths = ['file:///pnfs/knu.ac.kr/data/cms/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210001/b71df8a5-d018-4039-a52e-e29e43bf3d0c.root',
              'file:///pnfs/knu.ac.kr/data/cms/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210002/7b430187-1c89-432f-a854-46173c2d939f.root',
              'file:///pnfs/knu.ac.kr/data/cms/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210002/eeb61afc-9f37-4546-8888-8096016cbd40.root'
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
    fileName = cms.string("data_2022E_3.root")
)

process.Analysis = cms.EDAnalyzer('Analysis',
    muons = cms.InputTag("muons"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    vertices = cms.InputTag("offlinePrimaryVertices"),  
    hltPath = cms.string("HLT_IsoMu24_v")
)

process.p = cms.Path(process.Analysis)

process.options.wantSummary = True
