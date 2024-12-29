import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("AnalysisMC")

# Configure the MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")

# Reduce the verbosity of the MessageLogger
process.MessageLogger.cerr.FwkReport.reportEvery = 10000  # Report every 10000 events

# Optional: Suppress all other messages except for errors
#process.MessageLogger.cerr.threshold = 'ERROR'
#process.MessageLogger.cerr.default.limit = 10

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

directory_path = '/eos/cms/store/mc/Run3Summer22DRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v1/70000/'
file_list = os.listdir(directory_path)
#file_paths = [f'file://{directory_path}{filename}' for filename in file_list]
#file_paths =  ['file:///eos/cms/store/mc/Run3Summer22DRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v1/70000/00775a3b-7386-44ac-8ade-e3aafbdc0261.root',
#               'file:///eos/cms/store/mc/Run3Summer22DRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v1/70000/00656c16-fb80-415d-b50b-15fbe6824f22.root',
#               'file:///eos/cms/store/mc/Run3Summer22DRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v1/70000/0066b911-997b-4e7a-b322-442102f0a30e.root']
#file_paths for QCD
file_paths = ['file:///eos/cms/store/mc/Run3Winter24Reco/QCD_PT-15to7000_TuneCP5_Flat2022_13p6TeV_pythia8/AODSIM/EpsilonPU_133X_mcRun3_2024_realistic_v9-v2/2550000/fe84c63e-1bf1-44bb-b647-68c0fa9c89c4.root',
              'file:///eos/cms/store/mc/Run3Winter24Reco/QCD_PT-15to7000_TuneCP5_Flat2022_13p6TeV_pythia8/AODSIM/EpsilonPU_133X_mcRun3_2024_realistic_v9-v2/2550000/fc55c85e-8a47-4161-84ff-6504088cd3e7.root',
              'file:///eos/cms/store/mc/Run3Winter24Reco/QCD_PT-15to7000_TuneCP5_Flat2022_13p6TeV_pythia8/AODSIM/EpsilonPU_133X_mcRun3_2024_realistic_v9-v2/2550000/f93c98ae-51de-49d3-8b94-be36d0438a96.root']

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*file_paths)
)

# Enable multithreading
process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(8),  # Number of threads
    numberOfStreams = cms.untracked.uint32(0),  # Number of streams
#    sizeOfStackForThreadsInKB = cms.untracked.uint32(10*1024)  # Stack size per thread in KB
)

# Use TFileService to manage output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("mc_qcd_data.root")
)

process.Analysis = cms.EDAnalyzer('AnalysisMC',
    muons = cms.InputTag("muons"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    vertices = cms.InputTag("offlinePrimaryVertices"),  
    hltPath = cms.string("HLT_IsoMu24_v"),
    genEventInfo = cms.InputTag("generator"),
    genMuons = cms.InputTag("genParticles"),
    lheEvent = cms.InputTag("externalLHEProducer")
)

process.p = cms.Path(process.Analysis)

process.options.wantSummary = True
