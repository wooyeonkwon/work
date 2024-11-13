import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("skim")

# Configure the MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000  # Report every 10000 events

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
file_paths = ['file:///eos/cms/store/mc/Run3Summer22EEDRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v4/2810005/d995c0ec-3c59-4c6c-9349-ac0e6c0f59f8.root',
              'file:///eos/cms/store/mc/Run3Summer22EEDRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v4/2810002/5a18b409-966b-4ecf-9f5b-6864bad80de8.root',
              'file:///eos/cms/store/mc/Run3Summer22EEDRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v4/70003/7e974abe-f7ee-44a3-82fd-1f88c3cd6ac3.root'
              ]
#file_paths =  ['file:///eos/cms/store/mc/Run3Summer22DRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v1/70000/00775a3b-7386-44ac-8ade-e3aafbdc0261.root']
#directory_path= '/eos/cms/store/mc/Run3Summer22EEDRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v4/2810005/'
#directory_path = '/eos/cms/store/mc/Run3Summer22DRPremix/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v1/70000/'
#file_list = os.listdir(directory_path)
#file_paths = [f'file://{directory_path}{filename}' for filename in file_list]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*file_paths)
)

# Enable multithreading 
process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(16),
    numberOfStreams = cms.untracked.uint32(0)
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:skimmed_data_22EEmc.root'),
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
        ,'keep recoGenParticles_genParticles__*'
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