import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022B-v1/2580000/0e19b0b1-5b63-4730-98c1-7b416e1d620f.root',
        'root://210.117.210.34//home/dndus0107/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022B-v1/2580000/0e19b0b1-5b63-4730-98c1-7b416e1d620f.root',
        'root://210.117.210.34//home/dndus0107/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022B-v1/2580000/43d9dc02-7182-4559-98da-934918643eae.root',
        'root://210.117.210.34//home/dndus0107/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022B-v1/2580000/7c61e272-d81c-4e56-acdf-0d492284a957.root'
        
    )
)

# Enable multithreading
process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),  # Number of threads
    numberOfStreams = cms.untracked.uint32(0)   # Number of streams
)

process.Analysis = cms.EDAnalyzer('Analysis',
    muons = cms.InputTag("muons")
)

process.p = cms.Path(process.Analysis)

