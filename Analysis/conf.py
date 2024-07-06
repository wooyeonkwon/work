import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022C-v1/2580000/14939243-8dba-4bcb-b6fc-2747d52953d2.root',
        "file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022C-v1/2580000/28ffedd5-5e5f-4369-87cb-e985331b09f5.root",
        "file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022C-v1/2580000/338fbd05-df4d-4fda-b9d8-b7a8d87b4f22.root",
        "file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022C-v1/2580000/3c9d32de-78e8-40c1-865f-6822c3ed025b.root",
        "file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022C-v1/2580000/4053989c-154d-4b5d-b6de-2fda2ad33236.root",
        "file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022C-v1/2580000/43fe20d4-beee-4809-9300-420701a15478.root",
        "file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022C-v1/2580000/5d784d2a-e07d-46f8-b6a3-d95fec88e8c9.root",
        "file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022C-v1/2580000/7f339389-b9c6-4bac-b218-eb312eba0c78.root",
        "file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022C-v1/2580000/8d5b1a99-f249-4fd4-9c1d-d28bf7288487.root",
        "file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022C-v1/2580000/92e31a2e-cddb-4f7f-85bd-7115d23ed9ec.root",
        "file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022C-v1/2580000/b54dda5c-cc07-4cd9-8237-ee5aabba6759.root",
        "file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022C-v1/2580000/b98a1122-0029-4bc9-9603-e69c8b21e973.root",
        "file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022C-v1/2580000/d21672c9-fcbb-45a0-a014-500e3ad49d92.root",
        "file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022C-v1/2580000/d38c009e-e4f9-4367-94ac-b6b371bf606a.root",
        "file:///pnfs/knu.ac.kr/data/cms/store/relval/CMSSW_14_0_0/DoubleMuon/AOD/140X_dataRun3_v3_STD_2022_Data_RelVal_2022C-v1/2580000/ed0ea63c-169c-4820-97e0-5af651d5dd7f.root"
        
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

