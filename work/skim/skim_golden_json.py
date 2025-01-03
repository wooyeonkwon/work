import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import os
process = cms.Process("ANALYSIS")

# Golden JSON 파일 로드
directory_path = '/data1/users/dndus0107/AnalysisResults/Muon/crab_Run2022E-27Jun2023-v1/241231_154236/0000/'
file_list = os.listdir(directory_path)
file_paths = [f'file://{directory_path}{filename}' for filename in file_list]


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000 

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

golden_json_path = "/home/dndus0107/CMSSW_14_0_19_patch2/src/work/skim/Cert_Collisions2022_355100_362760_Golden.json"
lumi_list = LumiList.LumiList(filename=golden_json_path).getCMSSWString().split(',')

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(0),
    numberOfStreams = cms.untracked.uint32(0)
)

# Process Source 설정
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(*file_paths),
    lumisToProcess = cms.untracked.VLuminosityBlockRange(*lumi_list)
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:/data1/users/dndus0107/AnalysisResults/Muon/crab_Run2022E-27Jun2023-v1/241231_154236/1000/skimmed_data.root'),
)

process.outpath = cms.EndPath(process.out)

process.options.wantSummary = True