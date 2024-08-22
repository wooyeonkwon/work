from CRABClient.UserUtilities import config

config = config()

config.Site.blacklist = ['T2_US_MIT']

config.General.requestName = 'MuonAnalysis_Run2022E'
config.General.workArea = 'crab_projects' 

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'conf.py'

config.Data.inputDataset = '/Muon/Run2022E-27Jun2023-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting    = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.lumiMask = 'Cert_Collisions2022_355100_362760_Golden.json'
config.Data.outLFNDirBase = '/store/user/wkwon/AnalysisResults/'
config.Data.publication = False

config.Site.storageSite = 'T3_KR_KNU'
