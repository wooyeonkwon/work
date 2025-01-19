from CRABClient.UserUtilities import config

config = config()

#config.Site.blacklist = ['T2_US_MIT']

config.General.requestName = 'Run2022G-PromptReco-v1'
config.General.workArea = 'crab_projects_2' 

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'conf.py'
config.JobType.numCores = 1
config.JobType.maxMemoryMB = 2000


config.Data.inputDataset = '/Muon/Run2022G-PromptReco-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting    = 'Automatic'
config.Data.unitsPerJob = 1800
config.Data.lumiMask = 'Cert_Collisions2022_355100_362760_Golden.json'
config.Data.outLFNDirBase = '/store/user/wkwon/AnalysisResults/'
config.Data.publication = False

config.Site.storageSite = 'T3_KR_KNU'
config.JobType.outputFiles = ['skimmed_data.root']
