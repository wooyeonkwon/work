from CRABClient.UserUtilities import config

config = config()

#config.Site.blacklist = ['T2_US_MIT']

config.General.requestName = 'MuonAnalysis_TnP_Run2024G1'
config.General.workArea = 'crab_projects' 

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'conf.py'
config.JobType.numCores = 1
config.JobType.maxMemoryMB = 2000


config.Data.inputDataset = '/Muon1/Run2024G-PromptReco-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting    = 'FileBased'
config.Data.unitsPerJob = 100
config.Data.lumiMask = 'Cert_Collisions2024_378981_385194_Golden.json'
config.Data.outLFNDirBase = '/store/user/wkwon/AnalysisResults/'
config.Data.publication = False

config.Site.storageSite = 'T3_KR_KNU'
config.JobType.outputFiles = ['data.root']
