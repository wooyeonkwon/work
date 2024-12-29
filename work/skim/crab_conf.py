from CRABClient.UserUtilities import config

config = config()

#config.Site.blacklist = ['T2_US_MIT']

config.General.requestName = 'Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v4'
config.General.workArea = 'crab_projects' 

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'conf.py'
config.JobType.numCores = 1
config.JobType.maxMemoryMB = 2000


config.Data.inputDataset = '/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v4/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting    = 'Automatic'
config.Data.unitsPerJob = 1800
#config.Data.lumiMask = 'Cert_Collisions2024_378981_385194_Golden.json'
config.Data.outLFNDirBase = '/store/user/wkwon/AnalysisResults/'
config.Data.publication = False

config.Site.storageSite = 'T3_KR_KNU'
config.JobType.outputFiles = ['mc_skimmed_data.root']
