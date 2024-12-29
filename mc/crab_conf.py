from CRABClient.UserUtilities import config

config = config()

#config.Site.blacklist = ['T2_US_MIT']

config.General.requestName = 'MuonAnalysis_DYJetsToLL_Run3Summer22EEDR'
config.General.workArea = 'crab_projects' 

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'conf.py'
config.JobType.numCores = 1
config.JobType.maxMemoryMB = 2000

#['/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v2/AODSIM',
#'/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM',
#'/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer22EEDRPremix-forPOG_124X_mcRun3_2022_realistic_postEE_v1-v3/AODSIM',
#'/QCD_PT-15to7000_TuneCP5_Flat2022_13p6TeV_pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12-v3/AODSIM',
#'/QCD_PT-15to7000_TuneCP5_Flat2022_13p6TeV_pythia8/Run3Summer22EEDRPremix-124X_mcRun3_2022_realistic_postEE_v1-v3/AODSIM',
#'/QCD_PT-15to7000_TuneCP5_13p6TeV_pythia8/Run3Summer22EEDRPremix-castor_124X_mcRun3_2022_realistic_postEE_v1-v2/AODSIM'
#]
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer22EEDRPremix-forPOG_124X_mcRun3_2022_realistic_postEE_v1-v3/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting    = 'FileBased'
config.Data.unitsPerJob = 100
config.Data.outLFNDirBase = '/store/user/wkwon/AnalysisResults/'
config.Data.publication = False

config.Site.storageSite = 'T3_KR_KNU'
config.JobType.outputFiles = ['mc.root']