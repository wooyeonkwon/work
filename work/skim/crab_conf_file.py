from CRABClient.UserUtilities import config

config = config()

#config.Site.blacklist = ['T2_US_MIT']

config.General.requestName = 'MuonSkimming_Run2022E_sup'
config.General.workArea = 'crab_projects' 

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'conf.py'
config.JobType.numCores = 1
config.JobType.maxMemoryMB = 2000


config.Data.userInputFiles = [
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/320000/89e4b866-1c43-4fa7-85db-19bb80e070db.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/25410000/fa9ce585-64f1-4cd4-a6ac-b3e1f2b81f1f.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60000/65b493ad-30d6-4cb0-b1d9-717c735c7c26.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/25410000/2cea81ee-cf9a-4cd0-8910-340aff43e8b7.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/320000/5d0df6f4-27a3-41df-a495-a7ed6ff9f099.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/25410000/6410f66c-e8b9-42cb-8e26-faa3d57772ca.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/25410000/c4002ecc-360c-4ab3-bcdc-cbbfa3795dd1.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/25410000/73d8e4c3-8bfe-472c-8583-94f2c4251ef0.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60000/a61e3bf5-0675-4664-b1e8-310b004939ec.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/25410000/c5fecf8d-8524-48db-8935-12d0e3f72636.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/320000/d359cc9c-d8af-4724-9e57-ce17bf4e069f.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/25410000/d7e7e941-6272-489e-b3a0-bea407fcb54e.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/320001/c494bd7b-cc18-4675-bd96-d36b96d4d50a.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/320001/1205ffa0-2254-474f-a90d-6cb0109b088d.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60000/09cb0da5-a038-47f0-a2c3-e55ede200224.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/25410000/9c701c5e-f5c9-4eb0-a135-776c6ca8f609.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/25410000/bccb0287-0155-47e7-9ef7-a2f0f858b54c.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/25410000/dd007f25-c56d-4e0d-883d-53a2c911443a.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/320000/1d8d809d-fae5-41dd-bc3f-79d20301d157.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/25410000/19a12ecd-9627-461e-9667-474f0952ca7b.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60005/ce03d8dd-55cb-4899-b5a4-09b7a580a5b9.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60005/dfcc9f72-eb00-46d9-84f1-2ec82fd528ab.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/320000/282810ab-1fb8-46d5-a587-af6d78b956dd.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/320000/2acf0621-b923-416a-8513-110f9546a80c.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60005/4a157ca0-47c7-49ff-90e0-85638bb581bc.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60005/2204d32e-eac5-4542-868e-6ec2c0efb43f.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/320000/1dc98250-0252-4934-b23a-701e580e541a.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60005/e283bd6d-6cef-48c2-a15b-cfab6820f579.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60005/a62f6318-2964-4003-87da-1630e79caaab.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60005/1691e8b1-5bf4-4c30-8d4a-04048d58f8a5.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60005/816cd1ef-b101-4d3d-b7ac-a3664d9b4599.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60005/9d72c014-09fb-42d7-8043-2fc20122492a.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60005/8fa336e1-2c7a-4840-994c-279f1cb53ca4.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60005/51216214-b3f9-40f5-9794-83a60c4757bd.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60005/e20327ab-576b-4f93-bebd-8aac8ec3381f.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210001/fb60350a-6a06-427f-b65c-33439d675e75.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210000/f4f03f9d-ec7c-4c67-af4e-a8bf5c3c4eaa.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210001/ad86f8c5-33cf-4cf4-80e2-748dd70c1efd.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210002/eeb61afc-9f37-4546-8888-8096016cbd40.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210000/c7669ae0-be24-483e-af16-5d0f409f7ee8.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60003/5d435217-23bc-4a7d-a53d-da7b659baf14.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210000/a6ddd57c-101a-4154-832b-da2ee0835ea8.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210003/28a92ef5-4908-4ebf-bca6-5c597e9c7d62.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210002/5ff35525-fb19-48f9-a597-d82a1e18b360.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210002/d6102ec8-691d-4d02-bb70-8f04f09bde9e.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210001/73639574-9473-486a-bd5d-cf5f636ace97.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60003/146b74c1-63de-444a-ab44-018f6bbff969.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210001/28225e9c-4b4a-439a-9e3a-ccfe6e515646.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210002/7b430187-1c89-432f-a854-46173c2d939f.root'
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210002/7b430187-1c89-432f-a854-46173c2d939f.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210000/7a658316-6997-42ad-848c-adb1560d2ab2.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210001/85714d30-405d-436b-99db-b2f0d35332e7.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210000/d3846a79-a860-4a12-a856-371610e74ae5.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210004/2b5348bc-73ee-4f60-8c9f-c691c40110bb.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210000/61a07451-30d1-4514-82a4-898eae7e60b6.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210000/789bf0d6-c0ab-4f0d-bfa9-f6c841124d33.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210001/d126bbf1-be4c-4b35-b5ae-6190a0966b0d.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210003/37a41580-4c02-412e-8e91-ce33b321d369.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210003/96847165-5a1a-4e79-8397-9c8794982d7a.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210002/a44e0527-9d9e-4c6e-a0e7-f7f738267492.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210002/84536a0e-c567-4c9e-ad71-9ca87e5ae9d3.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210001/b71df8a5-d018-4039-a52e-e29e43bf3d0c.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210001/b21c35de-1a5a-46be-9776-95320c1b10fe.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210003/65bf6e85-b42f-4a63-ae7d-f989ef9f7663.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210002/cc189f1b-d870-49b9-a9e1-78adb657107a.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210003/67efdee3-2d55-4f98-9168-b91aa5cbff60.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210001/fb60350a-6a06-427f-b65c-33439d675e75.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210000/f4f03f9d-ec7c-4c67-af4e-a8bf5c3c4eaa.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210001/ad86f8c5-33cf-4cf4-80e2-748dd70c1efd.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210002/eeb61afc-9f37-4546-8888-8096016cbd40.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210000/c7669ae0-be24-483e-af16-5d0f409f7ee8.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60003/5d435217-23bc-4a7d-a53d-da7b659baf14.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210000/a6ddd57c-101a-4154-832b-da2ee0835ea8.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210003/28a92ef5-4908-4ebf-bca6-5c597e9c7d62.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210002/5ff35525-fb19-48f9-a597-d82a1e18b360.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210002/d6102ec8-691d-4d02-bb70-8f04f09bde9e.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210001/73639574-9473-486a-bd5d-cf5f636ace97.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/60003/146b74c1-63de-444a-ab44-018f6bbff969.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210001/28225e9c-4b4a-439a-9e3a-ccfe6e515646.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210002/7b430187-1c89-432f-a854-46173c2d939f.root',
    '/store/data/Run2022E/Muon/AOD/27Jun2023-v1/28210000/7a658316-6997-42ad-848c-adb1560d2ab2.root'
]
config.Data.inputDBS = 'global'
config.Data.splitting    = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.lumiMask = 'Cert_Collisions2022_355100_362760_Golden.json'
#config.Data.lumiMask = 'Cert_Collisions2024_378981_385194_Golden.json'
config.Data.outLFNDirBase = '/store/user/wkwon/AnalysisResults/'
config.Data.publication = False

config.Site.storageSite = 'T3_KR_KNU'
config.JobType.outputFiles = ['skimmed_data.root']
