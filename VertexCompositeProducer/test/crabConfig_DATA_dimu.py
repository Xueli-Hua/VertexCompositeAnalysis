# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from WMCore.Configuration import Configuration
#from CRABClient.UserUtilities import getUsername

config = Configuration()

config.section_("General")
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'Cert_Collisions2023HI_374288_375823_Muon.json'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True
config.Data.allowNonValidInputDataset = True

config.section_('Site')
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T1_US_*','T2_US_*','T1_FR_*','T2_FR_*','T2_CH_CERN','T2_BE_IIHE']
config.Site.storageSite = 'T2_CH_CERN'


dataMap = {
            "HIPhysicsRawPrime0": "/HIPhysicsRawPrime0/HIRun2023A-PromptReco-v2/MINIAOD",
            "HIPhysicsRawPrime1": "/HIPhysicsRawPrime1/HIRun2023A-PromptReco-v2/MINIAOD",
            "HIPhysicsRawPrime2": "/HIPhysicsRawPrime2/HIRun2023A-PromptReco-v2/MINIAOD",
          }

## Submit the muon PDs
config.General.requestName = 'Jpsi_HIPhysicsRawPrime16_HIRun2023A-PromptRec_20240617'
config.Data.inputDataset = '/HIPhysicsRawPrime16/HIRun2023A-PromptReco-v2/MINIAOD'
config.Data.unitsPerJob = 20
#config.Data.totalUnits = 200
config.JobType.maxMemoryMB = 2500
config.JobType.maxJobRuntimeMin = 2100
config.JobType.psetName = 'PbPbSkimAndTree2023_DiMuContBoth_ZDC_MiniAOD_cfg.py'
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/group/phys_heavyions/xueli/HIPhysicsRawPrime/' 
