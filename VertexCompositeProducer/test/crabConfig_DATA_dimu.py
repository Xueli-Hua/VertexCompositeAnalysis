# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
#from WMCore.Configuration import Configuration
#config = Configuration()

from datetime import datetime
from CRABAPI.RawCommand import crabCommand
from CRABClient.UserUtilities import config
config = config()
date = datetime.now().strftime('%y%m%d')
date_time = datetime.now().strftime('%y%m%d_%H%M%S')

config.section_("General")
config.General.workArea = 'crab_projects/EPPF'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'Cert_Collisions2023HI_374288_375823_Good_ZDC_Muon.json'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True
config.Data.allowNonValidInputDataset = True

config.section_('Site')
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T1_US_*','T2_US_*','T1_FR_*','T2_FR_*','T2_CH_CERN','T2_BE_IIHE']
config.Site.storageSite = 'T2_CH_CERN'

## Submit the muon PDs
#config.General.requestName = 'Jpsi_HIPhysicsRawPrime31_HIRun2023A-PromptRec_Cen40_0812'
#config.Data.inputDataset = '/HIPhysicsRawPrime31/HIRun2023A-PromptReco-v2/MINIAOD'
config.Data.unitsPerJob = 10
#config.Data.totalUnits = 10
config.JobType.maxMemoryMB = 2500
config.JobType.maxJobRuntimeMin = 2100
config.JobType.psetName = 'PbPbSkimAndTree2023_DiMuContBoth_ZDC_TrksvPFEP_MiniAOD.py'
#config.Data.outputDatasetTag = config.General.requestName 
config.Data.outLFNDirBase = '/store/group/phys_heavyions/xueli/HIPhysicsRawPrimeEPPF/EtCut0p01to30/'

#config.Data.runRange = '374961'
#config.Site.storageSite = 'T3_CH_CERNBOX'

## Submit PDs ###############################################################################
for i in range(0, 1):
    config.General.requestName = f'PFEt0p01to30_HIPhysicsRawPrime{i}_HIRun2023A-PromptRec_Cen40_'+ date_time
    config.Data.inputDataset = f'/HIPhysicsRawPrime{i}/HIRun2023A-PromptReco-v2/MINIAOD'
    config.Data.outputDatasetTag = config.General.requestName

    crabCommand('submit', config = config, dryrun=False)

print('='*50)
print('All jobs submitted.')
print(f'Output directory: {config.Data.outLFNDirBase}')
print('='*50)
