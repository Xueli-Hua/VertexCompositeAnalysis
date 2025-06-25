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
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True

## Submit
config.General.requestName = 'coh_jpsi_dimu_DIGIRAW_2025_04_14_Lay2e3'

config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.inputDataset = '/STARLIGHT_5p36TeV_2024Run3/phys_heavyions-xueli-coh_jpsi_dimu_GENSIM_2025_04_11-f8c014688e587c1253dff6fa6fdf1d8c/USER' 
#config.Data.outputPrimaryDataset = 'STARLIGHT_5p36TeV_2024Run3'
config.Data.splitting = 'FileBased'
config.Data.publication = False
config.Data.allowNonValidInputDataset = True
config.Data.unitsPerJob = 1
#config.Data.totalUnits = 3000
#config.JobType.maxMemoryMB = 2500
#config.JobType.maxJobRuntimeMin = 2100
config.JobType.psetName = 'step1_DIGI_L1_DIGI2RAW_HLT_L1out.py'
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/group/phys_heavyions/xueli/MC/'

config.section_('Site')
#config.Site.storageSite = 'T3_CH_CERNBOX'
config.Site.storageSite = 'T2_CH_CERN'

