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
config.JobType.pluginName = 'PrivateMC'
config.JobType.allowUndistributedCMSSW = True

## Submit
config.General.requestName = 'coh_jpsi_dimu_GENSIM_2025_04_11'

config.section_('Data')
config.Data.outputPrimaryDataset = 'STARLIGHT_5p36TeV_2024Run3'
config.Data.splitting = 'EventBased'
config.Data.publication = True
config.Data.allowNonValidInputDataset = True
config.Data.unitsPerJob = 5000
config.Data.totalUnits = 9000000
config.JobType.maxMemoryMB = 2500
config.JobType.maxJobRuntimeMin = 2100
config.JobType.psetName = 'STARLIGHT_coh_jpsi_dimu_LHE_GEN_SIM.py'
config.JobType.inputFiles = ['/eos/user/x/xueli/MC/CMSSW_15_0_2/src/starlight_coherent_jpsi_dimuon_el8_amd64_gcc12_CMSSW_15_0_2_tarball.tgz']
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/group/phys_heavyions/xueli/MC/'

config.section_('Site')
#config.Site.storageSite = 'T3_CH_CERNBOX'
config.Site.storageSite = 'T2_CH_CERN'

