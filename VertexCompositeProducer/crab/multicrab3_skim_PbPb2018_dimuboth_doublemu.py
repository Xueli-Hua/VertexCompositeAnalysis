if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    from CRABClient.UserUtilities import config, getUsername #FromSiteDB
    config = config()
    print "fortest"
    config.General.workArea = '/eos/cms/store/group/phys_heavyions/xueli/test/VertexCompositeAna'
    config.General.transferOutputs = True
    config.General.transferLogs = False
    config.JobType.pluginName = 'Analysis'
#    config.JobType.maxMemoryMB = 4000
#    config.JobType.maxJobRuntimeMin = 2750
    config.JobType.psetName = '../test/PbPbSkimAndTree2018_DiMuContBoth_cfg.py'
    config.JobType.inputFiles=['../test/HeavyIonRPRcd_PbPb2018_offline.db']
    config.Data.unitsPerJob = 100    # 20->5 One could also change the number of luminosity sections to analyze per job (Data.unitsPerJob); e.g. one could decrease it so that to have shorter jobs.
    config.Data.totalUnits = 100   # 100->5 totalUnits = NJobs*unitsPerJob
    config.Data.splitting = 'LumiBased'
#    config.Data.splitting = 'Automatic' (no unitsPerJob and totalUnits)
    config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON_MuonPhys.txt'
    #config.Data.lumiMask = 'Cert_326381-327489_HI_PromptReco_Collisions18_JSON_MuonPhys.txt'
#    config.Data.lumiMask = 'jsondiff_327560_327489.txt'
#    config.Data.lumiMask = 'json_DCSONLY_HI_327327.txt'
#    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsername())#FromSiteDB())
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/xueli'
    config.Data.publication = False
    config.Site.storageSite = 'T2_US_Vanderbilt'
#    config.Site.storageSite = 'T2_US_MIT'
#    config.Site.storageSite = 'T3_US_Rice'
#    config.Site.storageSite = 'T2_CH_CERN'
    print "fortest"
    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
            #print("Failed submitting task: %s" % hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)
            #print("Failed submitting task: %s" % cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

#    config.General.requestName = 'PbPbDiMu_v1_MP327564_v5'
#    config.Data.inputDataset = '/HIDoubleMuon/HIRun2018A-PromptReco-v1/AOD'
#    config.Data.outputDatasetTag = '2018Skimv1_DiMuCont_MuonPhysics_HLTMass7toInf_v5'
#    config.JobType.psetName = '../test/PbPbSkimAndTree2018_DiMuContBoth_cfg.py'
#    submit(config)
    print "fortest"
    config.General.requestName = 'PbPbDiMu_v2_MP327564'
    #config.Data.inputDataset = '/HIDoubleMuon/HIRun2018A-PromptReco-v2/AOD'
    config.Data.inputDataset = '/HIDoubleMuon/HIRun2018A-04Apr2019-v1/AOD'
    #config.Data.inputDBS = 'global'
    config.Data.outputDatasetTag = '2018Skimv2_DiMuCont_MuonPhysics_HLTMass2p5toInf_v2' # 7->2p5
    config.JobType.psetName = '../test/PbPbSkimAndTree2018_DiMuContBoth_cfg.py'
    submit(config)
    print "fortest"
#    config.General.requestName = 'PbPbDiMuPeri_v2_MP327560_v4'
#    config.Data.inputDataset = '/HIDoubleMuonPsiPeri/HIRun2018A-PromptReco-v2/AOD'
#    config.Data.outputDatasetTag = '2018Skimv2_DiMuCont_MuonPhysics_HLT40100OS_v4'
#    config.JobType.psetName = '../test/PbPbSkimAndTree2018_DiMuContBoth_peripheral_cfg.py'
#    submit(config)

#    config.General.requestName = 'PbPbDiMuUPC_v1_MP327564_v4'
#    config.Data.inputDataset = '/HIForward/HIRun2018A-PromptReco-v1/AOD'
#    config.Data.outputDatasetTag = '2018Skimv1_DiMuCont_MuonPhysics_DiMuOpen_v4'
#    config.JobType.psetName = '../test/PbPbSkimAndTree2018_DiMuContBoth_UPC_cfg.py'
#    submit(config)

#    config.General.requestName = 'PbPbDiMuUPC_v2_MP327564_v4'
#    config.Data.inputDataset = '/HIForward/HIRun2018A-PromptReco-v2/AOD'
#    config.Data.outputDatasetTag = '2018Skimv2_DiMuCont_MuonPhysics_DiMuOpen_v4'
#    submit(config)
