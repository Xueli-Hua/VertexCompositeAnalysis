import FWCore.ParameterSet.Config as cms
process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 2

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
#'/store/user/davidlw/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/nonprompt_pt1p2_y2p4_hi1031p1_RECO_v2/181112_180307/0000/step3_RAW2DIGI_L1Reco_RECO_955.root'
'/store/user/mnguyen/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_2018/crab_Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/181018_065404/0000/step3_pre5DIGI_RAW2DIGI_L1Reco_RECO_1.root',
'/store/user/mnguyen/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_2018/crab_Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/181018_065404/0000/step3_pre5DIGI_RAW2DIGI_L1Reco_RECO_2.root',
'/store/user/mnguyen/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_2018/crab_Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/181018_065404/0000/step3_pre5DIGI_RAW2DIGI_L1Reco_RECO_3.root',
'/store/user/mnguyen/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_2018/crab_Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/181018_065404/0000/step3_pre5DIGI_RAW2DIGI_L1Reco_RECO_4.root',
'/store/user/mnguyen/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_2018/crab_Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/181018_065404/0000/step3_pre5DIGI_RAW2DIGI_L1Reco_RECO_5.root',
'/store/user/mnguyen/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_2018/crab_Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/181018_065404/0000/step3_pre5DIGI_RAW2DIGI_L1Reco_RECO_6.root',
'/store/user/mnguyen/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_2018/crab_Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/181018_065404/0000/step3_pre5DIGI_RAW2DIGI_L1Reco_RECO_7.root',
'/store/user/mnguyen/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_2018/crab_Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/181018_065404/0000/step3_pre5DIGI_RAW2DIGI_L1Reco_RECO_8.root',
'/store/user/mnguyen/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_2018/crab_Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/181018_065404/0000/step3_pre5DIGI_RAW2DIGI_L1Reco_RECO_9.root',
'/store/user/mnguyen/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_2018/crab_Hydjet_Quenched_Drum5Ev8_PbPbMinBias_5020GeV_RECO_103X_upgrade2018_realistic_HI_v4/181018_065404/0000/step3_pre5DIGI_RAW2DIGI_L1Reco_RECO_10.root',
)
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(800))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = '103X_upgrade2018_realistic_HI_v6'

# =============== Import Sequences =====================
#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
#process.hltHM.HLTPaths = ['HLT_PAFullTracks_Multiplicity280_v*']
process.hltHM.HLTPaths = ['HLT_PAFullTracks_Multiplicity280_v*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 50 && position.Rho <= 5 && tracksSize >= 2"),
#    cut = cms.string("!isFake && abs(z) <= 1 && position.Rho <= 2 && tracksSize >= 5"),
    filter = cms.bool(True),   # otherwise it won't filter the events
)

#Reject beam scraping events standard pp configuration
process.NoScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

process.PAcollisionEventSelection = cms.Sequence(
                                         process.hfCoincFilter * 
                                         process.PAprimaryVertexFilter *
                                         process.NoScraping
                                         )

process.eventFilter_HM = cms.Sequence( 
#    process.hltHM *
    process.PAcollisionEventSelection
)

process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

#process.dEdx_step = cms.Path( process.eventFilter_HM * process.produceEnergyLoss )

########## D0 candidate rereco ###############################################################
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalD0Candidates_cff")
process.generalD0CandidatesNew = process.generalD0Candidates.clone()
process.generalD0CandidatesNew.tkPtSumCut = cms.double(2.1)
process.generalD0CandidatesNew.tkEtaDiffCut = cms.double(1.0)
process.generalD0CandidatesNew.tkNhitsCut = cms.int32(11)
process.generalD0CandidatesNew.tkPtErrCut = cms.double(0.1)
process.generalD0CandidatesNew.tkPtCut = cms.double(0.7)
process.generalD0CandidatesNew.alphaCut = cms.double(1.0)
process.generalD0CandidatesNew.alpha2DCut = cms.double(1.0)

process.generalD0CandidatesNewWrongSign = process.generalD0Candidates.clone(isWrongSign = cms.bool(True))

process.d0rereco_step = cms.Path( process.eventFilter_HM * process.generalD0CandidatesNew * process.generalD0CandidatesNewWrongSign )

###############################################################################################

#process.load("RiceHIG.Skim2013.ppanalysisSkimContentFull_cff")
#process.load("RiceHIG.Skim2013.ppanalysisSkimContentSlim_cff")
process.load("VertexCompositeAnalysis.VertexCompositeProducer.ppanalysisSkimContentD0_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('/eos/cms/store/group/phys_heavyions/flowcorr/PbPb_MC.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_HM_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('AOD')
    )
)

process.output_HM_step = cms.EndPath(process.output_HM)

process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.d0rereco_step,
    process.output_HM_step
)