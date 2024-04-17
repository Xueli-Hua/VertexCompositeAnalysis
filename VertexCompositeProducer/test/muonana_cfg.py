import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM', eras.Run3_2023)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring("file:/eos/cms/store/group/phys_heavyions/anstahll/CERN/PbPb2023/SKIM/SKIM_AOD_HIForward0_HIRun2023A_20231009/HIForward0/SKIM_AOD_HIForward0_HIRun2023A_20231009/231009_081732/0000/reco_RAW2DIGI_L1Reco_RECO_89.root"),
    fileNames = cms.untracked.vstring("root://cmsxrootd.fnal.gov//store/hidata/HIRun2023A/HIMinimumBias0/AOD/PromptReco-v1/000/373/870/00000/82016a82-f1a3-45d7-8090-26ce3f361113.root"),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('132X_dataRun3_Prompt_v4')

# Set ZDC information
process.es_pool = cms.ESSource("PoolDBESSource",
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(cms.PSet(record = cms.string("HcalElectronicsMapRcd"), tag = cms.string("HcalElectronicsMap_2021_v2.0_data"))),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    authenticationMethod = cms.untracked.uint32(1)
)
process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')
process.es_ascii = cms.ESSource('HcalTextCalibrations',
    input = cms.VPSet(cms.PSet(object = cms.string('ElectronicsMap'), file = cms.FileInPath("/eos/cms/store/group/phys_heavyions/xueli/test/VertexComposite2024/CMSSW_13_2_6_patch2/src/VertexCompositeAnalysis/VertexCompositeProducer/test/emap_2023_newZDC_v3.txt")))
)

# Add PbPb centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289"),
        connect = cms.string("sqlite_file:CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run3v1302x04_offline_374289.db"),
        label = cms.untracked.string("HFtowers")
        )
    ]
)
process.cent_seq = cms.Sequence(process.centralityBin)



process.eventinfoana = cms.EDAnalyzer('EventInfoTreeProducer',
  vtxInputTag = cms.untracked.InputTag("offlinePrimaryVertices"),
  trkInputTag = cms.untracked.InputTag("generalTracks"),
  muInputTag = cms.untracked.InputTag("muons"),
  caloTowerInputTag = cms.untracked.InputTag("towerMaker"),

  isCentrality = cms.bool(True),
  centBinLabelTag = cms.InputTag("centralityBin","HFtowers"),
  centSrcTag = cms.InputTag("hiCentrality"),
)

# Add PbPb collision event selection
process.load('JpsiAnalysis.muAnalyzer.collisionEventSelection_cff')
process.load('JpsiAnalysis.muAnalyzer.clusterCompatibilityFilter_cfi')
process.load('JpsiAnalysis.muAnalyzer.hfCoincFilter_cff')
process.colEvtSel = cms.Sequence(process.hfCoincFilter2Th4 * process.primaryVertexFilter*process.clusterCompatibilityFilter)

# Define the analysis steps
process.pcentandep_step = cms.Path(process.cent_seq)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('test.root'))
process.p = cms.EndPath(process.eventinfoana)

# Define the process schedule
process.schedule = cms.Schedule(
    process.pcentandep_step,
    process.p
)

# Add the event selection filters
#process.Flag_colEvtSel = cms.Path(process.colEvtSel)
#process.Flag_hfCoincFilter2Th4 = cms.Path(process.hfCoincFilter2Th4)
#process.Flag_primaryVertexFilter = cms.Path(process.primaryVertexFilter)
#process.Flag_clusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
#eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter2Th4 , process.Flag_primaryVertexFilter , process.Flag_clusterCompatibilityFilter ]
#for P in eventFilterPaths:
#    process.schedule.insert(0, P)
