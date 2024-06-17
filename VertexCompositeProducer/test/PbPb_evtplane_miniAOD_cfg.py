import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM',eras.Run3_pp_on_PbPb_2023)

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
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/hidata/HIRun2023A/HIPhysicsRawPrime0/MINIAOD/PromptReco-v2/000/374/668/00000/f9afd210-86f4-46c0-8a2a-05f48d59cc64.root'
        #'root://cms-xrd-global.cern.ch//store/hidata/HIRun2023A/HIMinimumBias0/MINIAOD/PromptReco-v2/000/374/668/00000/99ceb28c-542e-4b2c-ac5d-bf44588e288d.root',
        #"root://cmsxrootd.fnal.gov//store/hidata/HIRun2023A/HIMinimumBias0/AOD/PromptReco-v1/000/373/870/00000/82016a82-f1a3-45d7-8090-26ce3f361113.root",
        #'root://cms-xrd-global.cern.ch//store/hidata/HIRun2023A/HIZeroBias0/MINIAOD/PromptReco-v2/000/375/695/00000/fdce424a-4085-45cd-bd49-f21bebe054dc.root',
        #'root://cms-xrd-global.cern.ch//store/hidata/HIRun2023A/HIMinimumBias0/MINIAOD/PromptReco-v2/000/374/668/00000/99ceb28c-542e-4b2c-ac5d-bf44588e288d.root',
    ),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('132X_dataRun3_Prompt_v7')

# Set ZDC information
process.load('VertexCompositeAnalysis.VertexCompositeProducer.QWZDC2018Producer_cfi')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.QWZDC2018RecHit_cfi')
process.zdcdigi.SOI = cms.untracked.int32(2)
process.es_pool = cms.ESSource("PoolDBESSource",
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(cms.PSet(record = cms.string("HcalElectronicsMapRcd"), tag = cms.string("HcalElectronicsMap_2021_v2.0_data"))),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    authenticationMethod = cms.untracked.uint32(1)
)
process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')
process.es_ascii = cms.ESSource('HcalTextCalibrations',
    input = cms.VPSet(cms.PSet(object = cms.string('ElectronicsMap'), file = cms.FileInPath("VertexCompositeAnalysis/VertexCompositeProducer/data/emap_2023_newZDC_v3.txt")))
)


# Add PbPb centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")
process.cent_seq = cms.Sequence(process.centralityBin)

process.pcentandep_step = cms.Path(process.cent_seq)


process.EvtPlane = cms.EDAnalyzer('PATEventPlane',
  doRecoNtuple = cms.untracked.bool(True),
  beamSpotSrc = cms.untracked.InputTag("offlineBeamSpot"),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  MuonCollection = cms.untracked.InputTag("muons"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  caloTowerInputTag = cms.untracked.InputTag("towerMaker"),

  isCentrality = cms.bool(True),
  centralityBinLabel = cms.InputTag("centralityBin","HFtowers"),
  centralitySrc = cms.InputTag("hiCentrality"),

  saveTree = cms.untracked.bool(True),
  saveHistogram = cms.untracked.bool(True)
)

process.p = cms.EndPath(process.EvtPlane)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('evtplane.root'))

# Define the process schedule
process.schedule = cms.Schedule(
    process.pcentandep_step,
    process.p
)

process.load('VertexCompositeAnalysis.VertexCompositeProducer.unpackedTracksAndVertices_cfi')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.unpackedMuons_cfi')
process.PATInput = cms.Sequence(process.unpackedTracksAndVertices*process.unpackedMuons)
process.PATInput_path = cms.Path(process.PATInput)

process.schedule.insert(0, process.PATInput_path)


from FWCore.ParameterSet.MassReplace import massReplaceInputTag as MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","unpackedTracksAndVertices")
process = MassReplaceInputTag(process,"generalTracks","unpackedTracksAndVertices")
process = MassReplaceInputTag(process,"muons","unpackedMuons")

