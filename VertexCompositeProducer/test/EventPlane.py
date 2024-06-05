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
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/hidata/HIRun2023A/HIPhysicsRawPrime0/MINIAOD/PromptReco-v2/000/374/668/00000/f9afd210-86f4-46c0-8a2a-05f48d59cc64.root'),
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


process.load('VertexCompositeAnalysis.VertexCompositeProducer.unpackedTracksAndVertices_cfi')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.unpackedMuons_cfi')
process.eventinfoana = cms.EDAnalyzer('ForEventPlane',
  vtxInputTag = cms.untracked.InputTag("offlinePrimaryVertices"),
  trkInputTag = cms.untracked.InputTag("generalTracks"),
  muInputTag = cms.untracked.InputTag("muons"),
  caloTowerInputTag = cms.untracked.InputTag("towerMaker"),

  isCentrality = cms.bool(True),
  centBinLabelTag = cms.InputTag("centralityBin","HFtowers"),
  centSrcTag = cms.InputTag("hiCentrality"),
)
'''
# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.clusterCompatibilityFilter_cfi')
process.colEvtSel = cms.Sequence(process.primaryVertexFilter * process.clusterCompatibilityFilter)
'''

# Define the analysis steps
process.pcentandep_step = cms.Path(process.cent_seq)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('testMBMINIAOD.root'))
process.p = cms.EndPath(process.eventinfoana)


# Define the process schedule
process.schedule = cms.Schedule(
    process.pcentandep_step,
    process.p
)
'''
# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_hfCoincFilter2Th4 = cms.Path(process.eventFilter_HM * process.hfCoincFilter2Th4)
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter_HM * process.primaryVertexFilter)
process.Flag_clusterCompatibilityFilter = cms.Path(process.eventFilter_HM * process.clusterCompatibilityFilter)
eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_primaryVertexFilter, process.Flag_clusterCompatibilityFilter ]

for P in eventFilterPaths:
    process.schedule.insert(0, P)
'''
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD
changeToMiniAOD(process)

