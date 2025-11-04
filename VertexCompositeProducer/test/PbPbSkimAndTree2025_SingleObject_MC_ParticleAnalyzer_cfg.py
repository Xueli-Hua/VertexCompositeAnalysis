import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM', eras.Run3_pp_on_PbPb_2025)

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
    fileNames = cms.untracked.vstring("root://xrootd-cms.infn.it//store/user/fdamas/HLT/PbPb2025/RunPrepMC/GammaGammaToDimuon_Starlight_1510pre6/RECO_151X_mcRun3_2025_realistic_HI_v1/250922_191258/0000/step3_RAW2DIGI_L1Reco_RECO_1.root")
)
#from GGTMM_List import inputFileNames
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFileNames))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('151X_mcRun3_2025_realistic_HI_v1')

#* Set ZDC information
#process.es_pool = cms.ESSource("PoolDBESSource",
#    timetype = cms.string('runnumber'),
#    toGet = cms.VPSet(cms.PSet(record = cms.string("HcalElectronicsMapRcd"), tag = cms.string("HcalElectronicsMap_2021_v2.0_data"))),
#    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
#    authenticationMethod = cms.untracked.uint32(1)
#)
#process.es_prefer = cms.ESPrefer('HcalTextCalibrations', 'es_ascii')
#process.es_ascii = cms.ESSource('HcalTextCalibrations',
#    input = cms.VPSet(cms.PSet(object = cms.string('ElectronicsMap'), file = cms.FileInPath("emap_2023_newZDC_v3.txt")))
#)

# Add the Particle producer
from VertexCompositeAnalysis.VertexCompositeProducer.generalParticles_cff import generalParticles

process.muons = generalParticles.clone(
    pdgId = cms.uint32(13),
    muons = cms.InputTag('patMuons'),
)

process.electrons = generalParticles.clone(
    pdgId = cms.uint32(11),
    electrons = cms.InputTag('patElectrons')
)

process.lowPtElectrons = generalParticles.clone(
    pdgId = cms.uint32(11),
    electrons = cms.InputTag('patLowPtElectrons')
)

process.photons = generalParticles.clone(
    pdgId = cms.uint32(22),
    photons = cms.InputTag('patPhotons')
)

process.convertedPhotons = generalParticles.clone(
    pdgId = cms.uint32(22),
    conversions = cms.InputTag('allConversions')
)

process.tracks = generalParticles.clone(
    tracks = cms.InputTag('generalTracks'),
    dEdxInputs = cms.vstring('dedxHarmonic2', 'dedxPixelHarmonic2')
)

process.pixelTracks = generalParticles.clone(
    tracks = cms.InputTag('hiConformalPixelTracks')
)

process.pfCandidates = generalParticles.clone(
    pfParticles = cms.InputTag('particleFlow'),
    tracks = cms.InputTag('')
)

# Add PAT objects
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import doPATMuons, doPATElectrons, doPATPhotons
doPATMuons(process)
doPATElectrons(process)
doPATPhotons(process)

process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.colEvtSel = cms.Sequence(process.hiClusterCompatibility)

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.colEvtSel *
    process.primaryVertexFilter
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

event_filter = cms.untracked.vstring(
        "Flag_colEvtSel",
        "Flag_clusterCompatibilityFilter",
        "Flag_primaryVertexFilter",
        "Flag_hfPosFilterNTh7",
        "Flag_hfPosFilterNTh7p3",
        "Flag_hfPosFilterNTh8",
        "Flag_hfPosFilterNTh10",
        "Flag_hfNegFilterNTh7",
        "Flag_hfNegFilterNTh7p6",
        "Flag_hfNegFilterNTh8",
        "Flag_hfNegFilterNTh10",
)

trig_info = cms.untracked.VPSet([
    # UPC triggers
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuCosmic_NotMBHF2AND_MaxPixelCluster1000_v*'), filter = cms.string('hltL1sSingleMuCosmicNotMBHF2AND'), minN = cms.int32(1)),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuCosmic_NotMBHF2OR_MaxPixelCluster1000_v*'), filter = cms.string('hltL1sSingleMuCosmicNotMBHF2OR'), minN = cms.int32(1)),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuOpen_NotMBHF2AND_MaxPixelCluster1000_v*'), filter = cms.string('hltL1sSingleMuOpenNotMBHF2AND'), minN = cms.int32(1)),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuOpen_NotMBHF2OR_MaxPixelCluster1000_v*'), filter = cms.string('hltL1sSingleMuOpenNotMBHF2OR'), minN = cms.int32(1)),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleEG3_NotMBHF2AND_SinglePixelTrack_MaxPixelTrack_v*'), filter = cms.string('hltL1sSingleEG3NotHF2AND'), minN = cms.int32(1)),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleEG3_NotMBHF2OR_SinglePixelTrack_MaxPixelTrack_v*'), filter = cms.string('hltL1sSingleEG3NotHF2OR'), minN = cms.int32(1)),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleEG5_NotMBHF2AND_SinglePixelTrack_MaxPixelTrack_v*'), filter = cms.string('hltL1sSingleEG5NotHF2AND'), minN = cms.int32(1)),
  ])

"""
trig_info = cms.untracked.VPSet([
    # UPC triggers
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuCosmic_NotMBHF2AND_MaxPixelCluster1000_v*'), minN = cms.int32(1)),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuCosmic_NotMBHF2OR_MaxPixelCluster1000_v*'), minN = cms.int32(1)),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuOpen_NotMBHF2AND_MaxPixelCluster1000_v*'), minN = cms.int32(1)),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleMuOpen_NotMBHF2OR_MaxPixelCluster1000_v*'), minN = cms.int32(1)),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleEG3_NotMBHF2AND_SinglePixelTrack_MaxPixelTrack_v*'), minN = cms.int32(1)),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleEG3_NotMBHF2OR_SinglePixelTrack_MaxPixelTrack_v*'), minN = cms.int32(1)),
    cms.PSet(path = cms.string('HLT_HIUPC_SingleEG5_NotMBHF2AND_SinglePixelTrack_MaxPixelTrack_v*'), minN = cms.int32(1)),
  ])
"""

# Define the analysis steps
process.muon_step = cms.Path(process.eventFilter_HM * process.patMuonSequence * process.muons)
process.electron_step = cms.Path(process.eventFilter_HM * process.patElectronSequence * process.electrons * process.lowPtElectrons)
process.photon_step = cms.Path(process.eventFilter_HM * process.patPhotonSequence * process.photons * process.convertedPhotons)
process.track_step = cms.Path(process.eventFilter_HM * process.tracks  * process.pixelTracks * process.pfCandidates)

# Add the Particle tree
from VertexCompositeAnalysis.VertexCompositeAnalyzer.particle_tree_cff import particleAna_mc
#from TrackingTools.TrackPropagator.default_cff import *

process.hltESPSteppingHelixPropagatorAlong = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ApplyRadX0Correction = cms.bool(True),
    AssumeNoMaterial = cms.bool(False),
    ComponentName = cms.string('hltESPSteppingHelixPropagatorAlong'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('alongMomentum'),
    SetVBFPointer = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    debug = cms.bool(False),
    endcapShiftInZNeg = cms.double(0.0),
    endcapShiftInZPos = cms.double(0.0),
    returnTangentPlane = cms.bool(True),
    sendLogWarning = cms.bool(False),
    useEndcapShiftsInZ = cms.bool(False),
    useInTeslaFromMagField = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    useMagVolumes = cms.bool(True),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(False)
)

process.muonAna = particleAna_mc.clone(
  recoParticles = cms.InputTag("muons"),
  selectEvents = cms.string(""),
  eventFilterNames = event_filter,
  addTrgObj = cms.untracked.bool(True),
  triggerInfo = trig_info,
  maxGenDeltaR = cms.untracked.double(0.03),
  maxGenDeltaPtRel = cms.untracked.double(0.5),
  #propToMuon = cms.untracked.bool(True),
  useTrack = cms.string('tracker'),
  useState = cms.string('atVertex'),
  useSimpleGeometry = cms.bool(True),
  useStation2 = cms.bool(True),
  fallbackToME1 = cms.bool(True),
  useMB2InOverlap = cms.bool(True),
  cosmicPropagationHypothesis = cms.bool(False),
  propagatorAlong = cms.ESInputTag('', 'hltESPSteppingHelixPropagatorAlong'),
  propagatorAny = cms.ESInputTag('', 'SteppingHelixPropagatorAny'),
  propagatorOpposite = cms.ESInputTag('', 'hltESPSteppingHelixPropagatorOpposite'),

)

# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('obj_ana_mcmu.root'))
process.p = cms.EndPath( process.muonAna)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.muon_step,
    process.electron_step,
    process.photon_step,
    process.track_step,
    process.p
)

process.Flag_colEvtSel = cms.Path(process.colEvtSel)
process.Flag_clusterCompatibilityFilter = cms.Path(process.eventFilter_HM * process.hiClusterCompatibility)
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter_HM * process.primaryVertexFilter)
process.Flag_hfPosFilterNTh7 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh7_seq)
process.Flag_hfPosFilterNTh7p3 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh7p3_seq)
process.Flag_hfPosFilterNTh8 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh8_seq)
process.Flag_hfPosFilterNTh10 = cms.Path(process.eventFilter_HM * process.hfPosFilterNTh10_seq)
process.Flag_hfNegFilterNTh7 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh7_seq)
process.Flag_hfNegFilterNTh7p6 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh7p6_seq)
process.Flag_hfNegFilterNTh8 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh8_seq)
process.Flag_hfNegFilterNTh10 = cms.Path(process.eventFilter_HM * process.hfNegFilterNTh10_seq)

eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_clusterCompatibilityFilter , process.Flag_primaryVertexFilter , process.Flag_hfPosFilterNTh7 , process.Flag_hfPosFilterNTh7p3 , process.Flag_hfPosFilterNTh8 , process.Flag_hfPosFilterNTh10 , process.Flag_hfNegFilterNTh7 , process.Flag_hfNegFilterNTh7p6 , process.Flag_hfNegFilterNTh8 , process.Flag_hfNegFilterNTh10 ]

for P in eventFilterPaths:
    process.schedule.insert(0, P)
