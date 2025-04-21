import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM',eras.Run3_2024_UPC)

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
        'root://cms-xrd-global.cern.ch//store/group/phys_heavyions/xueli/MC/STARLIGHT_5p36TeV_2024Run3/coh_jpsi_dimu_RECO_2025_04_11/250411_104338/0000/AODSIM_1.root',
        'root://cms-xrd-global.cern.ch//store/group/phys_heavyions/xueli/MC/STARLIGHT_5p36TeV_2024Run3/coh_jpsi_dimu_RECO_2025_04_11/250411_104338/0000/AODSIM_2.root',
        'root://cms-xrd-global.cern.ch//store/group/phys_heavyions/xueli/MC/STARLIGHT_5p36TeV_2024Run3/coh_jpsi_dimu_RECO_2025_04_11/250411_104338/0000/AODSIM_3.root',
        'root://cms-xrd-global.cern.ch//store/group/phys_heavyions/xueli/MC/STARLIGHT_5p36TeV_2024Run3/coh_jpsi_dimu_RECO_2025_04_11/250411_104338/0000/AODSIM_4.root',
        'root://cms-xrd-global.cern.ch//store/group/phys_heavyions/xueli/MC/STARLIGHT_5p36TeV_2024Run3/coh_jpsi_dimu_RECO_2025_04_11/250411_104338/0000/AODSIM_5.root',
    )
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('141X_mcRun3_2024_realistic_HI_v14')

# Add PbPb centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")
process.cent_seq = cms.Sequence(process.centralityBin)

# Add the VertexComposite producer
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cff")
process.generalMuMuMassMin0CandidatesWrongSign = process.generalMuMuMassMin0Candidates.clone(isWrongSign = cms.bool(True))
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import doPATMuons
doPATMuons(process, True)

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    # Double muon triggers
    'HLT_HIL2DoubleMu0_M1p5to6_Open_v*', # 
    'HLT_HIL2DoubleMu0_M7to15_Open_v*', # 
    'HLT_HIL3DoubleMu0_M2to4p5_Open_v*', # 
    'HLT_HIL3DoubleMu0_M7to15_Open_v*', # 
    'HLT_HIL3DoubleMu0_Quarkonia_Open_v*', # 
    # Single muon triggers
    'HLT_HIL3SingleMu12_v*', # 
    'HLT_HIUPC_SingleMuCosmic_NotMBHF2AND_v*', #
    'HLT_HIUPC_SingleMuCosmic_NotMBHF2OR_v*', #
    'HLT_HIUPC_SingleMuOpen_NotMBHF2AND_v*', #
    'HLT_HIUPC_SingleMuOpen_NotMBHF2OR_v*', #
    # Minimum bias triggers
    'HLT_HIMinimumBiasHF1ANDZDC1nOR_v*', # 
    'HLT_HIMinimumBiasHF1AND_v*', # 
    # Zero bias triggers
    'HLT_HIUPC_ZeroBias_SinglePixelTrack_MaxPixelTrack_v*', # 
    'HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v*', # 
    'HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v*', # 
    'HLT_HIZeroBias_HighRate_v*', # 
    # UPC
    'HLT_HIUPC_Random_HighRate_v*',
    ]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.clusterCompatibilityFilter_cfi')
process.colEvtSel = cms.Sequence(process.primaryVertexFilter * process.clusterCompatibilityFilter)

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.hltFilter
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

# Define the analysis steps
process.pcentandep_step = cms.Path(process.eventFilter_HM * process.cent_seq)# * process.evtplane_seq)
process.dimurereco_step = cms.Path(process.eventFilter_HM * process.patMuonSequence * process.generalMuMuMassMin0Candidates)
process.dimurerecowrongsign_step = cms.Path(process.eventFilter_HM * process.patMuonSequence * process.generalMuMuMassMin0CandidatesWrongSign)

# Add the VertexComposite tree
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_tree_cff")
process.dimucontana_mc.selectEvents = cms.untracked.string("eventFilter_HM_step")
process.dimucontana_mc.VertexCompositeCollection = cms.untracked.InputTag("generalMuMuMassMin0Candidates:DiMu")
process.dimucontana_wrongsign_mc = process.dimucontana_mc.clone(VertexCompositeCollection = cms.untracked.InputTag("generalMuMuMassMin0CandidatesWrongSign:DiMu"))


# Define the output
process.TFileService = cms.Service("TFileService", fileName = cms.string('dimuanamc_Lay4e3.root'))
process.p = cms.EndPath(process.dimucontana_mc * process.dimucontana_wrongsign_mc)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.pcentandep_step,
    process.dimurereco_step,
    process.dimurerecowrongsign_step,
    process.p
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.eventFilter_HM * process.colEvtSel)
process.Flag_hfCoincFilter2Th4 = cms.Path(process.eventFilter_HM * process.hfCoincFilter2Th4)
process.Flag_primaryVertexFilter = cms.Path(process.eventFilter_HM * process.primaryVertexFilter)
process.Flag_clusterCompatibilityFilter = cms.Path(process.eventFilter_HM * process.clusterCompatibilityFilter)
eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_primaryVertexFilter, process.Flag_clusterCompatibilityFilter ]

for P in eventFilterPaths:
    process.schedule.insert(0, P)

#from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD
#changeToMiniAOD(process)
