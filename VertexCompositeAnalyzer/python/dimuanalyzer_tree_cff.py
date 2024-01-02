import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeAnalyzer.dimuanalyzer_tree_cfi import *

dimucontana = dimuana.clone(
    VertexCompositeCollection = cms.untracked.InputTag("generalMuMuMassMin2p5Candidates:DiMu")
)
dimucontana_mc = dimuana_mc.clone(
    VertexCompositeCollection = cms.untracked.InputTag("generalDiMuCandidates:DiMu")
)

dimucontana_wrongsign = dimuana.clone(
    VertexCompositeCollection = cms.untracked.InputTag("generalMuMuMassMin2p5CandidatesWrongSign:DiMu")
)
dimucontana_wrongsign_mc = dimuana.clone(
    VertexCompositeCollection = cms.untracked.InputTag("generalDiMuCandidatesWrongSign:DiMu")
)
