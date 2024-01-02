# VertexCompositeAnalysis

Example of setting up and running gamma+gamma to dimuon tree

cmsrel CMSSW_10_3_3_patch1

cd CMSSW_10_3_3_patch1/src

cmsenv

git clone -b 10_3_X https://github.com/davidlw/VertexCompositeAnalysis

cd VertexCompositeAnalysis

scram b -j8

cd VertexCompositeProducer/test

cmsRun PbPbSkimAndTree2018_DiMuContBothGammaGamma_mc_cfg.py 

## Changes
1. generalMuMuMassMin7Candidates -> generalMuMuMassMin2p5Candidates
   
   files:
   
          VetexCompositeProducer/test/PbPbSkimAndTree2018_DiMuContBoth_cfg.py
   
          VertexCompositeAnalyzer/python/dimuanalyzer_tree_cff.py
   
          VertexCompositeProducer/python/generalDiMuCandidates_cff.py
