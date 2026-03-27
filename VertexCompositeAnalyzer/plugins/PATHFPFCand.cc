// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMath.h>

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDefs.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include <HepMC/PdfInfo.h>

#include <iostream>
#include <fstream>

//
// constants, enums and typedefs
//

#define PI 3.1416
#define MAXCAN 50000
#define MAXDAU 3
#define MAXGDAU 2
#define MAXTRG 1024
#define MAXSEL 100
#define MAXPFCAN 50000

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;
typedef ROOT::Math::SVector<double, 6> SVector6;


//
// class decleration
//

class PATHFPFCand : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  explicit PATHFPFCand(const edm::ParameterSet&);
  ~PATHFPFCand();


private:
  virtual void beginJob();
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void endRun(const edm::Run&, const edm::EventSetup&) {};
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void fillpfCand(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void initTree();
  virtual void initHistogram();
  // ----------member data ---------------------------

  edm::Service<TFileService> fs;

  TTree* PATPFCandNtuple;

  bool   saveTree_;
  bool   saveHistogram_;

  //options
  bool doPFCandNtuple_;

  //cut variables

  //tree branches
  //event info
  uint  runNb;
  uint  eventNb;
  uint  lsNb;

  short centrality;
  int   Ntrkoffline;
  int   NtrkHP;

  //int nmuons = 0;
  bool isCentrality_;

  //PF candidate info
  float hiHFEPlus_pf[MAXPFCAN];
  float hiHFEMinus_pf[MAXPFCAN];
  float hiHFPlus_pfle;
  float hiHFMinus_pfle;
  uint nCountsHF_pf;
  uint nCountsHFPlus_pf;
  uint nCountsHFMinus_pf;

  //token
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandidateTag_;
  edm::EDGetTokenT<int> tok_centBinLabel_;
  edm::EDGetTokenT<reco::Centrality> tok_centSrc_;
  edm::EDGetTokenT<reco::TrackCollection> tok_tracks_;

};

//
// static data member definitions
//

//
// constructors and destructor
//

PATHFPFCand::PATHFPFCand(const edm::ParameterSet& iConfig) :
  pfCandidateTag_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidateSrc")))
{
  //options
  doPFCandNtuple_ = iConfig.getUntrackedParameter<bool>("doPFCandNtuple");
  saveTree_ = iConfig.getUntrackedParameter<bool>("saveTree");
  saveHistogram_ = iConfig.getUntrackedParameter<bool>("saveHistogram");

  isCentrality_ = (iConfig.exists("isCentrality") ? iConfig.getParameter<bool>("isCentrality") : false);
  if(isCentrality_)
  {
    tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
    tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
  }

  tok_tracks_ = consumes<reco::TrackCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("TrackCollection")));
}


PATHFPFCand::~PATHFPFCand()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PATHFPFCand::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //check event
  if(doPFCandNtuple_) fillpfCand(iEvent,iSetup);
  //if(saveTree_&&centrality>=80) PATPFCandNtuple->Fill();
  if(saveTree_) PATPFCandNtuple->Fill();
}


void
PATHFPFCand::fillpfCand(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  runNb = iEvent.id().run();
  eventNb = iEvent.id().event();
  lsNb = iEvent.luminosityBlock();

  centrality = -1;
  if(isCentrality_)
  {
    const auto& cent = iEvent.getHandle(tok_centSrc_);
    Ntrkoffline = (cent.isValid() ? cent->Ntracks() : -1);
    edm::Handle<int> cbin;
    iEvent.getByToken(tok_centBinLabel_, cbin);
    centrality = (cbin.isValid() ? *cbin : -1);

    edm::Handle<pat::PackedCandidateCollection> pfCandidates;
    iEvent.getByToken(pfCandidateTag_, pfCandidates);

    hiHFPlus_pfle = 0; hiHFMinus_pfle = 0;
    nCountsHF_pf = 0; nCountsHFPlus_pf = 0; nCountsHFMinus_pf = 0;

    for (const auto& pfcand : *pfCandidates) {
      if (pfcand.pdgId() == 1 || pfcand.pdgId() == 2){
        const bool eta_plus = (pfcand.eta() > 3.0) && (pfcand.eta() < 6.0);
        const bool eta_minus = (pfcand.eta() < -3.0) && (pfcand.eta() > -6.0);
        if (pfcand.et() < 0.0) continue;
        if (eta_plus || eta_minus)
        {   
          nCountsHF_pf++;
          if(eta_plus){
            hiHFEPlus_pf[nCountsHFPlus_pf] = pfcand.energy();
            if(pfcand.energy() >= hiHFPlus_pfle) hiHFPlus_pfle = pfcand.energy();
            nCountsHFPlus_pf++;
          }
          else if(eta_minus){
            hiHFEMinus_pf[nCountsHFMinus_pf] = pfcand.energy();
            if(pfcand.energy() >= hiHFMinus_pfle) hiHFMinus_pfle = pfcand.energy();
            nCountsHFMinus_pf++;
          }
        }
      }
    }
  }

  NtrkHP = -1;
  const auto& trackColl = iEvent.getHandle(tok_tracks_);
  if(trackColl.isValid()) 
  {
    NtrkHP = 0;
    for (const auto& trk : *trackColl) { if (trk.quality(reco::TrackBase::highPurity)) NtrkHP++; }
  }
}


// ------------ method called once each job just before starting event
//loop  ------------
void
PATHFPFCand::beginJob()
{
  TH1D::SetDefaultSumw2();

  // Check inputs
  if(!doPFCandNtuple_) throw cms::Exception("PATCompositeAnalyzer") << "No output for RECO Fix config!!" << std::endl;
  if(saveTree_) initTree();
  if(saveHistogram_) initHistogram();
  
}


void 
PATHFPFCand::initTree()
{ 
  PATPFCandNtuple = fs->make< TTree>("hiHF_PFCandidate","hiHF_PFCandidate");

  if(doPFCandNtuple_)
  {
    // Event info
    
    PATPFCandNtuple->Branch("RunNb",&runNb,"RunNb/i");
    PATPFCandNtuple->Branch("LSNb",&lsNb,"LSNb/i");
    PATPFCandNtuple->Branch("EventNb",&eventNb,"EventNb/i");
    
    if(isCentrality_) 
    {
      PATPFCandNtuple->Branch("centrality",&centrality,"centrality/S");
      PATPFCandNtuple->Branch("Ntrkoffline",&Ntrkoffline,"Ntrkoffline/I");
      PATPFCandNtuple->Branch("NtrkHP",&NtrkHP,"NtrkHP/I");
      PATPFCandNtuple->Branch("nCountsHF_pf",&nCountsHF_pf,"nCountsHF_pf/I");
      PATPFCandNtuple->Branch("nCountsHFPlus_pf",&nCountsHFPlus_pf,"nCountsHFPlus_pf/I");
      PATPFCandNtuple->Branch("nCountsHFMinus_pf",&nCountsHFMinus_pf,"nCountsHFMinus_pf/I");
      PATPFCandNtuple->Branch("hiHFEPlus_pf",hiHFEPlus_pf,"hiHFEPlus_pf[nCountsHFPlus_pf]/F");
      PATPFCandNtuple->Branch("hiHFEMinus_pf",hiHFEMinus_pf,"hiHFEMinus_pf[nCountsHFMinus_pf]/F");
      PATPFCandNtuple->Branch("hiHFPlus_pfle",&hiHFPlus_pfle,"hiHFPlus_pfle/F");
      PATPFCandNtuple->Branch("hiHFMinus_pfle",&hiHFMinus_pfle,"hiHFMinus_pfle/F");
    }

  } // doPFCandNtuple_

}

void
PATHFPFCand::initHistogram()
{
}


//--------------------------------------------------------------------------------------------------
void 
PATHFPFCand::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
}


// ------------ method called once each job just after ending the event
//loop  ------------
void 
PATHFPFCand::endJob()
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATHFPFCand);
