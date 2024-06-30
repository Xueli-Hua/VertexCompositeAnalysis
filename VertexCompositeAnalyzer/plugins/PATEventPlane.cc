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

//
// constants, enums and typedefs
//

#define PI 3.1416
#define MAXCAN 50000
#define MAXDAU 3
#define MAXGDAU 2
#define MAXTRG 1024
#define MAXSEL 100

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;
typedef ROOT::Math::SVector<double, 6> SVector6;

//
// class decleration
//

class PATEventPlane : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  explicit PATEventPlane(const edm::ParameterSet&);
  ~PATEventPlane();


private:
  virtual void beginJob();
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void endRun(const edm::Run&, const edm::EventSetup&) {};
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void fillRECO(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void initTree();
  virtual void initHistogram();
  // ----------member data ---------------------------

  edm::Service<TFileService> fs;

  TTree* PATCompositeNtuple;
  TH1D* htrkpt;
  TH2D* hEtavsPt_DauTrk;
  TH2D* houtmu;
  TH2D* htrk;
  TH2D* hMassvsPt_dimu;
  TH2D* hEtavsPt_mu1;
  TH2D* hEtavsPt_mu2;
  TH2D* htrkd1;
  TH2D* htrkd2;
  TH1D* htwEt;
  TH1D* htwqx;
  TH1D* htwqy;

  bool   saveTree_;
  bool   saveHistogram_;

  //options
  bool doRecoNtuple_;

  //cut variables

  //tree branches
  //event info
  uint  runNb;
  uint  eventNb;
  uint  lsNb;
  short centrality;
  int   Ntrkoffline;
  int   NtrkHP;
  short nPV;
  float bestvx;
  float bestvy;
  float bestvz;
  float bestvxError;
  float bestvyError;
  float bestvzError;

  double trkQx;
  double trkQy;
  double twQx;
  double twQy;

  double all_trkQx;
  double all_trkQy;


  bool isCentrality_;

  //token
  edm::EDGetTokenT<reco::BeamSpot> tok_offlineBS_;
  edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;

  edm::EDGetTokenT<pat::MuonCollection> tok_muoncol_;
  //edm::EDGetTokenT<reco::MuonCollection> tok_muoncol_;
  edm::EDGetTokenT<reco::TrackCollection> tok_tracks_;

  edm::EDGetTokenT<int> tok_centBinLabel_;
  edm::EDGetTokenT<reco::Centrality> tok_centSrc_;
  edm::EDGetTokenT<CaloTowerCollection> caloTowerToken_;
 
  std::vector<double> d1Eta;
  std::vector<double> d2Eta;
  std::vector<double> d1Phi;
  std::vector<double> d2Phi; 

};

//
// static data member definitions
//

//
// constructors and destructor
//

PATEventPlane::PATEventPlane(const edm::ParameterSet& iConfig)
{
  //options
  doRecoNtuple_ = iConfig.getUntrackedParameter<bool>("doRecoNtuple");
  saveTree_ = iConfig.getUntrackedParameter<bool>("saveTree");
  saveHistogram_ = iConfig.getUntrackedParameter<bool>("saveHistogram");

  //input tokens
  tok_offlineBS_ = consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotSrc"));
  tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCollection"));
  tok_muoncol_ = consumes<pat::MuonCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("MuonCollection")));
  //tok_muoncol_ = consumes<reco::MuonCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("MuonCollection")));
  tok_tracks_ = consumes<reco::TrackCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("TrackCollection")));
  caloTowerToken_ = consumes<CaloTowerCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("caloTowerInputTag")));


  isCentrality_ = (iConfig.exists("isCentrality") ? iConfig.getParameter<bool>("isCentrality") : false);
  if(isCentrality_)
  {
    tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
    tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
  }

}


PATEventPlane::~PATEventPlane()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PATEventPlane::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //check event
  if(doRecoNtuple_) fillRECO(iEvent,iSetup);
  if(saveTree_) PATCompositeNtuple->Fill();
}


void
PATEventPlane::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get collection
  edm::Handle<reco::BeamSpot> beamspot;
  iEvent.getByToken(tok_offlineBS_, beamspot);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(tok_offlinePV_, vertices);
  if(!vertices.isValid()) throw cms::Exception("PATEventPlane") << "Primary vertices  collection not found!" << std::endl;


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
  }
  
  NtrkHP = -1;
  const auto& trackColl = iEvent.getHandle(tok_tracks_);
  if(trackColl.isValid()) 
  {
    NtrkHP = 0;
    for (const auto& trk : *trackColl) { if (trk.quality(reco::TrackBase::highPurity)) NtrkHP++; }
  }

  nPV = vertices->size();
  //best vertex
  const auto& vtxPrimary = (vertices->size()>0 ? (*vertices)[0] : reco::Vertex());
  const bool& isPV = (!vtxPrimary.isFake() && vtxPrimary.tracksSize()>=2);
  const auto& bs = (!isPV ? reco::Vertex(beamspot->position(), beamspot->covariance3D()) : reco::Vertex());
  const reco::Vertex& vtx = (isPV ? vtxPrimary : bs);
  bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
  const math::XYZPoint bestvtx(bestvx, bestvy, bestvz);
  bestvzError = vtx.zError(), bestvxError = vtx.xError(), bestvyError = vtx.yError();

 
  //edm::Handle<reco::MuonCollection> recoMuons;
  edm::Handle<pat::MuonCollection> recoMuons;
  iEvent.getByToken(tok_muoncol_, recoMuons);

  //reco::MuonCollection muonColl;
  pat::MuonCollection muonColl;
  for (const auto& muon : *recoMuons) {
    const reco::TrackRef& trackRef = muon.track();
    if(trackRef.isNull()) continue;
    muonColl.push_back(muon);
  }

  //auto out = std::make_unique<std::vector<reco::Muon>>();
  auto out = std::make_unique<std::vector<pat::Muon>>();
  for (uint ic = 0; ic < muonColl.size(); ic++) {
    const pat::Muon& cand1 = muonColl[ic];

    for(uint fc = ic+1; fc < muonColl.size(); fc++) {
       const pat::Muon& cand2 = muonColl[fc];

       const auto cand1P4 = math::PtEtaPhiMLorentzVector(cand1.pt(), cand1.eta(), cand1.phi(), 0.10565837);
       const auto cand2P4 = math::PtEtaPhiMLorentzVector(cand2.pt(), cand2.eta(), cand2.phi(), 0.10565837);

       const double& mass = (cand1P4 + cand2P4).mass();
       const double& pt = (cand1P4 + cand2P4).pt();

       //if (abs(mass-MASS_.at(443)) < WIDTH_.at(443) && pt < 0.2) {
       if (mass>2.9 && mass<3.2 && pt < 0.2) {
         hMassvsPt_dimu->Fill(mass,pt);
         hEtavsPt_mu1->Fill(cand1.eta(),cand1.pt());
         hEtavsPt_mu2->Fill(cand2.eta(),cand2.pt());
         out->push_back(cand1);out->push_back(cand2);
	 d1Eta.push_back(cand1P4.eta());
         d1Phi.push_back(cand1P4.phi());
         d2Eta.push_back(cand2P4.eta());
         d2Phi.push_back(cand2P4.phi());
      }
    }
  }

  //track info
  double trkqx = 0;
  double trkqy = 0;
  double trkPt = 0;
  trkQx = -1;
  trkQy = -1;
  bool DauTrk = false;

  double all_trkqx= 0;
  double all_trkqy = 0;
  double all_trkPt = 0;
  all_trkQx = -1;
  all_trkQy = -1;

  for(unsigned it=0; it<trackColl->size(); ++it){
	DauTrk = false;
	reco::TrackRef track(trackColl, it);
	double pt  = track->pt();
	double phi = track->phi();
	double eta = track->eta();
	htrk->Fill(track->eta(),track->pt());

    	all_trkqx += pt*cos(2*phi);
    	all_trkqy += pt*sin(2*phi);
    	all_trkPt += pt;

    	for (unsigned i=0; i<d1Eta.size(); ++i)
    	{
            if( fabs(eta-d1Eta[i]) <0.001 && fabs(phi-d1Phi[i]) <0.001 ) htrkd1->Fill(track->eta(),track->pt());
      	    if( fabs(eta-d2Eta[i]) <0.001 && fabs(phi-d2Phi[i]) <0.001) htrkd2->Fill(track->eta(),track->pt());
   	}

	reco::TrackRef muonTrack;
    	for (std::vector<pat::Muon>::const_iterator muon = out->begin(); muon < out->end(); muon++) {
            muonTrack = muon->innerTrack();
            //if (muonTrack == track) DauTrk = true;
	    if (track->charge() == muonTrack->charge() && std::abs(muonTrack->eta() - track->eta()) < 1.E-4 && std::abs(reco::deltaPhi(muonTrack->phi(), track->phi())) < 1.E-4 && std::abs(muonTrack->pt() - track->pt()) / muonTrack->pt() < 1.E-4) DauTrk = true; 
    	}

    	if (DauTrk == true) {
      	    hEtavsPt_DauTrk->Fill(track->eta(),track->pt());
      	    continue;
    	}
    	htrkpt->Fill(pt);
    
    	trkqx += pt*cos(2*phi);
    	trkqy += pt*sin(2*phi);
    	trkPt += pt;
  }
  trkQx = trkqx/trkPt;
  trkQy = trkqy/trkPt;
  all_trkQx = all_trkqx/all_trkPt;
  all_trkQy = all_trkqy/all_trkPt;

  //Calo tower info
  //
  edm::Handle<CaloTowerCollection> towers;
  iEvent.getByToken(caloTowerToken_, towers);
  if(!towers.isValid()) return;
  
  double twqx = 0;
  double twqy = 0;
  double twEt = 0;
  twQx = -1;
  twQy = -1;
  for(unsigned itw = 0; itw < towers->size(); ++itw){

    const CaloTower & hit= (*towers)[itw];

    double et = hit.et(bestvz);
    double caloPhi = hit.phi();

    twqx += et*cos(2*caloPhi);
    twqy += et*sin(2*caloPhi);
    twEt += et;

    htwqx->Fill(twqx);
    htwqy->Fill(twqy);
    htwEt->Fill(twEt);

  }
  twQx = twqx/twEt;
  twQy = twqy/twEt;

}


// ------------ method called once each job just before starting event
//loop  ------------
void
PATEventPlane::beginJob()
{
  TH1D::SetDefaultSumw2();

  // Check inputs
  if(!doRecoNtuple_) throw cms::Exception("PATCompositeAnalyzer") << "No output for RECO Fix config!!" << std::endl;
  if(saveTree_) initTree();
  if(saveHistogram_) initHistogram();
  
}


void 
PATEventPlane::initTree()
{ 
  PATCompositeNtuple = fs->make< TTree>("EventPlane","EventPlane");

  if(doRecoNtuple_)
  {
    // Event info
    
    PATCompositeNtuple->Branch("RunNb",&runNb,"RunNb/i");
    PATCompositeNtuple->Branch("LSNb",&lsNb,"LSNb/i");
    PATCompositeNtuple->Branch("EventNb",&eventNb,"EventNb/i");
    PATCompositeNtuple->Branch("nPV",&nPV,"nPV/S");
    PATCompositeNtuple->Branch("bestvtxX",&bestvx,"bestvtxX/F");
    PATCompositeNtuple->Branch("bestvtxY",&bestvy,"bestvtxY/F");
    PATCompositeNtuple->Branch("bestvtxZ",&bestvz,"bestvtxZ/F");
    
    if(isCentrality_) 
    {
      PATCompositeNtuple->Branch("centrality",&centrality,"centrality/S");
      PATCompositeNtuple->Branch("Ntrkoffline",&Ntrkoffline,"Ntrkoffline/I");
      PATCompositeNtuple->Branch("NtrkHP",&NtrkHP,"NtrkHP/I");
    }
    
    PATCompositeNtuple->Branch("trkQx",&trkQx,"trkQx/D");
    PATCompositeNtuple->Branch("trkQy",&trkQy,"trkQy/D");
    PATCompositeNtuple->Branch("twQx",&twQx,"twQx/D");
    PATCompositeNtuple->Branch("twQy",&twQy,"twQy/D");
    PATCompositeNtuple->Branch("all_trkQx",&all_trkQx,"all_trkQx/D");
    PATCompositeNtuple->Branch("all_trkQy",&all_trkQy,"all_trkQy/D");

  } // doRecoNtuple_

}

void
PATEventPlane::initHistogram()
{
  htrkpt = fs->make<TH1D>("hTrk",";pT",100,0,10);
  hEtavsPt_DauTrk = fs->make<TH2D>("hEtavsPt_DauTrk",";Eta;Pt",200,-10,10,100,0,10);
  houtmu = fs->make<TH2D>("houtmu",";Eta;Pt",200,-10,10,100,0,10);
  htrk = fs->make<TH2D>("htrk",";Eta;Pt",200,-10,10,100,0,10);
  hMassvsPt_dimu = fs->make<TH2D>("hMassvsPt_dimu",";Mass;Pt",20,2.6,4.2,100,0,10);
  hEtavsPt_mu1 = fs->make<TH2D>("hEtavsPt_mu1",";Eta;Pt",500,-10,10,500,0,10);
  hEtavsPt_mu2 = fs->make<TH2D>("hEtavsPt_mu2",";Eta;Pt",500,-10,10,500,0,10);
  htrkd1 = fs->make<TH2D>("htrkd1",";Eta;Pt",500,-10,10,500,0,10);
  htrkd2 = fs->make<TH2D>("htrkd2",";Eta;Pt",500,-10,10,500,0,10);

  htwEt=fs->make<TH1D>("htwEt",";Et",100,0,100);
  htwqx=fs->make<TH1D>("htwqx",";qx",100,0,100);
  htwqy=fs->make<TH1D>("htwqy",";qy",100,0,100);
}


//--------------------------------------------------------------------------------------------------
void 
PATEventPlane::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
}


// ------------ method called once each job just after ending the event
//loop  ------------
void 
PATEventPlane::endJob()
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATEventPlane);
