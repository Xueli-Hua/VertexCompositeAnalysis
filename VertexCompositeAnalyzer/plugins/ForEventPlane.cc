// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//from chenyan
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
//#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

//
// constants, enums and typedefs
//

#define PI 3.1416
#define MAXTRG 1024
#define MAXSEL 100


//
// class decleration
//

#define PI 3.1416
#define MAXCAN 10000
#define MAXTRG 1024
#define MAXSEL 100

using namespace std;

class ForEventPlane : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  explicit ForEventPlane(const edm::ParameterSet&);
  ~ForEventPlane();

private:
  virtual void beginJob();
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void endRun(const edm::Run&, const edm::EventSetup&) {};
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void fillRECO(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  virtual void initTree();

  // ----------member data ---------------------------

  edm::Service<TFileService> fs;

  TTree* EventInfoNtuple;

  TH1D* htrkpt;
  TH2D* hEtavsPt_DauTrk;
  TH2D* hMassvsPt_dimu;
  TH2D* hEtavsPt_mu1;
  TH2D* hEtavsPt_mu2;

  //tree branches
  //event info
  uint  runNb;
  uint  eventNb;
  uint  lsNb;
  short centrality;
  int   Ntrkoffline;
  int   NtrkHP;
  uint candSize;

  double trkQx;
  double trkQy;
  double twQx;
  double twQy;

  double all_trkQx;
  double all_trkQy;
  
  bool isCentrality_;
    
  //tokens
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<reco::TrackCollection> generalTrkToken_;
  edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
  edm::EDGetTokenT<CaloTowerCollection> caloTowerToken_;
  edm::EDGetTokenT<reco::Centrality> centSrcToken_;
  edm::EDGetTokenT<int> centBinLabelToken_;

  const std::map<uint, double> MASS_ = {
  {11, 0.000511}, {13, 0.10565837}, {15, 1.77686}, // leptons
  {22, 0}, {23, 91.188}, {24, 80.38}, // bosons
  {113, 0.77526}, {211, 0.13957018}, {310, 0.497614}, {321, 0.493677}, {333, 1.019445}, // light and strange mesons
  {411, 1.86962}, {421, 1.86484}, {431, 1.96847}, {511, 5.27929}, // charmed and bottom mesons
  {443, 3.096900}, {100443, 3.68609}, {553, 9.4603}, {100553, 10.0233}, {200553, 10.3552}, // quarkonia
  {2212, 0.938272013}, {3122, 1.115683}, {3312, 1.32171}, {3334, 1.67245}, {4122, 2.28646}, // baryons
  {3872, 3.87169} // exotic hadrons
};
  const std::map<uint, float> WIDTH_ = {{211, 3.5E-7f}, {321, 1.6E-5f}, {2212, 1.6E-5f}, {443, 4.51623e-02}};

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

ForEventPlane::ForEventPlane(const edm::ParameterSet& ps)
{
  //input tokens
  vtxToken_ = consumes<reco::VertexCollection>(ps.getUntrackedParameter<edm::InputTag>("vtxInputTag"));
  generalTrkToken_ = consumes<reco::TrackCollection>(ps.getUntrackedParameter<edm::InputTag>("trkInputTag"));
  muonsToken_ = consumes<reco::MuonCollection>(ps.getUntrackedParameter<edm::InputTag>("muInputTag"));
  caloTowerToken_ = consumes<CaloTowerCollection>(ps.getUntrackedParameter<edm::InputTag>("caloTowerInputTag"));

  isCentrality_ = (ps.exists("isCentrality") ? ps.getParameter<bool>("isCentrality") : false);
  if(isCentrality_)
  {
    centBinLabelToken_ = consumes<int>(ps.getParameter<edm::InputTag>("centBinLabelTag"));
    centSrcToken_ = consumes<reco::Centrality>(ps.getParameter<edm::InputTag>("centSrcTag"));
  }

}


ForEventPlane::~ForEventPlane()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ForEventPlane::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
    
  fillRECO(iEvent,iSetup);
  EventInfoNtuple->Fill();
}

void
ForEventPlane::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get collections
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_,vertices);
  if(!vertices.isValid()) throw cms::Exception("ForEventPlane") << "Primary vertices  collection not found!" << std::endl;

  //best vertex
    double bestvz=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z();
 
  edm::Handle<CaloTowerCollection> towers;
  iEvent.getByToken(caloTowerToken_, towers);
  if(!towers.isValid()) return;

  runNb = iEvent.id().run();
  eventNb = iEvent.id().event();
  lsNb = iEvent.luminosityBlock();

  centrality = -1;
  if(isCentrality_)
  {
    edm::Handle<reco::Centrality> cent;
    iEvent.getByToken(centSrcToken_, cent);
    Ntrkoffline = (cent.isValid() ? cent->Ntracks() : -1);

    edm::Handle<int> cbin;
    iEvent.getByToken(centBinLabelToken_, cbin);
    centrality = (cbin.isValid() ? *cbin : -1);
  }

  NtrkHP = -1;
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(generalTrkToken_, tracks);
  if(tracks.isValid()) 
  {
    NtrkHP = 0;
    for (const auto& trk : *tracks) { if (trk.quality(reco::TrackBase::highPurity)) NtrkHP++; }
  }
    
  //RECO Candidate info; muon info
  edm::Handle<reco::MuonCollection> recoMuons;
  iEvent.getByToken(muonsToken_, recoMuons);

  reco::MuonCollection muonColl;
  for (const auto& muon : *recoMuons) {
    const reco::TrackRef& trackRef = muon.track();
    if(trackRef.isNull()) continue;
    muonColl.push_back(muon);
  }

  auto out = std::make_unique<std::vector<reco::Muon>>();
  for (uint ic = 0; ic < muonColl.size(); ic++) {
    const reco::Muon& cand1 = muonColl[ic];

    for(uint fc = ic+1; fc < muonColl.size(); fc++) {
       const reco::Muon& cand2 = muonColl[fc];

       //auto cand1P4 = cand1.polarP4(); cand1P4.SetM(MASS_.at(13));
       //auto cand2P4 = cand2.polarP4(); cand2P4.SetM(MASS_.at(13));
       const auto cand1P4 = math::PtEtaPhiMLorentzVector(cand1.pt(), cand1.eta(), cand1.phi(), MASS_.at(13));
       const auto cand2P4 = math::PtEtaPhiMLorentzVector(cand2.pt(), cand2.eta(), cand2.phi(), MASS_.at(13));
      
       const double& mass = (cand1P4 + cand2P4).mass();
       const double& pt = (cand1P4 + cand2P4).pt();

       //if (abs(mass-MASS_.at(443)) < WIDTH_.at(443) && pt < 0.2) {
       if (mass>2.9 && mass<3.2 && pt < 0.2) {
        d1Eta.push_back(cand1P4.eta());
        d1Phi.push_back(cand1P4.phi());
        d2Eta.push_back(cand2P4.eta());
        d2Phi.push_back(cand2P4.phi());
        out->push_back(cand1);out->push_back(cand2);
        hMassvsPt_dimu->Fill(mass,pt);
        hEtavsPt_mu1->Fill(cand1.eta(),cand1.pt());
        hEtavsPt_mu2->Fill(cand2.eta(),cand2.pt());
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
  
  for(unsigned it=0; it<tracks->size(); ++it){
    DauTrk = false;
    reco::TrackRef track(tracks, it);
    double pt  = track->pt();
    double phi = track->phi();
    for (std::vector<reco::Muon>::const_iterator muon = out->begin(); muon < out->end(); muon++) {
        if (muon->innerTrack() == track) DauTrk = true;
    }

    /*double eta = track->eta();
    for(unsigned i=0; i<d1Eta.size(); ++i)
    {
      if( fabs(eta-d1Eta[i]) <0.03 && fabs(phi-d1Phi[i]) <0.03 ) DauTrk = true;
      if( fabs(eta-d2Eta[i]) <0.03 && fabs(phi-d2Phi[i]) <0.03) DauTrk = true;
    }*/

    all_trkqx += pt*cos(2*phi);
    all_trkqy += pt*sin(2*phi);
    all_trkPt += pt;
    
    if(DauTrk == true) {
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
        
  }
  twQx = twqx/twEt;
  twQy = twqy/twEt;
}

// ------------ method called once each job just before starting event
//loop  ------------
void
ForEventPlane::beginJob()
{
    TH1D::SetDefaultSumw2();

    initTree();

    htrkpt = fs->make<TH1D>("hTrk",";pT",100,0,10);
    hEtavsPt_DauTrk = fs->make<TH2D>("hEtavsPt_DauTrk",";Eta;Pt",200,-10,10,100,0,10);
    hMassvsPt_dimu = fs->make<TH2D>("hMassvsPt_dimu",";Mass;Pt",20,2.6,4.2,100,0,10);
    hEtavsPt_mu1 = fs->make<TH2D>("hEtavsPt_mu1",";Eta;Pt",100,-10,10,100,0,10);
    hEtavsPt_mu2 = fs->make<TH2D>("hEtavsPt_mu2",";Eta;Pt",100,-10,10,100,0,10);
}

void 
ForEventPlane::initTree()
{ 
  EventInfoNtuple = fs->make< TTree>("EventInfoNtuple","EventInfoNtuple");

  // Event info
  EventInfoNtuple->Branch("RunNb",&runNb,"RunNb/i");
  EventInfoNtuple->Branch("LSNb",&lsNb,"LSNb/i");
  EventInfoNtuple->Branch("EventNb",&eventNb,"EventNb/i");
  if(isCentrality_) 
  {
    EventInfoNtuple->Branch("centrality",&centrality,"centrality/S");
    EventInfoNtuple->Branch("Ntrkoffline",&Ntrkoffline,"Ntrkoffline/I");
    EventInfoNtuple->Branch("NtrkHP",&NtrkHP,"NtrkHP/I");
  }
  EventInfoNtuple->Branch("trkQx",&trkQx,"trkQx/D");
  EventInfoNtuple->Branch("trkQy",&trkQy,"trkQy/D");
  EventInfoNtuple->Branch("twQx",&twQx,"twQx/D");
  EventInfoNtuple->Branch("twQy",&twQy,"twQy/D");
  EventInfoNtuple->Branch("all_trkQx",&all_trkQx,"all_trkQx/D");
  EventInfoNtuple->Branch("all_trkQy",&all_trkQy,"all_trkQy/D");
}

//--------------------------------------------------------------------------------------------------
void 
ForEventPlane::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
}

// ------------ method called once each job just after ending the event
//loop  ------------
void
ForEventPlane::endJob() {
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(ForEventPlane);
