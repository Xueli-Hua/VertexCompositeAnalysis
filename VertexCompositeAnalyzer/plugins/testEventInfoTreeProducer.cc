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

#include "DataFormats/TrackReco/interface/Track.h"
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

class testEventInfoTreeProducer : public edm::EDAnalyzer {
public:
  explicit testEventInfoTreeProducer(const edm::ParameterSet&);
  ~testEventInfoTreeProducer();


private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void fillRECO(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  virtual void initTree();

  // ----------member data ---------------------------

  edm::Service<TFileService> fs;

  TTree* EventInfoNtuple;

  //tree branches
  //event info
  uint  runNb;
  uint  eventNb;
  uint  lsNb;
  short centrality;
  int   NtrkHP;
  uint candSize;
  
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
  const std::map<uint, float> WIDTH_ = {{211, 3.5E-7f}, {321, 1.6E-5f}, {2212, 1.6E-5f}, {443, 6E-6f}};

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

testEventInfoTreeProducer::testEventInfoTreeProducer(const edm::ParameterSet& ps)
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


testEventInfoTreeProducer::~testEventInfoTreeProducer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
testEventInfoTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
    
  fillRECO(iEvent,iSetup);
  EventInfoNtuple->Fill();
}

void
testEventInfoTreeProducer::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get collections
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_,vertices);
  if(!vertices.isValid()) throw cms::Exception("testEventInfoTreeProducer") << "Primary vertices  collection not found!" << std::endl;

  //best vertex
    double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(generalTrkToken_, tracks);
    
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

    edm::Handle<int> cbin;
    iEvent.getByToken(centBinLabelToken_, cbin);
    centrality = (cbin.isValid() ? *cbin : -1);
  }

  NtrkHP = -1;
  if(tracks.isValid()) 
  {
    NtrkHP = 0;
    for (const auto& trk : *tracks) { if (trk.quality(reco::TrackBase::highPurity)) NtrkHP++; }
  }
    
  //RECO Candidate info; muon info
  edm::Handle<reco::MuonCollection> recoMuons;
  iEvent.getByToken(muonsToken_, recoMuons);

  for (const auto& mu1 : *recoMuons) {
    if (!(mu1.isPFMuon() || mu1.isGlobalMuon() || mu1.isTrackerMuon()))
      continue;

    for (const auto& mu2 : *recoMuons) {
      if (!(mu2.isPFMuon() || mu2.isGlobalMuon() || mu2.isTrackerMuon()))
        continue;
      double dimuCharge = mu1.charge()+mu2.charge();
      double dimuMass = mu1.mass()+mu2.mass();

      if (dimuCharge==0 && abs(dimuMass-MASS_.at(443)) < WIDTH_.at(443)) {
        d1Eta.push_back(mu1.eta());
        d1Phi.push_back(mu1.phi());
        d2Eta.push_back(mu2.eta());
        d2Phi.push_back(mu2.phi());
      }
    }
  }

  //track info and Ntrkoffline
  //int Ntrkoffline = 0;
  double trkqx = 0;
  double trkqy = 0;
  double trkPt = 0;
  double trkQx = 0;
  double trkQy = 0;
  //bool DauTrk = false;
  

  for(unsigned it=0; it<tracks->size(); ++it){
        
    const reco::Track & trk = (*tracks)[it];
    double eta = trk.eta();
    double pt  = trk.pt();
    double phi = trk.phi();

    //for(unsigned i=0; i<d1Eta.size(); ++i)
   // {
    //  if( fabs(eta-d1Eta[i]) <0.03 && fabs(phi-d1Phi[i]) <0.03 ) DauTrk = true;
     // if( fabs(eta-d2Eta[i]) <0.03 && fabs(phi-d2Phi[i]) <0.03) DauTrk = true;
    //}

    //if(DauTrk == true) continue;

    trkqx += pt*cos(2*phi);
    trkqy += pt*sin(2*phi);
    trkPt += pt;
  }
  trkQx = trkqx/trkPt;
  trkQy = trkqy/trkPt;
    
  //Calo tower info
  double twqx = 0;
  double twqy = 0;
  double twEt = 0;
  for(unsigned itw = 0; itw < towers->size(); ++itw){
        
    const CaloTower & hit= (*towers)[itw];
        
    //double caloEta = hit.eta();
    double et = hit.et(bestvz);
    double caloPhi = hit.phi();

    twqx += et*cos(2*caloPhi);
    twqy += et*sin(2*caloPhi);
    twEt += et;
        
  }
  double twQx = twqx/twEt;
  double twQy = twqy/twEt;
}

// ------------ method called once each job just before starting event
//loop  ------------
void
testEventInfoTreeProducer::beginJob()
{
    TH1D::SetDefaultSumw2();

    initTree();
}

void 
testEventInfoTreeProducer::initTree()
{ 
  EventInfoNtuple = fs->make< TTree>("EventInfoNtuple","EventInfoNtuple");

  // Event info
  EventInfoNtuple->Branch("RunNb",&runNb,"RunNb/i");
  EventInfoNtuple->Branch("LSNb",&lsNb,"LSNb/i");
  EventInfoNtuple->Branch("EventNb",&eventNb,"EventNb/i");
  if(isCentrality_) 
  {
    EventInfoNtuple->Branch("centrality",&centrality,"centrality/S");
    EventInfoNtuple->Branch("NtrkHP",&NtrkHP,"NtrkHP/I");
  }
  EventInfoNtuple->Branch("trkQx",&trkQx,"trkQx/D");
  EventInfoNtuple->Branch("trkQy",&trkQy,"trkQy/D");
  EventInfoNtuple->Branch("twQx",&twQx,"twQx/D");
  EventInfoNtuple->Branch("twQy",&twQy,"twQy/D");
}

// ------------ method called once each job just after ending the event
//loop  ------------
void
testEventInfoTreeProducer::endJob() {
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(testEventInfoTreeProducer);
