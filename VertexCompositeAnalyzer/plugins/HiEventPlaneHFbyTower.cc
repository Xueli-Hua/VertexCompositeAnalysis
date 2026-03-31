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

class HiEventPlaneHFbyTower : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  explicit HiEventPlaneHFbyTower(const edm::ParameterSet&);
  ~HiEventPlaneHFbyTower();

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

  //tree branches
  //event info
  uint  runNb;
  uint  eventNb;
  uint  lsNb;
  short centrality;
  int   Ntrkoffline;
  int   NtrkHP;
  uint candSize;

  double twQx;
  double twQy;
  double twQx_forw;
  double twQy_forw;
  double twQx_afterw;
  double twQy_afterw;

  bool isCentrality_;
    
  //tokens
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<reco::TrackCollection> generalTrkToken_;
  edm::EDGetTokenT<CaloTowerCollection> caloTowerToken_;
  edm::EDGetTokenT<reco::Centrality> centSrcToken_;
  edm::EDGetTokenT<int> centBinLabelToken_;

};

//
// static data member definitions
//

//
// constructors and destructor
//

HiEventPlaneHFbyTower::HiEventPlaneHFbyTower(const edm::ParameterSet& ps)
{
  //input tokens
  vtxToken_ = consumes<reco::VertexCollection>(ps.getUntrackedParameter<edm::InputTag>("vtxInputTag"));
  generalTrkToken_ = consumes<reco::TrackCollection>(ps.getUntrackedParameter<edm::InputTag>("trkInputTag"));
  caloTowerToken_ = consumes<CaloTowerCollection>(ps.getUntrackedParameter<edm::InputTag>("caloTowerInputTag"));

  isCentrality_ = (ps.exists("isCentrality") ? ps.getParameter<bool>("isCentrality") : false);
  if(isCentrality_)
  {
    centBinLabelToken_ = consumes<int>(ps.getParameter<edm::InputTag>("centBinLabelTag"));
    centSrcToken_ = consumes<reco::Centrality>(ps.getParameter<edm::InputTag>("centSrcTag"));
  }

}

HiEventPlaneHFbyTower::~HiEventPlaneHFbyTower()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HiEventPlaneHFbyTower::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
    
  fillRECO(iEvent,iSetup);
  EventInfoNtuple->Fill();
}

void
HiEventPlaneHFbyTower::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get collections
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_,vertices);
  if(!vertices.isValid()) throw cms::Exception("ForEventPlaneHF") << "Primary vertices  collection not found!" << std::endl;

  //best vertex
    double bestvz=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z();
 
  edm::Handle<CaloTowerCollection> towers;
  iEvent.getByToken(caloTowerToken_, towers);
  //if(!towers.isValid()) return;
  if(!towers.isValid()) throw cms::Exception("ForEventPlaneHF") << "HF tower collection not found!" << std::endl;

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
    
  //Calo tower info
  double twqx = 0;
  double twqy = 0;
  double twEt = 0;
  twQx = -1;
  twQy = -1;

  double twqx_forw = 0;
  double twqy_forw = 0;
  double twEt_forw = 0;
  twQx_forw = -1;
  twQy_forw = -1;
  double twqx_afterw = 0;
  double twqy_afterw = 0;
  double twEt_afterw = 0;
  twQx_afterw = -1;
  twQy_afterw = -1;

  for(unsigned itw = 0; itw < towers->size(); ++itw){
        
    const CaloTower & hit= (*towers)[itw];
    
    double et = hit.et(bestvz);
    double caloPhi = hit.phi();
    double caloEta = hit.eta();

    //if(et<0.05) continue;
    if(abs(caloEta)>=5 || abs(caloEta)<=3) continue;
    twqx += et*cos(2*caloPhi);
    twqy += et*sin(2*caloPhi);
    twEt += et;

    if(caloEta>3 && caloEta<5) {
      twqx_forw += et*cos(2*caloPhi);
      twqy_forw += et*sin(2*caloPhi);
      twEt_forw += et;
    }
    if(caloEta>-5 && caloEta<-3) {
      twqx_afterw += et*cos(2*caloPhi);
      twqy_afterw += et*sin(2*caloPhi);
      twEt_afterw += et;
    }        
  }
  twQx = twqx/twEt;
  twQy = twqy/twEt;

  twQx_forw = twqx_forw/twEt_forw;
  twQy_forw = twqy_forw/twEt_forw;
  twQx_afterw = twqx_afterw/twEt_afterw;
  twQy_afterw = twqy_afterw/twEt_afterw;
}

// ------------ method called once each job just before starting event
//loop  ------------
void
HiEventPlaneHFbyTower::beginJob()
{
    TH1D::SetDefaultSumw2();

    initTree();

    htrkpt = fs->make<TH1D>("hTrk",";pT",100,0,10);
}

void 
HiEventPlaneHFbyTower::initTree()
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

  EventInfoNtuple->Branch("twQx",&twQx,"twQx/D");
  EventInfoNtuple->Branch("twQy",&twQy,"twQy/D");
  EventInfoNtuple->Branch("twQx_forw",&twQx_forw,"twQx_forw/D");
  EventInfoNtuple->Branch("twQy_forw",&twQx_forw,"twQy_forw/D");
  EventInfoNtuple->Branch("twQx_afterw",&twQx_afterw,"twQx_afterw/D");
  EventInfoNtuple->Branch("twQy_afterw",&twQy_afterw,"twQy_afterw/D"); 
}

//--------------------------------------------------------------------------------------------------
void 
HiEventPlaneHFbyTower::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
}

// ------------ method called once each job just after ending the event
//loop  ------------
void
HiEventPlaneHFbyTower::endJob() {
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiEventPlaneHFbyTower);
