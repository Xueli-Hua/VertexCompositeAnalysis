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

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;
typedef ROOT::Math::SVector<double, 6> SVector6;


//
// class decleration
//

class HiEventPlaneTrackMC : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  explicit HiEventPlaneTrackMC(const edm::ParameterSet&);
  ~HiEventPlaneTrackMC();


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

  TTree* PATEventPlaneNtuple;
  TH1D* hdR;
  TH1D* htrketa_forw;
  TH1D* htrketa_afterw;
  TH1D* htrkpt;
  TH1D* htrketa;

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

  Double_t trkQx;
  Double_t trkQy;
  Double_t all_trkQx;
  Double_t all_trkQy;
  Double_t trkQx_forw;
  Double_t trkQy_forw;
  Double_t trkQx_afterw;
  Double_t trkQy_afterw;

  Double_t trkQx_v3;
  Double_t trkQy_v3;
  Double_t trkQx_v3_forw;
  Double_t trkQy_v3_forw;
  Double_t trkQx_v3_afterw;
  Double_t trkQy_v3_afterw;

  //int nmuons = 0;
  bool isCentrality_;
  double maxDR_;
  std::vector<int> pdgIds_;
  //std::vector<int> statusIds_;

  //token
  edm::EDGetTokenT<reco::BeamSpot> tok_offlineBS_;
  edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
  edm::EDGetTokenT<reco::TrackCollection> tok_tracks_;
  edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;

  edm::EDGetTokenT<int> tok_centBinLabel_;
  edm::EDGetTokenT<reco::Centrality> tok_centSrc_;
 
};

//
// static data member definitions
//

//
// constructors and destructor
//

HiEventPlaneTrackMC::HiEventPlaneTrackMC(const edm::ParameterSet& iConfig)
{
  //options
  doRecoNtuple_ = iConfig.getUntrackedParameter<bool>("doRecoNtuple");
  saveTree_ = iConfig.getUntrackedParameter<bool>("saveTree");
  saveHistogram_ = iConfig.getUntrackedParameter<bool>("saveHistogram");

  maxDR_ = iConfig.getParameter<double>("maxDR");
  pdgIds_ = iConfig.getParameter<std::vector<int>>("pdgIds");
  //statusIds_ = iConfig.getParameter<std::vector<int>>("statusIds");

  //input tokens
  tok_offlineBS_ = consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotSrc"));
  tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCollection"));
  tok_tracks_ = consumes<reco::TrackCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("TrackCollection")));
  tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection")));
  
  isCentrality_ = (iConfig.exists("isCentrality") ? iConfig.getParameter<bool>("isCentrality") : false);
  if(isCentrality_)
  {
    tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
    tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
  }

}


HiEventPlaneTrackMC::~HiEventPlaneTrackMC()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HiEventPlaneTrackMC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //check event
  if(doRecoNtuple_) fillRECO(iEvent,iSetup);
  if(saveTree_) PATEventPlaneNtuple->Fill();
}


void
HiEventPlaneTrackMC::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get collection
  edm::Handle<reco::BeamSpot> beamspot;
  iEvent.getByToken(tok_offlineBS_, beamspot);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(tok_offlinePV_, vertices);
  if(!vertices.isValid()) throw cms::Exception("PATEventPlaneTrack") << "Primary vertices  collection not found!" << std::endl;
  if(!vertices.isValid()) return;

  edm::Handle<reco::TrackCollection> trackColl;
  iEvent.getByToken(tok_tracks_, trackColl);
  //const auto& trackColl = iEvent.getHandle(tok_tracks_);

  edm::Handle<reco::GenParticleCollection> genpars;                                                                         iEvent.getByToken(tok_genParticle_, genpars);                                                                             if(!genpars.isValid()) throw cms::Exception("PATCompositeAnalyzer") << "Gen matching cannot be done without Gen collection!" << std::endl;
  if(!genpars.isValid()) return;

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

  //track info
  double trkqx = 0;
  double trkqy = 0;
  double trkPt = 0;
  trkQx = -1;
  trkQy = -1;
  double all_trkqx= 0;
  double all_trkqy = 0;
  double all_trkPt = 0;
  all_trkQx = -1;
  all_trkQy = -1;

  double trkPt_forw = 0;
  double trkPt_afterw = 0;
  double trkqx_forw = 0;
  double trkqy_forw = 0;
  double trkqx_afterw = 0;
  double trkqy_afterw = 0;
  trkQx_forw = -1;
  trkQy_forw = -1;
  trkQx_afterw = -1;
  trkQy_afterw = -1;

  double trkqx_v3 = 0;
  double trkqy_v3 = 0;
  double trkqx_v3_forw = 0;
  double trkqy_v3_forw = 0;
  double trkqx_v3_afterw = 0;
  double trkqy_v3_afterw = 0;
  trkQx_v3 = -1;
  trkQy_v3 = -1;
  trkQx_v3_forw = -1;
  trkQy_v3_forw = -1;
  trkQx_v3_afterw = -1;
  trkQy_v3_afterw = -1;

  for(unsigned it=0; it<trackColl->size(); ++it){
    bool isGenMuon = false;
    reco::TrackRef track(trackColl, it);

    double dzvtx = track->dz(bestvtx);
    double dxyvtx = track->dxy(bestvtx);
    double dzerror = sqrt(track->dzError()*track->dzError()+bestvzError*bestvzError);
    double dxyerror = sqrt(track->d0Error()*track->d0Error()+bestvxError*bestvyError);
        
    if(!track->quality(reco::TrackBase::highPurity)) continue;
    if(fabs(track->ptError())/track->pt()>0.10) continue;
    if(fabs(dzvtx/dzerror) > 3) continue;
    if(fabs(dxyvtx/dxyerror) > 3) continue;
    if(track->pt()<=0.3 || track->pt()>=3.0) continue;
    if(abs(track->eta())>=2.4) continue;

    htrkpt->Fill(track->pt());
    htrketa->Fill(track->eta());
          
    double pt  = track->pt();
    double phi = track->phi();
    double eta = track->eta();

    all_trkqx += pt*cos(2*phi);
    all_trkqy += pt*sin(2*phi);
    all_trkPt += pt;

    for(unsigned idx=0; idx<genpars->size(); ++idx){
      const auto& gentrk = reco::GenParticleRef(genpars, idx);
      std::vector<reco::GenParticleRef> dauVec;
      if (gentrk.isNull()) continue; //check gen particle ref

      /*if (eventNb==67099007&&std::abs(gentrk->pdgId())==13) {
          cout << "!gentrk.isNull()" << endl;
          cout << "absId, stId: " << std::abs(gentrk->pdgId()) << ", " << gentrk->status() << endl;
          cout << "gen eta phi: " << gentrk->eta() << " " << gentrk->phi() << endl;
          cout << "trk eta phi: " << track->eta() << " " << track->phi() << endl;
      }*/

      int absId = std::abs(gentrk->pdgId());
      if (std::find(pdgIds_.begin(), pdgIds_.end(), absId) == pdgIds_.end()) continue;
      /*int stId = gentrk->status();
      if (std::find(statusIds_.begin(), statusIds_.end(), stId) == statusIds_.end()) continue;

      float dEta = track->eta() - gentrk->eta();
      float dPhi = reco::deltaPhi(track->phi(), gentrk->phi());
      float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
      */
      double dR = reco::deltaR(track->eta(), track->phi(), gentrk->eta(), gentrk->phi());
      hdR->Fill(dR);

      if (dR < maxDR_) {
        isGenMuon = true;
        break;
	    }
    }

    if (isGenMuon == true) continue;

    trkqx += pt*cos(2*phi);
    trkqy += pt*sin(2*phi);
    trkPt += pt;

    trkqx_v3 += pt*cos(3*phi);
    trkqy_v3 += pt*sin(3*phi);

    if (eta>0.05&&eta<2.4) {
      trkPt_forw += pt;
      trkqx_forw += pt*cos(2*phi);
      trkqy_forw += pt*sin(2*phi);
      trkqx_v3_forw += pt*cos(3*phi);
      trkqy_v3_forw += pt*sin(3*phi);
      htrketa_forw->Fill(eta);
	  }
    if (eta>-2.4&&eta<-0.05) {
      trkPt_afterw += pt;
      trkqx_afterw += pt*cos(2*phi);
      trkqy_afterw += pt*sin(2*phi);
      trkqx_v3_afterw += pt*cos(3*phi);
      trkqy_v3_afterw += pt*sin(3*phi);
     	htrketa_afterw->Fill(eta);
	  }
  }
  trkQx = trkqx/trkPt;
  trkQy = trkqy/trkPt;
  all_trkQx = all_trkqx/all_trkPt;
  all_trkQy = all_trkqy/all_trkPt;
  trkQx_v3 = trkqx_v3/trkPt;
  trkQy_v3 = trkqy_v3/trkPt;

  trkQx_forw = trkqx_forw/trkPt_forw;
  trkQy_forw = trkqy_forw/trkPt_forw;
  trkQx_v3_forw = trkqx_v3_forw/trkPt_forw;
  trkQy_v3_forw = trkqy_v3_forw/trkPt_forw;
  trkQx_afterw = trkqx_afterw/trkPt_afterw;
  trkQy_afterw = trkqy_afterw/trkPt_afterw;
  trkQx_v3_afterw = trkqx_v3_afterw/trkPt_afterw;
  trkQy_v3_afterw = trkqy_v3_afterw/trkPt_afterw;
}


// ------------ method called once each job just before starting event
//loop  ------------
void
HiEventPlaneTrackMC::beginJob()
{
  TH1D::SetDefaultSumw2();

  // Check inputs
  if(!doRecoNtuple_) throw cms::Exception("PATCompositeAnalyzer") << "No output for RECO Fix config!!" << std::endl;
  if(saveTree_) initTree();
  if(saveHistogram_) initHistogram();
  
}


void 
HiEventPlaneTrackMC::initTree()
{ 
  PATEventPlaneNtuple = fs->make< TTree>("EventPlane","EventPlane");

  if(doRecoNtuple_)
  {
    // Event info
    
    PATEventPlaneNtuple->Branch("RunNb",&runNb,"RunNb/i");
    PATEventPlaneNtuple->Branch("LSNb",&lsNb,"LSNb/i");
    PATEventPlaneNtuple->Branch("EventNb",&eventNb,"EventNb/i");
    PATEventPlaneNtuple->Branch("nPV",&nPV,"nPV/S");
    PATEventPlaneNtuple->Branch("bestvtxX",&bestvx,"bestvtxX/F");
    PATEventPlaneNtuple->Branch("bestvtxY",&bestvy,"bestvtxY/F");
    PATEventPlaneNtuple->Branch("bestvtxZ",&bestvz,"bestvtxZ/F");
    
    if(isCentrality_) 
    {
      PATEventPlaneNtuple->Branch("centrality",&centrality,"centrality/S");
      PATEventPlaneNtuple->Branch("Ntrkoffline",&Ntrkoffline,"Ntrkoffline/I");
      PATEventPlaneNtuple->Branch("NtrkHP",&NtrkHP,"NtrkHP/I");
    }
    
    PATEventPlaneNtuple->Branch("trkQx",&trkQx,"trkQx/D");
    PATEventPlaneNtuple->Branch("trkQy",&trkQy,"trkQy/D");
    PATEventPlaneNtuple->Branch("all_trkQx",&all_trkQx,"all_trkQx/D");
    PATEventPlaneNtuple->Branch("all_trkQy",&all_trkQy,"all_trkQy/D");
    PATEventPlaneNtuple->Branch("trkQx_forw",&trkQx_forw,"trkQx_forw/D");
    PATEventPlaneNtuple->Branch("trkQy_forw",&trkQy_forw,"trkQy_forw/D");
    PATEventPlaneNtuple->Branch("trkQx_afterw",&trkQx_afterw,"trkQx_afterw/D");
    PATEventPlaneNtuple->Branch("trkQy_afterw",&trkQy_afterw,"trkQy_afterw/D");
    PATEventPlaneNtuple->Branch("trkQx_v3",&trkQx_v3,"trkQx_v3/D");
    PATEventPlaneNtuple->Branch("trkQy_v3",&trkQy_v3,"trkQy_v3/D");
    PATEventPlaneNtuple->Branch("trkQx_v3_forw",&trkQx_v3_forw,"trkQx_v3_forw/D");
    PATEventPlaneNtuple->Branch("trkQy_v3_forw",&trkQy_v3_forw,"trkQy_v3_forw/D");
    PATEventPlaneNtuple->Branch("trkQx_v3_afterw",&trkQx_v3_afterw,"trkQx_v3_afterw/D");
    PATEventPlaneNtuple->Branch("trkQy_v3_afterw",&trkQy_v3_afterw,"trkQy_v3_afterw/D");

  } // doRecoNtuple_

}

void
HiEventPlaneTrackMC::initHistogram()
{
  hdR = fs->make<TH1D>("hdR",";dR",1000000,0,10);
  htrketa_forw = fs->make<TH1D>("htrketa_forw",";Track Eta in forward",500,-2.5,2.5);
  htrketa_afterw = fs->make<TH1D>("htrketa_afterw",";Track Eta in afterward",500,-2.5,2.5);
  htrkpt = fs->make<TH1D>("htrkpt", ";track pt", 100,0,5);
  htrketa = fs->make<TH1D>("htrketa",";track eta", 300, -3.,3.);
}


//--------------------------------------------------------------------------------------------------
void 
HiEventPlaneTrackMC::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
}


// ------------ method called once each job just after ending the event
//loop  ------------
void 
HiEventPlaneTrackMC::endJob()
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiEventPlaneTrackMC);
