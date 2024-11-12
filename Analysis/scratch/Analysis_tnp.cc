#include <memory>
#include <mutex>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdlib>

class Analysis : public edm::global::EDAnalyzer<> {
public:
  explicit Analysis(const edm::ParameterSet&);
  ~Analysis() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::StreamID, const edm::Event&, const edm::EventSetup&) const override;
  void endJob() override;

  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  std::string hltPath_;

  TTree* treeGlobal;
  TTree* treeTracker;
  TTree* treeStandAlone;
  TTree* treeCalo;
  TTree* treePF;
  TTree* treeRPC;
  TTree* treeGEM;
  TTree* treeME0;

  mutable double zBosonMassGlobal;
  mutable double zBosonMassTracker;
  mutable double zBosonMassStandAlone;
  mutable double zBosonMassCalo;
  mutable double zBosonMassPF;
  mutable double zBosonMassRPC;
  mutable double zBosonMassGEM;  
  mutable double zBosonMassME0;

  mutable int eventNumber;
  mutable int runNumber;
  mutable int lumiSection;

  mutable int muonSizeGlobal;
  mutable int muonSizeTracker;
  mutable int muonSizeStandAlone;
  mutable int muonSizeCalo;
  mutable int muonSizePF;
  mutable int muonSizeRPC;
  mutable int muonSizeGEM;
  mutable int muonSizeME0;

  mutable std::mutex mtx_;
};

Analysis::Analysis(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
      vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      hltPath_(iConfig.getParameter<std::string>("hltPath")) {
}

Analysis::~Analysis() {
  std::lock_guard<std::mutex> lock(mtx_);
}

void Analysis::beginJob() {
  std::lock_guard<std::mutex> lock(mtx_);


  treeGlobal = new TTree("GlobalMuons", "Global Muons");
  if (!treeGlobal) {
    throw std::runtime_error("Failed to create treeGlobal");
  }
  treeGlobal = new TTree("GlobalMuons", "Global Muons");
  treeTracker = new TTree("TrackerMuons", "Tracker Muons");
  treeStandAlone = new TTree("StandAloneMuons", "StandAlone Muons");
  treeCalo = new TTree("CaloMuons", "Calo Muons");
  treePF = new TTree("PFMuons", "PF Muons");
  treeRPC = new TTree("RPCMuons", "RPC Muons");
  treeGEM = new TTree("GEMMuons", "GEM Muons");
  treeME0 = new TTree("ME0Muons", "ME0 Muons");

  // 
  auto initializeTree = [](TTree* tree, const std::string& name, const std::string& title) {
    tree->Branch("zBosonMass", &zBosonMassGlobal, "zBosonMass/D");
    tree->Branch("tagPt", &tagPt, "tagPt/D");
    tree->Branch("tagEta", &tagEta, "tagEta/D");
    tree->Branch("tagIso", &tagIso, "tagIso/D");
    tree->Branch("probePt", &probePt, "probePt/D");
    tree->Branch("probeEta", &probeEta, "probeEta/D");
    tree->Branch("probeIso", &probeIso, "probeIso/D");
    tree->Branch("passingProbe", &passingProbe, "passingProbe/O");
    tree->Branch("eventNumber", &eventNumber, "eventNumber/I");
    tree->Branch("runNumber", &runNumber, "runNumber/I");
    tree->Branch("lumiSection", &lumiSection, "lumiSection/I");
    tree->Branch("muonSize", &muonSizeGlobal, "muonSize/I");
  };

  initializeTree(treeGlobal, "GlobalMuons", "Global Muons");
  initializeTree(treeTracker, "TrackerMuons", "Tracker Muons");
  initializeTree(treeStandAlone, "StandAloneMuons", "StandAlone Muons");
  initializeTree(treeCalo, "CaloMuons", "Calo Muons");
  initializeTree(treePF, "PFMuons", "PF Muons");
  initializeTree(treeRPC, "RPCMuons", "RPC Muons");
  initializeTree(treeGEM, "GEMMuons", "GEM Muons");
  initializeTree(treeME0, "ME0Muons", "ME0 Muons");
}

void Analysis::analyze(const edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerToken_, triggerResults);

  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerResults);
  bool passTrigger = false;

  // Check if the event passes the HLT_IsoMu24_v trigger
  for (unsigned int i = 0; i < triggerResults->size(); ++i) {
    if (triggerNames.triggerName(i).find(hltPath_) != std::string::npos) {
      if (triggerResults->accept(i)) {
        passTrigger = true;
        break;
      }
    }
  }

  if (!passTrigger) return;  // Skip event if it doesn't pass the trigger

  reco::Vertex primaryVertex;
  if (!vertices->empty()) {
    primaryVertex = vertices->front();  // get first good vertex
  }
  std::vector<reco::Muon> tagMuons;
  // muons below are probe muons
  std::vector<reco::Muon> globalMuons;
  std::vector<reco::Muon> trackerMuons;
  std::vector<reco::Muon> standAloneMuons;
  std::vector<reco::Muon> caloMuons;
  std::vector<reco::Muon> pfMuons;
  std::vector<reco::Muon> rpcMuons;
  std::vector<reco::Muon> gemMuons;
  std::vector<reco::Muon> me0Muons;
  
  eventNumber = iEvent.id().event();
  runNumber = iEvent.id().run();
  lumiSection = iEvent.luminosityBlock();


  //muon selection
  for (const auto& muon : *muons) {
    if (muon.pt() > 24 && fabs(muon.eta()) < 2.4 && muon::isTightMuon(muon, primaryVertex)) {
      float isolation = (muon.pfIsolationR04().sumChargedHadronPt + 
                         ((muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) > 0 ? 
                          (muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) : 0)) / muon.pt();
      if (isolation < 0.15) {
        tagMuons.push_back(muon);
      }
    }
    if (muon.pt() > 15 && fabs(muon.eta()) < 2.4){
      if (muon.isGlobalMuon()){
        globalMuons.push_back(muon);
        if (muon.isTrackerMuon()) {
          trackerMuons.push_back(muon);
        }
        if (muon.isStandAloneMuon()) {
          standAloneMuons.push_back(muon);
          }
        if (muon.isCaloMuon()) {
          caloMuons.push_back(muon);
          }
        if (muon.isPFMuon()) {
          pfMuons.push_back(muon);
          }
        if (muon.isRPCMuon()) {
          rpcMuons.push_back(muon); 
          }
        if (muon.isGEMMuon()) {
          gemMuons.push_back(muon);
          }
        if (muon.isME0Muon()) {
          me0Muons.push_back(muon);  
        }
      }            
    }
  }




  auto reconstructZBoson = [&](const std::vector<reco::Muon>& tagMuons, const std::vector<reco::Muon>& probeMuons) -> double {
    if (tagMuons.empty() || probeMuons.empty()) return 0.0;
    math::XYZTLorentzVector goodZBoson;
    double goodMass = 0.0;
    
    for (const auto& tagMuon : tagMuons) {
      for (const auto& probeMuon : probeMuons) {
        if ((tagMuon.charge() != probeMuon.charge()) && (fabs(tagMuon.vz() - probeMuon.vz()) < 0.5)) {
          math::XYZTLorentzVector zBoson = tagMuon.p4() + probeMuon.p4();
          double mass = zBoson.M();
          if (fabs(mass - 91.1876) < fabs(goodMass - 91.1876)) {
            goodMass = mass;
            goodZBoson = zBoson;
            tagPt = tagMuon.pt();
            tagEta = tagMuon.eta();
            tagIso = (tagMuon.pfIsolationR04().sumChargedHadronPt + 
                      ((tagMuon.pfIsolationR04().sumNeutralHadronEt + tagMuon.pfIsolationR04().sumPhotonEt - 0.5 * tagMuon.pfIsolationR04().sumPUPt) > 0 ? 
                       (tagMuon.pfIsolationR04().sumNeutralHadronEt + tagMuon.pfIsolationR04().sumPhotonEt - 0.5 * tagMuon.pfIsolationR04().sumPUPt) : 0)) / tagMuon.pt();
            probePt = probeMuon.pt();
            probeEta = probeMuon.eta();
            probeIso = (probeMuon.pfIsolationR04().sumChargedHadronPt + 
                        ((probeMuon.pfIsolationR04().sumNeutralHadronEt + probeMuon.pfIsolationR04().sumPhotonEt - 0.5 * probeMuon.pfIsolationR04().sumPUPt) > 0 ? 
                         (probeMuon.pfIsolationR04().sumNeutralHadronEt + probeMuon.pfIsolationR04().sumPhotonEt - 0.5 * probeMuon.pfIsolationR04().sumPUPt) : 0)) / probeMuon.pt();
            passingProbe = 1;
          }
        }
      }
    }
    return goodMass;
  };

  zBosonMassGlobal = reconstructZBoson(tagMuons,globalMuons);
  zBosonMassTracker = reconstructZBoson(tagMuons,trackerMuons);
  zBosonMassStandAlone = reconstructZBoson(tagMuons,standAloneMuons);
  zBosonMassCalo = reconstructZBoson(tagMuons,caloMuons); 
  zBosonMassPF = reconstructZBoson(tagMuons,pfMuons);
  zBosonMassRPC = reconstructZBoson(tagMuons,rpcMuons);
  zBosonMassGEM = reconstructZBoson(tagMuons,gemMuons);
  zBosonMassME0 = reconstructZBoson(tagMuons,me0Muons);

  muonSizeGlobal = globalMuons.size();
  muonSizeTracker = trackerMuons.size();
  muonSizeStandAlone = standAloneMuons.size();
  muonSizeCalo = caloMuons.size();
  muonSizePF = pfMuons.size();
  muonSizeRPC = rpcMuons.size();
  muonSizeGEM = gemMuons.size();
  muonSizeME0 = me0Muons.size();

  if (zBosonMassGlobal > 0.0) {
    treeGlobal->Fill();
  }
  
  if (zBosonMassTracker > 0.0) {
    treeTracker->Fill();
  }
  
  if (zBosonMassStandAlone > 0.0) {
    treeStandAlone->Fill();
  }
  
  if (zBosonMassCalo > 0.0) {
    treeCalo->Fill();
  }
  
  if (zBosonMassPF > 0.0) {
    treePF->Fill();
  }
  
  if (zBosonMassRPC > 0.0) {
    treeRPC->Fill();
  }
  if (zBosonMassGEM > 0.0) {
    treeGEM->Fill();
  }

  if (zBosonMassME0 > 0.0) {
    treeME0->Fill();
  }
}

void Analysis::endJob() {
  std::lock_guard<std::mutex> lock(mtx_);
  std::cout << "endJob" << std::endl;
}

void Analysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muons", edm::InputTag("muons"));
  desc.add<edm::InputTag>("triggerResults", edm::InputTag("TriggerResults"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<std::string>("hltPath", "HLT_IsoMu24_v*");
  descriptions.add("analysis", desc);
}

DEFINE_FWK_MODULE(Analysis);


 
