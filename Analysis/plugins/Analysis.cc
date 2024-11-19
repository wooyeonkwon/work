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
  TTree* treeRPC;
  TTree* treeGEM;
  TTree* treeReco;

  mutable double zBosonMass;
  mutable int eventNumber;
  mutable int runNumber;
  mutable int lumiSection;


  // Reco muon variables
  mutable std::vector<double> muonPt;
  mutable std::vector<double> muonEta;
  mutable std::vector<double> muonPhi;
  mutable std::vector<double> muonIso;
  mutable int muonSize;

  mutable double vertexdz;
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
  treeTracker = new TTree("TrackerMuons", "Tracker Muons");
  treeStandAlone = new TTree("StandAloneMuons", "StandAlone Muons");
  treeRPC = new TTree("RPCMuons", "RPC Muons");
  treeGEM = new TTree("GEMMuons", "GEM Muons");
  treeReco = new TTree("RecoMuons", "Reco Muons");

  auto initializeTree = [this](TTree* tree) {
    tree->Branch("zBosonMass", &zBosonMass, "zBosonMass/D");
    tree->Branch("vertexdz", &vertexdz, "vertexdz/D");
    tree->Branch("eventNumber", &eventNumber, "eventNumber/I");
    tree->Branch("runNumber", &runNumber, "runNumber/I");
    tree->Branch("lumiSection", &lumiSection, "lumiSection/I");  
  };

  initializeTree(treeGlobal);
  initializeTree(treeTracker);
  initializeTree(treeStandAlone);
  initializeTree(treeRPC);
  initializeTree(treeGEM);

  // Reco branches - now using vectors
  treeReco->Branch("muonSize", &muonSize, "muonSize/I");
  treeReco->Branch("muonPt", &muonPt);
  treeReco->Branch("muonEta", &muonEta);
  treeReco->Branch("muonPhi", &muonPhi);
  treeReco->Branch("muonIso", &muonIso);
}

void Analysis::analyze(const edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const {
    std::lock_guard<std::mutex> lock(mtx_);

    edm::Handle<reco::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vertexToken_, vertices);



    // Check if vertices exist and are non-empty

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

    if (!passTrigger) return;

  std::vector<reco::Muon> globalMuons;
  std::vector<reco::Muon> trackerMuons;
  std::vector<reco::Muon> standAloneMuons;
  std::vector<reco::Muon> rpcMuons;
  std::vector<reco::Muon> gemMuons;
  
  eventNumber = iEvent.id().event();
  runNumber = iEvent.id().run();
  lumiSection = iEvent.luminosityBlock();

  // Clear vectors for new event
  muonPt.clear();
  muonEta.clear();
  muonPhi.clear();
  muonIso.clear();

  // Muon selection
  for (const auto& muon : *muons) {
    if (muon.pt() > 24 && fabs(muon.eta()) < 2.4) {
      double iso = (muon.pfIsolationR04().sumChargedHadronPt + 
                   std::max(0.0, muon.pfIsolationR04().sumNeutralHadronEt + 
                           muon.pfIsolationR04().sumPhotonEt - 
                           0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt();
      
      if (iso < 0.15 && muon.isGlobalMuon() && muon.isRPCMuon()) {  // Added isTightMuon check
        // Fill Reco muon information
        muonPt.push_back(muon.pt());
        muonEta.push_back(muon.eta());
        muonPhi.push_back(muon.phi());
        muonIso.push_back(iso);
        globalMuons.push_back(muon);
        if (muon.isTrackerMuon()) {
          trackerMuons.push_back(muon);
        }
        if (muon.isStandAloneMuon()) {
          standAloneMuons.push_back(muon);
        }
        if (muon.isRPCMuon()) {
          rpcMuons.push_back(muon);
        }
        if (muon.isGEMMuon()) {
          gemMuons.push_back(muon);
        }
      }
    }
  }

  muonSize = muonPt.size();
  if (muonSize > 0) {
    treeReco->Fill();
  }

  // Z reconstruction  
  auto reconstructZBoson = [](const std::vector<reco::Muon>& muons) -> std::pair<double, double> { // Change return type to pair
    if (muons.size() < 2) return {0.0, 0.0}; // Return pair with zero values

    math::XYZTLorentzVector goodZBoson;
    double goodMass = 0.0;
    double dvz = 0.0;
    double minDvz = std::numeric_limits<double>::max(); // Initialize minDvz

    for (size_t i = 0; i < muons.size() - 1; ++i) {
      for (size_t j = i + 1; j < muons.size(); ++j) {
        dvz = fabs(muons[i].vz() - muons[j].vz());
        if ((muons[i].charge() + muons[j].charge() == 0) && (dvz < 0.5)) {
          math::XYZTLorentzVector zBoson = muons[i].p4() + muons[j].p4();
          double mass = zBoson.M();
          if (fabs(mass - 91.1876) < fabs(goodMass - 91.1876)) {
            goodMass = mass;
            goodZBoson = zBoson;
            minDvz = dvz; // Store the minimum dvz
          }
        }
      }
    }
    
    return {goodMass, minDvz}; // Return pair of goodMass and minDvz
  };


  {
    auto [zBosonMassResult, vertexdzResult] = reconstructZBoson(globalMuons);
    zBosonMass = zBosonMassResult;
    vertexdz = vertexdzResult;
    if (zBosonMass > 0.0) {
    treeGlobal->Fill();    
    }   
  }

  {
    auto [zBosonMassResult, vertexdzResult] = reconstructZBoson(trackerMuons);
    zBosonMass = zBosonMassResult;
    vertexdz = vertexdzResult;
    if (zBosonMass > 0.0) {
      treeTracker->Fill();    
    }
  }

  {
    auto [zBosonMassResult, vertexdzResult] = reconstructZBoson(standAloneMuons);
    zBosonMass = zBosonMassResult;
    vertexdz = vertexdzResult;
    if (zBosonMass > 0.0) {
      treeStandAlone->Fill();    
    }
  }

  {
    auto [zBosonMassResult, vertexdzResult] = reconstructZBoson(rpcMuons);
    zBosonMass = zBosonMassResult;
    vertexdz = vertexdzResult; 
    if (zBosonMass > 0.0) {
      treeRPC->Fill();    
    }
  }

  {
    auto [zBosonMassResult, vertexdzResult] = reconstructZBoson(gemMuons);
    zBosonMass = zBosonMassResult;
    vertexdz = vertexdzResult;
    if (zBosonMass > 0.0) {
      treeGEM->Fill();    
    }
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