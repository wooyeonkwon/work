#include <memory>
#include <mutex>
#include <thread>
#include <array>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <set>
#include <iostream>

class Analysis : public edm::one::EDAnalyzer<> {
public:
  explicit Analysis(const edm::ParameterSet&);
  ~Analysis() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  std::string hltPath_;
  
  TTree* tree;

  mutable std::mutex mtx;

  static thread_local int eventNumber;
  static thread_local int runNumber;
  static thread_local int lumiSection;
  static thread_local std::array<double, 6> zBosonMass;
  static thread_local std::array<double, 6> vertexdz;

  static thread_local std::vector<reco::Muon> recoMuons;
  static thread_local std::vector<reco::Muon> recoTightMuons;
  static thread_local std::vector<reco::Muon> rpcMuons;
  static thread_local std::vector<reco::Muon> rpcTightMuons;  
  static thread_local std::vector<reco::Muon> gemMuons;
  static thread_local std::vector<reco::Muon> gemTightMuons;  
  static thread_local std::vector<reco::Muon> muonLists[];

  std::set<int> eventNumbers;

};
thread_local int Analysis::eventNumber = -1;
thread_local int Analysis::runNumber = -1;
thread_local int Analysis::lumiSection = -1;
thread_local std::array<double, 6> Analysis::zBosonMass;
thread_local std::array<double, 6> Analysis::vertexdz;
thread_local std::vector<reco::Muon> Analysis::recoMuons;
thread_local std::vector<reco::Muon> Analysis::recoTightMuons;
thread_local std::vector<reco::Muon> Analysis::rpcMuons;
thread_local std::vector<reco::Muon> Analysis::rpcTightMuons;
thread_local std::vector<reco::Muon> Analysis::gemMuons;
thread_local std::vector<reco::Muon> Analysis::gemTightMuons; 
thread_local std::vector<reco::Muon> muonLists[6];

Analysis::Analysis(const edm::ParameterSet& iConfig)
  : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
    vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    hltPath_(iConfig.getParameter<std::string>("hltPath")) {
}

Analysis::~Analysis() {
}

void Analysis::beginJob() {
  edm::Service<TFileService> fs;
  //TFile &file = fs->file();
  //file.SetCompressionLevel(3);
  tree = fs->make<TTree>("Analysis", "skimed Muon and reconstructed Z");
  tree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  tree->Branch("runNumber", &runNumber, "runNumber/I");
  tree->Branch("lumiSection", &lumiSection, "lumiSection/I");
  tree->Branch("zBosonMass", &zBosonMass);
  tree->Branch("vertexdz", &vertexdz);
  tree->Branch("recoMuons", &recoMuons);

  tree->SetAutoFlush(10000);
  tree->SetCacheSize(10000000);
  tree->SetMaxVirtualSize(1000000);
  std::cout << "beginJob" << std::endl;

}

void Analysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Clear vectors for new event
  zBosonMass = {{-1.0,-1.0,-1.0,-1.0,-1.0,-1.0}};
  vertexdz = {{-1.0,-1.0,-1.0,-1.0,-1.0,-1.0}};
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

  recoMuons.clear();
  recoTightMuons.clear();
  rpcMuons.clear();
  rpcTightMuons.clear();
  gemMuons.clear();
  gemTightMuons.clear();
  std::vector<reco::Muon>* muonLists[] = {&recoMuons, &recoTightMuons, &rpcMuons, &rpcTightMuons, &gemMuons, &gemTightMuons};

  // Muon selection
  for (const auto& muon : *muons) {
    if (muon.pt() > 24 && fabs(muon.eta()) < 2.4) {
      muonLists[0]->push_back(muon);  // Reco Muons
      if (muon.passed(reco::Muon::PFIsoTight) && muon.passed(reco::Muon::CutBasedIdTight)) {
        muonLists[1]->push_back(muon);  // Tight RPC
      }
      if (muon.isRPCMuon()) { // RPC Muons
        muonLists[2]->push_back(muon);
        if (muon.passed(reco::Muon::PFIsoTight) && muon.passed(reco::Muon::CutBasedIdTight)) {
          muonLists[3]->push_back(muon);  // Tight RPC
        }
      }
      if (muon.isGEMMuon()) { // GEM Muons
        muonLists[4]->push_back(muon);
        if (muon.passed(reco::Muon::PFIsoTight) && muon.passed(reco::Muon::CutBasedIdTight)) {
          muonLists[5]->push_back(muon);  // Tight GEM
        }  
      }            
    }
  }

  // Z Mass reconstruction
  auto reconstructZBoson = [](const std::vector<reco::Muon>& muons) -> std::pair<double, double> {
    if (muons.size() < 2) return {-1.0, -1.0};
    double bestMass = std::numeric_limits<double>::max();
    double minDvz = std::numeric_limits<double>::max();
    for (size_t i = 0; i < muons.size() - 1; ++i) {
      for (size_t j = i + 1; j < muons.size(); ++j) {
        if (muons[i].charge() + muons[j].charge() != 0) continue;
        double dvz = fabs(muons[i].vz() - muons[j].vz());
        math::XYZTLorentzVector zBoson = muons[i].p4() + muons[j].p4();
        double mass = zBoson.M();
        if (fabs(mass - 91.1876) < fabs(bestMass - 91.1876)) {
          bestMass = mass;
          minDvz = dvz;
        }
      }
    }
    if ((bestMass == std::numeric_limits<double>::max()) || (minDvz == std::numeric_limits<double>::max())) return {-1.0, -1.0};
    return {bestMass, minDvz};
  };
  
  
  for (size_t i = 0; i < 6; ++i) {
    std::tie(zBosonMass[i], vertexdz[i]) = reconstructZBoson(*muonLists[i]);
  }
  
  eventNumber = iEvent.id().event();
  runNumber = iEvent.id().run();
  lumiSection = iEvent.luminosityBlock();
  tree->Fill();
//
//  {  
//    bool i = false;
//    int attempt = 0;
//    while (!i) {
//      std::unique_lock<std::mutex> lock(mtx, std::defer_lock);
//      auto start = std::chrono::steady_clock::now();
//      while (true) {
//        if (lock.try_lock()) {
//          std::this_thread::sleep_for(std::chrono::milliseconds(10));
//          if (i == false) tree->Fill();
//          std::this_thread::sleep_for(std::chrono::milliseconds(10));
//          lock.unlock();
//          i = true;
//          break;
//        }
//        ++attempt;
//        // check timeout
//        if (std::chrono::steady_clock::now() - start > std::chrono::milliseconds(100)) {
//          break;
//        }
//      }
//      std::this_thread::sleep_for(std::chrono::milliseconds(10));
//      break;
//    }
//
//
//  }
//

    
}

void Analysis::endJob() {
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
