#include <memory>
#include <mutex>
#include <thread>
#include <array>
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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

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


  TTree* trees[5];
  std::array<std::string, 5> treeNames = {{"RecoMuons", "TightRPCMuons", "RPCMuons", "TightGEMMuons", "GEMMuons"}};
  mutable std::array<std::mutex, 5> mtx;
  
  static thread_local int eventNumber;
  static thread_local int runNumber;
  static thread_local int lumiSection;
  static thread_local int muonSize;
  static thread_local std::vector<double> muonPt;
  static thread_local std::vector<double> muonEta;
  static thread_local std::vector<double> muonPhi;
  static thread_local std::vector<double> muonIso;
  static thread_local double zBosonMass;
  static thread_local double vertexdz;


  static thread_local std::vector<reco::Muon> recoMuons;
  static thread_local std::vector<reco::Muon> rpcTightMuons;
  static thread_local std::vector<reco::Muon> rpcMuons;
  static thread_local std::vector<reco::Muon> gemTightMuons;
  static thread_local std::vector<reco::Muon> gemMuons;
  static thread_local std::vector<reco::Muon> muonLists[];

};


thread_local std::vector<reco::Muon> Analysis::recoMuons;
thread_local std::vector<reco::Muon> Analysis::rpcTightMuons;
thread_local std::vector<reco::Muon> Analysis::rpcMuons;
thread_local std::vector<reco::Muon> Analysis::gemTightMuons;
thread_local std::vector<reco::Muon> Analysis::gemMuons;
thread_local std::vector<reco::Muon> muonLists[5];
thread_local int Analysis::eventNumber = -1;
thread_local int Analysis::runNumber = -1;
thread_local int Analysis::lumiSection = -1;
thread_local int Analysis::muonSize = -1;
thread_local std::vector<double> Analysis::muonPt;
thread_local std::vector<double> Analysis::muonEta;
thread_local std::vector<double> Analysis::muonPhi;
thread_local std::vector<double> Analysis::muonIso;
thread_local double Analysis::zBosonMass = -1.0;
thread_local double Analysis::vertexdz = -1.0;

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
  for (size_t i = 0; i < treeNames.size(); ++i) {
    trees[i] = fs->make<TTree>(treeNames[i].c_str(), treeNames[i].c_str());
    trees[i]->Branch("zBosonMass", &zBosonMass, "zBosonMass/D");
    trees[i]->Branch("vertexdz", &vertexdz, "vertexdz/D");
    trees[i]->Branch("eventNumber", &eventNumber, "eventNumber/I");
    trees[i]->Branch("runNumber", &runNumber, "runNumber/I");
    trees[i]->Branch("lumiSection", &lumiSection, "lumiSection/I");
    trees[i]->Branch("muonSize", &muonSize, "muonSize/I");
    trees[i]->Branch("muonPt", &muonPt);
    trees[i]->Branch("muonEta", &muonEta);
    trees[i]->Branch("muonPhi", &muonPhi);
    trees[i]->Branch("muonIso", &muonIso);
  }

  std::cout << "beginJob" << std::endl;
}

void Analysis::analyze(const edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  // Clear vectors for new event
  muonPt.clear();
  muonEta.clear();
  muonPhi.clear();
  muonIso.clear();
  zBosonMass = -1.0;
  vertexdz = 99.0;
  eventNumber = -1;
  runNumber = -1;
  lumiSection = -1;
  muonSize = 0;
  
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

  std::vector<reco::Muon> recoMuons;
  std::vector<reco::Muon> rpcTightMuons;
  std::vector<reco::Muon> rpcMuons;
  std::vector<reco::Muon> gemTightMuons;
  std::vector<reco::Muon> gemMuons;
  std::vector<reco::Muon> muonLists[] =  {recoMuons, rpcTightMuons, rpcMuons, gemTightMuons, gemMuons};

  
  eventNumber = iEvent.id().event();
  runNumber = iEvent.id().run();
  lumiSection = iEvent.luminosityBlock();

  // Muon selection
  for (const auto& muon : *muons) {
    if (muon.pt() > 24 && fabs(muon.eta()) < 2.4) {
        double iso = (muon.pfIsolationR04().sumChargedHadronPt +
                      std::max(0.0, muon.pfIsolationR04().sumNeutralHadronEt +
                      muon.pfIsolationR04().sumPhotonEt -
                      0.5 * muon.pfIsolationR04().sumPUPt)) /
                      muon.pt();
      muonPt.push_back(muon.pt());
      muonEta.push_back(muon.eta());
      muonPhi.push_back(muon.phi());
      muonIso.push_back(iso);

      muonLists[0].push_back(muon);  // Reco Muons
      if (muon.isRPCMuon()) { // RPC Muons
        muonLists[2].push_back(muon);
        if (iso < 0.15 && muon.passed(reco::Muon::CutBasedIdTight)) {
          muonLists[1].push_back(muon);  // Tight RPC
        }
      }
      if (muon.isGEMMuon()) { // GEM Muons
        muonLists[4].push_back(muon);
        if (iso < 0.15 && muon.passed(reco::Muon::CutBasedIdTight)) {
          muonLists[3].push_back(muon);  // Tight GEM
        }  
      }            
    }
  }

  // Z Mass reconstruction
  auto reconstructZBoson = [](const std::vector<reco::Muon>& muons) -> std::pair<double, double> {
    if (muons.size() < 2) return {0.0, 0.0};
    double bestMass = -1.0, minDvz = std::numeric_limits<double>::max();
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
    return {bestMass, minDvz};
  };
  
  for (size_t i = 0; i < 5; ++i) {
    auto [zMass, dz] = reconstructZBoson(muonLists[i]);
    zBosonMass = zMass; 
    vertexdz = dz;
    std::lock_guard<std::mutex> lock(mtx[i]);
    muonSize = muonLists[i].size();
    trees[i]->Fill();
  }

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
