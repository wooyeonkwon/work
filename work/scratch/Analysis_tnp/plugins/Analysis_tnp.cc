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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

class AnalysisTnp : public edm::global::EDAnalyzer<> {
public:
  explicit AnalysisTnp(const edm::ParameterSet&);
  ~AnalysisTnp() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::StreamID, const edm::Event&, const edm::EventSetup&) const override;
  void endJob() override;

  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  std::string hltPath_;

  TTree* treeRPCPass;
  TTree* treeRPCFail;
  TTree* treeRPCEff;
  TTree* treeGEMPass;
  TTree* treeGEMFail;
  TTree* treeGEMEff;

  mutable double zBosonMass;
  mutable double efficiency;

  mutable int muonSizeRPC;
  mutable int muonSizeGEM;

  mutable std::mutex mtx_;
};

AnalysisTnp::AnalysisTnp(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
      vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      hltPath_(iConfig.getParameter<std::string>("hltPath")) {}

AnalysisTnp::~AnalysisTnp() {
  std::lock_guard<std::mutex> lock(mtx_);
}

void AnalysisTnp::beginJob() {
  std::lock_guard<std::mutex> lock(mtx_);
  edm::Service<TFileService> fs;

  treeRPCPass = fs->make<TTree>("RPCPassMuons", "RPC Pass Muons");
  treeRPCPass->Branch("zBosonMass", &zBosonMass, "zBosonMass/D");

  treeRPCFail = fs->make<TTree>("RPCFailMuons", "RPC Fail Muons");
  treeRPCFail->Branch("zBosonMass", &zBosonMass, "zBosonMass/D");

  treeRPCEff = fs->make<TTree>("RPCEffMuons", "RPC Eff Muons");
  treeRPCEff->Branch("efficiency", &efficiency, "efficiency/D");

  treeGEMPass = fs->make<TTree>("GEMPassMuons", "GEM Pass Muons");
  treeGEMPass->Branch("zBosonMass", &zBosonMass, "zBosonMass/D");

  treeGEMFail = fs->make<TTree>("GEMFailMuons", "GEM Fail Muons");
  treeGEMFail->Branch("zBosonMass", &zBosonMass, "zBosonMass/D");

  treeGEMEff = fs->make<TTree>("GEMEffMuons", "GEM Eff Muons");
  treeGEMEff->Branch("efficiency", &efficiency, "efficiency/D");

}

void AnalysisTnp::analyze(const edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerToken_, triggerResults);

  const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);
  bool passTrigger = false;

  for (unsigned int i = 0; i < triggerResults->size(); ++i) {
    if (triggerNames.triggerName(i).find(hltPath_) != std::string::npos) {
      if (triggerResults->accept(i)) {
        passTrigger = true;
        break;
      }
    }
  }

  if (!passTrigger) return;

  std::vector<const reco::Muon*> tagMuons;
  std::vector<const reco::Muon*> passingProbeRPCMuons;
  std::vector<const reco::Muon*> failingProbeRPCMuons;
  std::vector<const reco::Muon*> passingProbeGEMMuons;
  std::vector<const reco::Muon*> failingProbeGEMMuons;

  for (const auto& muon : *muons) {
    if (muon.pt() > 24 && fabs(muon.eta()) < 2.4 && muon.passed(reco::Muon::CutBasedIdTight) && muon.passed(reco::Muon::PFIsoTight)) {
      tagMuons.push_back(&muon);
    } else {
      if (muon.pt() > 15 && fabs(muon.eta()) < 2.4) {
        if (muon.isRPCMuon()){
          passingProbeRPCMuons.push_back(&muon);
        } else {
          failingProbeRPCMuons.push_back(&muon);
        }
      }
      if (muon.isGEMMuon()) {
        if(muon.pt() > 15 && fabs(muon.eta()) < 2.4){
          passingProbeGEMMuons.push_back(&muon);
        } else {
          failingProbeGEMMuons.push_back(&muon);
        }
      }
    }
  }

  auto reconstructZBoson = [](const reco::Muon* tagMuon, const reco::Muon* probeMuon) {
    if ((tagMuon->charge() != probeMuon->charge()) && (fabs(tagMuon->vz() - probeMuon->vz()) < 0.5 )) {
      auto zBoson = tagMuon->p4() + probeMuon->p4();
      return zBoson.M();
    }
    return -1.0;
  };


  int passCountRPC = 0;
  int failCountRPC = 0;
  int passCountGEM = 0;
  int failCountGEM = 0;

  for (const auto& tagMuon : tagMuons) {
    
    for (const auto& probeMuon : passingProbeRPCMuons) {
      double mass = reconstructZBoson(tagMuon, probeMuon);
      if (mass > 0) {
        zBosonMass = mass;
        treeRPCPass->Fill();
        if (mass > 60 && mass < 120) passCountRPC++;
      }
    }
    
    for (const auto& probeMuon : failingProbeRPCMuons) {
      double mass = reconstructZBoson(tagMuon, probeMuon);
      if (mass > 0) {
        zBosonMass = mass;
        treeRPCFail->Fill();
        if (mass > 60 && mass < 120) failCountRPC++;
      }
    }
    
    for (const auto& probeMuon : passingProbeGEMMuons) {
      double mass = reconstructZBoson(tagMuon, probeMuon);
      if (mass > 0) {
        zBosonMass = mass;
        treeGEMPass->Fill();
        if (mass > 60 && mass < 120) passCountGEM++;
      }
    }
    
    for (const auto& probeMuon : failingProbeGEMMuons) {
      double mass = reconstructZBoson(tagMuon, probeMuon);
      if (mass > 0) {
        zBosonMass = mass;
        treeGEMFail->Fill();
        if (mass > 60 && mass < 120) failCountGEM++;
      }
    }
  }

  if (!(passCountRPC == 0 && failCountRPC == 0)){
    efficiency = passCountRPC / static_cast<double>(passCountRPC + failCountRPC);
  } else {
    efficiency = -1;
  }
  treeRPCEff->Fill(); 

  if (!(passCountGEM == 0 && failCountGEM == 0)){
    efficiency = (passCountGEM) / static_cast<double>(passCountGEM + failCountGEM);
    treeGEMEff->Fill();
  }

}

void AnalysisTnp::endJob() {
  std::lock_guard<std::mutex> lock(mtx_);
  std::cout << "endJob" << std::endl;
}

void AnalysisTnp::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muons", edm::InputTag("muons"));
  desc.add<edm::InputTag>("triggerResults", edm::InputTag("TriggerResults"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<std::string>("hltPath", "HLT_IsoMu24_v*");
  descriptions.add("analysisTnp", desc);
}

DEFINE_FWK_MODULE(AnalysisTnp);
