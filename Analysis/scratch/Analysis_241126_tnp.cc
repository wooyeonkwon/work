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

  TTree* treeRPC;
  TTree* treeGEM;


  mutable double zBosonMassPass;
  mutable double zBosonMassFail;


  mutable double efficiency;


  mutable int muonSizeRPC;
  mutable int muonSizeGEM;

  mutable std::mutex mtx_;
};

Analysis::Analysis(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
      vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      hltPath_(iConfig.getParameter<std::string>("hltPath")),
      isTagToken_(produces<edm::ValueMap<bool>>("isTag")),
      isPassingProbeRPCToken_(produces<edm::ValueMap<bool>>("isPassingProbeRPC")),
      isPassingProbeGEMToken_(produces<edm::ValueMap<bool>>("isPassingProbeGEM")) {
        
      }


Analysis::~Analysis() {
  std::lock_guard<std::mutex> lock(mtx_);
}

void Analysis::beginJob() {
  std::lock_guard<std::mutex> lock(mtx_);

  treeRPC = new TTree("RPCMuons", "RPC Muons");
  treeGEM = new TTree("GEMMuons", "GEM Muons");

  auto initializeTree = [this](TTree* tree) {
    tree->Branch("zBosonMassPass", &zBosonMassPass, "zBosonMassPass/D");
    tree->Branch("zBosonMassFail", &zBosonMassFail, "zBosonMassFail/D");
    tree->Branch("efficiency", &efficiency,"efficiency/D");
  };

  initializeTree(treeRPC);
  initializeTree(treeGEM);

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


  std::vector<bool> isTag, isPassingProbeRPC, isPassingProbeGEM;


  for (const auto& muon : *muons) {
    bool tag = muon.pt() > 24 && fabs(muon.eta()) < 2.4 &&
               muon.passed(reco::Muon::CutBasedIdTight) &&
               muon.passed(reco::Muon::PFIsoTight);
    isTag.push_back(tag);

    if (!tag) {
      isPassingProbeRPC.push_back(muon.isRPCMuon() && muon.pt() > 15 && fabs(muon.eta()) < 2.4);
      isPassingProbeGEM.push_back(muon.isGEMMuon() && muon.pt() > 15 && fabs(muon.eta()) < 2.4);
  }

  // Make ValueMap 
  edm::EDPutTokenT<edm::ValueMap<bool>> isTagToken_;
  edm::EDPutTokenT<edm::ValueMap<bool>> isPassingProbeRPCToken_;
  edm::EDPutTokenT<edm::ValueMap<bool>> isPassingProbeGEMToken_;


  auto createValueMap = [&](const std::vector<bool>& flags, edm::EDPutTokenT<edm::ValueMap<bool>> token) {
      auto valueMap = std::make_unique<edm::ValueMap<bool>>();
      edm::ValueMap<bool>::Filler filler(*valueMap);
      filler.insert(muons, flags.begin(), flags.end());
      filler.fill();
      iEvent.put(token, std::move(valueMap));
  };

  createValueMap(isTag, "isTag");
  createValueMap(isPassingProbeRPC, "isPassingProbeRPC");
  createValueMap(isPassingProbeGEM, "isPassingProbeGEM");

  // Z reconstruction
  auto reconstructZBoson = [&](const std::vector<const reco::Muon*>& tagMuon, const std::vector<const reco::Muon*>& probeMuon,) {
    double zMass = 0.0;
    if (tagMuon->charge() != probeMuon->charge()) {
      auto zBoson = tagMuon->p4() + probeMuon->p4();
      zMass = zBoson.M();
    }
    return Zmass;
  };

  for (const auto& muon : *muons) {
      std::vector<const reco::Muon*> tagMuons;
      std::vector<const reco::Muon*> passingProbeRPCMuons;
      std::vector<const reco::Muon*> failingProbeRPCMuons;
      std::vector<const reco::Muon*> passingProbeGEMMuons;
      std::vector<const reco::Muon*> failingProbeGEMMuons;

      // 태그된 뮤온과 프로브를 분류합니다.
      for (size_t i = 0; i < muons->size(); ++i) {
          const auto& muon = (*muons)[i];

          if (isTag[i]) {
              tagMuons.push_back(&muon);
          } else {
              if (isPassingProbeRPC[i]) {
                  passingProbeRPCMuons.push_back(&muon);
              } else {
                  failingProbeRPCMuons.push_back(&muon);
              }

              if (isPassingProbeGEM[i]) {
                  passingProbeGEMMuons.push_back(&muon);
              } else {
                  failingProbeGEMMuons.push_back(&muon);
              }
          }
      }

      // 가능한 태그-프로브 조합으로 Z 보존을 재구성합니다.
      auto fillTree = [&](TTree* tree, double& branchMass, double& branchEfficiency, 
                          const std::vector<const reco::Muon*>& probes) {
          if (tagMuons.empty() || probes.empty()) return;

          std::vector<double> masses;
          for (const auto* tagMuon : tagMuons) {
              for (const auto* probeMuon : probes) {
                  double zMass = reconstructZBoson(tagMuon, probeMuon);
                  if (zMass > 0) {
                      masses.push_back(zMass);
                  }
              }
          }

          if (!masses.empty()) {
              // 평균 Z 보존 질량 계산
              branchMass = std::accumulate(masses.begin(), masses.end(), 0.0) / masses.size();
          }

          // 효율 계산 (passingProbe / 전체 probe)
          branchEfficiency = static_cast<double>(passingProbeRPCMuons.size()) / 
                            (passingProbeRPCMuons.size() + failingProbeRPCMuons.size());

          tree->Fill();
      };

      fillTree(treeRPC, zBosonMassPass, efficiency, passingProbeRPCMuons);
      fillTree(treeRPC, zBosonMassFail, efficiency, failingProbeRPCMuons);
      fillTree(treeGEM, zBosonMassPass, efficiency, passingProbeGEMMuons);
      fillTree(treeGEM, zBosonMassFail, efficiency, failingProbeGEMMuons);
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


 
