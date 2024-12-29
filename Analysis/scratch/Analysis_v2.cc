#include <memory>
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

  TFile* outputFile;
  TTree* treeGlobal;
  TTree* treeRPC;
  
  // New TTrees for dynamic pt/isolation combinations
  std::map<std::string, TTree*> globalMuonsTrees;
  std::map<std::string, TTree*> rpcMuonsTrees;

  mutable double zBosonMass;
  mutable int eventNumber;
  mutable int runNumber;
  mutable int lumiSection;
  mutable double efficiency; // New efficiency

  // Declare here the thresholds for pt and isolation
  std::vector<int> ptThresholds = {20, 30, 40, 50, 60, 70};
  std::vector<double> isolationThresholds = {0.15, 0.20, 0.25};

  // Z boson reconstruction function
  double reconstructZBoson(const std::vector<reco::Muon>& tagMuons, const std::vector<reco::Muon>& probeMuons) const;
};

Analysis::Analysis(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
      vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      hltPath_(iConfig.getParameter<std::string>("hltPath")) {
}

Analysis::~Analysis() {
}

// In the analyze function
void Analysis::analyze(const edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerToken_, triggerResults);

  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerResults);
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

  reco::Vertex primaryVertex;
  if (!vertices->empty()) {
    primaryVertex = vertices->front();
  }

  // Loop over muons with various pt and isolation cuts
  for (int ptCut : ptThresholds) {
    for (double isolationCut : isolationThresholds) {
      std::string muonBranchName = "muon_" + std::to_string(ptCut) + "_" + std::to_string(static_cast<int>(isolationCut * 100));

      std::vector<reco::Muon> selectedTagMuonsGlobal, selectedProbeMuonsGlobal;
      std::vector<reco::Muon> selectedTagMuonsRPC, selectedProbeMuonsRPC;

      for (const auto& muon : *muons) {
        float isolation = (muon.pfIsolationR04().sumChargedHadronPt + 
                         ((muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) > 0 ? 
                          (muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) : 0)) / muon.pt();

        if (muon.pt() > ptCut && isolation < isolationCut && muon::isTightMuon(muon, primaryVertex) && fabs(muon.eta()) < 2.4) {
          if (muon.isGlobalMuon()) {
            selectedTagMuonsGlobal.push_back(muon);
            if (muon.isRPCMuon()) {
              selectedTagMuonsRPC.push_back(muon);
            }
          }
        } else if (muon.pt() > 15 && fabs(muon.eta()) < 2.4) {
          if (muon.isGlobalMuon()) {
            selectedProbeMuonsGlobal.push_back(muon);
            if (muon.isRPCMuon()) {
              selectedProbeMuonsRPC.push_back(muon);
            }
          }
        }
      }

      // Z boson mass reconstruction
      double zBosonMassGlobal = reconstructZBoson(selectedTagMuonsGlobal, selectedProbeMuonsGlobal);
      double zBosonMassRPC = reconstructZBoson(selectedTagMuonsRPC, selectedProbeMuonsRPC);

      // Calculate efficiency (e.g., ratio of passing probe muons to all probe muons)
      efficiency = selectedProbeMuonsGlobal.size() > 0 ? static_cast<double>(selectedTagMuonsGlobal.size()) / selectedProbeMuonsGlobal.size() : 0.0;

      // Fill the trees for Global and RPC muons
      if (zBosonMassGlobal > 0) {
        globalMuonsTrees.at(muonBranchName)->Fill();
      }
      if (zBosonMassRPC > 0) {
        rpcMuonsTrees.at(muonBranchName)->Fill();
      }
    }
  }
}

// Function to reconstruct Z boson mass from tag and probe muons
double Analysis::reconstructZBoson(const std::vector<reco::Muon>& tagMuons, const std::vector<reco::Muon>& probeMuons) const {
  if (tagMuons.empty() || probeMuons.empty()) return 0.0;

  math::XYZTLorentzVector goodZBoson;
  double goodMass = 0.0;

  for (const auto& tagMuon : tagMuons) {
    for (const auto& probeMuon : probeMuons) {
      if ((tagMuon.charge() != probeMuon.charge()) && (fabs(tagMuon.vz()-probeMuon.vz()) < 0.5)) {
        math::XYZTLorentzVector zBoson = tagMuon.p4() + probeMuon.p4();
        double mass = zBoson.M();
        if (fabs(mass - 91.2) < fabs(goodMass - 91.2)) {
          goodMass = mass;
          goodZBoson = zBoson;
        }
      }
    }
  }

  return goodMass;
}

// In beginJob, setup the tree structure
void Analysis::beginJob() {
  const char* filename = std::getenv("CRAB_OUTPUT_FILENAME");
  if (filename == nullptr) {
    filename = "data_D.root";
  }
  outputFile = new TFile(filename, "RECREATE");

  for (int ptCut : ptThresholds) {
    for (double isolationCut : isolationThresholds) {
      std::string branchName = "muon_" + std::to_string(ptCut) + "_" + std::to_string(static_cast<int>(isolationCut * 100));

      globalMuonsTrees[branchName] = new TTree(("GlobalMuons_" + branchName).c_str(), ("Global Muons with pt > " + std::to_string(ptCut) + " and isolation < " + std::to_string(isolationCut)).c_str());
      globalMuonsTrees[branchName]->Branch("zBosonMass", &zBosonMass, "zBosonMass/D");
      globalMuonsTrees[branchName]->Branch("eventNumber", &eventNumber, "eventNumber/I");
      globalMuonsTrees[branchName]->Branch("runNumber", &runNumber, "runNumber/I");
      globalMuonsTrees[branchName]->Branch("lumiSection", &lumiSection, "lumiSection/I");
      globalMuonsTrees[branchName]->Branch("efficiency", &efficiency, "efficiency/D");

      rpcMuonsTrees[branchName] = new TTree(("RPCMuons_" + branchName).c_str(), ("RPC Muons with pt > " + std::to_string(ptCut) + " and isolation < " + std::to_string(isolationCut)).c_str());
      rpcMuonsTrees[branchName]->Branch("zBosonMass", &zBosonMass, "zBosonMass/D");
      rpcMuonsTrees[branchName]->Branch("eventNumber", &eventNumber, "eventNumber/I");
      rpcMuonsTrees[branchName]->Branch("runNumber", &runNumber, "runNumber/I");
      rpcMuonsTrees[branchName]->Branch("lumiSection", &lumiSection, "lumiSection/I");
      rpcMuonsTrees[branchName]->Branch("efficiency", &efficiency, "efficiency/D");
    }
  }
}

void Analysis::endJob() {
  outputFile->Write();
  for (auto& tree : globalMuonsTrees) {
    delete tree.second;
  }
  for (auto& tree : rpcMuonsTrees) {
    delete tree.second;
  }
  outputFile->Close();
  delete outputFile;
}


void Analysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(Analysis);



