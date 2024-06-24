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
#include "TFile.h"
#include "TTree.h"

class Analysis : public edm::global::EDAnalyzer<> {
public:
  explicit Analysis(const edm::ParameterSet&);
  ~Analysis() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  
  TFile* outputFile;
  TTree* treeGlobal;
  TTree* treeRPC;
  
  std::vector<reco::Muon> globalMuons;
  std::vector<reco::Muon> rpcMuons;
  double zBosonMassGlobal;
  double zBosonMassRPC;
};

Analysis::Analysis(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))) {
    usesResource("TFileService");
    outputFile = new TFile("data.root", "RECREATE");
    treeGlobal = new TTree("GlobalMuons", "Global Muons");
    treeRPC = new TTree("RPCMuons", "RPC Muons");
    treeGlobal->Branch("zBosonMass", &zBosonMassGlobal, "zBosonMass/D");
    treeRPC->Branch("zBosonMass", &zBosonMassRPC, "zBosonMass/D");
}

Analysis::~Analysis() {
    outputFile->Close();
    delete outputFile;
}

void Analysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  globalMuons.clear();
  rpcMuons.clear();
  zBosonMassGlobal = -1;
  zBosonMassRPC = -1;

  for (const auto& muon : *muons) {
    if (muon.isGlobalMuon() && muon.pt() > 20 && fabs(muon.eta()) < 2.4) {
      globalMuons.push_back(muon);
      if (muon.isRPCMuon()) {
        rpcMuons.push_back(muon);
      }
    }
  }

  auto reconstructZBoson = [](const std::vector<reco::Muon>& muons) {
    math::XYZTLorentzVector maxZBoson;
    double maxMass = 0.0;
  
    for (size_t i = 0; i < muons.size() - 1; ++i) {
      for (size_t j = i + 1; j < muons.size(); ++j) {
        if (muons[i].charge() + muons[j].charge() == 0) {
          math::XYZTLorentzVector zBoson = muons[i].p4() + muons[j].p4();
          double mass = zBoson.M();
          if (mass > maxMass) {
            maxMass = mass;
            maxZBoson = zBoson;
          }
        }
      }
    }
    
    return maxMass;
  };

  zBosonMassGlobal = reconstructZBoson(globalMuons);
  zBosonMassRPC = reconstructZBoson(rpcMuons);

  if (zBosonMassGlobal > 0.0) {
    treeGlobal->Fill();
  }

  if (zBosonMassRPC > 0.0) {
    treeRPC->Fill();
  }
}

void Analysis::beginJob() {}

void Analysis::endJob() {
    outputFile->Write();
}

void Analysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(Analysis);
