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
  void analyze(const edm::StreamID, const edm::Event&, const edm::EventSetup&) const override;
  void endJob() override;

  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  
  TFile* outputFile;
  TTree* treeGlobal;
  TTree* treeRPC;
  TTree* treenotRPC;

  mutable double zBosonMassGlobal;
  mutable double zBosonMassRPC;
  mutable double zBosonMassnotRPC;

  mutable int eventNumber;
  mutable int runNumber;
  mutable int lumiSection;

  mutable int muonSizeGlobal;
  mutable int muonSizeRPC;
  mutable int muonSizeNotRPC;
};

Analysis::Analysis(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))) {
}

Analysis::~Analysis() {
}

void Analysis::analyze(const edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  std::vector<reco::Muon> globalMuons;
  std::vector<reco::Muon> rpcMuons;
  std::vector<reco::Muon> notrpcMuons;

  eventNumber = iEvent.id().event();
  runNumber = iEvent.id().run();
  lumiSection = iEvent.luminosityBlock();

//muon selection
  for (const auto& muon : *muons) {
    if (muon.pt() > 20 && fabs(muon.eta()) < 2.4) {
      float isolation = (muon.pfIsolationR04().sumChargedHadronPt + 
                         ((muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) > 0 ? 
                          (muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) : 0)) / muon.pt();
      if (isolation < 0.15) {
        if (muon.isGlobalMuon()) {
          globalMuons.push_back(muon);
          if (muon.isRPCMuon()) {
            rpcMuons.push_back(muon);
          } else {
            notrpcMuons.push_back(muon);
        }
          }
      }
    }
  }

//Z reconstruction  
  auto reconstructZBoson = [](const std::vector<reco::Muon>& muons) -> double {
    if (muons.size() < 2) return 0.0;

    math::XYZTLorentzVector goodZBoson;
    double goodMass = 0.0;
    
    for (size_t i = 0; i < muons.size() - 1; ++i) {
      for (size_t j = i + 1; j < muons.size(); ++j) {
        if (muons[i].charge() + muons[j].charge() == 0) {
          math::XYZTLorentzVector zBoson = muons[i].p4() + muons[j].p4();
          double mass = zBoson.M();
          if (fabs(mass - 91.2) < fabs(goodMass - 91.2)) {
            goodMass = mass;
            goodZBoson = zBoson;
          }
        }
      }
    }
    
    return goodMass;
  };

  zBosonMassGlobal = reconstructZBoson(globalMuons);
  zBosonMassRPC = reconstructZBoson(rpcMuons);
  zBosonMassnotRPC = reconstructZBoson(notrpcMuons);

  muonSizeGlobal = globalMuons.size();
  muonSizeRPC = rpcMuons.size();
  muonSizeNotRPC = notrpcMuons.size();

  if (zBosonMassGlobal > 0.0) {
    treeGlobal->Fill();
  }

  if (zBosonMassRPC > 0.0) {
    treeRPC->Fill();
  }
  
  if (zBosonMassnotRPC > 0.0) {
    treenotRPC->Fill();
  }  
}

void Analysis::beginJob() {
  outputFile = new TFile("data3.root", "RECREATE");
  treeGlobal = new TTree("GlobalMuons", "Global Muons");
  treeRPC = new TTree("RPCMuons", "RPC Muons");
  treenotRPC = new TTree("notRPCMuons", "Global & notRPC Muons");

  treeGlobal->Branch("zBosonMass", &zBosonMassGlobal, "zBosonMass/D");
  treeGlobal->Branch("eventNumber", &eventNumber, "eventNumber/I");
  treeGlobal->Branch("runNumber", &runNumber, "runNumber/I");
  treeGlobal->Branch("lumiSection", &lumiSection, "lumiSection/I");
  treeGlobal->Branch("muonSizeGlobal", &muonSizeGlobal, "muonSizeGlobal/I");

  treeRPC->Branch("zBosonMass", &zBosonMassRPC, "zBosonMass/D");
  treeRPC->Branch("eventNumber", &eventNumber, "eventNumber/I");
  treeRPC->Branch("runNumber", &runNumber, "runNumber/I");
  treeRPC->Branch("lumiSection", &lumiSection, "lumiSection/I");
  treeRPC->Branch("muonSizeRPC", &muonSizeRPC, "muonSizeRPC/I");

  treenotRPC->Branch("zBosonMass", &zBosonMassnotRPC, "zBosonMass/D");
  treenotRPC->Branch("eventNumber", &eventNumber, "eventNumber/I");
  treenotRPC->Branch("runNumber", &runNumber, "runNumber/I");
  treenotRPC->Branch("lumiSection", &lumiSection, "lumiSection/I");
  treenotRPC->Branch("muonSizeNotRPC", &muonSizeNotRPC, "muonSizeNotRPC/I");
}

void Analysis::endJob() {
  outputFile->Write();
  delete treeGlobal;
  delete treeRPC;
  delete treenotRPC;
  outputFile->Close();
  delete outputFile;
}

void Analysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(Analysis);
