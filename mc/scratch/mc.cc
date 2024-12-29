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
#include "FWCore/Common/interface/TriggerNames.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdlib>

class mc : public edm::global::EDAnalyzer<> {
public:
  explicit mc(const edm::ParameterSet&);
  ~mc() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::StreamID, const edm::Event&, const edm::EventSetup&) const override;
  void endJob() override;

  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  std::string hltPath_;

  TFile* outputFile;
  TTree* treeGlobal;
  TTree* treeRPC;

  mutable double zBosonMassGlobal;
  mutable double zBosonMassRPC;

  mutable int eventNumber;
  mutable int runNumber;
  mutable int lumiSection;

  mutable int tagMuonSizeGlobal;
  mutable int probeMuonSizeGlobal;
  mutable int tagMuonSizeRPC;
  mutable int probeMuonSizeRPC;
};

mc::mc(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
      hltPath_(iConfig.getParameter<std::string>("hltPath")) {
}

mc::~mc() {
}

void mc::analyze(const edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerToken_, triggerResults);

  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerResults);
  bool passTrigger = false;

  // Check if the event passes the specified HLT trigger
  for (unsigned int i = 0; i < triggerResults->size(); ++i) {
    if (triggerNames.triggerName(i).find(hltPath_) != std::string::npos) {
      if (triggerResults->accept(i)) {
        passTrigger = true;
        break;
      }
    }
  }

  if (!passTrigger) return;  // Skip event if it doesn't pass the trigger

  std::vector<reco::Muon> globalTagMuons;
  std::vector<reco::Muon> globalProbeMuons;
  std::vector<reco::Muon> rpcTagMuons;
  std::vector<reco::Muon> rpcProbeMuons;

  eventNumber = iEvent.id().event();
  runNumber = iEvent.id().run();
  lumiSection = iEvent.luminosityBlock();

  // Muon selection
  for (const auto& muon : *muons) {
    float isolation = (muon.pfIsolationR04().sumChargedHadronPt + 
                         ((muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) > 0 ? 
                          (muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) : 0)) / muon.pt();
    if (muon.pt() > 26 && fabs(muon.eta()) < 2.4 && isolation < 0.2) {
      if (muon.isGlobalMuon()) {
        globalTagMuons.push_back(muon);
        if (muon.isRPCMuon()) {
          rpcTagMuons.push_back(muon);
        }
      }
    } else if (muon.pt() > 15 && fabs(muon.eta()) < 2.4) {
      if (muon.isGlobalMuon()) {
        globalProbeMuons.push_back(muon);
        if (muon.isRPCMuon()) {
          rpcProbeMuons.push_back(muon);
        }
      }
    }
  }

  // Z reconstruction function
  auto reconstructZBoson = [](const std::vector<reco::Muon>& tagMuons, const std::vector<reco::Muon>& probeMuons) -> double {
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
  };

  zBosonMassGlobal = reconstructZBoson(globalTagMuons, globalProbeMuons);
  zBosonMassRPC = reconstructZBoson(rpcTagMuons, rpcProbeMuons);

  tagMuonSizeGlobal = globalTagMuons.size();
  probeMuonSizeGlobal = globalProbeMuons.size();
  tagMuonSizeRPC = rpcTagMuons.size();
  probeMuonSizeRPC = rpcProbeMuons.size();

  if (zBosonMassGlobal > 0.0) {
    treeGlobal->Fill();
  }

  if (zBosonMassRPC > 0.0) {
    treeRPC->Fill();
  }
}

void mc::beginJob() {
  
  outputFile = new TFile("mc_data.root", "RECREATE");
    
  treeGlobal = new TTree("GlobalMuons", "Global Muons (Tag and Probe)");
  treeRPC = new TTree("RPCMuons", "RPC Muons (Tag and Probe)");

  treeGlobal->Branch("zBosonMass", &zBosonMassGlobal, "zBosonMass/D");
  treeGlobal->Branch("eventNumber", &eventNumber, "eventNumber/I");
  treeGlobal->Branch("runNumber", &runNumber, "runNumber/I");
  treeGlobal->Branch("lumiSection", &lumiSection, "lumiSection/I");
  treeGlobal->Branch("tagMuonSize", &tagMuonSizeGlobal, "tagMuonSize/I");
  treeGlobal->Branch("probeMuonSize", &probeMuonSizeGlobal, "probeMuonSize/I");

  treeRPC->Branch("zBosonMass", &zBosonMassRPC, "zBosonMass/D");
  treeRPC->Branch("eventNumber", &eventNumber, "eventNumber/I");
  treeRPC->Branch("runNumber", &runNumber, "runNumber/I");
  treeRPC->Branch("lumiSection", &lumiSection, "lumiSection/I");
  treeRPC->Branch("tagMuonSize", &tagMuonSizeRPC, "tagMuonSize/I");
  treeRPC->Branch("probeMuonSize", &probeMuonSizeRPC, "probeMuonSize/I");
}

void mc::endJob() {
  outputFile->Write();
  delete treeGlobal;
  delete treeRPC;
  outputFile->Close();
  delete outputFile;
}

void mc::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(mc);
