#include <memory>
#include <thread>
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


  TTree* treeRPC;
  TTree* treeTightRPC;
  TTree* treeGEM;
  TTree* treeTightGEM;
  TTree* treeReco;

  const tree* Ttrees[];
  const char* treeNames[];

  static thread_local double zBosonMass;
  static thread_local int eventNumber;
  static thread_local int runNumber;
  static thread_local int lumiSection;


  // Reco muon variables
  static thread_local std::vector<double> muonPt;
  static thread_local std::vector<double> muonEta;
  static thread_local std::vector<double> muonPhi;
  static thread_local std::vector<double> muonIso;
  static thread_local int muonSize;

  static thread_local double vertexdz;

};



thread_local double Analysis::zBosonMass = 0.0;
thread_local int Analysis::eventNumber = 0;
thread_local int Analysis::runNumber = 0;
thread_local int Analysis::lumiSection = 0;
thread_local std::vector<double> Analysis::muonPt;
thread_local std::vector<double> Analysis::muonEta;
thread_local std::vector<double> Analysis::muonPhi;
thread_local std::vector<double> Analysis::muonIso;
thread_local int Analysis::muonSize = 0;
thread_local double Analysis::vertexdz = 0.0;


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
  Ttrees[] = {treeReco, treeRPC, treeTightRPC, treeGEM, treeTightGEM};
  treeNames[] = {"RecoMuons", "TightRPCMuons", "RPCMuons", "TightGEMMuons", "GEMMuons",};

  for (size_t i = 0; i < 5; ++i) {
    Ttrees[i] = fs->make<TTree>(treeNames[i], treeNames[i]);    
  }

  auto initializeTree = [this](TTree* tree) {
    tree->Branch("zBosonMass", &zBosonMass, "zBosonMass/D");
    tree->Branch("vertexdz", &vertexdz, "vertexdz/D");
    tree->Branch("eventNumber", &eventNumber, "eventNumber/I");
    tree->Branch("runNumber", &runNumber, "runNumber/I");
    tree->Branch("lumiSection", &lumiSection, "lumiSection/I");
    tree->Branch("muonSize", &muonSize, "muonSize/I");
    tree->Branch("muonPt", &muonPt);
    tree->Branch("muonEta", &muonEta);
    tree->Branch("muonPhi", &muonPhi);
    tree->Branch("muonIso", &muonIso);

  };

  for (size_t i = 0; i < 5; ++i) {
    initializeTree(Ttrees[i]);  
  }

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

  const std::vector<reco::Muon>* MuonList = [recoMuons, rpcTightMuons, rpcMuons, gemTightMuons, gemMuons];

  
  eventNumber = iEvent.id().event();
  runNumber = iEvent.id().run();
  lumiSection = iEvent.luminosityBlock();

  // Muon selection
  for (const auto& muon : *muons) {
    if (muon.pt() > 24 && fabs(muon.eta()) < 2.4) {
      double iso = (muon.pfIsolationR04().sumChargedHadronPt + 
                   std::max(0.0, muon.pfIsolationR04().sumNeutralHadronEt + 
                           muon.pfIsolationR04().sumPhotonEt - 
                           0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt();
      muonPt.push_back(muon.pt());
      muonEta.push_back(muon.eta());
      muonPhi.push_back(muon.phi());
      muonIso.push_back(iso);
      if (muon.isRPCMuon()) {
        rpcMuons.push_back(muon);
      }
      if (muon.isGEMMuon()) {
        gemMuons.push_back(muon);
      }

      
      if (iso < 0.15 && muon.passed(reco::Muon::CutBasedIdTight)) {  // Added isTightMuon check
        // Fill Reco muon information
        recoMuons.push_back(muon);
        if (muon.isRPCMuon()) {
          rpcTightMuons.push_back(muon);
        }
        if (muon.isGEMMuon()) {
          gemTightMuons.push_back(muon);
        }
      }
    }
  }




  


  // Z reconstruction  
  auto reconstructZBoson = [](const std::vector<reco::Muon>& muons) -> std::pair<double, double> { // Change return type to pair
    if (muons.size() < 2) return {0.0, 0.0}; // Return pair with zero values

    math::XYZTLorentzVector goodZBoson;
    double goodMass = -1.0;
    double dvz = 99.0;
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

  for (size_t i = 0; i < 5; ++i) {
    auto [zBosonMassResult, vertexdzResult] = reconstructZBoson(MuonList[i]);
    zBosonMass = zBosonMassResult;
    vertexdz = vertexdzResult;
    muonSize = MuonList[i].size();
    if (muonSize > 0) {
      Ttrees[i]->Fill();
    }  
    
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
