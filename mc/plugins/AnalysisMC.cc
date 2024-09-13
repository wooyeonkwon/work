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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h" // For generator-level info
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"     // For PDF weights
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdlib>

class AnalysisMC : public edm::global::EDAnalyzer<> {
public:
  explicit AnalysisMC(const edm::ParameterSet&);
  ~AnalysisMC() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::StreamID, const edm::Event&, const edm::EventSetup&) const override;
  void endJob() override;

  // Tokens for the event data
  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genEventToken_;       // For MC weight
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genMuonToken_; // For gen-level muons
  edm::EDGetTokenT<LHEEventProduct> lheEventToken_;           // For PDF weights
  
  std::string hltPath_;

  // Output ROOT file and trees
  TFile* outputFile;
  TTree* treeReco;
  TTree* treeGen;

  // Variables for tree branches
  mutable double zBosonMassReco;
  mutable double zBosonMassGen;
  mutable int eventNumber;
  mutable int runNumber;
  mutable int lumiSection;
  mutable double mcWeight;
  mutable double pdfWeight;

  mutable int recoMuonSize;
  mutable int genMuonSize;
};

AnalysisMC::AnalysisMC(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
      vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      genEventToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),
      genMuonToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genMuons"))),
      lheEventToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEvent"))),
      hltPath_(iConfig.getParameter<std::string>("hltPath")) {
}

AnalysisMC::~AnalysisMC() {}

void AnalysisMC::analyze(const edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  // Get the muons and event information
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerToken_, triggerResults);

  edm::Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByToken(genEventToken_, genEventInfo);

  edm::Handle<edm::View<reco::GenParticle>> genMuons;
  iEvent.getByToken(genMuonToken_, genMuons);

  edm::Handle<LHEEventProduct> lheEventInfo;
  iEvent.getByToken(lheEventToken_, lheEventInfo);

  // Check if the event passes the trigger
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
    primaryVertex = vertices->front();  // Get first good vertex
  }

  // Store reco muons
  std::vector<reco::Muon> recoMuons;
  for (const auto& muon : *muons) {
    if (muon.pt() > 24 && fabs(muon.eta()) < 2.4) {
      recoMuons.push_back(muon);
    }
  }

  recoMuonSize = recoMuons.size();

  // Store gen-level muons
  std::vector<reco::GenParticle> genMuonCandidates;
  for (const auto& genMuon : *genMuons) {
    if (fabs(genMuon.pdgId()) == 13 && genMuon.status() == 1 && genMuon.pt() > 24 && fabs(genMuon.eta()) < 2.4) {
      genMuonCandidates.push_back(genMuon);
    }
  }

  genMuonSize = genMuonCandidates.size();

  // Reco-level Z boson mass
  auto reconstructZBoson = [](const std::vector<reco::Muon>& muons) -> double {
    if (muons.size() < 2) return 0.0;
    math::XYZTLorentzVector zBoson = muons[0].p4() + muons[1].p4();
    return zBoson.M();
  };

  zBosonMassReco = recoMuons.size() >= 2 ? reconstructZBoson(recoMuons) : 0;

  // Gen-level Z boson mass
  auto reconstructGenZBoson = [](const std::vector<reco::GenParticle>& genMuons) -> double {
    if (genMuons.size() < 2) return 0.0;
    math::XYZTLorentzVector zBoson = genMuons[0].p4() + genMuons[1].p4();
    return zBoson.M();
  };

  zBosonMassGen = genMuonCandidates.size() >= 2 ? reconstructGenZBoson(genMuonCandidates) : 0;

  // Get event info and MC weights
  eventNumber = iEvent.id().event();
  runNumber = iEvent.id().run();
  lumiSection = iEvent.luminosityBlock();
  mcWeight = genEventInfo->weight();

  // Get PDF weights
  if (lheEventInfo.isValid() && lheEventInfo->weights().size() > 0) {
    pdfWeight = lheEventInfo->weights()[0].wgt;
  }

  // Fill trees
  treeReco->Fill();
  treeGen->Fill();
}

void AnalysisMC::beginJob() {
  // Setup output file and trees
//  const char* filename = std::getenv("CRAB_OUTPUT_FILENAME"); // For CRAB jobs
//  if (filename == nullptr) {
//    filename = "mc_data.root";
//  }
  
  outputFile = new TFile("mc_data.root", "RECREATE");

  // Trees for Reco-level and Gen-level information
  treeReco = new TTree("RecoTree", "Reco-level Muons and Z Boson");
  treeGen = new TTree("GenTree", "Gen-level Muons and Z Boson");

  treeReco->Branch("zBosonMassReco", &zBosonMassReco, "zBosonMassReco/D");
  treeReco->Branch("eventNumber", &eventNumber, "eventNumber/I");
  treeReco->Branch("runNumber", &runNumber, "runNumber/I");
  treeReco->Branch("lumiSection", &lumiSection, "lumiSection/I");
  treeReco->Branch("mcWeight", &mcWeight, "mcWeight/D");
  treeReco->Branch("recoMuonSize", &recoMuonSize, "recoMuonSize/I");

  treeGen->Branch("zBosonMassGen", &zBosonMassGen, "zBosonMassGen/D");
  treeGen->Branch("genMuonSize", &genMuonSize, "genMuonSize/I");
  treeGen->Branch("pdfWeight", &pdfWeight, "pdfWeight/D");
}

void AnalysisMC::endJob() {
  outputFile->Write();
  outputFile->Close();
}

void AnalysisMC::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// Define the plugin
DEFINE_FWK_MODULE(AnalysisMC);
