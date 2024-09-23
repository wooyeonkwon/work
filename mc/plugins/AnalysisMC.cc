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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"    
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdlib>
#include "DataFormats/Math/interface/deltaR.h"
#include <mutex>

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
  edm::EDGetTokenT<GenEventInfoProduct> genEventToken_;       
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genMuonToken_; 
  edm::EDGetTokenT<LHEEventProduct> lheEventToken_;           
  
  std::string hltPath_;

  // Output ROOT file and trees
  TFile* outputFile = nullptr;
  TTree* treeReco = nullptr;
  TTree* treeGen = nullptr;

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
  mutable double scaleFactor;

  // Mutex for thread safety
  mutable std::mutex mtx_;
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

AnalysisMC::~AnalysisMC() {
  if (outputFile) {
    outputFile->Close();
    delete outputFile;
  }
}

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

  if (!muons.isValid() || !vertices.isValid() || !triggerResults.isValid() || !genEventInfo.isValid() || !genMuons.isValid() || !lheEventInfo.isValid()) {
    return;
  }

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
    primaryVertex = vertices->front();  
  }

  // Store reco muons
  std::vector<reco::Muon> recoMuons;
  for (const auto& muon : *muons) {
    if (muon.pt() > 24 && fabs(muon.eta()) < 2.4) {
      float isolation = (muon.pfIsolationR04().sumChargedHadronPt + 
                         ((muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) > 0 ? 
                          (muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) : 0)) / muon.pt();
      if (isolation < 0.15) {
        recoMuons.push_back(muon);
      }
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

  // MC Truth Matching - dR
  double deltaRCut = 0.15;
  std::vector<reco::Muon> matchedRecoMuons;
  std::vector<reco::GenParticle> matchedGenMuons;

  for (const auto& recoMuon : recoMuons) {
    for (const auto& genMuon : genMuonCandidates) {
      double deltaR = reco::deltaR(recoMuon.eta(), recoMuon.phi(), genMuon.eta(), genMuon.phi());
      if (deltaR < deltaRCut) {
        matchedRecoMuons.push_back(recoMuon);
        matchedGenMuons.push_back(genMuon);
        break;  
      }
    }
  }

  // Scale Factor 
  // ...?

  // Reco-level Z boson mass
  auto reconstructZBoson = [](const std::vector<reco::Muon>& muons) -> double {
    if (muons.size() < 2) return 0.0;

    math::XYZTLorentzVector goodZBoson;
    double goodMass = 0.0;
    
    for (size_t i = 0; i < muons.size() - 1; ++i) {
      for (size_t j = i + 1; j < muons.size(); ++j) {
        if ((muons[i].charge() + muons[j].charge() == 0) && (fabs(muons[i].vz()-muons[j].vz()) < 0.5)) {
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

  zBosonMassReco = recoMuons.size() >= 2 ? reconstructZBoson(recoMuons) : 0;

  // Gen-level Z boson mass
  auto reconstructGenZBoson = [](const std::vector<reco::GenParticle>& genMuons) -> double {
    if (genMuons.size() < 2) return 0.0;

    math::XYZTLorentzVector goodZBoson;
    double goodMass = 0.0;
    
    for (size_t i = 0; i < genMuons.size() - 1; ++i) {
      for (size_t j = i + 1; j < genMuons.size(); ++j) {
        if ((genMuons[i].charge() + genMuons[j].charge() == 0) && (fabs(genMuons[i].vz()-genMuons[j].vz()) < 0.5)) {
          math::XYZTLorentzVector zBoson = genMuons[i].p4() + genMuons[j].p4();
          double mass = zBoson.M();
          if (fabs(mass - 91.1876) < fabs(goodMass - 91.1876)) {
            goodMass = mass;
            goodZBoson = zBoson;
          }
        }
      }
    }
    
    return goodMass;
  };

  zBosonMassGen = genMuonCandidates.size() >= 2 ? reconstructGenZBoson(genMuonCandidates) : 0;

  // Store weights
  mcWeight = genEventInfo->weight();

  // Write to tree
  {
    std::lock_guard<std::mutex> lock(mtx_);
    treeReco->Fill();
    treeGen->Fill();
  }
}

void AnalysisMC::beginJob() {
  std::lock_guard<std::mutex> lock(mtx_);
  outputFile = new TFile("z_analysis_output.root", "RECREATE");

  treeReco = new TTree("treeReco", "Reconstructed objects");
  treeReco->Branch("zBosonMassReco", &zBosonMassReco, "zBosonMassReco/D");
  treeReco->Branch("recoMuonSize", &recoMuonSize, "recoMuonSize/I");
  treeReco->Branch("scaleFactor", &scaleFactor, "scaleFactor/D");

  treeGen = new TTree("treeGen", "Generated objects");
  treeGen->Branch("zBosonMassGen", &zBosonMassGen, "zBosonMassGen/D");
  treeGen->Branch("genMuonSize", &genMuonSize, "genMuonSize/I");
}

void AnalysisMC::endJob() {
  std::lock_guard<std::mutex> lock(mtx_);
  if (outputFile) {
    outputFile->Write();
    outputFile->Close();
    delete outputFile;
    outputFile = nullptr;
  }
}

void AnalysisMC::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(AnalysisMC);
