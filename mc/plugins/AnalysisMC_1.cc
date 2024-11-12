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
#include <iostream>




class AnalysisMC : public edm::global::EDAnalyzer<> {
public:
  explicit AnalysisMC(const edm::ParameterSet&);
  ~AnalysisMC() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  // Mutable variables for tree branches

  
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
  TTree* treeReco;
  TTree* treeGen;
  TTree* treeGlobal;
  TTree* treeTracker;
  TTree* treeStandAlone;
  TTree* treeCalo; 
  TTree* treePF;
  TTree* treeRPC;
  TTree* treeGEM;
  TTree* treeME0;

  mutable double zBosonMassReco;
  mutable double zBosonMassGlobal;
  mutable double zBosonMassTracker;
  mutable double zBosonMassStandAlone;
  mutable double zBosonMassCalo;
  mutable double zBosonMassPF;
  mutable double zBosonMassRPC;
  mutable double zBosonMassGEM;
  mutable double zBosonMassME0;

  mutable double efficiencyReco;
  mutable double efficiencyGlobal;
  mutable double efficiencyTracker;
  mutable double efficiencyStandAlone;
  mutable double efficiencyCalo;
  mutable double efficiencyPF;
  mutable double efficiencyRPC;
  mutable double efficiencyGEM;
  mutable double efficiencyME0;

  mutable double zBosonMassGen;
  mutable int genMuonSize;
  mutable double mcWeight;
  mutable double pdfWeight;

  mutable double scaleFactor;

  mutable double isZBosonMatched;
  mutable int muonSize;

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
  std::lock_guard<std::mutex> lock(mtx_);
}


void AnalysisMC::analyze(const edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  // Get the muons and event information
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  std::cout << "Muons retrieved: " << muons->size() << " muons found." << std::endl;

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);
  std::cout << "Vertices retrieved: " << vertices->size() << " vertices found." << std::endl;

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerToken_, triggerResults);
  std::cout << "Trigger results retrieved." << std::endl;

  edm::Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByToken(genEventToken_, genEventInfo);
  std::cout << "Gen event info retrieved." << std::endl;

  edm::Handle<edm::View<reco::GenParticle>> genMuons;
  iEvent.getByToken(genMuonToken_, genMuons);
  std::cout << "Gen muons retrieved: " << genMuons->size() << " gen muons found." << std::endl;

  edm::Handle<LHEEventProduct> lheEventInfo;
  iEvent.getByToken(lheEventToken_, lheEventInfo);
  std::cout << "LHE event info retrieved." << std::endl;

  if (!muons.isValid()) {
    std::cout << "Muon handle is invalid." << std::endl;
  }
  if (!vertices.isValid()) {
    std::cout << "Vertex handle is invalid." << std::endl;
  }
  if (!triggerResults.isValid()) {
    std::cout << "Trigger results handle is invalid." << std::endl;
  }
  if (!genEventInfo.isValid()) {
    std::cout << "Gen event info handle is invalid." << std::endl;
  }
  if (!genMuons.isValid()) {
    std::cout << "Gen muons handle is invalid." << std::endl;
  }
  if (!lheEventInfo.isValid()) {
    std::cout << "LHE event info handle is invalid." << std::endl;
  }

  if (!muons.isValid() || !vertices.isValid() || !triggerResults.isValid() || !genEventInfo.isValid() || !genMuons.isValid() || !lheEventInfo.isValid()) {
    std::cout << "One or more handles are invalid. Exiting analyze." << std::endl;
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

  std::vector<reco::Muon> globalMuons;
  std::vector<reco::Muon> trackerMuons;
  std::vector<reco::Muon> standAloneMuons;
  std::vector<reco::Muon> caloMuons;
  std::vector<reco::Muon> pfMuons;
  std::vector<reco::Muon> rpcMuons;
  std::vector<reco::Muon> gemMuons;
  std::vector<reco::Muon> me0Muons;

  // Store reco muons
  std::vector<reco::Muon> recoMuons;
  for (const auto& muon : *muons) {
    if (muon.pt() > 24 && fabs(muon.eta()) < 2.4 && muon::isTightMuon(muon, primaryVertex)) {
      float isolation = (muon.pfIsolationR04().sumChargedHadronPt + 
                         ((muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) > 0 ? 
                          (muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) : 0)) / muon.pt();
      if (isolation < 0.15) {
        recoMuons.push_back(muon);
        if (muon.isGlobalMuon()) {
          globalMuons.push_back(muon);
          if (muon.isTrackerMuon()) {
            trackerMuons.push_back(muon);
          }
          if (muon.isStandAloneMuon()) {
            standAloneMuons.push_back(muon);
          }
          if (muon.isCaloMuon()) {
            caloMuons.push_back(muon);
          }
          if (muon.isPFMuon()) {
            pfMuons.push_back(muon);
          }
          if (muon.isRPCMuon()) {
            rpcMuons.push_back(muon); 
          }
          if (muon.isGEMMuon()) {
            gemMuons.push_back(muon);
          }
          if (muon.isME0Muon()) {
            me0Muons.push_back(muon);
          }
        }
      }
    }
  }

  std::cout << "Reconstructed muons: " << recoMuons.size() << " muons found." << std::endl;

  muonSize = recoMuons.size();

  // Store gen-level muons
  std::vector<reco::GenParticle> genMuonCandidates;
  for (const auto& genMuon : *genMuons) {
    if (fabs(genMuon.pdgId()) == 13 && genMuon.status() == 1 && genMuon.pt() > 24 && fabs(genMuon.eta()) < 2.4) {
      genMuonCandidates.push_back(genMuon);
    }
  }

  std::cout << "Gen-level muons: " << genMuonCandidates.size() << " gen muons found." << std::endl;

  genMuonSize = genMuonCandidates.size();

  // MC Truth Matching - dR, charge, and PDG ID matching
  double deltaRCut = 0.15;
  auto matchMuons = [&](const std::vector<reco::Muon>& recoMuons, const std::vector<reco::GenParticle>& genMuons) -> int {
    int matchedCount = 0;
    std::vector<bool> genMuonMatched(genMuons.size(), false);

    for (const auto& recoMuon : recoMuons) {
      double bestDeltaR = deltaRCut;
      int bestGenMuonIndex = -1;

      for (size_t i = 0; i < genMuons.size(); ++i) {
        if (genMuonMatched[i]) continue;

        double deltaR = reco::deltaR(recoMuon.eta(), recoMuon.phi(), genMuons[i].eta(), genMuons[i].phi());
        if (deltaR < bestDeltaR && recoMuon.charge() == genMuons[i].charge() && recoMuon.pdgId() == genMuons[i].pdgId()) {
          bestDeltaR = deltaR;
          bestGenMuonIndex = i;
        }
      }

      if (bestGenMuonIndex >= 0) {
        matchedCount++;
        genMuonMatched[bestGenMuonIndex] = true;
      }
    }

    return matchedCount;
  };
  int matchedReco = matchMuons(recoMuons, genMuonCandidates);
  int matchedGlobal = matchMuons(globalMuons, genMuonCandidates);
  int matchedTracker = matchMuons(trackerMuons, genMuonCandidates);
  int matchedStandAlone = matchMuons(standAloneMuons, genMuonCandidates);
  int matchedCalo = matchMuons(caloMuons, genMuonCandidates);
  int matchedPF = matchMuons(pfMuons, genMuonCandidates);
  int matchedRPC = matchMuons(rpcMuons, genMuonCandidates);
  int matchedGEM = matchMuons(gemMuons, genMuonCandidates);
  int matchedME0 = matchMuons(me0Muons, genMuonCandidates);

  efficiencyReco = genMuonSize > 0 ? static_cast<double>(matchedReco) / genMuonSize : 0.0;
  efficiencyGlobal = genMuonSize > 0 ? static_cast<double>(matchedGlobal) / genMuonSize : 0.0;
  efficiencyTracker = genMuonSize > 0 ? static_cast<double>(matchedTracker) / genMuonSize : 0.0;
  efficiencyStandAlone = genMuonSize > 0 ? static_cast<double>(matchedStandAlone) / genMuonSize : 0.0;
  efficiencyCalo = genMuonSize > 0 ? static_cast<double>(matchedCalo) / genMuonSize : 0.0;
  efficiencyPF = genMuonSize > 0 ? static_cast<double>(matchedPF) / genMuonSize : 0.0;
  efficiencyRPC = genMuonSize > 0 ? static_cast<double>(matchedRPC) / genMuonSize : 0.0;
  efficiencyGEM = genMuonSize > 0 ? static_cast<double>(matchedGEM) / genMuonSize : 0.0;
  efficiencyME0 = genMuonSize > 0 ? static_cast<double>(matchedME0) / genMuonSize : 0.0;

  // Reco-level Z boson mass
  auto reconstructZBoson = [](const std::vector<reco::Muon>& muons) -> double {
    if (muons.size() < 2) return -1;

    math::XYZTLorentzVector goodZBoson;
    double goodMass = -1;
    
    for (size_t i = 0; i < muons.size() - 1; ++i) {
      for (size_t j = i + 1; j < muons.size(); ++j) {
        if ((muons[i].charge() + muons[j].charge() == 0) && (fabs(muons[i].vz()-muons[j].vz()) < 0.5)) {
          math::XYZTLorentzVector zBoson = muons[i].p4() + muons[j].p4();
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

  zBosonMassReco = recoMuons.size() >= 2 ? reconstructZBoson(recoMuons) : -1;
  zBosonMassGlobal = globalMuons.size() >= 2 ? reconstructZBoson(globalMuons) : -1;
  zBosonMassTracker = trackerMuons.size() >= 2 ? reconstructZBoson(trackerMuons) : -1;
  zBosonMassStandAlone = standAloneMuons.size() >= 2 ? reconstructZBoson(standAloneMuons) : -1;
  zBosonMassCalo = caloMuons.size() >= 2 ? reconstructZBoson(caloMuons) : -1;
  zBosonMassPF = pfMuons.size() >= 2 ? reconstructZBoson(pfMuons) : -1;
  zBosonMassRPC = rpcMuons.size() >= 2 ? reconstructZBoson(rpcMuons) : -1;
  zBosonMassGEM = gemMuons.size() >= 2 ? reconstructZBoson(gemMuons) : -1;
  zBosonMassME0 = me0Muons.size() >= 2 ? reconstructZBoson(me0Muons) : -1;

  // Gen-level Z boson mass
  auto reconstructGenZBoson = [](const std::vector<reco::GenParticle>& genMuons) -> double {
    if (genMuons.size() < 2) return -1;

    math::XYZTLorentzVector goodZBoson;
    double goodMass = -1;
    
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

  zBosonMassGen = genMuonCandidates.size() >= 2 ? reconstructGenZBoson(genMuonCandidates) : -1;

  // Composite object matching for Z boson
  auto matchZBoson = [&](const std::vector<reco::Muon>& recoMuons, const std::vector<reco::GenParticle>& genMuons) -> bool {
    if (recoMuons.empty() || genMuons.empty()) return false;
    return reco::deltaR(recoMuons[0].eta(), recoMuons[0].phi(), genMuons[0].eta(), genMuons[0].phi()) < deltaRCut;
  };

  {
    std::lock_guard<std::mutex> lock(mtx_);
    isZBosonMatched = matchZBoson(recoMuons, genMuonCandidates);
    treeReco->Fill();

    isZBosonMatched = matchZBoson(globalMuons, genMuonCandidates);
    treeGlobal->Fill();

    isZBosonMatched = matchZBoson(trackerMuons, genMuonCandidates);
    treeTracker->Fill();

    isZBosonMatched = matchZBoson(standAloneMuons, genMuonCandidates);
    treeStandAlone->Fill();

    isZBosonMatched = matchZBoson(caloMuons, genMuonCandidates);
    treeCalo->Fill();

    isZBosonMatched = matchZBoson(pfMuons, genMuonCandidates);
    treePF->Fill();

    isZBosonMatched = matchZBoson(rpcMuons, genMuonCandidates);
    treeRPC->Fill();

    isZBosonMatched = matchZBoson(gemMuons, genMuonCandidates);
    treeGEM->Fill();

    isZBosonMatched = matchZBoson(me0Muons, genMuonCandidates);
    treeME0->Fill();
  }

  // Store weights
  mcWeight = genEventInfo->weight();

  // Write to tree
  {
    std::lock_guard<std::mutex> lock(mtx_);
    treeReco->Fill();
    treeGen->Fill();
    treeGlobal->Fill();
    treeTracker->Fill();
    treeStandAlone->Fill();
    treeCalo->Fill();
    treePF->Fill();
    treeRPC->Fill();
    treeGEM->Fill();
    treeME0->Fill();
  }
}

void AnalysisMC::beginJob() {
  std::lock_guard<std::mutex> lock(mtx_);

  treeGen = new TTree("GenMuons", "Generated Muons");
  treeGen->Branch("zBosonMass", &zBosonMassGen, "zBosonMass/D");
  treeGen->Branch("muonSize", &genMuonSize, "muonSize/I");

  treeReco = new TTree("RecoMuons", "Reconstructed Muons");
  treeReco->Branch("zBosonMass", &zBosonMassReco, "zBosonMass/D");
  treeReco->Branch("isZBosonMatched", &isZBosonMatched, "isZBosonMatched/O");
  treeReco->Branch("efficiency", &efficiencyReco, "efficiency/D");
  treeReco->Branch("muonSize", &muonSize, "muonSize/I");

  treeGlobal = new TTree("GlobalMuons", "Global Muons");
  treeGlobal->Branch("zBosonMass", &zBosonMassGlobal, "zBosonMass/D");
  treeGlobal->Branch("isZBosonMatched", &isZBosonMatched, "isZBosonMatched/O");
  treeGlobal->Branch("efficiency", &efficiencyGlobal, "efficiency/D");
  treeGlobal->Branch("muonSize", &muonSize, "muonSize/I");

  treeTracker = new TTree("TrackerMuons", "Tracker Muons");
  treeTracker->Branch("zBosonMass", &zBosonMassTracker, "zBosonMass/D");
  treeTracker->Branch("isZBosonMatched", &isZBosonMatched, "isZBosonMatched/O");
  treeTracker->Branch("efficiency", &efficiencyTracker, "efficiency/D");
  treeTracker->Branch("muonSize", &muonSize, "muonSize/I");

  treeStandAlone = new TTree("StandAloneMuons", "StandAlone Muons");
  treeStandAlone->Branch("zBosonMass", &zBosonMassStandAlone, "zBosonMass/D");
  treeStandAlone->Branch("isZBosonMatched", &isZBosonMatched, "isZBosonMatched/O");
  treeStandAlone->Branch("efficiency", &efficiencyStandAlone, "efficiency/D");
  treeStandAlone->Branch("muonSize", &muonSize, "muonSize/I");

  treeCalo = new TTree("CaloMuons", "Calo Muons");
  treeCalo->Branch("zBosonMass", &zBosonMassCalo, "zBosonMass/D");
  treeCalo->Branch("isZBosonMatched", &isZBosonMatched, "isZBosonMatched/O");
  treeCalo->Branch("efficiency", &efficiencyCalo, "efficiency/D");
  treeCalo->Branch("muonSize", &muonSize, "muonSize/I");

  treePF = new TTree("PFMuons", "PF Muons");
  treePF->Branch("zBosonMass", &zBosonMassPF, "zBosonMass/D");
  treePF->Branch("isZBosonMatched", &isZBosonMatched, "isZBosonMatched/O");
  treePF->Branch("efficiency", &efficiencyPF, "efficiency/D");
  treePF->Branch("muonSize", &muonSize, "muonSize/I");

  treeRPC = new TTree("RPCMuons", "RPC Muons");
  treeRPC->Branch("zBosonMass", &zBosonMassRPC, "zBosonMass/D");
  treeRPC->Branch("isZBosonMatched", &isZBosonMatched, "isZBosonMatched/O");
  treeRPC->Branch("efficiency", &efficiencyRPC, "efficiency/D");
  treeRPC->Branch("muonSize", &muonSize, "muonSize/I");

  treeGEM = new TTree("GEMMuons", "GEM Muons");
  treeGEM->Branch("zBosonMass", &zBosonMassGEM, "zBosonMass/D");
  treeGEM->Branch("isZBosonMatched", &isZBosonMatched, "isZBosonMatched/O");
  treeGEM->Branch("efficiency", &efficiencyGEM, "efficiency/D");
  treeGEM->Branch("muonSize", &muonSize, "muonSize/I");

  treeME0 = new TTree("ME0Muons", "ME0 Muons");
  treeME0->Branch("zBosonMass", &zBosonMassME0, "zBosonMass/D");
  treeME0->Branch("isZBosonMatched", &isZBosonMatched, "isZBosonMatched/O");
  treeME0->Branch("efficiency", &efficiencyME0, "efficiency/D");
  treeME0->Branch("muonSize", &muonSize, "muonSize/I");
}

void AnalysisMC::endJob() {
  std::lock_guard<std::mutex> lock(mtx_);
}

void AnalysisMC::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(AnalysisMC);