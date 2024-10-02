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
  mutable double zBosonMassGlobal;
  mutable double zBosonMassTracker;
  mutable double zBosonMassStandAlone;
  mutable double zBosonMassCalo;
  mutable double zBosonMassPF;
  mutable double zBosonMassRPC;
  mutable double zBosonMassGEM;
  mutable double zBosonMassME0;

  mutable int eventNumber;
  mutable int runNumber;
  mutable int lumiSection;
  mutable double mcWeight;
  mutable double pdfWeight;

  mutable int recoMuonSize;
  mutable int genMuonSize;
  mutable double scaleFactor;

  mutable int globalMuonSize;
  mutable int trackerMuonSize;
  mutable int standAloneMuonSize;
  mutable int caloMuonSize;
  mutable int pfMuonSize;
  mutable int rpcMuonSize;
  mutable int gemMuonSize;
  mutable int me0MuonSize;

  mutable double efficiencyReco;
  mutable double efficiencyGlobal;
  mutable double efficiencyTracker;
  mutable double efficiencyStandAlone;
  mutable double efficiencyCalo;
  mutable double efficiencyPF;
  mutable double efficiencyRPC;
  mutable double efficiencyGEM;
  mutable double efficiencyME0;

  mutable bool isZBosonRecoMatched;
  mutable bool isZBosonGlobalMatched;
  mutable bool isZBosonTrackerMatched;
  mutable bool isZBosonStandAloneMatched;
  mutable bool isZBosonCaloMatched;
  mutable bool isZBosonPFMatched;
  mutable bool isZBosonRPCMatched;
  mutable bool isZBosonGEMMatched;
  mutable bool isZBosonME0Matched;

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

  recoMuonSize = recoMuons.size();
  globalMuonSize = globalMuons.size();
  trackerMuonSize = trackerMuons.size();
  standAloneMuonSize = standAloneMuons.size();
  caloMuonSize = caloMuons.size();
  pfMuonSize = pfMuons.size();
  rpcMuonSize = rpcMuons.size();
  gemMuonSize = gemMuons.size();
  me0MuonSize = me0Muons.size();

  // Store gen-level muons
  std::vector<reco::GenParticle> genMuonCandidates;
  for (const auto& genMuon : *genMuons) {
    if (fabs(genMuon.pdgId()) == 13 && genMuon.status() == 1 && genMuon.pt() > 24 && fabs(genMuon.eta()) < 2.4) {
      genMuonCandidates.push_back(genMuon);
    }
  }

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

  int matchedGlobal = matchMuons(globalMuons, genMuonCandidates);
  int matchedTracker = matchMuons(trackerMuons, genMuonCandidates);
  int matchedStandAlone = matchMuons(standAloneMuons, genMuonCandidates);
  int matchedCalo = matchMuons(caloMuons, genMuonCandidates);
  int matchedPF = matchMuons(pfMuons, genMuonCandidates);
  int matchedRPC = matchMuons(rpcMuons, genMuonCandidates);
  int matchedGEM = matchMuons(gemMuons, genMuonCandidates);
  int matchedME0 = matchMuons(me0Muons, genMuonCandidates);

  
  efficiencyReco = genMuonSize > 0 ? static_cast<double>(matchedGlobal) / genMuonSize : 0.0;
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

  zBosonMassReco = recoMuons.size() >= 2 ? reconstructZBoson(recoMuons) : -1  ;
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
 
  auto matchZBoson = [&](double zMassReco, double zMassGen) -> bool {
    if (recoMuons.empty() || genMuonCandidates.empty()) return false;
    return reco::deltaR(recoMuons[0].eta(), recoMuons[0].phi(), genMuonCandidates[0].eta(), genMuonCandidates[0].phi()) < deltaRCut;
  };

  isZBosonRecoMatched = matchZBoson(zBosonMassReco, zBosonMassGen);
  isZBosonGlobalMatched = matchZBoson(zBosonMassGlobal, zBosonMassGen);
  isZBosonTrackerMatched = matchZBoson(zBosonMassTracker, zBosonMassGen);
  isZBosonStandAloneMatched = matchZBoson(zBosonMassStandAlone, zBosonMassGen);
  isZBosonCaloMatched = matchZBoson(zBosonMassCalo, zBosonMassGen);
  isZBosonPFMatched = matchZBoson(zBosonMassPF, zBosonMassGen);
  isZBosonRPCMatched = matchZBoson(zBosonMassRPC, zBosonMassGen);
  isZBosonGEMMatched = matchZBoson(zBosonMassGEM, zBosonMassGen);
  isZBosonME0Matched = matchZBoson(zBosonMassME0, zBosonMassGen);


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
  treeReco = new TTree("treeReco", "Reconstructed objects");
  treeReco->Branch("zBosonMassReco", &zBosonMassReco, "zBosonMassReco/D");
  treeReco->Branch("zBosonMassGlobal", &zBosonMassGlobal, "zBosonMassGlobal/D");
  treeReco->Branch("zBosonMassTracker", &zBosonMassTracker, "zBosonMassTracker/D");
  treeReco->Branch("zBosonMassStandAlone", &zBosonMassStandAlone, "zBosonMassStandAlone/D");
  treeReco->Branch("zBosonMassCalo", &zBosonMassCalo, "zBosonMassCalo/D");
  treeReco->Branch("zBosonMassPF", &zBosonMassPF, "zBosonMassPF/D");
  treeReco->Branch("zBosonMassRPC", &zBosonMassRPC, "zBosonMassRPC/D");
  treeReco->Branch("zBosonMassGEM", &zBosonMassGEM, "zBosonMassGEM/D");
  treeReco->Branch("zBosonMassME0", &zBosonMassME0, "zBosonMassME0/D");

  treeReco->Branch("isZBosonRecoMatched", &isZBosonRecoMatched, "isZBosonRecoMatched/O");
  treeReco->Branch("isZBosonGlobalMatched", &isZBosonGlobalMatched, "isZBosonGlobalMatched/O");
  treeReco->Branch("isZBosonTrackerMatched", &isZBosonTrackerMatched, "isZBosonTrackerMatched/O");
  treeReco->Branch("isZBosonStandAloneMatched", &isZBosonStandAloneMatched, "isZBosonStandAloneMatched/O");
  treeReco->Branch("isZBosonCaloMatched", &isZBosonCaloMatched, "isZBosonCaloMatched/O");
  treeReco->Branch("isZBosonPFMatched", &isZBosonPFMatched, "isZBosonPFMatched/O");
  treeReco->Branch("isZBosonRPCMatched", &isZBosonRPCMatched, "isZBosonRPCMatched/O");
  treeReco->Branch("isZBosonGEMMatched", &isZBosonGEMMatched, "isZBosonGEMMatched/O");
  treeReco->Branch("isZBosonME0Matched", &isZBosonME0Matched, "isZBosonME0Matched/O");

  treeReco->Branch("efficiencyReco", &efficiencyReco, "efficiencyReco/D");
  treeReco->Branch("efficiencyGlobal", &efficiencyGlobal, "efficiencyGlobal/D");
  treeReco->Branch("efficiencyTracker", &efficiencyTracker, "efficiencyTracker/D");
  treeReco->Branch("efficiencyStandAlone", &efficiencyStandAlone, "efficiencyStandAlone/D");
  treeReco->Branch("efficiencyCalo", &efficiencyCalo, "efficiencyCalo/D");
  treeReco->Branch("efficiencyPF", &efficiencyPF, "efficiencyPF/D");
  treeReco->Branch("efficiencyRPC", &efficiencyRPC, "efficiencyRPC/D");
  treeReco->Branch("efficiencyGEM", &efficiencyGEM, "efficiencyGEM/D");
  treeReco->Branch("efficiencyME0", &efficiencyME0, "efficiencyME0/D");

  treeGen = new TTree("treeGen", "Generated objects");
  treeGen->Branch("zBosonMassGen", &zBosonMassGen, "zBosonMassGen/D");
  treeGen->Branch("genMuonSize", &genMuonSize, "genMuonSize/I");
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