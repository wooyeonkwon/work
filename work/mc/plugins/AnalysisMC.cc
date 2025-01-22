#include <memory>
#include <mutex>
#include <array>
#include <iostream>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h" 
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdlib>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TSystem.h"


class AnalysisMC : public edm::one::EDAnalyzer<> {
public:
  explicit AnalysisMC(const edm::ParameterSet&);
  ~AnalysisMC() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> muonMCMatchToken_;

  TTree* tree;

  // Event information
  unsigned long long eventNumber;
  unsigned int runNumber;
  unsigned int lumiSection;
  float genWeight;
  int muonMCMatchCount;
  // Z boson information - arrays for different categories
  static const int nCategories = 11; // +1 for Gen Z
  std::array<double, nCategories> zMass;
  std::array<double, nCategories> zDvz;

  // Reconstructed Muon information
  std::vector<int> muon_charge;
  std::vector<double> muon_pt;
  std::vector<double> muon_eta;
  std::vector<double> muon_phi;
  std::vector<double> muon_iso;
  std::vector<double> muon_vz;

  // Reconstructed Muon flags
  std::vector<bool> muon_isGlobal;
  std::vector<bool> muon_isTracker;
  std::vector<bool> muon_isRPC;
  std::vector<bool> muon_isGEM;
  std::vector<bool> muon_isTight;
  
  // Z boson flags for reconstructed muons
  std::vector<bool> muon_isRecoZ;
  std::vector<bool> muon_isGlobalZ;
  std::vector<bool> muon_isTightGlobalZ;
  std::vector<bool> muon_isTrackerZ;
  std::vector<bool> muon_isRPCZ;
  std::vector<bool> muon_isTunedRPCZ;
  std::vector<bool> muon_isTightRPCZ;
  std::vector<bool> muon_isGEMZ;
  std::vector<bool> muon_isTunedGEMZ;
  std::vector<bool> muon_isTightGEMZ;

  // Generator level muon information
  std::vector<int> genMuon_charge;
  std::vector<double> genMuon_pt;
  std::vector<double> genMuon_eta;
  std::vector<double> genMuon_phi;
  std::vector<double> genMuon_vz;
  std::vector<bool> genMuon_isGenZ;

  // muon MC Matching information
  std::vector<bool> muonMCMatched;
  bool muonMCMatchPass;
};

AnalysisMC::AnalysisMC(const edm::ParameterSet& iConfig)
  : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    genParticleToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
    muonMCMatchToken_(consumes<edm::Association<reco::GenParticleCollection>>(iConfig.getParameter<edm::InputTag>("muonMCMatch"))) {
}

AnalysisMC::~AnalysisMC() {
}

void AnalysisMC::beginJob() {
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("Analysis", "skimmed Muon and reconstructed Z");

  // Event information branches
  tree->Branch("eventNumber", &eventNumber, "eventNumber/l");
  tree->Branch("runNumber", &runNumber, "runNumber/i");
  tree->Branch("lumiSection", &lumiSection, "lumiSection/i");
  tree->Branch("genWeight", &genWeight, "genWeight/F");

  // Z boson information branches
  tree->Branch("zMass_Reco", &zMass[0]);
  tree->Branch("zMass_Global", &zMass[1]);
  tree->Branch("zMass_TightGlobal", &zMass[2]);
  tree->Branch("zMass_Tracker", &zMass[3]);
  tree->Branch("zMass_RPC", &zMass[4]);
  tree->Branch("zMass_TunedRPC", &zMass[5]);
  tree->Branch("zMass_TightRPC", &zMass[6]);
  tree->Branch("zMass_GEM", &zMass[7]);
  tree->Branch("zMass_TunedGEM", &zMass[8]);
  tree->Branch("zMass_TightGEM", &zMass[9]);
  tree->Branch("zMass_Gen", &zMass[10]);

  tree->Branch("zDvz_Reco", &zDvz[0]);
  tree->Branch("zDvz_Global", &zDvz[1]);
  tree->Branch("zDvz_TightGlobal", &zDvz[2]);
  tree->Branch("zDvz_Tracker", &zDvz[3]);
  tree->Branch("zDvz_RPC", &zDvz[4]);
  tree->Branch("zDvz_TunedRPC", &zDvz[5]);
  tree->Branch("zDvz_TightRPC", &zDvz[6]);
  tree->Branch("zDvz_GEM", &zDvz[7]);
  tree->Branch("zDvz_TunedGEM", &zDvz[8]);
  tree->Branch("zDvz_TightGEM", &zDvz[9]);
  tree->Branch("zDvz_Gen", &zDvz[10]);

  // Reconstructed muon branches
  tree->Branch("muon_charge", &muon_charge);
  tree->Branch("muon_pt", &muon_pt);
  tree->Branch("muon_eta", &muon_eta);
  tree->Branch("muon_phi", &muon_phi);
  tree->Branch("muon_iso", &muon_iso);
  tree->Branch("muon_vz", &muon_vz);
  
  tree->Branch("muon_isGlobal", &muon_isGlobal);
  tree->Branch("muon_isTracker", &muon_isTracker);
  tree->Branch("muon_isRPC", &muon_isRPC);
  tree->Branch("muon_isGEM", &muon_isGEM);
  tree->Branch("muon_isTight", &muon_isTight);

  tree->Branch("muon_isRecoZ", &muon_isRecoZ);
  tree->Branch("muon_isGlobalZ", &muon_isGlobalZ);
  tree->Branch("muon_isTightGlobalZ", &muon_isTightGlobalZ);
  tree->Branch("muon_isTrackerZ", &muon_isTrackerZ);
  tree->Branch("muon_isRPCZ", &muon_isRPCZ);
  tree->Branch("muon_isTunedRPCZ", &muon_isTunedRPCZ);
  tree->Branch("muon_isTightRPCZ", &muon_isTightRPCZ);
  tree->Branch("muon_isGEMZ", &muon_isGEMZ);
  tree->Branch("muon_isTunedGEMZ", &muon_isTunedGEMZ);
  tree->Branch("muon_isTightGEMZ", &muon_isTightGEMZ);

  // Generator level muon branches
  tree->Branch("genMuon_charge", &genMuon_charge);
  tree->Branch("genMuon_pt", &genMuon_pt);
  tree->Branch("genMuon_eta", &genMuon_eta);
  tree->Branch("genMuon_phi", &genMuon_phi);
  tree->Branch("genMuon_vz", &genMuon_vz);
  tree->Branch("genMuon_isGenZ", &genMuon_isGenZ);

  tree->Branch("muonMCMatched", &muonMCMatched);
  tree->Branch("muonMCMatchPass", &muonMCMatchPass);

  tree->SetCacheSize(10000000);
  tree->SetMaxVirtualSize(1000000);
}

void AnalysisMC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<reco::MuonCollection> muons;
  edm::Handle<reco::GenParticleCollection> genParticles;
  edm::Handle<GenEventInfoProduct> genInfo;
  edm::Handle<edm::Association<reco::GenParticleCollection>> muonMCMatch;

  iEvent.getByToken(muonToken_, muons);
  iEvent.getByToken(genParticleToken_, genParticles);
  iEvent.getByToken(generatorToken_, genInfo);
  iEvent.getByToken(muonMCMatchToken_, muonMCMatch);

  // Clear all vectors
  muon_charge.clear();
  muon_pt.clear();
  muon_eta.clear();
  muon_phi.clear();
  muon_iso.clear();
  muon_vz.clear();
  muon_isGlobal.clear();
  muon_isTracker.clear();
  muon_isRPC.clear();
  muon_isGEM.clear();
  muon_isTight.clear();
  muon_isRecoZ.clear();
  muon_isGlobalZ.clear();
  muon_isTightGlobalZ.clear();
  muon_isTrackerZ.clear();
  muon_isRPCZ.clear();
  muon_isTunedRPCZ.clear();
  muon_isTightRPCZ.clear();
  muon_isGEMZ.clear();
  muon_isTunedGEMZ.clear();
  muon_isTightGEMZ.clear();
  genMuon_charge.clear();
  genMuon_pt.clear();
  genMuon_eta.clear();
  genMuon_phi.clear();
  genMuon_vz.clear();
  genMuon_isGenZ.clear();
  muonMCMatched.clear();
  muonMCMatchPass = false;
  // Initialize Z boson variables
  std::fill(zMass.begin(), zMass.end(), -1.0);
  std::fill(zDvz.begin(), zDvz.end(), -1.0);

  muonMCMatchCount = 0;
  // Reconstructed muon information
  for (const auto& muon : *muons) {
    //
    muon_isGlobal.push_back(muon.isGlobalMuon());
    muon_isTracker.push_back(muon.isTrackerMuon());
    muon_isRPC.push_back(muon.isRPCMuon());
    muon_isGEM.push_back(muon.isGEMMuon());

    bool isTight = muon.passed(reco::Muon::PFIsoTight) && muon.passed(reco::Muon::CutBasedIdTight);
    muon_isTight.push_back(isTight);

    // Get the reference to the current muon in the collection
    edm::Ref<reco::MuonCollection> muonRef(muons, &muon - &(*muons)[0]);
    // Check if there is a match in the muonMCMatch association
    if (muonMCMatch.isValid()) {
      // Get the index of the genParticle matched with this muon
      const edm::Ref<reco::GenParticleCollection> genMatch = (*muonMCMatch)[muonRef];
      muonMCMatched.push_back(genMatch.isAvailable());
      if (genMatch.isAvailable()) muonMCMatchCount++;
    }

    double iso = (muon.pfIsolationR04().sumChargedHadronPt +
                  std::max(0.0, muon.pfIsolationR04().sumNeutralHadronEt +
                  muon.pfIsolationR04().sumPhotonEt -
                  0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt();

    muon_charge.push_back(muon.charge());
    muon_pt.push_back(muon.pt());
    muon_eta.push_back(muon.eta());
    muon_phi.push_back(muon.phi());
    muon_iso.push_back(iso);
    muon_vz.push_back(muon.vz());

    muon_isRecoZ.push_back(false);
    muon_isGlobalZ.push_back(false);
    muon_isTightGlobalZ.push_back(false);
    muon_isTrackerZ.push_back(false);
    muon_isRPCZ.push_back(false);
    muon_isTunedRPCZ.push_back(false);
    muon_isTightRPCZ.push_back(false);
    muon_isGEMZ.push_back(false);
    muon_isTunedGEMZ.push_back(false);
    muon_isTightGEMZ.push_back(false);
  }
  if (muonMCMatchCount > 1) muonMCMatchPass = true;
  // Generator level muon information
  for (const auto& genParticle : *genParticles) {
    if (abs(genParticle.pdgId()) == 13) {
      genMuon_charge.push_back(genParticle.charge());
      genMuon_pt.push_back(genParticle.pt());
      genMuon_eta.push_back(genParticle.eta());
      genMuon_phi.push_back(genParticle.phi());
      genMuon_vz.push_back(genParticle.vz());
      genMuon_isGenZ.push_back(false);
    }
  }
  // Generator level Z reconstruction
  auto reconstructGenZBoson = [this]() {
    double bestMass = std::numeric_limits<double>::max();
    double minDvz = std::numeric_limits<double>::max();
    int bestMuon1 = -1;
    int bestMuon2 = -1;

    for (size_t i = 0; i < genMuon_pt.size(); ++i) {
      if (genMuon_pt[i] > 24 && fabs(genMuon_eta[i]) < 2.4) continue;
      for (size_t j = i + 1; j < genMuon_pt.size(); ++j) {
        if (genMuon_charge[i] + genMuon_charge[j] != 0) continue;
        if (genMuon_pt[j] > 24 && fabs(genMuon_eta[j]) < 2.4) continue;

        double dvz = fabs(genMuon_vz[i] - genMuon_vz[j]);
        math::PtEtaPhiMLorentzVector p4_1(genMuon_pt[i], genMuon_eta[i], genMuon_phi[i], 0.10566);
        math::PtEtaPhiMLorentzVector p4_2(genMuon_pt[j], genMuon_eta[j], genMuon_phi[j], 0.10566);
        double mass = (p4_1 + p4_2).mass();

        if (fabs(mass - 91.1876) < fabs(bestMass - 91.1876)) {
          bestMass = mass;
          minDvz = dvz;
          bestMuon1 = i;
          bestMuon2 = j;
        }
      }
    }

    if (bestMuon1 != -1) {
      genMuon_isGenZ[bestMuon1] = genMuon_isGenZ[bestMuon2] = true;
      zMass[10] = bestMass;
      zDvz[10] = minDvz;
    }
  };

  // Reconstructed Z reconstruction helper
  auto reconstructZBoson = [this](int categoryIndex, std::function<bool(size_t)> filter) {
    double bestMass = std::numeric_limits<double>::max();
    double minDvz = std::numeric_limits<double>::max();
    int bestMuon1 = -1;
    int bestMuon2 = -1;

    for (size_t i = 0; i < muon_pt.size(); ++i) {
      if (!filter(i)) continue;
      for (size_t j = i + 1; j < muon_pt.size(); ++j) {
        if (!filter(j)) continue;
        if (muon_charge[i] + muon_charge[j] != 0) continue;

        double dvz = fabs(muon_vz[i] - muon_vz[j]);
        math::PtEtaPhiMLorentzVector p4_1(muon_pt[i], muon_eta[i], muon_phi[i], 0.10566);
        math::PtEtaPhiMLorentzVector p4_2(muon_pt[j], muon_eta[j], muon_phi[j], 0.10566);
        double mass = (p4_1 + p4_2).mass();

        if (fabs(mass - 91.1876) < fabs(bestMass - 91.1876)) {
          bestMass = mass;
          minDvz = dvz;
          bestMuon1 = i;
          bestMuon2 = j;
        }
      }
    }




    if (bestMuon1 != -1) {
      switch(categoryIndex) {
        case 0: 
          muon_isRecoZ[bestMuon1] = muon_isRecoZ[bestMuon2] = true; 
          break;
        case 1: 
          muon_isGlobalZ[bestMuon1] = muon_isGlobalZ[bestMuon2] = true; 
          break;
        case 2: 
          muon_isTightGlobalZ[bestMuon1] = muon_isTightGlobalZ[bestMuon2] = true; 
          break;
        case 3: 
          muon_isTrackerZ[bestMuon1] = muon_isTrackerZ[bestMuon2] = true; 
          break;
        case 4: 
          muon_isRPCZ[bestMuon1] = muon_isRPCZ[bestMuon2] = true; 
          break;
        case 5: 
          muon_isTunedRPCZ[bestMuon1] = muon_isTunedRPCZ[bestMuon2] = true; 
          break;
        case 6: 
          muon_isTightRPCZ[bestMuon1] = muon_isTightRPCZ[bestMuon2] = true; 
          break;
        case 7: 
          muon_isGEMZ[bestMuon1] = muon_isGEMZ[bestMuon2] = true; 
          break;
        case 8: 
          muon_isTunedGEMZ[bestMuon1] = muon_isTunedGEMZ[bestMuon2] = true; 
          break;
        case 9: 
          muon_isTightGEMZ[bestMuon1] = muon_isTightGEMZ[bestMuon2] = true; 
          break;
      }
      
      zMass[categoryIndex] = bestMass;
      zDvz[categoryIndex] = minDvz;
    }
  };

  // Apply Z boson reconstruction for both gen and reco categories
  reconstructGenZBoson();
  reconstructZBoson(0, [this](size_t i) { return (muonMCMatched[i] && muon_pt[i] > 24 && fabs(muon_eta[i]) < 2.4); });
  reconstructZBoson(1, [this](size_t i) { return (muonMCMatched[i] && muon_pt[i] > 24 && fabs(muon_eta[i]) < 2.4) && muon_isGlobal[i]; });
  reconstructZBoson(2, [this](size_t i) { return (muonMCMatched[i] && muon_pt[i] > 24 && fabs(muon_eta[i]) < 2.4) && muon_isGlobal[i] && muon_isTight[i]; });
  reconstructZBoson(3, [this](size_t i) { return (muonMCMatched[i] && muon_pt[i] > 24 && fabs(muon_eta[i]) < 2.4) && muon_isTracker[i]; });
  reconstructZBoson(4, [this](size_t i) { return (muonMCMatched[i] && muon_pt[i] > 24 && fabs(muon_eta[i]) < 2.4) && muon_isRPC[i]; });
  reconstructZBoson(5, [this](size_t i) { return (muonMCMatched[i] && muon_pt[i] > 24 && fabs(muon_eta[i]) < 2.4) && muon_isRPC[i] && muon_isGlobal[i]; });
  reconstructZBoson(6, [this](size_t i) { return (muonMCMatched[i] && muon_pt[i] > 24 && fabs(muon_eta[i]) < 2.4) && muon_isRPC[i] && muon_isTight[i]; });
  reconstructZBoson(7, [this](size_t i) { return (muonMCMatched[i] && muon_pt[i] > 24 && fabs(muon_eta[i]) < 2.4) && muon_isGEM[i]; });
  reconstructZBoson(8, [this](size_t i) { return (muonMCMatched[i] && muon_pt[i] > 24 && fabs(muon_eta[i]) < 2.4) && muon_isGEM[i] && muon_isGlobal[i]; });
  reconstructZBoson(9, [this](size_t i) { return (muonMCMatched[i] && muon_pt[i] > 24 && fabs(muon_eta[i]) < 2.4) && muon_isGEM[i] && muon_isTight[i]; });

  // Event data
  eventNumber = iEvent.id().event();
  runNumber = iEvent.id().run();
  lumiSection = iEvent.luminosityBlock();
  genWeight = genInfo->weight();

  tree->Fill();
}

void AnalysisMC::endJob() {
  std::cout << "endJob" << std::endl;
}

void AnalysisMC::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muons", edm::InputTag("muons"));
  desc.add<edm::InputTag>("genParticles", edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("generator", edm::InputTag("generator"));
  desc.add<edm::InputTag>("muonMCMatch", edm::InputTag("muonMCMatch"));
  descriptions.add("analysisMC", desc);
}

DEFINE_FWK_MODULE(AnalysisMC);
