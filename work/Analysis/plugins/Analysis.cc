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
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdlib>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TSystem.h"



class Analysis : public edm::one::EDAnalyzer<> {
public:
  explicit Analysis(const edm::ParameterSet&);
  ~Analysis() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  
  TTree* tree;

  // Event information
  unsigned long long eventNumber;
  unsigned int runNumber;
  unsigned int lumiSection;

  // Z boson information - arrays for different categories
  static const int nCategories = 6;
  std::array<double, nCategories> zMass;
  std::array<double, nCategories> zDvz;

  // Muon information
  std::vector<int> muon_charge;
  std::vector<double> muon_pt;
  std::vector<double> muon_eta;
  std::vector<double> muon_phi;
  std::vector<double> muon_iso;
  std::vector<double> muon_vz;
  std::vector<int> muon_size;
  
  // Muon flags
  std::vector<bool> muon_isReco;
  std::vector<bool> muon_isTightReco;
  std::vector<bool> muon_isRPC;
  std::vector<bool> muon_isTightRPC;
  std::vector<bool> muon_isGEM;
  std::vector<bool> muon_isTightGEM;
  
  // Z boson flags for each muon
  std::vector<bool> muon_isRecoZ;
  std::vector<bool> muon_isTightRecoZ;
  std::vector<bool> muon_isRPCZ;
  std::vector<bool> muon_isTightRPCZ;
  std::vector<bool> muon_isGEMZ;
  std::vector<bool> muon_isTightGEMZ;
};

Analysis::Analysis(const edm::ParameterSet& iConfig)
  : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))){
}

Analysis::~Analysis() {
}


void Analysis::beginJob() {
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("Analysis", "skimmed Muon and reconstructed Z");

  // Event information branches
  tree->Branch("eventNumber", &eventNumber, "eventNumber/l");
  tree->Branch("runNumber", &runNumber, "runNumber/i");
  tree->Branch("lumiSection", &lumiSection, "lumiSection/i");

  // Z boson information branches
  tree->Branch("zMass_Reco", &zMass[0]);
  tree->Branch("zMass_TightReco", &zMass[1]);
  tree->Branch("zMass_RPC", &zMass[2]);
  tree->Branch("zMass_TightRPC", &zMass[3]);
  tree->Branch("zMass_GEM", &zMass[4]);
  tree->Branch("zMass_TightGEM", &zMass[5]);

  tree->Branch("zDvz_Reco", &zDvz[0]);
  tree->Branch("zDvz_TightReco", &zDvz[1]);
  tree->Branch("zDvz_RPC", &zDvz[2]);
  tree->Branch("zDvz_TightRPC", &zDvz[3]);
  tree->Branch("zDvz_GEM", &zDvz[4]);
  tree->Branch("zDvz_TightGEM", &zDvz[5]);

  // Muon information branches
  tree->Branch("muon_charge", &muon_charge);
  tree->Branch("muon_pt", &muon_pt);
  tree->Branch("muon_eta", &muon_eta);
  tree->Branch("muon_phi", &muon_phi);
  tree->Branch("muon_iso", &muon_iso);
  tree->Branch("muon_vz", &muon_vz);
  tree->Branch("muon_size", &muon_size);
  
  // Muon flags branches
  tree->Branch("muon_isReco", &muon_isReco);
  tree->Branch("muon_isTightReco", &muon_isTightReco);
  tree->Branch("muon_isRPC", &muon_isRPC);
  tree->Branch("muon_isTightRPC", &muon_isTightRPC);
  tree->Branch("muon_isGEM", &muon_isGEM);
  tree->Branch("muon_isTightGEM", &muon_isTightGEM);
  
  // Z boson flag branches
  tree->Branch("muon_isRecoZ", &muon_isRecoZ);
  tree->Branch("muon_isTightRecoZ", &muon_isTightRecoZ);
  tree->Branch("muon_isRPCZ", &muon_isRPCZ);
  tree->Branch("muon_isTightRPCZ", &muon_isTightRPCZ);
  tree->Branch("muon_isGEMZ", &muon_isGEMZ);
  tree->Branch("muon_isTightGEMZ", &muon_isTightGEMZ);

  tree->SetCacheSize(10000000);
  tree->SetMaxVirtualSize(1000000);
}

void Analysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  // Clear all vectors
  muon_charge.clear();
  muon_pt.clear();
  muon_eta.clear();
  muon_phi.clear();
  muon_iso.clear();
  muon_vz.clear();
  muon_size.clear();
  muon_isReco.clear();
  muon_isTightReco.clear();
  muon_isRPC.clear();
  muon_isTightRPC.clear();
  muon_isGEM.clear();
  muon_isTightGEM.clear();
  muon_isRecoZ.clear();
  muon_isTightRecoZ.clear();
  muon_isRPCZ.clear();
  muon_isTightRPCZ.clear();
  muon_isGEMZ.clear();
  muon_isTightGEMZ.clear();

  // Initialize Z boson variables
  std::fill(zMass.begin(), zMass.end(), -1.0);
  std::fill(zDvz.begin(), zDvz.end(), -1.0);

  // Muon selection
  for (const auto& muon : *muons) {
    if (muon.pt() > 24 && fabs(muon.eta()) < 2.4) {
      double iso = (muon.pfIsolationR04().sumChargedHadronPt +
                   std::max(0.0, muon.pfIsolationR04().sumNeutralHadronEt +
                           muon.pfIsolationR04().sumPhotonEt -
                           0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt();

      // Store muon basic information
      muon_charge.push_back(muon.charge());
      muon_pt.push_back(muon.pt());
      muon_eta.push_back(muon.eta());
      muon_phi.push_back(muon.phi());
      muon_iso.push_back(iso);
      muon_vz.push_back(muon.vz());
      

      // Store muon flags
      bool isTight = muon.passed(reco::Muon::PFIsoTight) && muon.passed(reco::Muon::CutBasedIdTight);
      muon_isReco.push_back(true);
      muon_isTightReco.push_back(isTight);
      muon_isRPC.push_back(muon.isRPCMuon());
      muon_isTightRPC.push_back(muon.isRPCMuon() && isTight);
      muon_isGEM.push_back(muon.isGEMMuon());
      muon_isTightGEM.push_back(muon.isGEMMuon() && isTight);

      // Initialize Z flags
      muon_isRecoZ.push_back(false);
      muon_isTightRecoZ.push_back(false);
      muon_isRPCZ.push_back(false);
      muon_isTightRPCZ.push_back(false);
      muon_isGEMZ.push_back(false);
      muon_isTightGEMZ.push_back(false);
    }
  }

  // Z Mass reconstruction (helper function)
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
          muon_isTightRecoZ[bestMuon1] = muon_isTightRecoZ[bestMuon2] = true; 
          break;
        case 2: 
          muon_isRPCZ[bestMuon1] = muon_isRPCZ[bestMuon2] = true; 
          break;
        case 3: 
          muon_isTightRPCZ[bestMuon1] = muon_isTightRPCZ[bestMuon2] = true; 
          break;
        case 4: 
          muon_isGEMZ[bestMuon1] = muon_isGEMZ[bestMuon2] = true; 
          break;
        case 5: 
          muon_isTightGEMZ[bestMuon1] = muon_isTightGEMZ[bestMuon2] = true; 
          break;
      }
      
      zMass[categoryIndex] = bestMass;
      zDvz[categoryIndex] = minDvz;
    }
  };

  // Apply Z boson reconstruction for each category
  reconstructZBoson(0, [this](size_t i) { return muon_isReco[i]; });
  reconstructZBoson(1, [this](size_t i) { return muon_isTightReco[i]; });
  reconstructZBoson(2, [this](size_t i) { return muon_isRPC[i]; });
  reconstructZBoson(3, [this](size_t i) { return muon_isTightRPC[i]; });
  reconstructZBoson(4, [this](size_t i) { return muon_isGEM[i]; });
  reconstructZBoson(5, [this](size_t i) { return muon_isTightGEM[i]; });

  // Event data
  eventNumber = iEvent.id().event();
  runNumber = iEvent.id().run();
  lumiSection = iEvent.luminosityBlock();

  tree->Fill();
}

void Analysis::endJob() {
  std::cout << "endJob" << std::endl;
}

void Analysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muons", edm::InputTag("muons"));
  descriptions.add("analysis", desc);
}

DEFINE_FWK_MODULE(Analysis);
