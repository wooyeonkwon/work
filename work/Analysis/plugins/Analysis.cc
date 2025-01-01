#include <memory>
#include <mutex>
#include <thread>
#include <array>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdlib>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <set>
#include <iostream>
#include "ZBosonInfo.h"
#include "SelectedMuon.h"


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

  mutable std::mutex mtx;
  


  static thread_local int eventNumber;
  static thread_local int runNumber;
  static thread_local int lumiSection;
  static thread_local std::vector<ZBosonInfo> zBoson;
  static thread_local std::vector<SelectedMuon> selectedMuons;
 


  std::set<int> eventNumbers;

};
thread_local int Analysis::eventNumber = -1;
thread_local int Analysis::runNumber = -1;
thread_local int Analysis::lumiSection = -1;
thread_local std::vector<ZBosonInfo> Analysis::zBoson;
thread_local std::vector<SelectedMuon> Analysis::selectedMuons;



Analysis::Analysis(const edm::ParameterSet& iConfig)
  : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))){
}

Analysis::~Analysis() {
}

void Analysis::beginJob() {
  edm::Service<TFileService> fs;
  //TFile &file = fs->file();
  //file.SetCompressionLevel(3);
  tree = fs->make<TTree>("Analysis", "skimed Muon and reconstructed Z");
  tree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  tree->Branch("runNumber", &runNumber, "runNumber/I");
  tree->Branch("lumiSection", &lumiSection, "lumiSection/I");
  tree->Branch("zBoson", &zBoson);
  tree->Branch("muons", &selectedMuons);

  tree->SetAutoFlush(10000);
  tree->SetCacheSize(10000000);
  tree->SetMaxVirtualSize(1000000);
  std::cout << "beginJob" << std::endl;

}

void Analysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Clear vectors for new event
  zBoson.clear();
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  selectedMuons.clear();

  // Muon selection
  for (const auto& muon : *muons) {
    if (muon.pt() > 24 && fabs(muon.eta()) < 2.4) {
      double iso = (muon.pfIsolationR04().sumChargedHadronPt +
                    std::max(0.0, muon.pfIsolationR04().sumNeutralHadronEt +
                            muon.pfIsolationR04().sumPhotonEt -
                            0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt();

      SelectedMuon selectedMuon = {
        .charge = muon.charge(),
        .pt = muon.pt(),
        .eta = muon.eta(),
        .phi = muon.phi(),
        .iso = iso,
        .vz = muon.vz(),
        .p4 = muon.p4(),
        .isReco = true,
        .isTightReco = muon.passed(reco::Muon::PFIsoTight) && muon.passed(reco::Muon::CutBasedIdTight),
        .isRPC = muon.isRPCMuon(),
        .isTightRPC = muon.isRPCMuon() && muon.passed(reco::Muon::PFIsoTight) && muon.passed(reco::Muon::CutBasedIdTight),
        .isGEM = muon.isGEMMuon(),
        .isTightGEM = muon.isGEMMuon() && muon.passed(reco::Muon::PFIsoTight) && muon.passed(reco::Muon::CutBasedIdTight),
        .isRecoZ = false,
        .isTightRecoZ = false,
        .isRPCZ = false,
        .isTightRPCZ = false,
        .isGEMZ = false,
        .isTightGEMZ = false,
      };

      selectedMuons.push_back(selectedMuon);
    }
  }

  // Z Mass reconstruction
  auto reconstructZBoson = [](std::vector<SelectedMuon>& muons, std::function<bool(const SelectedMuon&)> filter) -> ZBosonInfo {
    double bestMass = std::numeric_limits<double>::max();
    double minDvz = std::numeric_limits<double>::max();
    for (size_t i = 0; i < muons.size() - 1; ++i) {
      if (!filter(muons[i])) continue;
      for (size_t j = i + 1; j < muons.size(); ++j) {
        if (!filter(muons[j])) continue;
        if (muons[i].charge + muons[j].charge != 0) continue; // Charge matching

        double dvz = fabs(muons[i].vz - muons[j].vz);
        math::XYZTLorentzVector zBoson = muons[i].p4 + muons[j].p4;
        double mass = zBoson.M();

          if (fabs(mass - 91.1876) < fabs(bestMass - 91.1876)) {
            bestMass = mass;
            minDvz = dvz;
            muons[i].setZProperty(filter); // Mark as part of Z boson
            muons[j].setZProperty(filter);
          }
        }
      }
      return {bestMass == std::numeric_limits<double>::max() ? -1.0 : bestMass,
              minDvz == std::numeric_limits<double>::max() ? -1.0 : minDvz};
  };

  // Apply Z boson reconstruction for each category
  zBoson[0] = reconstructZBoson(selectedMuons, [](const SelectedMuon& muon) { return muon.isReco; });
  zBoson[1] = reconstructZBoson(selectedMuons, [](const SelectedMuon& muon) { return muon.isTightReco; });
  zBoson[2] = reconstructZBoson(selectedMuons, [](const SelectedMuon& muon) { return muon.isRPC; });
  zBoson[3] = reconstructZBoson(selectedMuons, [](const SelectedMuon& muon) { return muon.isTightRPC; });
  zBoson[4] = reconstructZBoson(selectedMuons, [](const SelectedMuon& muon) { return muon.isGEM; });
  zBoson[5] = reconstructZBoson(selectedMuons, [](const SelectedMuon& muon) { return muon.isTightGEM; });

  // Event data
  eventNumber = iEvent.id().event();
  runNumber = iEvent.id().run();
  lumiSection = iEvent.luminosityBlock();

  tree->Fill();
//
//  {  
//    bool i = false;
//    int attempt = 0;
//    while (!i) {
//      std::unique_lock<std::mutex> lock(mtx, std::defer_lock);
//      auto start = std::chrono::steady_clock::now();
//      while (true) {
//        if (lock.try_lock()) {
//          std::this_thread::sleep_for(std::chrono::milliseconds(10));
//          if (i == false) tree->Fill();
//          std::this_thread::sleep_for(std::chrono::milliseconds(10));
//          lock.unlock();
//          i = true;
//          break;
//        }
//        ++attempt;
//        // check timeout
//        if (std::chrono::steady_clock::now() - start > std::chrono::milliseconds(100)) {
//          break;
//        }
//      }
//      std::this_thread::sleep_for(std::chrono::milliseconds(10));
//      break;
//    }
//
//
//  }
//

    
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
