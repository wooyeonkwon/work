#include <memory>
#include <mutex>
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
#include "work/Analysis/interface/AnalysisClasses.h"
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

  static thread_local int eventNumber;
  static thread_local int runNumber;
  static thread_local int lumiSection;
  static thread_local std::vector<ZBosonInfo> zBoson;
  static thread_local std::vector<SelectedMuon> selectedMuons;
};

thread_local int Analysis::eventNumber;
thread_local int Analysis::runNumber;
thread_local int Analysis::lumiSection;
thread_local std::vector<ZBosonInfo> Analysis::zBoson;
thread_local std::vector<SelectedMuon> Analysis::selectedMuons;

Analysis::Analysis(const edm::ParameterSet& iConfig)
  : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))){
}

Analysis::~Analysis() {
}

void Analysis::beginJob() {
  gSystem->Load("libAnalysisClasses.so");
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
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  zBoson.clear();
  zBoson.resize(6);
  selectedMuons.clear();

  // Muon selection
  for (const auto& muon : *muons) {
    if (muon.pt() > 24 && fabs(muon.eta()) < 2.4) {
      double iso = (muon.pfIsolationR04().sumChargedHadronPt +
                    std::max(0.0, muon.pfIsolationR04().sumNeutralHadronEt +
                            muon.pfIsolationR04().sumPhotonEt -
                            0.5 * muon.pfIsolationR04().sumPUPt)) / muon.pt();

      SelectedMuon selectedMuon{
          muon.charge(), muon.pt(), muon.eta(), muon.phi(), iso, muon.vz(),
          muon.p4(), true,
          muon.passed(reco::Muon::PFIsoTight) && muon.passed(reco::Muon::CutBasedIdTight),
          muon.isRPCMuon(),
          muon.isRPCMuon() && muon.passed(reco::Muon::PFIsoTight) && muon.passed(reco::Muon::CutBasedIdTight),
          muon.isGEMMuon(),
          muon.isGEMMuon() && muon.passed(reco::Muon::PFIsoTight) && muon.passed(reco::Muon::CutBasedIdTight),
          false, false, false, false, false, false
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
