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
#include "interface/AnalysisClasses.h"
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
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;  // Added for generator weight
  
  TTree* tree;

  unsigned long long int eventNumber;
  unsigned int runNumber;
  unsigned int lumiSection;
  std::vector<ZBosonInfo> zBoson;
  std::vector<SelectedMuon> selectedMuons;
  
  // Gen Muon information
  std::vector<GenMuonInfo> genMuons;
};

AnalysisMC::AnalysisMC(const edm::ParameterSet& iConfig)
  : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    genParticleToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))) {
}

AnalysisMC::~AnalysisMC() {
}

void AnalysisMC::beginJob() {
  gSystem->Load("libAnalysisClasses.so");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("Analysis", "skimed Muon and reconstructed Z");
  tree->Branch("eventNumber", &eventNumber, "eventNumber/I");
  tree->Branch("runNumber", &runNumber, "runNumber/I");
  tree->Branch("lumiSection", &lumiSection, "lumiSection/I");
  tree->Branch("zBoson", &zBoson);
  tree->Branch("muons", &selectedMuons);
  tree->Branch("genMuons", &genMuons);

  tree->SetCacheSize(10000000);
  tree->SetMaxVirtualSize(1000000);
  std::cout << "beginJob" << std::endl;
}

void AnalysisMC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Clear vectors for new event
  edm::Handle<reco::MuonCollection> muons;
  edm::Handle<reco::GenParticleCollection> genParticles;
  edm::Handle<GenEventInfoProduct> genInfo;
  
  iEvent.getByToken(muonToken_, muons);
  iEvent.getByToken(genParticleToken_, genParticles);
  iEvent.getByToken(generatorToken_, genInfo);
  
  zBoson.clear();
  zBoson.resize(6);
  selectedMuons.clear();
  genMuons.clear();

  // Get generator weight
  float genWeight = genInfo->weight();

  // Gen Muon selection
  for (const auto& genParticle : *genParticles) {
    // Select only muons (PDG ID = ±13)
    if (abs(genParticle.pdgId()) == 13) {
      GenMuonInfo genMuon;
      genMuon.pt = genParticle.pt();
      genMuon.eta = genParticle.eta();
      genMuon.phi = genParticle.phi();
      genMuon.charge = genParticle.charge();
      genMuon.genWeight = genWeight;
      
      genMuons.push_back(genMuon);
    }
  }

  // Rest of the analysis code remains the same
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

  // Z Mass reconstruction logic remains the same
  auto reconstructZBoson = [](std::vector<SelectedMuon>& muons, std::function<bool(const SelectedMuon&)> filter) -> ZBosonInfo {
    double bestMass = std::numeric_limits<double>::max();
    double minDvz = std::numeric_limits<double>::max();
    for (size_t i = 0; i < muons.size(); ++i) {
      if (!filter(muons[i])) continue;
      for (size_t j = i + 1; j < muons.size(); ++j) {
        if (!filter(muons[j])) continue;
        if (muons[i].charge + muons[j].charge != 0) continue;

        double dvz = fabs(muons[i].vz - muons[j].vz);
        math::XYZTLorentzVector zBoson = muons[i].p4 + muons[j].p4;
        double mass = zBoson.M();

        if (fabs(mass - 91.1876) < fabs(bestMass - 91.1876)) {
          bestMass = mass;
          minDvz = dvz;
          muons[i].setZProperty(filter);
          muons[j].setZProperty(filter);
        }
      }
    }
    return {bestMass == std::numeric_limits<double>::max() ? -1.0 : bestMass,
            minDvz == std::numeric_limits<double>::max() ? -1.0 : minDvz};
  };

  zBoson[0] = reconstructZBoson(selectedMuons, [](const SelectedMuon& muon) { return muon.isReco; });
  zBoson[1] = reconstructZBoson(selectedMuons, [](const SelectedMuon& muon) { return muon.isTightReco; });
  zBoson[2] = reconstructZBoson(selectedMuons, [](const SelectedMuon& muon) { return muon.isRPC; });
  zBoson[3] = reconstructZBoson(selectedMuons, [](const SelectedMuon& muon) { return muon.isTightRPC; });
  zBoson[4] = reconstructZBoson(selectedMuons, [](const SelectedMuon& muon) { return muon.isGEM; });
  zBoson[5] = reconstructZBoson(selectedMuons, [](const SelectedMuon& muon) { return muon.isTightGEM; });

  eventNumber = iEvent.id().event();
  runNumber = iEvent.id().run();
  lumiSection = iEvent.luminosityBlock();

  tree->Fill();
}

void AnalysisMC::endJob() {
  std::cout << "endJob" << std::endl;
}

void AnalysisMC::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muons", edm::InputTag("muons"));
  desc.add<edm::InputTag>("genParticles", edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("generator", edm::InputTag("generator"));  // Added for generator weight
  descriptions.add("analysisMC", desc);
}

DEFINE_FWK_MODULE(AnalysisMC);