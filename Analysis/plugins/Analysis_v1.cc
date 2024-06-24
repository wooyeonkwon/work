#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"


#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"

#include <fstream>

class Analysis : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit Analysis(const edm::ParameterSet&);
  ~Analysis() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  
  TH1F* h_zBosonsGlobal;
  TH1F* h_zBosonsRPC;
  TFile* outputFile;
};

Analysis::Analysis(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))) {
    usesResource("TFileService");
    outputFile = new TFile("histograms.root", "RECREATE");
    h_zBosonsGlobal = new TH1F("zBosonsGlobal", "Z Bosons (Global);Mass (GeV);Events", 180, 0, 180);
    h_zBosonsRPC = new TH1F("zBosonsRPC", "Z Bosons (RPC);Mass (GeV);Events", 180, 0, 180);
}

Analysis::~Analysis() {
    delete h_zBosonsGlobal;
    delete h_zBosonsRPC;
    outputFile->Close();
    delete outputFile;
}

void Analysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  std::vector<reco::Muon> globalMuons;
  std::vector<reco::Muon> rpcMuons;
  std::vector<reco::Muon> notrpcMuons;

  for (const auto& muon : *muons) {
    if (muon.isGlobalMuon() && muon.pt() > 20 && muon.eta() < 2.4) { // muon_pt > 20 Gev, abs(muon_eta) < 2.4  
      globalMuons.push_back(muon);
      if (muon.isRPCMuon()) {
        rpcMuons.push_back(muon); //.. add muon size
      }
      else {
        notrpcMuons.push_back(muon);
      }
    }
  }
  auto reconstructZBoson = [](const std::vector<reco::Muon>& muons) {
    math::XYZTLorentzVector maxZBoson;
    double goodMass = 0.0;
  
    for (size_t i = 0; i < muons.size() - 1; ++i) {
      for (size_t j = i + 1; j < muons.size(); ++j) {
        if (muons[i].charge() + muons[j].charge() == 0) {  // opposite charge muons
          math::XYZTLorentzVector zBoson = muons[i].p4() + muons[j].p4();
          double mass = zBoson.M(); //min(abs(mass- 91.2)) 
          if (fabs(mass-91.2) < fabs(goodMass-91.2)) {
            goodMass = mass;
            goodZBoson = zBoson;
          }
        }
      }
    }
    
    std::vector<math::XYZTLorentzVector> zBosons;
    if (goodMass > 0.0) {
      zBosons.push_back(goodZBoson);  
    }
    
    return zBosons;
  };

  auto zBosonsGlobal = reconstructZBoson(globalMuons);
  auto zBosonsRPC = reconstructZBoson(rpcMuons);

  for (const auto& zBoson : zBosonsGlobal) {
    h_zBosonsGlobal->Fill(zBoson.M());
  }

  for (const auto& zBoson : zBosonsRPC) {
    h_zBosonsRPC->Fill(zBoson.M());
  }
}

void Analysis::beginJob() {}
void Analysis::endJob() {
    outputFile->cd();
    h_zBosonsGlobal->Write();
    h_zBosonsRPC->Write();

    // Gaussian fitting for the range 40-140
    TF1* fitFuncGlobal = new TF1("fitFuncGlobal", "gaus", 40, 140);
    h_zBosonsGlobal->Fit(fitFuncGlobal, "R");

    TF1* fitFuncRPC = new TF1("fitFuncRPC", "gaus", 40, 140);
    h_zBosonsRPC->Fit(fitFuncRPC, "R");

    // Write the fitted histograms
    h_zBosonsGlobal->Write();
    h_zBosonsRPC->Write();

    // Open a file to log the fit results
    std::ofstream logFile("fit_results.txt");
    if (logFile.is_open()) {
        logFile << "Global Fit Parameters:\n"
                << "Mean: " << fitFuncGlobal->GetParameter(1) << "\n"
                << "Width: " << fitFuncGlobal->GetParameter(2) << "\n";
        logFile << "RPC Fit Parameters:\n"
                << "Mean: " << fitFuncRPC->GetParameter(1) << "\n"
                << "Width: " << fitFuncRPC->GetParameter(2) << "\n";
        logFile.close();
    } else {
        edm::LogError("ZBosonAnalyzer") << "Unable to open log file for writing fit results.";
    }

    delete fitFuncGlobal;
    delete fitFuncRPC;
}

void Analysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(Analysis);
