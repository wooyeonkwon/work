#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Math/interface/PtEtaPhiM4D.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdlib>


class Analysis : public edm::global::EDAnalyzer<> {
public:
  explicit Analysis(const edm::ParameterSet&);
  ~Analysis() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::StreamID, const edm::Event&, const edm::EventSetup&) const override;
  void endJob() override;

  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  std::string hltPath_;

  TFile* outputFile;
  TTree* treeGlobal;
  TTree* treeTracker;
  TTree* treeStandAlone;
  TTree* treeCalo;
  TTree* treePF;
  TTree* treeRPC;
  TTree* treeGEM;
  TTree* treeME0;

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

  mutable int muonSizeGlobal;
  mutable int muonSizeTracker;
  mutable int muonSizeStandAlone;
  mutable int muonSizeCalo;
  mutable int muonSizePF;
  mutable int muonSizeRPC;
  mutable int muonSizeGEM;
  mutable int muonSizeME0;
  
};




Analysis::Analysis(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
      vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      hltPath_(iConfig.getParameter<std::string>("hltPath")) {
}

Analysis::~Analysis() {
}

void Analysis::analyze(const edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerToken_, triggerResults);

  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerResults);
  bool passTrigger = false;

  // Check if the event passes the HLT_IsoMu24_v trigger
  for (unsigned int i = 0; i < triggerResults->size(); ++i) {
    if (triggerNames.triggerName(i).find(hltPath_) != std::string::npos) {
      if (triggerResults->accept(i)) {
        passTrigger = true;
        break;
      }
    }
  }

  if (!passTrigger) return;  // Skip event if it doesn't pass the trigger

  reco::Vertex primaryVertex;
  if (!vertices->empty()) {
    primaryVertex = vertices->front();  // get first good vertex
  }
  
  std::vector<reco::Muon> globalMuons;
  std::vector<reco::Muon> trackerMuons;
  std::vector<reco::Muon> standAloneMuons;
  std::vector<reco::Muon> caloMuons;
  std::vector<reco::Muon> pfMuons;
  std::vector<reco::Muon> rpcMuons;
  std::vector<reco::Muon> gemMuons;
  std::vector<reco::Muon> me0Muons;
  
  eventNumber = iEvent.id().event();
  runNumber = iEvent.id().run();
  lumiSection = iEvent.luminosityBlock();

//muon selection
  for (const auto& muon : *muons) {
    if (muon.pt() > 24 && fabs(muon.eta()) < 2.4) {
      float isolation = (muon.pfIsolationR04().sumChargedHadronPt + 
                         ((muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) > 0 ? 
                          (muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) : 0)) / muon.pt();
      if (isolation < 0.15) {
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



//Z reconstruction  
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


  zBosonMassGlobal = reconstructZBoson(globalMuons);
  zBosonMassTracker = reconstructZBoson(trackerMuons);
  zBosonMassStandAlone = reconstructZBoson(standAloneMuons);
  zBosonMassCalo = reconstructZBoson(caloMuons);
  zBosonMassPF = reconstructZBoson(pfMuons);
  zBosonMassRPC = reconstructZBoson(rpcMuons);
  zBosonMassGEM = reconstructZBoson(gemMuons);
  zBosonMassME0 = reconstructZBoson(me0Muons);

  muonSizeGlobal = globalMuons.size();
  muonSizeTracker = trackerMuons.size();
  muonSizeStandAlone = standAloneMuons.size();
  muonSizeCalo = caloMuons.size();
  muonSizePF = pfMuons.size();
  muonSizeRPC = rpcMuons.size();
  muonSizeGEM = gemMuonss.size();
  muonSizeME0 = me0Muons.size();


  if (zBosonMassGlobal > 0.0) {
    treeGlobal->Fill();
  }
  
  if (zBosonMassTracker > 0.0) {
    treeTracker->Fill();
  }
  
  if (zBosonMassStandAlone > 0.0) {
    treeStandAlone->Fill();
  }
  
  if (zBosonMassCalo > 0.0) {
    treeCalo->Fill();
  }
  
  if (zBosonMassPF > 0.0) {
    treePF->Fill();
  }
  
  if (zBosonMassRPC > 0.0) {
    treeRPC->Fill();
  }
  
  if (zBosonMassGEM > 0.0) {
    treeGEM->Fill();
  }
  
  if (zBosonMassME0 > 0.0) {
    treeME0->Fill();
  }


void Analysis::beginJob() {
  const char* filename = std::getenv("CRAB_OUTPUT_FILENAME"); //use this with crab job
  if (filename == nullptr) {
    filename = "data_D.root";
  }
  outputFile = new TFile(filename, "RECREATE");


  treeGlobal = new TTree("GlobalMuons", "Global Muons");
  treeTracker = new TTree("TrackerMuons", "Tracker Muons");
  treeStandAlone = new TTree("StandAloneMuons", "StandAlone Muons");
  treeCalo = new TTree("CaloMuons", "Calo Muons");
  treePF = new TTree("PFMuons", "PF Muons");
  treeRPC = new TTree("RPCMuons", "RPC Muons");
  treeGEM = new TTree("GEMMuons", "GEM Muons");
  treeME0 = new TTree("ME0Muons", "ME0 Muons");


  treeGlobal->Branch("zBosonMass", &zBosonMassGlobal, "zBosonMass/D");
  treeGlobal->Branch("eventNumber", &eventNumber, "eventNumber/I");
  treeGlobal->Branch("runNumber", &runNumber, "runNumber/I");
  treeGlobal->Branch("lumiSection", &lumiSection, "lumiSection/I");
  treeGlobal->Branch("muonSize", &muonSizeGlobal, "muonSize/I");

  treeTracker->Branch("zBosonMass", &zBosonMassTracker, "zBosonMass/D");
  treeTracker->Branch("eventNumber", &eventNumber, "eventNumber/I");
  treeTracker->Branch("runNumber", &runNumber, "runNumber/I");
  treeTracker->Branch("lumiSection", &lumiSection, "lumiSection/I");
  treeTracker->Branch("muonSize", &muonSizeTracker, "muonSize/I");

  treeStandAlone->Branch("zBosonMass", &zBosonMassStandAlone, "zBosonMass/D");
  treeStandAlone->Branch("eventNumber", &eventNumber, "eventNumber/I");
  treeStandAlone->Branch("runNumber", &runNumber, "runNumber/I");
  treeStandAlone->Branch("lumiSection", &lumiSection, "lumiSection/I");
  treeStandAlone->Branch("muonSize", &muonSizeStandAlone, "muonSize/I");

  treeCalo->Branch("zBosonMass", &zBosonMassCalo, "zBosonMass/D");
  treeCalo->Branch("eventNumber", &eventNumber, "eventNumber/I");
  treeCalo->Branch("runNumber", &runNumber, "runNumber/I");
  treeCalo->Branch("lumiSection", &lumiSection, "lumiSection/I");
  treeCalo->Branch("muonSize", &muonSizeCalo, "muonSize/I");

  treePF->Branch("zBosonMass", &zBosonMassPF, "zBosonMass/D");
  treePF->Branch("eventNumber", &eventNumber, "eventNumber/I");
  treePF->Branch("runNumber", &runNumber, "runNumber/I");
  treePF->Branch("lumiSection", &lumiSection, "lumiSection/I");
  treePF->Branch("muonSize", &muonSizePF, "muonSize/I");

  treeRPC->Branch("zBosonMass", &zBosonMassRPC, "zBosonMass/D");
  treeRPC->Branch("eventNumber", &eventNumber, "eventNumber/I");
  treeRPC->Branch("runNumber", &runNumber, "runNumber/I");
  treeRPC->Branch("lumiSection", &lumiSection, "lumiSection/I");
  treeRPC->Branch("muonSize", &muonSizeRPC, "muonSize/I");

  treeGEM->Branch("zBosonMass", &zBosonMassGEM, "zBosonMass/D");
  treeGEM->Branch("eventNumber", &eventNumber, "eventNumber/I");
  treeGEM->Branch("runNumber", &runNumber, "runNumber/I");
  treeGEM->Branch("lumiSection", &lumiSection, "lumiSection/I");
  treeGEM->Branch("muonSize", &muonSizeGEM, "muonSize/I");

  treeME0->Branch("zBosonMass", &zBosonMassME0, "zBosonMass/D");
  treeME0->Branch("eventNumber", &eventNumber, "eventNumber/I");
  treeME0->Branch("runNumber", &runNumber, "runNumber/I");
  treeME0->Branch("lumiSection", &lumiSection, "lumiSection/I");
  treeME0->Branch("muonSize", &muonSizeME0, "muonSize/I");
}


void Analysis::endJob() {
  outputFile->Write();

  delete treeGlobal;
  delete treeTracker;
  delete treeStandAlone;
  delete treeCalo;
  delete treePF;
  delete treeRPC;
  delete treeGEM;
  delete treeME0;

  outputFile->Close();
  delete outputFile;
}


void Analysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(Analysis);


 