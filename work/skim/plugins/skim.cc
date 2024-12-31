#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"


class skim : public edm::stream::EDFilter<> {
public:
  explicit skim(const edm::ParameterSet&);
  ~skim() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  bool filter(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;
  std::string hltPath_;
  edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
};

skim::skim(const edm::ParameterSet& iConfig) 
  : hltPath_(iConfig.getParameter<std::string>("hltPath")),
    triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))){

}

skim::~skim() {
}

bool skim::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerToken_, triggerResults);
  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerResults);

  // Check if the event passes the HLT_IsoMu24_v trigger
  for (unsigned int i = 0; i < triggerResults->size(); ++i) {
    if (triggerNames.triggerName(i).find(hltPath_) != std::string::npos) {
      if (triggerResults->accept(i)) {
        return true;
      }
    }
  }

  return false;
}

void skim::beginStream(edm::StreamID) {
}

void skim::endStream() {
}
void skim::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("triggerResults", edm::InputTag("TriggerResults"));
  desc.add<std::string>("hltPath", "HLT_IsoMu24_v");
  descriptions.add("skim", desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(skim);
