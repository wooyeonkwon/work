#ifndef ANALYSIS_CLASSES_H
#define ANALYSIS_CLASSES_H

// classes.h
#include <vector>
#include "DataFormats/Math/interface/LorentzVector.h"
#include <functional>

struct ZBosonInfo {
  double mass;
  double dvz;

  ClassDef(ZBosonInfo, 1); 
};

struct SelectedMuon {
  int charge;
  double pt;
  double eta;
  double phi;
  double iso;
  double vz;
  math::XYZTLorentzVector p4;
  bool isReco;
  bool isTightReco;
  bool isRPC;
  bool isTightRPC;
  bool isGEM;
  bool isTightGEM;
  bool isRecoZ;
  bool isTightRecoZ;
  bool isRPCZ;
  bool isTightRPCZ;
  bool isGEMZ;
  bool isTightGEMZ;

  // Z 보손 플래그를 설정하는 메서드
  void setZProperty(const std::function<bool(const SelectedMuon&)>& filter) {

    if (!filter) {
      std::cerr << "Error: filter is null!" << std::endl;
      return;  // 함수 종료
    }

    // filter 함수 호출
    if (filter(*this)) {
      if (isReco) isRecoZ = true;
      if (isTightReco) isTightRecoZ = true;
      if (isRPC) isRPCZ = true;
      if (isTightRPC) isTightRPCZ = true;
      if (isGEM) isGEMZ = true;
      if (isTightGEM) isTightGEMZ = true;
    }
  }

  ClassDef(SelectedMuon, 1);
};

// Structure for GenMuon information
struct GenMuonInfo {
  int charge;
  double pt;
  double eta;
  double phi;
  double vz;
  math::XYZTLorentzVector p4;
  bool isGenZ;
  
  ClassDef(GenMuonInfo, 1); 
};

namespace {
  struct dictionary {
    std::vector<ZBosonInfo> dummy1;
    std::vector<SelectedMuon> dummy2;
    std::vector<GenMuonInfo> dummy3;
  };
}

#endif