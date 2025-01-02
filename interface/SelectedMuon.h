// SelectedMuon.h
#ifndef SELECTEDMUON_H
#define SELECTEDMUON_H

#include "Rtypes.h"
#include <vector>
#include "DataFormats/Math/interface/LorentzVector.h"
#include <functional>

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
    void setZProperty(std::function<bool(const SelectedMuon&)> filter) {
        if (filter(*this)) {
            if (isReco) isRecoZ = true;
            if (isTightReco) isTightRecoZ = true;
            if (isRPC) isRPCZ = true;
            if (isTightRPC) isTightRPCZ = true;
            if (isGEM) isGEMZ = true;
            if (isTightGEM) isTightGEMZ = true;
        }
    }

};

#endif