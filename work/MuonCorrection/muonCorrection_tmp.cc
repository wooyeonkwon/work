#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <memory>
#include <tuple>
#include <array>
#include "ROOT/RDataFrame.hxx"
#include <TSystem.h>

const char* path = "/data1/users/dndus0107/AnalysisResults/processed_data/";
const char* gen_name = "mini_DYto2Mu_MLL-50to120_22EEDR.root";
const char* reco_name = "mini_Analysis_Data_22EFG";



class MuonAnalyzer {
private:

/////////////////////////////
//                         //
//         Binning         //
//                         //
/////////////////////////////

    // Binning parameters
    std::vector<double> eta_bins;
    std::vector<double> phi_bins;
    std::vector<int> charge_bins;

    // Fine tuning factors
    std::map<std::tuple<int, int, int>, double> fine_tune_factors;
    std::map<std::tuple<int, int>, std::array<double, 2>> correction_factors;
    // Constants
    const double Z_MASS = 91.1876;
    const double CONVERGENCE_THRESHOLD = 0.001;


    int getEtaBin(double eta) {
        // eta의 절대값 사용
        double absEta = std::abs(eta);
        for (int i = 0; i < eta_bins.size() - 1; ++i) {
            if (absEta >= eta_bins[i] && absEta < eta_bins[i + 1]) {
                return i;
            }
        }
        return -1;
    }

    int getPhiBin(double phi) {
        // phi의 절대값 사용하고 2π로 정규화
        double absPhi = std::abs(phi);
        // phi가 -π에서 π 범위에 있다고 가정
        if (absPhi > M_PI) {
            absPhi = 2 * M_PI - absPhi;
        }
        
        for (int i = 0; i < phi_bins.size() - 1; ++i) {
            if (absPhi >= phi_bins[i] && absPhi < phi_bins[i + 1]) {
                return i;
            }
        }
        return -1;
    }

    int getChargeBin(int charge) {
        // charge_bins는 {-1, 1}로 초기화되어 있다고 가정
        for (int i = 0; i < charge_bins.size(); ++i) {
            if (charge == charge_bins[i]) {
                return i;  // charge 값에 맞는 bin 인덱스 반환
            }
        }
        return -1;  // charge 값이 범위에 없을 경우 -1 반환
    }

//////////////////////////////////////////////////////////////////////////////////////////


    // Calculate Z mass statistics for a given bin
    void calculateZMassStats(ROOT::RDF::RNode& frame) {
        std::cout << "Start calculating Z mass statistics..." << std::endl;

        try {
            // 새로운 데이터 프레임 정의
            auto df = frame
                // isTightRecoZ 조건을 만족하는 뮤온 선택
                .Define("selected_muons", [](const std::vector<SelectedMuon>& muons) {
                    std::vector<SelectedMuon> selected;
                    for (const auto& muon : muons) {
                        if (muon.isTightRecoZ) {
                            selected.push_back(muon);
                        }
                    }
                    return selected;
                }, {"muons"})

                // Z 질량 계산
                .Define("z_mass", [](const std::vector<SelectedMuon>& selected_muons) {
                    if (selected_muons.size() != 2) {
                        std::cerr << "Error: Exactly 2 muons with isTightRecoZ = true are required!" << std::endl;
                        return -1.0; // 오류 상황 처리
                    }
                    // 두 뮤온의 4-벡터 합산 후 질량 계산
                    auto p4 = selected_muons[0].p4 + selected_muons[1].p4;
                    return p4.M();
                }, {"selected_muons"})

                // zBoson[1].mass 업데이트
                .Define("updated_zBosons", [](const std::vector<ZBosonInfo>& zBosons, double z_mass) {
                    if (zBosons.size() > 1) {
                        // zBosons 크기가 2 이상일 때만 업데이트
                        std::vector<ZBosonInfo> updatedBosons = zBosons;
                        updatedBosons[1].mass = z_mass; // zBoson[1]의 질량 업데이트
                        return updatedBosons;
                    } else {
                        // 경계 조건에서 원래 zBosons 반환
                        std::cerr << "Warning: Not enough Z boson entries in the data! Returning original zBosons." << std::endl;
                        return std::vector<ZBosonInfo>(zBosons); // 항상 유효한 객체 반환
                    }
                }, {"zBoson", "z_mass"});
            // bin별 통계 수집
            std::map<std::tuple<int, int, int>, std::pair<double, int>> bin_stats;

            // 각 bin에 대한 통계 계산
            auto process_event = [this, &bin_stats](const std::vector<SelectedMuon>& muons, double z_mass) {
                if (z_mass <= 0 || muons.size() != 2) return false;

                for (const auto& muon : muons) {
                    int eta_bin = this->getEtaBin(muon.eta);
                    int phi_bin = this->getPhiBin(muon.phi);

                    if (eta_bin >= 0 && phi_bin >= 0) {
                        auto key = std::make_tuple(muon.charge, eta_bin, phi_bin);
                        auto& stats = bin_stats[key];
                        stats.first += z_mass;  // sum_mass
                        stats.second += 1;      // count
                    }
                }

                return true; // 반환값 추가
            };

            // 이벤트 처리
            *df.Define("dummy", process_event, {"muons", "z_mass"})
            .Count();

            // fine tuning factors 업데이트
            for (const auto& [key, stats] : bin_stats) {
                if (stats.second > 0) {
                    double avg_mass = stats.first / stats.second;
                    fine_tune_factors[key] = Z_MASS / avg_mass;
                    std::cout << "Bin " << std::get<0>(key) << " " << std::get<1>(key) << " " 
                             << std::get<2>(key) << ": avg_mass = " << avg_mass 
                             << ", factor = " << fine_tune_factors[key] << std::endl;
                }
            }    

        } catch (const std::exception& e) {
            std::cerr << "Error in calculateZMassStats: " << e.what() << std::endl;
            throw;
        }
    }

    void updateFineTuneFactors(TTree* tree, int eta_bin, int phi_bin) {
        std::vector<ZBosonInfo> zBoson;
        std::vector<SelectedMuon> muons;
        
        SetBranchAddress(tree, "muons", &muons);
        SetBranchAddress(tree, "zBoson", &zBoson);
        
        for (int charge : charge_bins) {
            double sum_Z_mass = 0;
            int nZ = 0;
            
            Long64_t nentries = tree->GetEntries();
            for (Long64_t i = 0; i < nentries; ++i) {
                tree->GetEntry(i);
                
                if (zBoson.empty() || zBoson[1].mass == -1) continue;

                for (const auto& muon : muons) {
                    if (muon.charge != charge) continue;

                    int muon_eta_bin = getEtaBin(muon.eta);
                    int muon_phi_bin = getPhiBin(muon.phi);

                    if (muon_eta_bin == eta_bin && muon_phi_bin == phi_bin) {
                        sum_Z_mass += zBoson[1].mass;
                        nZ++;
                    }
                }
            }
            
            if (sum_Z_mass > 0 && nZ > 0) {
                auto key = std::make_tuple(charge, eta_bin, phi_bin);
                fine_tune_factors[key] = 1 + (2 * (sum_Z_mass - nZ * Z_MASS) / sum_Z_mass);
            }
        }
    }

    ROOT::RDF::RNode ptCorrection(ROOT::RDataFrame& frame, const std::string& reco_name) {
        std::cout << "Start pt correction for all muons in the dataset" << std::endl;

        // 'muons' 컬럼 존재 여부 확인
        auto colNames = frame.GetColumnNames();
        bool hasMuons = false;
        for (const auto& name : colNames) {
            if (name == "muons") hasMuons = true;
        }

        if (!hasMuons) {
            std::cerr << "Error: 'muons' column not found in the DataFrame" << std::endl;
            return frame;
        }
        try {
            // 데이터프레임 복제
            ROOT::RDF::RNode corrected_frame = frame;
            // 'muons' 컬럼을 보정하여 새로운 컬럼 생성
            corrected_frame = corrected_frame.Define("corrected_muons", [this](const std::vector<SelectedMuon>& input_muons) {
                std::vector<SelectedMuon> corrected_muons;
                for (const auto& muon : input_muons) {
                    SelectedMuon corrected_muon = muon;  // 기존 값 복제

                    // η 및 φ에 해당하는 bin 찾기
                    int eta_bin = this->getEtaBin(muon.eta);
                    int phi_bin = this->getPhiBin(muon.phi);

                    if (eta_bin >= 0 && phi_bin >= 0) {
                        auto factor_it = this->correction_factors.find(std::make_tuple(eta_bin, phi_bin));
                        if (factor_it != this->correction_factors.end()) {
                            std::array<double, 2> factor = factor_it->second;
                            double A = factor[0];
                            double M = factor[1];

                            // pt 보정 계산
                            if (muon.pt > 0) {
                                corrected_muon.pt = 1.0 / ((M / muon.pt) + A * muon.charge);

                                // 보정된 pt를 사용해 새로운 4-벡터 생성
                                TLorentzVector new_p4;
                                new_p4.SetPtEtaPhiM(
                                    corrected_muon.pt,
                                    corrected_muon.p4.Eta(),
                                    corrected_muon.p4.Phi(),
                                    corrected_muon.p4.M()
                                );
                                corrected_muon.p4 = new_p4;
                            }
                        }
                    }
                    corrected_muons.push_back(corrected_muon);  // 보정된 뮤온 추가
                }
                return corrected_muons;  // 보정된 뮤온 리스트 반환
            }, {"muons"});  // 'muons' 열을 기반으로 계산
            // 결과를 새로운 ROOT 파일로 저장
            std::string output_filename = reco_name + "_corrected.root";
            //corrected_frame.Snapshot("Analysis", output_filename);
            std::cout << "Finished pt correction and saved to " << output_filename << std::endl;
            return corrected_frame;  // 보정된 데이터프레임 반환
        } catch (const std::exception& e) {
            std::cerr << "Error in ptCorrection: " << e.what() << std::endl;
            throw;
        }
    }



    double calculateInversePtAverage(ROOT::RDataFrame& frame, int charge, int eta_bin, int phi_bin, bool is_gen = false) {

        if (is_gen) {
            auto filterFunc = [this, charge, eta_bin, phi_bin] -> bool {
                return std::any_of(muons.begin(), muons.end(), 
                    [&](const GenMuonInfo& muon) {
                        return muon.isGenZ && 
                               muon.charge == charge && 
                               this->getEtaBin(muon.eta) == eta_bin && 
                               this->getPhiBin(muon.phi) == phi_bin;
                    });
            };

            auto calculateInvPt = [this, eta_bin, phi_bin](const std::vector<GenMuonInfo>& muons) -> double {
                double sum_inv_pt = 0.0;
                int count = 0;
                for (const auto& muon : muons) {
                    if (this->getEtaBin(muon.eta) == eta_bin && 
                        this->getPhiBin(muon.phi) == phi_bin && 
                        muon.pt > 0) {
                        sum_inv_pt += 1.0 / muon.pt;
                        count++;
                    }
                }
                return count > 0 ? sum_inv_pt / count : 0.0;
            };

            auto filtered_frame = frame.Filter(filterFunc, {"genMuon_charge", "genMuon_eta", "genMuon_phi"});  // genMuons 브랜치 사용
            auto inv_pt_frame = filtered_frame.Define("inv_pt", calculateInvPt, {"genMuons"});
            return *inv_pt_frame.Mean("inv_pt");
        } else {
            auto filterFunc = [this, charge, eta_bin, phi_bin](const std::vector<SelectedMuon>& muons) -> bool {
                for (const auto& muon : muons) {
                    if (muon.isTightRecoZ && muon.charge == charge) {
                        int current_eta_bin = this->getEtaBin(muon.eta);
                        int current_phi_bin = this->getPhiBin(muon.phi);
                        if (current_eta_bin == eta_bin && current_phi_bin == phi_bin) {
                            return true;
                        }
                    }
                }
                return false;
            };

            auto filtered_frame = frame.Filter(filterFunc, {"muons"});

            auto calculateInvPt = [this, eta_bin, phi_bin](const std::vector<SelectedMuon>& muons) -> double {
                double sum_inv_pt = 0.0;
                int nMuon = 0;
                for (const auto& muon : muons) {
                    int current_eta_bin = this->getEtaBin(muon.eta);
                    int current_phi_bin = this->getPhiBin(muon.phi);
                    if (current_eta_bin == eta_bin && current_phi_bin == phi_bin && muon.pt > 0) {
                        sum_inv_pt += 1.0 / muon.pt;
                        nMuon++;
                    }
                }
                return (nMuon > 0) ? sum_inv_pt / nMuon : 0.0;
            };

            auto inv_pt_frame = filtered_frame.Define("inv_pt", calculateInvPt, {"muons"});
            double average = *inv_pt_frame.Mean("inv_pt");
            return average;
        }
    }

    std::array<double, 2> calculateCorrectionFactors(ROOT::RDataFrame& gen_frame, ROOT::RDataFrame& reco_frame, int eta_bin, int phi_bin) {        
        double gen_inv_pt_m = calculateInversePtAverage(gen_frame, -1, eta_bin, phi_bin, true);
        double gen_inv_pt_p = calculateInversePtAverage(gen_frame, 1, eta_bin, phi_bin, true);

        double reco_inv_pt_avg_m = fine_tune_factors[std::make_tuple(-1, eta_bin, phi_bin)] * 
                                  calculateInversePtAverage(reco_frame, -1, eta_bin, phi_bin, false);
        double reco_inv_pt_avg_p = fine_tune_factors[std::make_tuple(1, eta_bin, phi_bin)] * 
                                  calculateInversePtAverage(reco_frame, 1, eta_bin, phi_bin, false);

        std::cout << "GEN inverse pt averages (minus/plus): " << gen_inv_pt_m << " / " << gen_inv_pt_p << std::endl;
        std::cout << "Data inverse pt averages (minus/plus): " << reco_inv_pt_avg_m << " / " << reco_inv_pt_avg_p << std::endl;

        double C_m = gen_inv_pt_m - reco_inv_pt_avg_m;
        double C_p = gen_inv_pt_p - reco_inv_pt_avg_p;

        double D_m = (C_p + C_m) / 2.0;
        double D_a = (C_p - C_m) / 2.0;

        double denominator = reco_inv_pt_avg_m + reco_inv_pt_avg_p;
        if (denominator == 0) {
            std::cerr << "Warning: denominator is 0 for eta_bin: " << eta_bin << ", phi_bin: " << phi_bin << std::endl;
            return {0.0, 1.0};  // 기본값 반환
        }

        double A = D_a - ((D_m * (reco_inv_pt_avg_m - reco_inv_pt_avg_p)) / denominator);
        double M = 1 + (2 * D_m / denominator);

        std::cout << "Calculated correction factors A: " << A << ", M: " << M << std::endl;
        return {A, M};
    }




/////////////////////////////
//                         //
//         public          //
//                         //
/////////////////////////////

public:
    
    std::vector<TTree*> splited_reco_trees;
    std::vector<TTree*> splited_gen_trees;
    std::vector<TFile*> output_files;

    TTree* GEN_TREE;
    TTree* RECO_TREE;
    
    void set_bin() {
        // Initialize binning
        eta_bins = {0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4};
        phi_bins.resize(11);
        double phi_step = M_PI / 10;
        for (int i = 0; i < 11; ++i) {
            phi_bins[i] = i * phi_step;
        }
        charge_bins = {-1, 1};

        for (size_t i = 0; i < eta_bins.size() - 1; ++i) {
            for (size_t j = 0; j < phi_bins.size() - 1; ++j) {
                for (int charge : charge_bins) {                    
                    fine_tune_factors[std::make_tuple(charge, i, j)] = 1.0;
                }
            }
        }
    }

    void factors_reset(){
        // reset fine tuning factors and correction factors
        for (size_t i = 0; i < eta_bins.size() - 1; ++i) {
            for (size_t j = 0; j < phi_bins.size() - 1; ++j) {
                correction_factors[std::make_tuple(i, j)] = std::array<double, 2>{0.0, 1.0};
            }
        }
    }


    void analyze(const std::string& gen_file_path, const std::string& reco_file_path) {
        std::unique_ptr<TFile> GEN_FILE(TFile::Open((gen_file_path).c_str(), "READ"));
        std::unique_ptr<TFile> RECO_FILE(TFile::Open((reco_file_path).c_str(), "READ"));

        if (!GEN_FILE || !RECO_FILE) {
            throw std::runtime_error("Failed to open input files");
        }

        std::unique_ptr<TTree> GEN_TREE(static_cast<TTree*>(GEN_FILE->Get("Analysis")));
        std::unique_ptr<TTree> RECO_TREE(static_cast<TTree*>(RECO_FILE->Get("Analysis")));

        if (!GEN_TREE || !RECO_TREE) {
            throw std::runtime_error("Failed to access trees");
        }

        ROOT::RDataFrame GEN_FRAME("Analysis", GEN_FILE.get());
        ROOT::RDataFrame RECO_FRAME("Analysis", RECO_FILE.get());
        
        try {
            factors_reset();

            for (const auto& [key, value] : correction_factors) {
                size_t eta_bin = std::get<0>(key);
                size_t phi_bin = std::get<1>(key);
                calculateCorrectionFactors(GEN_FRAME, RECO_FRAME, eta_bin, phi_bin);
            }
            
            // pt 보정
            auto corrected_frame = ptCorrection(RECO_FRAME, std::string(reco_name) + "_corrected" );

            // Z mass 재계산
            calculateZMassStats(corrected_frame);

            // fine_tune_factors의 제곱의 합 계산
            double sum_squared_diff = 0.0;
            int total_factors = 0;
            
            for (const auto& [key, value] : fine_tune_factors) {
                double diff_from_one = value - 1.0;
                sum_squared_diff += diff_from_one * diff_from_one;
                total_factors++;
            }
            
            double mean_squared_diff = sum_squared_diff / total_factors;
            std::cout << "Mean squared difference from 1.0: " << mean_squared_diff << std::endl;

            
        } catch (const std::exception& e) {
            std::cerr << "Exception occurred: " << e.what() << std::endl;
            throw;
        }
    }
};


int muonCorrection_tmp() {

    MuonAnalyzer analyzer;
    analyzer.set_bin();
    
    std::string gen_path = std::string(path) + gen_name;
    std::string reco_path = std::string(path) + reco_name +".root";
    
    // Check GEN file
    TFile* gen_file = TFile::Open(gen_path.c_str());
    if (!gen_file) {
        std::cerr << "Could not open gen file: " << gen_path << std::endl;
        return 1;
    }
    
    // List contents of GEN file
    gen_file->ls();
    
    // Check reco file
    TFile* reco_file = TFile::Open(reco_path.c_str());
    if (!reco_file) {
        std::cerr << "Could not open reco file: " << reco_path << std::endl;
        gen_file->Close();
        return 1;
    }
    
    // List contents of reco file
    reco_file->ls();

    // Clean up
    gen_file->Close();
    reco_file->Close();
    delete gen_file;
    delete reco_file;
    
    // Run the analysis
    analyzer.analyze(gen_path, reco_path);
    
    return 0;  // Return success
}