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

class MuonAnalyzer {
private:
    // Binning parameters
    std::vector<double> eta_bins;
    std::vector<double> phi_bins;
    std::vector<int> charge_bins;
    
    // Fine tuning factors
    std::map<std::tuple<int, int, int>, double> fine_tune_factors;
    
    // Constants
    const double Z_MASS = 91.1876;
    const double MAX_ITERATIONS = 10;
    const double CONVERGENCE_THRESHOLD = 0.001;

    // Helper function to get bin index
    int getBin(double value, const std::vector<double>& bins) {
        for (size_t i = 0; i < bins.size() - 1; ++i) {
            if (value >= bins[i] && value < bins[i + 1]) return i;
        }
        return -1;
    }

    // Calculate Z mass statistics for a given bin
    std::pair<double, double> calculateZMassStats(TTree* tree, int charge, int eta_bin, int phi_bin) {
        std::vector<double>* muon_pt = nullptr;
        std::vector<double>* muon_eta = nullptr;
        std::vector<double>* muon_phi = nullptr;
        std::vector<int>* muon_charge = nullptr;
        std::vector<bool>* muon_isTightReco = nullptr;
        std::vector<bool>* muon_isTightRecoZ = nullptr;
        std::vector<double>* z_mass = nullptr;

        tree->SetBranchAddress("muons.pt", &muon_pt);
        tree->SetBranchAddress("muons.eta", &muon_eta);
        tree->SetBranchAddress("muons.phi", &muon_phi);
        tree->SetBranchAddress("muons.charge", &muon_charge);
        tree->SetBranchAddress("muons.isTightReco", &muon_isTightReco);
        tree->SetBranchAddress("muons.isTightRecoZ", &muon_isTightRecoZ);
        tree->SetBranchAddress("zBoson.mass", &z_mass);

        double sum_mass = 0;
        double sum_dmass = 0;
        int count = 0;

        Long64_t nentries = tree->GetEntries();
        for (Long64_t i = 0; i < nentries; ++i) {
            tree->GetEntry(i);

            for (size_t j = 0; j < muon_pt->size(); ++j) {
                if (!(*muon_isTightReco)[j] || !(*muon_isTightRecoZ)[j]) continue;
                if ((*muon_charge)[j] != charge) continue;

                double eta = std::abs((*muon_eta)[j]);
                double phi = std::abs((*muon_phi)[j]);

                if (getBin(eta, eta_bins) == eta_bin && getBin(phi, phi_bins) == phi_bin) {
                    sum_mass += (*z_mass)[1];  // Using index 1 for tightReco Z mass
                    sum_dmass += (*z_mass)[1] - Z_MASS;
                    count++;
                }
            }
        }

        if (count == 0) return {0, 0};
        return {sum_mass / count, sum_dmass / count};
    }

    // Calculate inverse pt average
    double calculateInversePtAverage(TTree* tree, int charge, int eta_bin, int phi_bin, bool is_gen = false) {
        std::string muon_type = is_gen ? "genMuons" : "muons";
        std::vector<double>* muon_pt = nullptr;
        std::vector<double>* muon_eta = nullptr;
        std::vector<double>* muon_phi = nullptr;
        std::vector<int>* muon_charge = nullptr;
        std::vector<bool>* muon_isTightReco = nullptr;
        std::vector<bool>* muon_isTightRecoZ = nullptr;

        tree->SetBranchAddress((muon_type + ".pt").c_str(), &muon_pt);
        tree->SetBranchAddress((muon_type + ".eta").c_str(), &muon_eta);
        tree->SetBranchAddress((muon_type + ".phi").c_str(), &muon_phi);
        tree->SetBranchAddress((muon_type + ".charge").c_str(), &muon_charge);
        
        if (!is_gen) {
            tree->SetBranchAddress("muons.isTightReco", &muon_isTightReco);
            tree->SetBranchAddress("muons.isTightRecoZ", &muon_isTightRecoZ);
        }

        double sum_inv_pt = 0;
        int count = 0;

        Long64_t nentries = tree->GetEntries();
        for (Long64_t i = 0; i < nentries; ++i) {
            tree->GetEntry(i);

            for (size_t j = 0; j < muon_pt->size(); ++j) {
                if (!is_gen && (!(*muon_isTightReco)[j] || !(*muon_isTightRecoZ)[j])) continue;
                if ((*muon_charge)[j] != charge) continue;

                double eta = std::abs((*muon_eta)[j]);
                double phi = std::abs((*muon_phi)[j]);

                if (getBin(eta, eta_bins) == eta_bin && getBin(phi, phi_bins) == phi_bin) {
                    sum_inv_pt += 1.0 / (*muon_pt)[j];
                    count++;
                }
            }
        }

        if (count == 0) return 0;
        return sum_inv_pt / count;
    }

    // Calculate correction factors
    struct CorrectionFactors {
        std::map<std::tuple<int, int, int>, double> Cs;
        std::map<std::tuple<int, int>, double> D_ms;
        std::map<std::tuple<int, int>, double> D_as;
        std::map<std::tuple<int, int>, double> As;
        std::map<std::tuple<int, int>, double> Ms;
    };

    CorrectionFactors calculateCorrectionFactors(TTree* mc_tree, TTree* data_tree) {
        CorrectionFactors factors;
        std::map<std::tuple<int, int, int>, double> data_inv_pt_avg;

        // Calculate Cs factors
        for (size_t eta_bin = 0; eta_bin < eta_bins.size() - 1; ++eta_bin) {
            for (size_t phi_bin = 0; phi_bin < phi_bins.size() - 1; ++phi_bin) {
                for (int charge : charge_bins) {
                    auto key = std::make_tuple(charge, eta_bin, phi_bin);
                    double mc_inv_pt = calculateInversePtAverage(mc_tree, charge, eta_bin, phi_bin, true);
                    data_inv_pt_avg[key] = fine_tune_factors[key] * 
                                         calculateInversePtAverage(data_tree, charge, eta_bin, phi_bin);
                    factors.Cs[key] = mc_inv_pt - data_inv_pt_avg[key];
                }
            }
        }

        // Calculate D_ms, D_as, As, and Ms factors
        for (size_t eta_bin = 0; eta_bin < eta_bins.size() - 1; ++eta_bin) {
            for (size_t phi_bin = 0; phi_bin < phi_bins.size() - 1; ++phi_bin) {
                auto pos_key = std::make_tuple(1, eta_bin, phi_bin);
                auto neg_key = std::make_tuple(-1, eta_bin, phi_bin);
                auto bin_key = std::make_tuple(eta_bin, phi_bin);

                factors.D_ms[bin_key] = (factors.Cs[pos_key] + factors.Cs[neg_key]) / 2.0;
                factors.D_as[bin_key] = (factors.Cs[pos_key] - factors.Cs[neg_key]) / 2.0;

                double denominator = data_inv_pt_avg[pos_key] + data_inv_pt_avg[neg_key];
                if (denominator != 0) {
                    factors.As[bin_key] = factors.D_as[bin_key] - 
                        ((factors.D_ms[bin_key] * (data_inv_pt_avg[pos_key] - data_inv_pt_avg[neg_key])) / denominator);
                    factors.Ms[bin_key] = 1 + (2 * factors.D_ms[bin_key] / denominator);
                } else {
                    std::cout << "Warning: Zero denominator for eta_bin=" << eta_bin 
                              << ", phi_bin=" << phi_bin << std::endl;
                    factors.As[bin_key] = 0;
                    factors.Ms[bin_key] = 1;
                }
            }
        }

        return factors;
    }

    // Update fine tuning factors
    void updateFineTuneFactors(TTree* tree) {
        for (size_t eta_bin = 0; eta_bin < eta_bins.size() - 1; ++eta_bin) {
            for (size_t phi_bin = 0; phi_bin < phi_bins.size() - 1; ++phi_bin) {
                for (int charge : charge_bins) {
                    auto [zmass, dzmass] = calculateZMassStats(tree, charge, eta_bin, phi_bin);
                    if (zmass != 0) {
                        auto key = std::make_tuple(charge, eta_bin, phi_bin);
                        fine_tune_factors[key] = 1 + (2 * dzmass / zmass);
                    }
                }
            }
        }
    }

public:
    MuonAnalyzer() {
        // Initialize binning
        eta_bins = {0, 0.9, 1.2, 2.1, 2.4};
        phi_bins.resize(11);
        double phi_step = M_PI / 10;
        for (int i = 0; i < 11; ++i) {
            phi_bins[i] = i * phi_step;
        }
        charge_bins = {-1, 1};

        // Initialize fine tuning factors
        for (size_t i = 0; i < eta_bins.size() - 1; ++i) {
            for (size_t j = 0; j < phi_bins.size() - 1; ++j) {
                for (int charge : charge_bins) {
                    fine_tune_factors[std::make_tuple(charge, i, j)] = 1.0;
                }
            }
        }
    }

    void analyze(const std::string& mc_file_path, const std::string& data_file_path) {
        // Open files
        std::unique_ptr<TFile> mc_file(TFile::Open(mc_file_path.c_str(), "READ"));
        std::unique_ptr<TFile> data_file(TFile::Open(data_file_path.c_str(), "READ"));

        if (!mc_file || !data_file) {
            std::cerr << "Error opening files" << std::endl;
            return;
        }

        TTree* mc_tree = (TTree*)mc_file->Get("AnalysisMC/Analysis");
        TTree* data_tree = (TTree*)data_file->Get("AnalysisMC/Analysis");

        if (!mc_tree || !data_tree) {
            std::cerr << "Error accessing trees" << std::endl;
            return;
        }

        // Perform iterations
        int iteration = 0;
        bool convergence = false;

        try {
            while (!convergence && iteration < MAX_ITERATIONS) {
                std::cout << "Iteration " << (iteration + 1) << std::endl;

                // Calculate correction factors
                auto correction_factors = calculateCorrectionFactors(mc_tree, data_tree);

                // Store previous factors for convergence check
                auto previous_factors = fine_tune_factors;

                // Update fine tune factors
                updateFineTuneFactors(data_tree);

                // Check convergence
                double max_diff = 0.0;
                for (const auto& [key, value] : fine_tune_factors) {
                    double diff = std::abs(value - previous_factors[key]);
                    max_diff = std::max(max_diff, diff);
                }

                convergence = max_diff < CONVERGENCE_THRESHOLD;
                iteration++;

                std::cout << "Maximum difference: " << max_diff << std::endl;
            }
        }
        catch (const std::exception& e) {
            std::cerr << "Error during analysis: " << e.what() << std::endl;
        }

        // Save results
        std::unique_ptr<TFile> output_file(TFile::Open("data_corrected.root", "RECREATE"));
        if (output_file) {
            output_file->cd();
            // Here you would write the corrected data
            // This part would need to be implemented based on your specific needs
            std::cout << "Results saved to data_corrected.root" << std::endl;
        }
    }
};

int main() {
    MuonAnalyzer analyzer;
    analyzer.analyze(
        "/data1/users/dndus0107/AnalysisResults/processed_data/DYto2Mu_MLL-50to120_22EEDR.root",
        "/data1/users/dndus0107/AnalysisResults/processed_data/DYto2Mu_MLL-50to120_22EEDR.root"
    );
    return 0;
}