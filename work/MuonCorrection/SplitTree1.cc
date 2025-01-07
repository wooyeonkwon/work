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
#include <thread>
#include "interface/AnalysisClasses.h"
#include <mutex>
#include <TROOT.h>
#include <array>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/TTreeProcessorMT.hxx"

class MuonAnalyzer {
private:
    // Binning parameters
    std::vector<double> eta_bins;
    std::vector<double> phi_bins;
    std::vector<int> charge_bins;
    std::mutex treeMutex;



    int getEtaBin(double eta) {
        // eta_bins는 [0, 0.9, 1.2, 2.1, 2.4]로 초기화되어 있다고 가정
        for (int i = 0; i < eta_bins.size() - 1; ++i) {
            if (eta >= eta_bins[i] && eta < eta_bins[i + 1]) {
                return i;  // 해당 범위에 속하는 bin의 인덱스를 반환
            }
        }
        return -1;  // 범위를 벗어나는 경우 -1 반환
    }

    int getPhiBin(double phi) {
        // phi_bins는 [0, π/10, 2π/10, ..., 9π/10]으로 초기화되어 있다고 가정
        for (int i = 0; i < phi_bins.size() - 1; ++i) {
            if (phi >= phi_bins[i] && phi < phi_bins[i + 1]) {
                return i;  // 해당 범위에 속하는 bin의 인덱스를 반환
            }
        }
        return -1;  // 범위를 벗어나는 경우 -1 반환
    }


    template<typename MuonType>
    void processMuon(const MuonType& muon, bool isValid,
                    const std::vector<double>& eta_bins,
                    const std::vector<double>& phi_bins,
                    std::vector<TTree*>& splited_trees) {
        if (!isValid) return;

        double eta = std::abs(muon.eta);
        double phi = std::abs(muon.phi);

        for (size_t eta_bin = 0; eta_bin < eta_bins.size() - 1; ++eta_bin) {
            for (size_t phi_bin = 0; phi_bin < phi_bins.size() - 1; ++phi_bin) {
                if (eta >= eta_bins[eta_bin] && eta < eta_bins[eta_bin + 1] &&
                    phi >= phi_bins[phi_bin] && phi < phi_bins[phi_bin + 1]) {
                    size_t index = phi_bin + eta_bin * (phi_bins.size() - 1);
                    
                    std::lock_guard<std::mutex> lock(treeMutex);
                    splited_trees[index]->Fill();
                }
            }
        }
    }


public:
    

    TTree* MC_TREE;
    TTree* DATA_TREE;
    std::vector<std::thread> threads;
    

    void set_bin() {
        // Initialize binning
        eta_bins = {0, 0.9, 1.2, 2.1, 2.4};
        phi_bins.resize(11);
        double phi_step = M_PI / 10;
        for (int i = 0; i < 11; ++i) {
            phi_bins[i] = i * phi_step;
        }
        charge_bins = {-1, 1};
    }

    void processRange(TTree* tree, const std::vector<double>& eta_bins, 
                     const std::vector<double>& phi_bins, 
                     const std::vector<TFile*>& output_files,
                     std::vector<TTree*>& splited_trees, 
                     std::mutex& treeMutex, 
                     Long64_t start, Long64_t end, bool is_gen = false) {
        
        try {
            if (is_gen) {
                std::vector<GenMuonInfo>* genMuons = nullptr;
                {
                    std::lock_guard<std::mutex> lock(treeMutex);
                    tree->SetBranchAddress("genMuons", &genMuons);
                }
                for (Long64_t i = start; i < end; ++i) {
                    {
                        std::lock_guard<std::mutex> lock(treeMutex);
                        tree->GetEntry(i);
                    }

                    for (const auto& muon : *genMuons) {
                        processMuon(muon, muon.isGenZ, eta_bins, phi_bins, splited_trees);
                    }
                }
            } else {
                std::vector<SelectedMuon>* muons;
                {
                    std::lock_guard<std::mutex> lock(treeMutex);
                    tree->SetBranchAddress("muons", &muons);
                }
                for (Long64_t i = start; i < end; ++i) {
                    {
                        std::lock_guard<std::mutex> lock(treeMutex);
                        tree->GetEntry(i);
                    }

                    for (const auto& muon : *muons) {
                        processMuon(muon, muon.isTightRecoZ, eta_bins, phi_bins, splited_trees);
                    }
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Error in processRange: " << e.what() << std::endl;
            throw;
        }
    }




    void splitTreeMCAndSaveToFile(TTree* input_tree, const std::vector<double>& eta_bins, 
                                  const std::vector<double>& phi_bins) {
        std::string output_dir = "mc";
        if (mkdir(output_dir.c_str(), 0777) && errno != EEXIST) {
            throw std::runtime_error("Could not create output directory '" + output_dir + "'");
        }

        std::vector<TFile*> output_files;
        std::vector<TTree*> splited_trees;

        try {
            // Create output files and trees
            for (size_t eta_bin = 0; eta_bin < eta_bins.size() - 1; ++eta_bin) {
                for (size_t phi_bin = 0; phi_bin < phi_bins.size() - 1; ++phi_bin) {
                    std::string filename = output_dir + "/mc_tree_eta_" + 
                                        std::to_string(eta_bin) + "_phi_" + 
                                        std::to_string(phi_bin) + ".root";
                    
                    TFile* file = TFile::Open(filename.c_str(), "RECREATE");
                    if (!file) {
                        throw std::runtime_error("Failed to create output file: " + filename);
                    }
                    output_files.push_back(file);

                    std::string tree_name = "tree_eta_" + std::to_string(eta_bin) + 
                                          "_phi_" + std::to_string(phi_bin);
                    
                    TTree* new_tree = input_tree->CloneTree(0);
                    new_tree->SetName(tree_name.c_str());
                    new_tree->SetDirectory(file);
                    splited_trees.push_back(new_tree);
                }
            }

            // Process the tree in parallel
            Long64_t nentries = input_tree->GetEntries();
            unsigned int numThreads = std::max(1U, std::thread::hardware_concurrency());
            Long64_t entriesPerThread = nentries / numThreads;
            Long64_t remainder = nentries % numThreads;

            std::vector<std::thread> threads;
            for (unsigned int t = 0; t < numThreads; ++t) {
                Long64_t start = t * entriesPerThread;
                Long64_t end = (t == numThreads - 1) ? 
                              (start + entriesPerThread + remainder) : 
                              (start + entriesPerThread);

                threads.emplace_back(&MuonAnalyzer::processRange, this,
                                   input_tree, std::ref(eta_bins), std::ref(phi_bins),
                                   std::ref(output_files), std::ref(splited_trees),
                                   std::ref(treeMutex), start, end, false);
            }

            // Wait for all threads
            for (auto& thread : threads) {
                thread.join();
            }

            // Write and close files
            for (size_t i = 0; i < output_files.size(); ++i) {
                splited_trees[i]->Write();
                output_files[i]->Write();
                output_files[i]->Close();
            }

            // Cleanup
            for (auto file : output_files) {
                delete file; // This will also delete the associated trees
            }

        } catch (const std::exception& e) {
            // Cleanup in case of error
            for (auto file : output_files) {
                if (file) {
                    file->Close();
                    delete file;
                }
            }
            throw;
        }
    }
    void splitTreeDataAndSaveToFile(TTree* input_tree, const std::vector<double>& eta_bins, 
                                  const std::vector<double>& phi_bins) {
        std::string output_dir = "data";
        if (mkdir(output_dir.c_str(), 0777) && errno != EEXIST) {
            throw std::runtime_error("Could not create output directory '" + output_dir + "'");
        }

        std::vector<TFile*> output_files;
        std::vector<TTree*> splited_trees;

        try {
            // Create output files and trees
            for (size_t eta_bin = 0; eta_bin < eta_bins.size() - 1; ++eta_bin) {
                for (size_t phi_bin = 0; phi_bin < phi_bins.size() - 1; ++phi_bin) {
                    std::string filename = output_dir + "/data_tree_eta_" + 
                                        std::to_string(eta_bin) + "_phi_" + 
                                        std::to_string(phi_bin) + ".root";
                    
                    TFile* file = TFile::Open(filename.c_str(), "RECREATE");
                    if (!file) {
                        throw std::runtime_error("Failed to create output file: " + filename);
                    }
                    output_files.push_back(file);

                    std::string tree_name = "tree_eta_" + std::to_string(eta_bin) + 
                                          "_phi_" + std::to_string(phi_bin);
                    
                    TTree* new_tree = input_tree->CloneTree(0);
                    new_tree->SetName(tree_name.c_str());
                    new_tree->SetDirectory(file);
                    splited_trees.push_back(new_tree);
                }
            }

            // Process the tree in parallel
            Long64_t nentries = input_tree->GetEntries();
            unsigned int numThreads = std::max(1U, std::thread::hardware_concurrency());
            Long64_t entriesPerThread = nentries / numThreads;
            Long64_t remainder = nentries % numThreads;

            std::vector<std::thread> threads;
            for (unsigned int t = 0; t < numThreads; ++t) {
                Long64_t start = t * entriesPerThread;
                Long64_t end = (t == numThreads - 1) ? 
                              (start + entriesPerThread + remainder) : 
                              (start + entriesPerThread);

                threads.emplace_back(&MuonAnalyzer::processRange, this,
                                   input_tree, std::ref(eta_bins), std::ref(phi_bins),
                                   std::ref(output_files), std::ref(splited_trees),
                                   std::ref(treeMutex), start, end, false);
            }

            // Wait for all threads
            for (auto& thread : threads) {
                thread.join();
            }

            // Write and close files
            for (size_t i = 0; i < output_files.size(); ++i) {
                splited_trees[i]->Write();
                output_files[i]->Write();
                output_files[i]->Close();
            }

            // Cleanup
            for (auto file : output_files) {
                delete file; // This will also delete the associated trees
            }

        } catch (const std::exception& e) {
            // Cleanup in case of error
            for (auto file : output_files) {
                if (file) {
                    file->Close();
                    delete file;
                }
            }
            throw;
        }
    }

    // Similar changes for splitTreeMCAndSaveToFile...
    


    void analyze(const std::string& mc_file_path, const std::string& data_file_path) {
        // 입력 파일 열기
        std::unique_ptr<TFile> MC_FILE(TFile::Open(mc_file_path.c_str(), "READ"));
        std::unique_ptr<TFile> DATA_FILE(TFile::Open(data_file_path.c_str(), "READ"));

        if (!MC_FILE || !DATA_FILE) {
            throw std::runtime_error("Failed to open input files");
        }

        // MC와 데이터 트리 가져오기
        std::unique_ptr<TTree> MC_TREE(static_cast<TTree*>(MC_FILE->Get("AnalysisMC/Analysis")));
        std::unique_ptr<TTree> DATA_TREE(static_cast<TTree*>(DATA_FILE->Get("AnalysisMC/Analysis")));

        if (!MC_TREE || !DATA_TREE) {
            throw std::runtime_error("Failed to access trees");
        }

        // ROOT 멀티스레딩 비활성화
        //ROOT::DisableImplicitMT();
        // 데이터 트리 분할 및 저장

        // MC 트리 분할 및 저장
        std::cout << "Splitting and saving MC tree..." << std::endl;
        splitTreeMCAndSaveToFile(MC_TREE.get(), eta_bins, phi_bins);
        
        std::cout << "Splitting and saving DATA tree..." << std::endl;
        splitTreeDataAndSaveToFile(DATA_TREE.get(), eta_bins, phi_bins);



        std::cout << "Analysis complete!" << std::endl;
    }
};



int SplitTree1() {

    ROOT::EnableImplicitMT(56);

    MuonAnalyzer analyzer;
    analyzer.set_bin();

    // Load the library first
    if (gSystem->Load("libAnalysisClasses.so") < 0) {
        std::cerr << "Failed to load libAnalysisClasses.so" << std::endl;
        return 1;  // Return error code
    }

    // Check MC file
    const char* mc_path = "/data1/users/dndus0107/AnalysisResults/processed_data/DYto2Mu_MLL-50to120_22EEDR.root";
    TFile* mc_file = TFile::Open(mc_path);
    if (!mc_file) {
        std::cerr << "Could not open MC file: " << mc_path << std::endl;
        return 1;
    }
    
    // List contents of MC file
    mc_file->ls();
    
    // Check if AnalysisMC directory exists
    TDirectory* mc_dir = (TDirectory*)mc_file->Get("AnalysisMC");
    if (!mc_dir) {
        std::cerr << "Could not find AnalysisMC directory in MC file" << std::endl;
        mc_file->Close();
        return 1;
    }
    
    // List contents of AnalysisMC directory
    mc_dir->ls();
    
    // Check data file
    const char* data_path = "/data1/users/dndus0107/AnalysisResults/processed_data/DYto2Mu_MLL-50to120_22EEDR.root";
    TFile* data_file = TFile::Open(data_path);
    if (!data_file) {
        std::cerr << "Could not open data file: " << data_path << std::endl;
        mc_file->Close();
        return 1;
    }
    
    // List contents of data file
    data_file->ls();
    
    // Check if AnalysisMC directory exists in data file
    TDirectory* data_dir = (TDirectory*)data_file->Get("AnalysisMC");
    if (!data_dir) {
        std::cerr << "Could not find AnalysisMC directory in data file" << std::endl;
        mc_file->Close();
        data_file->Close();
        return 1;
    }
    
    // List contents of AnalysisMC directory in data file
    data_dir->ls();

    // Clean up
    mc_file->Close();
    data_file->Close();
    delete mc_file;
    delete data_file;
    
    // Run the analysis
    try {
        analyzer.analyze(mc_path,data_path);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;

}