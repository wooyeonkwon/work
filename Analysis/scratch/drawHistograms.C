#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"

void drawHistograms() {
    TFile *file = new TFile("data3.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    TTree *treeGlobal = nullptr;
    TTree *treeRPC = nullptr;
    TTree *treenotRPC = nullptr;
    file->GetObject("GlobalMuons", treeGlobal);
    file->GetObject("RPCMuons", treeRPC);
    file->GetObject("notRPCMuons", treenotRPC);

    if (!treeGlobal || !treeRPC || !treenotRPC) {
        std::cerr << "Error: Trees not found in file" << std::endl;
        file->Close();
        delete file;
        return;
    }

    std::cerr << "File read success" << std::endl;

    double zBosonMassGlobal;
    double zBosonMassRPC;
    double zBosonMassnotRPC;

    treeGlobal->SetBranchAddress("zBosonMass", &zBosonMassGlobal);
    treeRPC->SetBranchAddress("zBosonMass", &zBosonMassRPC);
    treenotRPC->SetBranchAddress("zBosonMass", &zBosonMassnotRPC);

    TH1F *h_zBosonsGlobal = new TH1F("zBosonsGlobal", "Z Bosons (Global);Mass (GeV);Events", 360, 0, 180);
    TH1F *h_zBosonsRPC = new TH1F("zBosonsRPC", "Z Bosons (RPC);Mass (GeV);Events", 360, 0, 180);
    TH1F *h_zBosonsnotRPC = new TH1F("zBosonsnotRPC", "Z Bosons (notRPC);Mass (GeV);Events", 360, 0, 180);
    
    Long64_t nEntriesGlobal = treeGlobal->GetEntries();
    for (Long64_t i = 0; i < nEntriesGlobal; ++i) {
        treeGlobal->GetEntry(i);
        h_zBosonsGlobal->Fill(zBosonMassGlobal);
    }

    Long64_t nEntriesRPC = treeRPC->GetEntries();
    for (Long64_t i = 0; i < nEntriesRPC; ++i) {
        treeRPC->GetEntry(i);
        h_zBosonsRPC->Fill(zBosonMassRPC);
    }
    
    Long64_t nEntriesnotRPC = treenotRPC->GetEntries();
    for (Long64_t i = 0; i < nEntriesnotRPC; ++i) {
        treenotRPC->GetEntry(i);
        h_zBosonsnotRPC->Fill(zBosonMassnotRPC);
    }
    
    TCanvas *c1 = new TCanvas("c1", "Histograms", 800, 600);
    h_zBosonsGlobal->Draw();
    TF1 *fitFuncGlobal = new TF1("fitFuncGlobal", "gaus(0)", 10, 170);
    
    // Set parameter limits for Global fit
    //fitFuncGlobal->SetParLimits(0, 500, 2000);  // Gaussian amplitude
    //fitFuncGlobal->SetParLimits(1, 80, 100);  // Gaussian mean
    //fitFuncGlobal->SetParLimits(2, 0, 5);    // Gaussian sigma
    //fitFuncGlobal->SetParLimits(3, 0, 10000);  // Exponential amplitude
    //fitFuncGlobal->SetParLimits(4, 0, 1);     // Exponential decay constant
    
    h_zBosonsGlobal->Fit(fitFuncGlobal, "R");
    c1->SaveAs("zBosonsGlobal.png");


    TCanvas *c2 = new TCanvas("c2", "Histograms", 800, 600);
    h_zBosonsRPC->Draw();
    TF1 *fitFuncRPC = new TF1("fitFuncRPC", "gaus(0)", 10, 170);
    
    // Set parameter limits for RPC fit
    //fitFuncRPC->SetParLimits(0, 500, 2000);  // Gaussian amplitude
    //fitFuncRPC->SetParLimits(1, 80, 100);  // Gaussian mean
    //fitFuncRPC->SetParLimits(2, 0, 5);    // Gaussian sigma
    //fitFuncRPC->SetParLimits(3, 0, 10000);  // Exponential amplitude
    //fitFuncRPC->SetParLimits(4, 0, 1);     // Exponential decay constant
    
    h_zBosonsRPC->Fit(fitFuncRPC, "R");
    c2->SaveAs("zBosonsRPC.png");


    TCanvas *c3 = new TCanvas("c3", "Histograms", 800, 600);
    h_zBosonsnotRPC->Draw();
    TF1 *fitFuncnotRPC = new TF1("fitFuncnotRPC", "gaus(0)", 10, 170);
    
    // Set parameter limits for notRPC fit
    //fitFuncRPC->SetParLimits(0, 0, 1000);  // Gaussian amplitude
    //fitFuncRPC->SetParLimits(1, 80, 100);  // Gaussian mean
    //fitFuncRPC->SetParLimits(2, 0, 5);    // Gaussian sigma
    //fitFuncRPC->SetParLimits(3, 0, 10000);  // Exponential amplitude
    //fitFuncRPC->SetParLimits(4, 0, 1);     // Exponential decay constant
    
    h_zBosonsnotRPC->Fit(fitFuncnotRPC, "R");
    c3->SaveAs("zBosonsnotRPC.png");


    std::ofstream logFile("fit_results.txt");
    if (logFile.is_open()) {
        logFile << "Global Fit Parameters:\n"
                << "Gaussian Amplitude: " << fitFuncGlobal->GetParameter(0) << "\n"
                << "Gaussian Mean: " << fitFuncGlobal->GetParameter(1) << "\n"
                << "Gaussian Sigma: " << fitFuncGlobal->GetParameter(2) << "\n";
//                << "Exponential Amplitude: " << fitFuncGlobal->GetParameter(3) << "\n"
//                << "Exponential Decay: " << fitFuncGlobal->GetParameter(4) << "\n";
        logFile << "RPC Fit Parameters:\n"
                << "Gaussian Amplitude: " << fitFuncRPC->GetParameter(0) << "\n"
                << "Gaussian Mean: " << fitFuncRPC->GetParameter(1) << "\n"
                << "Gaussian Sigma: " << fitFuncRPC->GetParameter(2) << "\n";
//                << "Exponential Amplitude: " << fitFuncRPC->GetParameter(3) << "\n"
//                << "Exponential Decay: " << fitFuncRPC->GetParameter(4) << "\n";
        logFile << "notRPC Fit Parameters:\n"
                << "Gaussian Amplitude: " << fitFuncnotRPC->GetParameter(0) << "\n"
                << "Gaussian Mean: " << fitFuncnotRPC->GetParameter(1) << "\n"
                << "Gaussian Sigma: " << fitFuncnotRPC->GetParameter(2) << "\n";
//                << "Exponential Amplitude: " << fitFuncnotRPC->GetParameter(3) << "\n"
//                << "Exponential Decay: " << fitFuncnotRPC->GetParameter(4) << "\n";
        logFile.close();
    } else {
        std::cerr << "Unable to open log file for writing fit results." << std::endl;
    }

    delete h_zBosonsGlobal;
    delete h_zBosonsRPC;
    delete h_zBosonsnotRPC;    
    delete fitFuncGlobal;
    delete fitFuncRPC;
    delete fitFuncnotRPC;
    delete c1;
    delete c2;
    delete c3;
    file->Close();
    delete file;
}
