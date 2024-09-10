#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"

double dscbf(double *x, double *par)
{ 
    double alpha_l = par[0]; 
    double alpha_h = par[1]; 
    double n_l     = par[2]; 
    double n_h     = par[3]; 
    double mean    = par[4]; 
    double sigma   = par[5];
    double N       = par[6];
    float t = (x[0]-mean)/sigma;
    double result;
    double fact1TLessMinosAlphaL = alpha_l/n_l;
    double fact2TLessMinosAlphaL = (n_l/alpha_l) - alpha_l -t;
    double fact1THihgerAlphaH = alpha_h/n_h;
    double fact2THigherAlphaH = (n_h/alpha_h) - alpha_h +t;
   
    if (-alpha_l <= t && alpha_h >= t)
    {
        result = exp(-0.5*t*t);
    }
    else if (t < -alpha_l)
    {
        result = exp(-0.5*alpha_l*alpha_l)*pow(fact1TLessMinosAlphaL*fact2TLessMinosAlphaL, -n_l);
    }
    else if (t > alpha_h)
    {
        result = exp(-0.5*alpha_h*alpha_h)*pow(fact1THihgerAlphaH*fact2THigherAlphaH, -n_h);
    }

    return N*result;
}

void drawHistograms_cbf(const char* filename) {
    TFile *file = new TFile(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }


    TTree *treeGlobal = nullptr;
    TTree *treeRPC = nullptr;
    file->GetObject("GlobalMuons", treeGlobal);
    file->GetObject("RPCMuons", treeRPC);

    if (!treeGlobal || !treeRPC) {
        std::cerr << "Error: Trees not found in file" << std::endl;
        file->Close();
        delete file;
        return;
    }

    std::cerr << "File read success" << std::endl;

    double zBosonMassGlobal;
    double zBosonMassRPC;

    treeGlobal->SetBranchAddress("zBosonMass", &zBosonMassGlobal);
    treeRPC->SetBranchAddress("zBosonMass", &zBosonMassRPC);

    TH1F *h_zBosonsGlobal = new TH1F("zBosonsGlobal", "Z Bosons (Global);Mass (GeV);Events", 360, 0, 180);
    TH1F *h_zBosonsRPC = new TH1F("zBosonsRPC", "Z Bosons (RPC);Mass (GeV);Events", 360, 0, 180);
    
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
    

    std::ofstream logFile("fit_results.txt");
    if (!logFile.is_open()) {
        std::cerr << "Unable to open log file for writing fit results." << std::endl;
        delete h_zBosonsGlobal;
        delete h_zBosonsRPC;
        file->Close();
        delete file;
        return;
    }

    // Fit and plot for Global Muons
    TCanvas *c1 = new TCanvas("c1", "Global Z mass", 800, 600);
    gStyle->SetOptStat(0);
    h_zBosonsGlobal->Draw();
    TF1 *fitFuncGlobal = new TF1("fitFuncGlobal", dscbf, 60, 120, 7);
    int NeventsGlobal = h_zBosonsGlobal->GetEntries();
    fitFuncGlobal->SetParNames("AlphaL", "AlphaH", "nL", "nH", "Mean", "Sigma", "Events");
    fitFuncGlobal->SetParameters(1.25, 1.44, 1.59, 1.8, 91.2, 2.0, NeventsGlobal);
    h_zBosonsGlobal->Fit(fitFuncGlobal, "R");

    // Add text to the canvas
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.70, 0.85, Form("Entries: %d", NeventsGlobal));
    latex.DrawLatex(0.70, 0.80, Form("Mean: %.2f", fitFuncGlobal->GetParameter(4)));
    latex.DrawLatex(0.70, 0.75, Form("Sigma: %.2f", fitFuncGlobal->GetParameter(5)));

    c1->SaveAs("zBosonsGlobal.png");

    logFile << "Global Fit Parameters:\n"
            << "AlphaL: " << fitFuncGlobal->GetParameter(0) << "\n"
            << "AlphaH: " << fitFuncGlobal->GetParameter(1) << "\n"
            << "nL: " << fitFuncGlobal->GetParameter(2) << "\n"
            << "nH: " << fitFuncGlobal->GetParameter(3) << "\n"
            << "Mean: " << fitFuncGlobal->GetParameter(4) << "\n"
            << "Sigma: " << fitFuncGlobal->GetParameter(5) << "\n"
            << "Amplitude: " << fitFuncGlobal->GetParameter(6) << "\n\n";

    // Fit and plot for RPC Muons
    TCanvas *c2 = new TCanvas("c2", "RPC Z mass", 800, 600);
    gStyle->SetOptStat(0);
    h_zBosonsRPC->Draw();
    TF1 *fitFuncRPC = new TF1("fitFuncRPC", dscbf, 60, 120, 7);
    int NeventsRPC = h_zBosonsRPC->GetEntries();
    fitFuncRPC->SetParNames("AlphaL", "AlphaH", "nL", "nH", "Mean", "Sigma", "Events");
    fitFuncRPC->SetParameters(1.25, 1.44, 1.59, 1.8, 91.2, 2.0, NeventsRPC);
    h_zBosonsRPC->Fit(fitFuncRPC, "R");

    // Add text to the canvas
    latex.DrawLatex(0.70, 0.85, Form("Entries: %d", NeventsRPC));
    latex.DrawLatex(0.70, 0.80, Form("Mean: %.2f", fitFuncRPC->GetParameter(4)));
    latex.DrawLatex(0.70, 0.75, Form("Sigma: %.2f", fitFuncRPC->GetParameter(5)));

    c2->SaveAs("zBosonsRPC.png");

    logFile << "RPC Fit Parameters:\n"
            << "AlphaL: " << fitFuncRPC->GetParameter(0) << "\n"
            << "AlphaH: " << fitFuncRPC->GetParameter(1) << "\n"
            << "nL: " << fitFuncRPC->GetParameter(2) << "\n"
            << "nH: " << fitFuncRPC->GetParameter(3) << "\n"
            << "Mean: " << fitFuncRPC->GetParameter(4) << "\n"
            << "Sigma: " << fitFuncRPC->GetParameter(5) << "\n"
            << "Amplitude: " << fitFuncRPC->GetParameter(6) << "\n\n";

    logFile.close();

    delete h_zBosonsGlobal;
    delete h_zBosonsRPC;
    delete fitFuncGlobal;
    delete fitFuncRPC;
    delete c1;
    delete c2;
    file->Close();
    delete file;
}
