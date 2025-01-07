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


void drawHistogramAndFit(TTree* tree, const char* branchName, const char* histTitle, const char* histFileName, const char* canvasTitle, const char* logFileName, double xMin = 60, double xMax = 120) {
    double zBosonMass;
    tree->SetBranchAddress(branchName, &zBosonMass);
    TH1F* hist = new TH1F(histTitle, Form("%s;M_{\\mu\\mu} (GeV);Counts /0.5GeV", histTitle), 360, 0, 180);



    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        hist->Fill(zBosonMass);
    }


    TCanvas* canvas = new TCanvas(canvasTitle, canvasTitle, 800, 600);
    gStyle->SetOptStat(0);
    hist->SetFillColor(kBlue);
    hist->Draw();
    TF1* fitFunc = new TF1(Form("fitFunc%s", canvasTitle), dscbf, xMin, xMax, 7);
    fitFunc->SetParNames("AlphaL", "AlphaH", "nL", "nH", "Mean", "Sigma", "Events");
    fitFunc->SetParameters(1.25, 1.44, 1.59, 1.8, 91.2, 2.0, nEntries);
    hist->Fit(fitFunc, "R");
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextFont(62); // bold font for "CMS"
    latex.DrawLatex(0.15, 0.92, "CMS");
    latex.SetTextFont(42); // normal font for "Preliminary"
    latex.DrawLatex(0.21, 0.92, "#it{In Progress}");
    latex.DrawLatex(0.60, 0.92, "#sqrt{s} = 13.6 TeV, L =  11.06/fb");
    latex.DrawLatex(0.70, 0.85, Form("nEntries: %d", (int)nEntries));
    latex.DrawLatex(0.70, 0.80, Form("Mean: %.2f", fitFunc->GetParameter(4)));
    latex.DrawLatex(0.70, 0.75, Form("Sigma: %.2f", fitFunc->GetParameter(5)));


    canvas->SaveAs(histFileName);


    std::ofstream logFile(logFileName, std::ios_base::app);
    if (logFile.is_open()) {
        logFile << canvasTitle << " Fit Parameters:\n"
                << "AlphaL: " << fitFunc->GetParameter(0) << "\n"
                << "AlphaH: " << fitFunc->GetParameter(1) << "\n"
                << "nL: " << fitFunc->GetParameter(2) << "\n"
                << "nH: " << fitFunc->GetParameter(3) << "\n"
                << "Mean: " << fitFunc->GetParameter(4) << "\n"
                << "Sigma: " << fitFunc->GetParameter(5) << "\n"
                << "Amplitude: " << fitFunc->GetParameter(6) << "\n\n";
        logFile.close();
    }

    delete hist;
    delete fitFunc;
    delete canvas;
}

void drawHistograms_cbf_2(const char* filename) {
    TFile* file = new TFile(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // define trees and get objects

    TTree *treeRPCPass = nullptr;
    TTree *treeRPCFail = nullptr;
    TTree *treeGEMPass = nullptr;
    TTree *treeGEMFail = nullptr;
    file->GetObject("RPCMuons", treeRPCPass);
    file->GetObject("RPCMuons", treeRPCFail);
    file->GetObject("GEMMuons", treeGEMPass);
    file->GetObject("GEMMuons", treeGEMFail);


    // draw hist for exist trees
    if (!treeRPCPass) std::cerr << "RPCMuons tree not found." << std::endl;
    else {
        drawHistogramAndFit(treeRPCPass, "zBosonMass", "RPC muon Pass", "fit_zBosonsRPCPass.png", "RPC Z mass Pass", "fit_results.txt");
    }

    if (!treeRPCFail) std::cerr << "RPCMuons tree not found." << std::endl;
    else {
        drawHistogramAndFit(treeRPCFail, "zBosonMass", "RPC muon Fail", "fit_zBosonsRPCFail.png", "RPC Z mass Fail", "fit_results.txt");
    }

    if (!treeGEMPass) std::cerr << "GEMMuons tree not found." << std::endl;
    else {
        drawHistogramAndFit(treeGEMPass, "zBosonMassPass", "GEM muon Pass", "fit_zBosonsGEMPass.png", "GEM Z mass Pass", "fit_results.txt");
    }

    if (!treeGEMFail) std::cerr << "GEMMuons tree not found." << std::endl;
    else {
        drawHistogramAndFit(treeGEMFail, "zBosonMassFail", "GEM muon Fail", "fit_zBosonsGEMFail.png", "GEM Z mass Fail", "fit_results.txt");
    }

    file->Close();
    delete file;
}