#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TMath.h>


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
    hist->Draw();
    TF1* fitFunc = new TF1(Form("fitFunc%s", canvasTitle), dscbf, xMin, xMax, 7);
    fitFunc->SetParNames("AlphaL", "AlphaH", "nL", "nH", "Mean", "Sigma", "Events");
    fitFunc->SetParameters(1.25, 1.44, 1.59, 1.8, 91.2, 2.0, nEntries);
    hist->Fit(fitFunc, "R");
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextFont(62); // bold font for "CMS"
    latex.DrawLatex(0.15, 0.95, "CMS");
    latex.SetTextFont(42); // normal font for "Preliminary"
    latex.DrawLatex(0.15, 0.92, "#it{Preliminary}");
    latex.DrawLatex(0.60, 0.92, "#sqrt{s} = 13.6 TeV, L = 120.57/fb");
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

// Efficiency drawing function
void drawEfficiency(TTree* tree, TH1F* histPt, TH1F* histEta, TH1F* histIso, int color) {
    double probePt, probeEta, probeIso;
    Bool_t passingProbe;

    tree->SetBranchAddress("probePt", &probePt);
    tree->SetBranchAddress("probeEta", &probeEta);
    tree->SetBranchAddress("probeIso", &probeIso);
    tree->SetBranchAddress("passingProbe", &passingProbe);

    // Temporary histograms to store total and passing events
    TH1F* histPtAll = (TH1F*)histPt->Clone("histPtAll");
    TH1F* histEtaAll = (TH1F*)histEta->Clone("histEtaAll");
    TH1F* histIsoAll = (TH1F*)histIso->Clone("histIsoAll");

    histPtAll->Reset();
    histEtaAll->Reset();
    histIsoAll->Reset();
    histPt->Reset();
    histEta->Reset();
    histIso->Reset();

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Fill the histograms with all events
        histPtAll->Fill(probePt);
        histEtaAll->Fill(probeEta);
        histIsoAll->Fill(probeIso);

        // Fill the passing histograms only if the probe passed
        if (passingProbe) {
            histPt->Fill(probePt);
            histEta->Fill(probeEta);
            histIso->Fill(probeIso);
        }
    }

    // Calculate efficiency by dividing the passing histogram by the total
    histPt->Divide(histPt, histPtAll);
    histEta->Divide(histEta, histEtaAll);
    histIso->Divide(histIso, histIsoAll);

    // Set histogram properties
    histPt->SetLineColor(color);
    histPt->SetMarkerColor(color);
    histPt->Draw();
    histEta->SetLineColor(color);
    histEta->SetMarkerColor(color);
    histEta->Draw();
    histIso->SetLineColor(color);
    histIso->SetMarkerColor(color);
    histEta->Draw();


    // Clean up temporary histograms
    delete histPtAll;
    delete histEtaAll;
    delete histIsoAll;
}


void drawHistograms(const char* filename) {
    TFile* file = new TFile(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Define trees and get objects
    TTree *treeGlobal = nullptr;
    TTree *treeTracker = nullptr;
    TTree *treeStandAlone = nullptr;
    TTree *treePF = nullptr;
    TTree *treeRPC = nullptr;
    TTree *treeGEM = nullptr;

    file->GetObject("GlobalMuons", treeGlobal);
    file->GetObject("TrackerMuons", treeTracker);
    file->GetObject("StandAloneMuons", treeStandAlone);
    file->GetObject("PFMuons", treePF);
    file->GetObject("RPCMuons", treeRPC);
    file->GetObject("GEMMuons", treeGEM);

    // Create histograms for each variable
    TH1F* histPt = new TH1F("histPt", "Muon Efficiency vs pT", 50, 0, 100);
    TH1F* histEta = new TH1F("histEta", "Muon Efficiency vs Eta", 50, -2.5, 2.5);
    TH1F* histIso = new TH1F("histIso", "Muon Efficiency vs Iso", 50, 0, 0.5);

    // Create canvases
    TCanvas* canvasPt = new TCanvas("canvas_pt", "Muon Efficiency vs pT", 800, 600);
    TCanvas* canvasEta = new TCanvas("canvas_eta", "Muon Efficiency vs eta", 800, 600);
    TCanvas* canvasIso = new TCanvas("canvas_iso", "Muon Efficiency vs Iso", 800, 600);

    // Create legends
    TLegend* legendPt = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend* legendEta = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend* legendIso = new TLegend(0.7, 0.7, 0.9, 0.9);

    // Draw histograms and efficiencies for existing trees
    if (!treeGlobal) std::cerr << "GlobalMuons tree not found." << std::endl;
    else {
        drawHistogramAndFit(treeGlobal, "zBosonMass", "Z (Global)", "zBosonsGlobal.png", "Global Z mass", "fit_results.txt");
        
        // Draw efficiency
        double probePt, probeEta, probeIso;
        Bool_t passingProbe;

        treeGlobal->SetBranchAddress("probePt", &probePt);
        treeGlobal->SetBranchAddress("probeEta", &probeEta);
        treeGlobal->SetBranchAddress("probeIso", &probeIso);
        treeGlobal->SetBranchAddress("passingProbe", &passingProbe);

        Long64_t nEntries = treeGlobal->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            treeGlobal->GetEntry(i);
            // Fill the histograms with all events
            histPt->Fill(probePt);
            histEta->Fill(probeEta);
            histIso->Fill(probeIso);

            // Fill the passing histograms only if the probe passed
            if (passingProbe) {
                histPt->Fill(probePt);
                histEta->Fill(probeEta);
                histIso->Fill(probeIso);
            }
        }
    }
    if (!treeTracker) std::cerr << "TrackerMuons tree not found." << std::endl;
    else {
        drawHistogramAndFit(treeTracker, "zBosonMass", "Z (Global&Tracker)", "zBosonsTracker.png", "Tracker Z mass", "fit_results.txt");
        // Efficiency logic as above
    }

    // Draw efficiency histograms
    canvasPt->cd();
    histPt->Draw();
    legendPt->AddEntry(histPt, "Muon Efficiency", "l");
    legendPt->Draw();
    canvasPt->SaveAs("muon_efficiency_pt.png");

    canvasEta->cd();
    histEta->Draw();
    legendEta->AddEntry(histEta, "Muon Efficiency", "l");
    legendEta->Draw();
    canvasEta->SaveAs("muon_efficiency_eta.png");

    canvasIso->cd();
    histIso->Draw();
    legendIso->AddEntry(histIso, "Muon Efficiency", "l");
    legendIso->Draw();
    canvasIso->SaveAs("muon_efficiency_iso.png");

    // Cleanup
    delete histPt;
    delete histEta;
    delete histIso;
    delete canvasPt;
    delete canvasEta;
    delete canvasIso;
    delete legendPt;
    delete legendEta;
    delete legendIso;

    file->Close();
    delete file;
}