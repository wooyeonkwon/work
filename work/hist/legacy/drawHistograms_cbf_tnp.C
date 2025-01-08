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
    latex.DrawLatex(0.10, 0.92, "CMS");
    latex.SetTextFont(42); // normal font for "Preliminary"
    latex.DrawLatex(0.16, 0.92, "#it{Preliminary}");
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
void drawEfficiency(TTree* tree, const char* treeName, TH1F* histPt, TH1F* histEta, TH1F* histIso, int color) {
    double probePt, probeEta, probeIso;
    Bool_t passingProbe;

    tree->SetBranchAddress("probePt", &probePt);
    tree->SetBranchAddress("probeEta", &probeEta);
    tree->SetBranchAddress("probeIso", &probeIso);
    tree->SetBranchAddress("passingProbe", &passingProbe);

    // Temporary histograms to store total and passing events
    TH1F* histPtAll = (TH1F*)histPt->Clone(Form("%s_PtAll", treeName));
    TH1F* histPtPass = (TH1F*)histPt->Clone(Form("%s_PtPass", treeName));
    TH1F* histEtaAll = (TH1F*)histEta->Clone(Form("%s_EtaAll", treeName));
    TH1F* histEtaPass = (TH1F*)histEta->Clone(Form("%s_EtaPass", treeName));
    TH1F* histIsoAll = (TH1F*)histIso->Clone(Form("%s_IsoAll", treeName));
    TH1F* histIsoPass = (TH1F*)histIso->Clone(Form("%s_IsoPass", treeName));

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Fill the histograms with all events
        histPtAll->Fill(probePt);
        histEtaAll->Fill(probeEta);
        histIsoAll->Fill(probeIso);

        // Fill the passing histograms only if the probe passed
        if (passingProbe) {
            histPtPass->Fill(probePt);
            histEtaPass->Fill(probeEta);
            histIsoPass->Fill(probeIso);
        }
    }

    // Calculate efficiency
    TGraphAsymmErrors* grPt = new TGraphAsymmErrors(histPtPass, histPtAll);
    TGraphAsymmErrors* grEta = new TGraphAsymmErrors(histEtaPass, histEtaAll);
    TGraphAsymmErrors* grIso = new TGraphAsymmErrors(histIsoPass, histIsoAll);

    // Set graph properties
    grPt->SetLineColor(color);
    grPt->SetMarkerColor(color);
    grPt->SetMarkerStyle(20);
    grEta->SetLineColor(color);
    grEta->SetMarkerColor(color);
    grEta->SetMarkerStyle(20);
    grIso->SetLineColor(color);
    grIso->SetMarkerColor(color);
    grIso->SetMarkerStyle(20);

    // Add graphs to the main histograms
    histPt->GetListOfFunctions()->Add(grPt);
    histEta->GetListOfFunctions()->Add(grEta);
    histIso->GetListOfFunctions()->Add(grIso);

    // Clean up temporary histograms
    delete histPtAll;
    delete histPtPass;
    delete histEtaAll;
    delete histEtaPass;
    delete histIsoAll;
    delete histIsoPass;
}

void drawHistograms_cbf_tnp(const char* filename) {
    TFile* file = new TFile(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Declare tree pointers
    TTree *treeGlobal = nullptr;
    TTree *treeTracker = nullptr;
    TTree *treeStandAlone = nullptr;
    TTree *treePF = nullptr;
    TTree *treeRPC = nullptr;
    TTree *treeGEM = nullptr;

    // Define trees and get objects
    std::vector<std::tuple<TTree**, const char*, const char*, const char*, const char*>> trees = {
        {&treeGlobal, "GlobalMuons", "Z (Global)", "zBosonsGlobal.png", "Global Z mass"},
        {&treeTracker, "TrackerMuons", "Z (Global&Tracker)", "zBosonsTracker.png", "Tracker Z mass"},
        {&treeStandAlone, "StandAloneMuons", "Z (Global&StandAlone)", "zBosonsStandAlone.png", "StandAlone Z mass"},
        {&treePF, "PFMuons", "Z (Global&PF)", "zBosonsPF.png", "PF Z mass"},
        {&treeRPC, "RPCMuons", "Z (Global&RPC)", "zBosonsRPC.png", "RPC Z mass"},
        {&treeGEM, "GEMMuons", "Z (Global&GEM)", "zBosonsGEM.png", "GEM Z mass"}
    };

    for (auto& tree : trees) {
        file->GetObject(std::get<1>(tree), *std::get<0>(tree));
        if (!*std::get<0>(tree)) std::cerr << std::get<1>(tree) << " tree not found." << std::endl;
    }

    // Create histograms for each variable
    TH1F* histPt = new TH1F("histPt", ";muon p_{T} (GeV);Efficiency", 40, 0, 80);
    TH1F* histEta = new TH1F("histEta", ";muon #eta;Efficiency", 50, -2.5, 2.5);
    TH1F* histIso = new TH1F("histIso", ";muon Iso;Efficiency", 50, 0, 0.5);

    // Set histogram style
    std::vector<TH1F*> histograms = {histPt, histEta, histIso};
    for (auto hist : histograms) {
        hist->SetStats(0);
        hist->GetXaxis()->SetTitleSize(0.05);
        hist->GetXaxis()->SetTitleOffset(0.9);
        hist->GetXaxis()->SetLabelSize(0.04);
        hist->GetYaxis()->SetTitleSize(0.05);
        hist->GetYaxis()->SetTitleOffset(0.9);
        hist->GetYaxis()->SetLabelSize(0.04);
        hist->SetMinimum(0);
        hist->SetMaximum(1.1);
    }

    // Create canvases
    TCanvas* canvasPt = new TCanvas("canvas_pt", "Muon Efficiency vs pT", 800, 600);
    TCanvas* canvasEta = new TCanvas("canvas_eta", "Muon Efficiency vs eta", 800, 600);
    TCanvas* canvasIso = new TCanvas("canvas_iso", "Muon Efficiency vs Iso", 800, 600);

    // Create legends
    TLegend* legendPt = new TLegend(0.7, 0.2, 0.9, 0.4);
    TLegend* legendEta = new TLegend(0.7, 0.2, 0.9, 0.4);
    TLegend* legendIso = new TLegend(0.7, 0.2, 0.9, 0.4);

    std::vector<TLegend*> legends = {legendPt, legendEta, legendIso};
    for (auto legend : legends) {
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
    }

    // Colors for different muon types
    int colors[] = {kRed, kBlue, kGreen+2, kMagenta, kOrange-3, kCyan+2};

    // Draw efficiency and histograms for each tree
    for (size_t i = 0; i < trees.size(); ++i) {
        if (*std::get<0>(trees[i])) {
            // Draw efficiency
            drawEfficiency(*std::get<0>(trees[i]), std::get<1>(trees[i]), histPt, histEta, histIso, colors[i]);
            
            // Add entry to legends (only muon type, no "Muons" suffix)
            std::string legendLabel = std::get<1>(trees[i]);
            legendLabel = legendLabel.substr(0, legendLabel.find("Muons"));
            legendPt->AddEntry(histPt->GetListOfFunctions()->Last(), legendLabel.c_str(), "lp");
            legendEta->AddEntry(histEta->GetListOfFunctions()->Last(), legendLabel.c_str(), "lp");
            legendIso->AddEntry(histIso->GetListOfFunctions()->Last(), legendLabel.c_str(), "lp");

            // Draw histogram and fit
            //drawHistogramAndFit(*std::get<0>(trees[i]), "zBosonMass", std::get<2>(trees[i]), std::get<3>(trees[i]), std::get<4>(trees[i]), "fit_results.txt");
        }
    }
    // Function to draw CMS labels
    auto drawCMSLabel = [](TCanvas* canvas) {
        canvas->cd();
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.05);
        latex.SetTextFont(62);
        latex.DrawLatex(0.10, 0.91, "CMS");
        latex.SetTextFont(52);
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.18, 0.91, "Preliminary");
        latex.SetTextFont(42);
        latex.SetTextSize(0.035);
        latex.DrawLatex(0.60, 0.91, "#sqrt{s} = 13.6 TeV, L = 120.57/fb");
    };
    
    // Draw and save efficiency histograms
    std::vector<std::pair<TCanvas*, TH1F*>> canvasHistPairs = {
        {canvasPt, histPt}, {canvasEta, histEta}, {canvasIso, histIso}
    };

    for (auto& pair : canvasHistPairs) {
        pair.first->cd();
        pair.second->Draw();  // "HIST"로 변경하여 히스토그램 스타일로 그리기
        gPad->SetTicks(1, 1);       // 눈금 표시
        if (pair.second == histPt) {
            pair.second->GetXaxis()->SetTitle("pT (GeV)");
            legendPt->Draw();
        } else if (pair.second == histEta) {
            pair.second->GetXaxis()->SetTitle("#eta");
            legendEta->Draw();
        } else {
            pair.second->GetXaxis()->SetTitle("Isolation");
            legendIso->Draw();
        }
        drawCMSLabel(pair.first);
        pair.first->SaveAs((std::string("efficiency_") + pair.second->GetName() + "_comparison.png").c_str());
    }
    // Cleanup
    delete canvasPt;
    delete canvasEta;
    delete canvasIso;
    delete histPt;
    delete histEta;
    delete histIso;
    delete legendPt;
    delete legendEta;
    delete legendIso;

    file->Close();
    delete file;
}