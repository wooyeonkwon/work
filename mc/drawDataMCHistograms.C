#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"  // 추가

void drawHistogramsForAllMuonTypes(const char* dataFilename, const char* mcFilename, const char* muonType) {
    const char* dataBranchName = "zBosonMass";
    const char* mcBranchName = "zBosonMass";
    
    const char* histTitle;
    if (strcmp(muonType, "GlobalMuons") == 0) {
        histTitle = "Z (GlobalMuons)";
    } else {
        histTitle = Form("Z (Global & %s)", muonType);
    }
    const char* histFileName = Form("zBosons%s.png", muonType);
    const char* canvasTitle = Form("%s Z mass", muonType);
    const char* logFileName = "fit_results.txt";
    double xMin = 60;
    double xMax = 120;

    TFile* dataFile = new TFile(dataFilename, "READ");
    TFile* mcFile = new TFile(mcFilename, "READ");

    if (!dataFile || dataFile->IsZombie()) {
        std::cerr << "Error opening data file: " << dataFilename << std::endl;
        return;
    }

    if (!mcFile || mcFile->IsZombie()) {
        std::cerr << "Error opening MC file: " << mcFilename << std::endl;
        return;
    }

    TTree* dataTree = nullptr;
    TTree* mcTree = nullptr;
    dataFile->GetObject(muonType, dataTree);
    mcFile->GetObject(muonType, mcTree);

    if (!dataTree) {
        std::cerr << muonType << " tree not found in data file." << std::endl;
        return;
    }

    if (!mcTree) {
        std::cerr << muonType << " tree not found in MC file." << std::endl;
        return;
    }

    TH1F* dataHist = new TH1F("dataHist", Form("%s;M_{\\mu\\mu} (GeV);Counts /0.5GeV", histTitle), 120, 60, 120);
    TH1F* mcHist = new TH1F("mcHist", Form("%s;M_{\\mu\\mu} (GeV);Counts /0.5GeV", histTitle), 120, 60, 120);

    double zBosonMass;
    dataTree->SetBranchAddress(dataBranchName, &zBosonMass);
    Long64_t nEntries = dataTree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        dataTree->GetEntry(i);
        dataHist->Fill(zBosonMass);
    }

    mcTree->SetBranchAddress(mcBranchName, &zBosonMass);
    nEntries = mcTree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        mcTree->GetEntry(i);
        mcHist->Fill(zBosonMass);
    }

    double dataEntries = dataHist->Integral();
    double mcEntries = mcHist->Integral();

    if (mcEntries > 0) {
        double scaleFactor = dataEntries / mcEntries;
        mcHist->Scale(scaleFactor);
    }

    TCanvas* canvas = new TCanvas(canvasTitle, canvasTitle, 800, 600);
    gStyle->SetOptStat(0);

    dataHist->SetMarkerStyle(20);
    dataHist->SetMarkerColor(kBlack);
    dataHist->SetLineColor(kBlack);
    mcHist->SetFillColor(kBlue);       
    mcHist->SetLineColor(kBlue);        


    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();

    
    mcHist->Draw("HIST");               
    dataHist->Draw("E same");           

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextFont(62); // bold font for "CMS"
    latex.DrawLatex(0.14, 0.91, "CMS");
    latex.SetTextFont(42); // normal font for "Preliminary"
    latex.DrawLatex(0.20, 0.91, "#it{Preliminary}");
    latex.DrawLatex(0.75, 0.91, "#sqrt{s} = 13.6 TeV, L = 120.57/fb");
    
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(dataHist, "Data", "lep");
    legend->AddEntry(mcHist, "MC", "f");
    legend->Draw();


    canvas->cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();

    
    TH1F* ratioHist = (TH1F*)dataHist->Clone("ratioHist");
    ratioHist->Divide(mcHist);
    ratioHist->SetTitle("");
    ratioHist->GetYaxis()->SetTitle("Data / MC");
    ratioHist->GetYaxis()->SetNdivisions(505);
    ratioHist->GetYaxis()->SetTitleSize(20);
    ratioHist->GetYaxis()->SetTitleFont(43);
    ratioHist->GetYaxis()->SetTitleOffset(1.55);
    ratioHist->GetYaxis()->SetLabelFont(43);
    ratioHist->GetYaxis()->SetLabelSize(15);
    ratioHist->GetXaxis()->SetTitleSize(20);
    ratioHist->GetXaxis()->SetTitleFont(43);
    ratioHist->GetXaxis()->SetTitleOffset(4.);
    ratioHist->GetXaxis()->SetLabelFont(43);
    ratioHist->GetXaxis()->SetLabelSize(15);
    ratioHist->Draw("ep");

    
    TLine *line = new TLine(60, 1, 120, 1);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->Draw();

    canvas->SaveAs(histFileName);

    delete dataHist;
    delete mcHist;
    delete legend;
    delete pad1;
    delete pad2;
    delete canvas;

    dataFile->Close();
    mcFile->Close();
    delete dataFile;
    delete mcFile;
}


void drawDataMCHistograms(const char* dataFilename, const char* mcFilename) {
    const char* muonTypes[] = {
        "GlobalMuons", "TrackerMuons", "StandAloneMuons", "CaloMuons", 
        "PFMuons", "RPCMuons", "GEMMuons", "ME0Muons"
    };

    for (const char* muonType : muonTypes) {
        drawHistogramsForAllMuonTypes(dataFilename, mcFilename, muonType);
    }
}