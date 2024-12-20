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
#include "TLine.h"


void drawHistogram(TTree* dataTree, TTree* mcTree, const char* treeName, const char* leafName, int bins, double xMin, double xMax, const char* histTitle, const char* xUnit, const char* yUnit) {
    
    dataTree->Draw(Form("%s>>dataHist(%d,%g,%g)", leafName, bins, xMin, xMax));
    TH1D* dataHist = (TH1D*)gDirectory->Get("dataHist");

    mcTree->Draw(Form("%s>>mcHist(%d,%g,%g)", leafName, bins, xMin, xMax));
    TH1D* mcHist = (TH1D*)gDirectory->Get("mcHist");


    

    if (mcHist->Integral() > 0) {
    double mcCrossSection = 6688.0; 
    double mcCrossSectionUnc = 83.99;
    double dataLumi = 5010.4;      
    double Nmc = 72909628;          
    double scaleFactor = mcCrossSection * dataLumi / Nmc;
    mcHist->Scale(scaleFactor);

        
        for (int bin = 1; bin <= mcHist->GetNbinsX(); ++bin) {
            double binContent = mcHist->GetBinContent(bin);
            double binStatUnc = mcHist->GetBinError(bin);
            double binScaleUnc = binContent * (mcCrossSectionUnc / mcCrossSection);
            double totalUnc = sqrt(binStatUnc * binStatUnc + binScaleUnc * binScaleUnc);
            mcHist->SetBinError(bin, totalUnc);
        }
    }


    TCanvas *canvas = new TCanvas(Form("canvas_%s", leafName), Form("Histogram of %s", leafName), 800, 600);
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
    latex.DrawLatex(0.20, 0.91, "#it{In Progress}");
    latex.DrawLatex(0.75, 0.91, "#sqrt{s} = 13.6 TeV, L = 11.06/fb");
    
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


    canvas->SaveAs(Form("dataMC_%s_%s_histogram.png",treeName, leafName));

    delete dataHist;
    delete mcHist;
    delete ratioHist;
    delete pad1;
    delete pad2;
    delete legend;
    delete line;



}


void drawDataMCHistograms_2(const char* dataFilename, const char* mcFilename) {
    TFile *datafile = TFile::Open(dataFilename, "READ");
    if (!datafile || datafile->IsZombie()) {
        printf("Error data File opening file!\n");
        return;
    }

    const char* dirName = "Analysis";
    TDirectoryFile *dataDir = (TDirectoryFile*)datafile->Get(dirName);
    if (!dataDir) {
        printf("Error: Directory '%s' not found in data file!\n", dirName);
        datafile->Close();
        delete datafile;
        return;
    }

    TFile *mcfile = TFile::Open(mcFilename, "READ");
    if (!mcfile || mcfile->IsZombie()) {
        printf("Error mc File opening file!\n");
        return;
    }

    TDirectoryFile *mcDir = (TDirectoryFile*)mcfile->Get(dirName);
    if (!mcDir) {
        printf("Error: Directory '%s' not found in mc file!\n", dirName);
        mcfile->Close();
        delete mcfile;
        return;
    }


    const char* treeNames[] = {"RPCMuons","GEMMuons", "RecoMuons"};
    const char* branches[] = {"zBosonMass", "vertexdz", "muonSize", "muonPt", "muonEta", "muonPhi", "muonIso"};
    const char* histTitles[] = {"diMuon Mass", "vertexdz", "muonSize", "muonPt", "muonEta", "muonPhi", "muonIso"};
    int bins[] = {60, 100, 5, 2000, 60, 70, 150};
    double xMin[] = {60.0, 0.0, 0.0, 0.0, -3.0, -3.5, 0};
    double xMax[] = {120.0, 0.1, 5.0, 200.0, 3.0, 3.5, 0.15};
    const char* xUnits[] = {"GeV", "cm", "", "GeV", "", "", ""};
    const char* yUnits[] = {"Events", "Events", "Multiplicity", "Multiplicity", "Multiplicity", "Multiplicity", "Multiplicity"};

    for (int i = 0; i < 3; ++i) {

        TTree *dataTree = (TTree*)dataDir->Get(treeNames[i]);
        if (!dataTree) {
            printf("Tree '%s' not found in data directory '%s'!\n", treeNames[i], dirName);
            continue;
        }
        TTree *mcTree = (TTree*)mcDir->Get(treeNames[i]);
        if (!mcTree) {
            printf("Tree '%s' not found in mc directory '%s'!\n", treeNames[i], dirName);
            continue;
        }

        for (int j = 0; j < 7; ++j ){
            drawHistogram(dataTree, mcTree, treeNames[i], branches[j], bins[j], xMin[j], xMax[j], histTitles[j], xUnits[j], yUnits[j]);
        }
        

    }

    datafile->Close();
    mcfile->Close();
    delete datafile;
    delete mcfile;
}
