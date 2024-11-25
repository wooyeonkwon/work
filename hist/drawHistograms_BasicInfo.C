#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <stdio.h>
#include <string.h>


void drawHistogram(TTree* tree, const char* leafName, int bins, double xMin, double xMax, const char* histTitle, const char* xUnit) {
    double binWidth = (xMax - xMin) / bins;
    const double tolerance = 1e-6;
    TH1D *hist = nullptr;

    if (fabs(binWidth - 1.0) < tolerance){
        hist = new TH1D(leafName, Form("%s;%s;Multiplicity %s", histTitle, xUnit, xUnit), bins, xMin, xMax);
    }
    else {
        hist = new TH1D(leafName, Form("%s;%s;Multiplicity/%g\n %s", histTitle, xUnit, binWidth, xUnit), bins, xMin, xMax);
    }
    
    // 히스토그램을 채우기 위해 leaf를 읽어들임
    tree->Draw(Form("%s>>%s", leafName, leafName));

    // 히스토그램 꾸미기 및 출력
    TCanvas *canvas = new TCanvas(Form("canvas_%s", leafName), Form("Histogram of %s", leafName), 800, 600);
    hist->SetStats(0);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetTitleSize(0.04);
    hist->GetYaxis()->SetTitleOffset(0.8);
    hist->GetYaxis()->SetLabelSize(0.03);
    hist->GetYaxis()->SetMaxDigits(2);
    hist->Draw();

    canvas->cd();
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.SetTextFont(62);
    latex.DrawLatex(0.17, 0.91, "CMS");
    latex.SetTextFont(52);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.25, 0.91, "In progress");
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.60, 0.91, "#sqrt{s} = 13.6 TeV, 11.06/fb");  //full data : L = 120.56 / 11.06/fb Nevent : 72909628


    canvas->SaveAs(Form("mc_%s_histogram.png", leafName));

    // 리소스 정리
    delete canvas;
    delete hist;
}


void drawHistograms_BasicInfo(const char* filename) {
    TFile *file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        printf("Error opening file!\n");
        return;
    }

    TTree *treeReco;
    TIter nextKey(file->GetListOfKeys());
    TKey *key;

    while ((key = (TKey*)nextKey())) {
        if (strcmp(key->GetClassName(), "TTree") == 0 && strstr(key->GetName(), "RecoMuons") != NULL) {
            treeReco = (TTree*)key->ReadObj();
            if (!treeReco) {
                printf("Error loading tree '%s'!\n", key->GetName());
                continue;
            }

            // 히스토그램을 그릴 leaf 이름과 범위 설정
            const char* leaves[] = {"muonSize", "muonPt", "muonEta", "muonPhi", "muonIso"};
            int bins[] = {10, 2000, 60, 70, 20};
            double xMin[] = {0, 0.0, -3.0, -3.5, 0.0};
            double xMax[] = {10, 200.0, 3.0, 3.5, 0.2};
            const char* histTitle[] = {"Muon Mutiplicity","Muon Pt", "Muon Eta", "Muon Phi", "Muon PFIso(dR = 0.4)"};
            const char* xUnit[] = {"","GeV","","",""};

            for (int i = 0; i < 5; ++i) {
                drawHistogram(treeReco, leaves[i], bins[i], xMin[i], xMax[i], histTitle[i], xUnit[i]);
            }
        }
    }

    file->Close();
    delete file;
}
