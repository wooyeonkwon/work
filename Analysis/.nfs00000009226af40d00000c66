#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <stdio.h>
#include <string.h>



void drawHistogram(TTree* tree, const char* leafName, int bins, double xMin, double xMax) {
    TH1F *hist = new TH1F(leafName, Form("%s Histogram;%s;Counts", leafName, leafName), bins, xMin, xMax);
    
    // 히스토그램을 채우기 위해 leaf를 읽어들임
    tree->Draw(Form("%s>>%s", leafName, leafName));

    // 히스토그램 꾸미기 및 출력
    TCanvas *canvas = new TCanvas(Form("canvas_%s", leafName), Form("Histogram of %s", leafName), 800, 600);
    hist->SetStats(0);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetLabelSize(0.04);

    canvas->cd();
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.SetTextFont(62);
    latex.DrawLatex(0.10, 0.91, "CMS");
    latex.SetTextFont(52);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.18, 0.91, "In progress");
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.60, 0.91, "#sqrt{s} = 13.6 TeV, L = 120.57/fb");

    hist->Draw();
    canvas->SaveAs(Form("%s_histogram.png", leafName));

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
            int bins[] = {100, 100, 100, 100, 100};
            double xMax[] = {200.0, 100.0, 3.0, 3.14, 1.0};

            for (int i = 0; i < 5; ++i) {
                drawHistogram(treeReco, leaves[i], bins[i], -xMax[i], xMax[i]);
            }
        }
    }

    file->Close();
    delete file;
}