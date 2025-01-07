#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <stdio.h>
#include <string.h>


void drawHistogram(TTree* tree,const char* treeName, const char* leafName, int bins, double xMin, double xMax, const char* histTitle, const char* xUnit, const char* yUnit) {
    double binWidth = (xMax - xMin) / bins;

    TH1D *hist = nullptr;

    hist = new TH1D(leafName, Form("%s;%s;%s/%g\n %s", histTitle, xUnit, yUnit, binWidth, xUnit), bins, xMin, xMax);
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
    latex.DrawLatex(0.60, 0.91, "#sqrt{s} = 13.6 TeV");  //full data : L = 120.56 / 11.06/fb Nevent : 72909628

    int bin1 = hist->FindBin(60);
    int bin2 = hist->FindBin(120);
    double entries = hist->Integral(bin1,bin2); 
    std::cout << "Number of entries: " << entries << std::endl;


    canvas->SaveAs(Form("%s_%s_histogram.png",treeName, leafName));

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

    const char* dirName = "Analysis";
    TDirectoryFile *dir = (TDirectoryFile*)file->Get(dirName);
    if (!dir) {
        printf("Error: Directory '%s' not found in file!\n", dirName);
        file->Close();
        delete file;
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

        TTree *tree = (TTree*)dir->Get(treeNames[i]);
        if (!tree) {
            printf("Tree '%s' not found in directory '%s'!\n", treeNames[i], dirName);
            continue;
        }
        for (int j = 0; j < 7; ++j ){
            drawHistogram(tree, treeNames[i], branches[j], bins[j], xMin[j], xMax[j], histTitles[j], xUnits[j], yUnits[j]);
        }
        // drawHistogram 함수 호출

    }

    file->Close();
    delete file;
}

