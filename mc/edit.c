#include <TFile.h>
#include <TTree.h>

int edit() {
    // mc_data.root 파일 열기
    TFile *file = TFile::Open("mc_data.root", "READ");
    TTree *treeReco = (TTree*)file->Get("treeReco");
    TTree *treeGen = (TTree*)file->Get("treeGen");

    // zBosonMassReco와 zBosonMassGen 변수 선언
    Double_t zBosonMassReco;
    Double_t zBosonMassGen;

    // 트리의 브랜치 설정
    treeReco->SetBranchAddress("zBosonMassReco", &zBosonMassReco);
    treeGen->SetBranchAddress("zBosonMassGen", &zBosonMassGen);

    // 새로운 ROOT 파일 생성
    TFile *new_file = new TFile("mc_data_modified.root", "RECREATE");
    TTree *new_treeReco = treeReco->CloneTree(0);
    TTree *new_treeGen = treeGen->CloneTree(0);

    // 조건에 따라 엔트리 복사
    Long64_t nentriesReco = treeReco->GetEntries();
    for (Long64_t i = 0; i < nentriesReco; i++) {
        treeReco->GetEntry(i);
        if (zBosonMassReco >= 2) {
            new_treeReco->Fill();
        }
    }

    Long64_t nentriesGen = treeGen->GetEntries();
    for (Long64_t i = 0; i < nentriesGen; i++) {
        treeGen->GetEntry(i);
        if (zBosonMassGen >= 2) {
            new_treeGen->Fill();
        }
    }

    // 새로운 파일에 트리 저장
    new_file->Write();
    new_file->Close();
    file->Close();

    return 0;
}
