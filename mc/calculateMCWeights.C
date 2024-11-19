void calculateMCWeights(const char* filename, const char* treeName) {
    // ROOT 파일 열기
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file." << std::endl;
        return;
    }

    // 트리 가져오기
    TTree* tree = (TTree*)file->Get(treeName);
    if (!tree) {
        std::cerr << "Error: Unable to find tree." << std::endl;
        return;
    }

    // 가중치 변수 설정 (weight branch 이름은 샘플에 따라 다를 수 있음)
    Float_t weight;
    tree->SetBranchAddress("eventWeight", &weight);

    // 가중치 합산 계산
    double totalWeight = 0.0;
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        totalWeight += weight; // positive + negative weights
    }

    std::cout << "Total weight (N_mc): " << totalWeight << std::endl;

    // 파일 닫기
    file->Close();
}
