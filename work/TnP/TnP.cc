
#include "TnPClasses.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// 전역 뮤텍스: 로그 기록 등 동기화에 사용 (필요시)
std::mutex logMutex;

// 두 뮤온의 invariant mass를 계산하는 함수 (muon mass = 0.105 GeV)
double calculate_z_mass(double pt1, double eta1, double phi1,
                        double pt2, double eta2, double phi2) {
    TLorentzVector vec1, vec2;
    vec1.SetPtEtaPhiM(pt1, eta1, phi1, 0.105);
    vec2.SetPtEtaPhiM(pt2, eta2, phi2, 0.105);
    return (vec1 + vec2).M();
}

int main(int argc, char** argv) {
    // 인자 체크: <root_file> [output_directory]
    ROOT::EnableThreadSafety();
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <root_file> [output_directory]" << std::endl;
        return 1;
    }
    std::string rootFileName = argv[1];
    std::string outputDir = ".";
    if (argc > 2) {
        outputDir = argv[2];
        gSystem->Exec(Form("mkdir -p %s", outputDir.c_str()));
    }

    gStyle->SetOptStat(0);
    gROOT->SetBatch(true);

    // 먼저 파일을 열어 전체 엔트리 수를 확인 (트리는 "Analysis/Analysis"에 존재한다고 가정)
    TFile* f = TFile::Open(rootFileName.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open file " << rootFileName << std::endl;
        return 1;
    }
    TDirectory* dir = (TDirectory*)f->Get("Analysis");
    if (!dir) {
        std::cerr << "Directory 'AnalysisMC' not found in file " << rootFileName << std::endl;
        return 1;
    }
    TTree* tree = (TTree*)dir->Get("Analysis");
    if (!tree) {
        std::cerr << "Tree 'Analysis' not found in directory 'AnalysisMC'" << std::endl;
        return 1;
    }
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Total entries: " << nEntries << std::endl;
    f->Close();

    // 사용 가능한 스레드 수 결정 (OpenMP 사용)
    int nThreads = 1;
#ifdef _OPENMP
    nThreads = omp_get_max_threads();
#endif

    #pragma omp parallel num_threads(nThreads)
    {
        int threadID = 0;
#ifdef _OPENMP
        threadID = omp_get_thread_num();
#endif

        // 각 스레드는 독립적으로 ROOT 파일과 트리를 연다.
        TFile* threadFile = TFile::Open(rootFileName.c_str(), "READ");
        if (!threadFile || threadFile->IsZombie()) {
            #pragma omp critical
            std::cerr << "Thread " << threadID << " cannot open file " << rootFileName << std::endl;
        } else {
            TDirectory* threadDir = (TDirectory*)threadFile->Get("Analysis");
            if (!threadDir) {
                #pragma omp critical
                std::cerr << "Thread " << threadID << " cannot find directory 'Analysis'" << std::endl;
            } else {
                TTree* threadTree = (TTree*)threadDir->Get("Analysis");
                if (!threadTree) {
                    #pragma omp critical
                    std::cerr << "Thread " << threadID << " cannot find tree 'Analysis'" << std::endl;
                } else {
                    // 스레드별 처리할 엔트리 범위 계산
                    Long64_t totalEntries = threadTree->GetEntries();
                    Long64_t eventsPerThread = totalEntries / nThreads;
                    Long64_t start = threadID * eventsPerThread;
                    Long64_t end = (threadID == nThreads - 1) ? totalEntries : start + eventsPerThread;

                    // TTreeReader를 사용하여 트리 읽기 준비
                    TTreeReader reader(threadTree);
                    
                    // Reconstructed Muon 정보
                    TTreeReaderArray<int> muon_charge(reader, "muon_charge");
                    TTreeReaderArray<double> muon_pt(reader, "muon_pt");
                    TTreeReaderArray<double> muon_eta(reader, "muon_eta");
                    TTreeReaderArray<double> muon_phi(reader, "muon_phi");
                    TTreeReaderArray<double> muon_iso(reader, "muon_iso");
                    TTreeReaderArray<double> muon_vz(reader, "muon_vz");
                    
                    // Reconstructed Muon Flags (vector<bool>는 TTreeReaderValue로 읽어 딕셔너리 문제 회피)
                    TTreeReaderValue<std::vector<bool>> muon_isGlobal(reader, "muon_isGlobal");
                    TTreeReaderValue<std::vector<bool>> muon_isTracker(reader, "muon_isTracker");
                    TTreeReaderValue<std::vector<bool>> muon_isRPC(reader, "muon_isRPC");
                    TTreeReaderValue<std::vector<bool>> muon_isGEM(reader, "muon_isGEM");
                    TTreeReaderValue<std::vector<bool>> muon_isTight(reader, "muon_isTight");
                    

                    // 출력 파일 및 TTree 생성 (스레드별 별도 파일)
                    std::string outFileName = outputDir + "/Analysis_TnP_thread" + std::to_string(threadID) + ".root";
                    TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");
                    TTree* outTree = new TTree("TnPTree", "Tag and Probe Tree");

                    std::vector<double>* zMassRPC = new std::vector<double>();
                    std::vector<double>* zDvzRPC  = new std::vector<double>();
                    std::vector<double>* zMassGEM = new std::vector<double>();
                    std::vector<double>* zDvzGEM  = new std::vector<double>();
                    std::vector<bool>* muon_isRPCPassProbe = new std::vector<bool>();
                    std::vector<bool>* muon_isGEMPassProbe = new std::vector<bool>();

                    std::vector<int>* probeMuon_charge = new std::vector<int>();
                    std::vector<double>* probeMuon_pt  = new std::vector<double>();
                    std::vector<double>* probeMuon_eta = new std::vector<double>();
                    std::vector<double>* probeMuon_phi = new std::vector<double>();
                    std::vector<double>* probeMuon_iso = new std::vector<double>();
                    std::vector<double>* probeMuon_vz = new std::vector<double>();

                    // Branch 생성 시, 두 번째 인자로 타입 이름을 명시하지 않고 포인터 주소만 넘김
                    outTree->Branch("zMassRPC", &zMassRPC);
                    outTree->Branch("zDvzRPC",  &zDvzRPC);
                    outTree->Branch("zMassGEM", &zMassGEM);
                    outTree->Branch("zDvzGEM",  &zDvzGEM);
                    outTree->Branch("muon_isRPCPassProbe", &muon_isRPCPassProbe);
                    outTree->Branch("muon_isGEMPassProbe", &muon_isGEMPassProbe);

                    outTree->Branch("probeMuon_charge", &probeMuon_charge);
                    outTree->Branch("probeMuon_pt",  &probeMuon_pt);
                    outTree->Branch("probeMuon_eta", &probeMuon_eta);
                    outTree->Branch("probeMuon_phi", &probeMuon_phi);
                    outTree->Branch("probeMuon_iso", &probeMuon_iso);
                    outTree->Branch("probeMuon_vz", &probeMuon_vz);

                    Long64_t entry = 0;
                    while (reader.Next()) {
                        if (entry < start) {
                            ++entry;
                            continue;
                        }
                        if (entry >= end)
                            break;

                        zMassRPC->clear();
                        zDvzRPC->clear();
                        zMassGEM->clear();
                        zDvzGEM->clear();
                        muon_isRPCPassProbe->clear();
                        muon_isGEMPassProbe->clear();
                        probeMuon_charge->clear();
                        probeMuon_pt->clear();
                        probeMuon_eta->clear();
                        probeMuon_phi->clear();
                        probeMuon_iso->clear();
                        probeMuon_vz->clear();
                        // vector<bool> 브랜치는 객체에서 참조를 얻음
                        const std::vector<bool>& isGlobal = *muon_isGlobal;
                        const std::vector<bool>& isTracker = *muon_isTracker;
                        const std::vector<bool>& isRPC = *muon_isRPC;
                        const std::vector<bool>& isGEM = *muon_isGEM;
                        const std::vector<bool>& isTight = *muon_isTight;

                        
                        size_t nMuons = muon_pt.GetSize();

                        // 태그 선택: pt > 24, |eta| < 2.4, isTight 및 isGlobal 조건 만족
                        std::vector<size_t> tagIndices;
                        for (size_t i = 0; i < nMuons; i++) {
                            if (muon_pt[i] > 24.0 &&
                                std::fabs(muon_eta[i]) < 2.4 &&
                                isTight[i] &&
                                isGlobal[i])
                            {
                                tagIndices.push_back(i);
                            }
                        }

                        // 프로브 선택: 태그와 다른 인덱스, 반대 부호, 그리고 (pt > 15, |eta| < 2.4, isTracker 조건)
                        for (auto tag : tagIndices) {
                            for (size_t probe = 0; probe < nMuons; probe++) {
                                if (tag == probe)
                                    continue;
                                if (muon_charge[tag] == muon_charge[probe])
                                    continue;
                                if (muon_pt[probe] < 15.0 ||
                                    std::fabs(muon_eta[probe]) >= 2.4 ||
                                    !isTracker[probe])
                                    continue;

                                // 태그와 프로브를 이용해 Z mass 및 dvz 계산
                                double mass = calculate_z_mass(muon_pt[tag], muon_eta[tag], muon_phi[tag],
                                                               muon_pt[probe], muon_eta[probe], muon_phi[probe]);
                                double dvz = std::fabs(muon_vz[tag] - muon_vz[probe]);

                                if (mass >= 60.0 && mass <= 120.0 && dvz < 0.5) {

                                    zMassRPC->push_back(mass);
                                    zDvzRPC->push_back(dvz);
                                    zMassGEM->push_back(mass);
                                    zDvzGEM->push_back(dvz);
                                    muon_isRPCPassProbe->push_back(isRPC[probe]);
                                    muon_isGEMPassProbe->push_back(isGEM[probe]);
                                    probeMuon_charge->push_back(muon_charge[probe]);
                                    probeMuon_pt->push_back(muon_pt[probe]);
                                    probeMuon_eta->push_back(muon_eta[probe]);
                                    probeMuon_phi->push_back(muon_phi[probe]);
                                    probeMuon_iso->push_back(muon_iso[probe]);
                                    probeMuon_vz->push_back(muon_vz[probe]);
                                    outTree->Fill();
                                }
                            }
                        }
                        ++entry;
                    } // end while

                    outFile->cd();
                    outTree->Write();
                    outFile->Close();
                    delete outFile;
                }
            }
            threadFile->Close();
        }
    } // end omp parallel

    return 0;
}
