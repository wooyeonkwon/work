import ROOT
from ROOT import RDataFrame
import numpy as np
import gc
from tqdm import tqdm
import concurrent.futures

# ROOT 설정
ROOT.EnableThreadSafety()
ROOT.EnableImplicitMT()

gc.collect()
ROOT.gSystem.Load("libAnalysisClasses.so")

# 파일 로드
mc_file = "../data/Analysis_Data_test.root"
data_file = "../data/Analysis_Data_test.root"

# RDataFrame 초기화
mc_df = RDataFrame("AnalysisMC/Analysis", mc_file)
data_df = RDataFrame("AnalysisMC/Analysis", data_file)

# Binning 설정
eta_bins = [0, 0.9, 1.2, 2.1, 2.4]
phi_bins = np.linspace(0, np.pi, 11)
charge_bins = [-1, 1]
fine_tune_factor = {}

# 초기 fine_tune_factor 설정
for eta_bin in range(len(eta_bins) - 1):
    for phi_bin in range(len(phi_bins) - 1):
        for charge in charge_bins:
            fine_tune_factor[(charge, eta_bin, phi_bin)] = 1.0

# Binning 함수 정의
def get_bin(value, bins):
    for i in range(len(bins) - 1):
        if bins[i] <= value < bins[i + 1]:
            return i
    return -1

# Z mass 평균 계산 함수 수정
def calculate_z_mass_stats(df, charge, eta_bin, phi_bin):
    filtered_df = df.Filter(
        f"muons.isTightReco && muons.isTightRecoZ &&"
        f"muons.charge == {charge} && "
        f"abs(muons.eta) >= {eta_bins[eta_bin]} && abs(muons.eta) < {eta_bins[eta_bin + 1]} && "
        f"abs(muons.phi) >= {phi_bins[phi_bin]} && abs(muons.phi) < {phi_bins[phi_bin + 1]}"
    )
    
    mean_mass = filtered_df.Mean("zBoson.mass[1]").GetValue()
    mean_dmass = filtered_df.Define("dMass", "zBoson.mass[1] - 91.1876").Mean("dMass").GetValue()
    
    return mean_mass, mean_dmass

# <1/pt> 평균값 계산 함수 수정
def calculate_inverse_pt_average(df, charge, eta_bin, phi_bin, is_gen=False):
    muon_type = "genMuons" if is_gen else "muons"
    # 조건 생성
    conditions = []
    if not is_gen:
        conditions.append(f"{muon_type}.isTightReco")
        conditions.append(f"{muon_type}.isTightRecoZ")
    
    conditions.append(f'All({muon_type}.charge == {charge})')
    conditions.append(f'All(abs({muon_type}.eta) >= {eta_bins[eta_bin]} && abs({muon_type}.eta) < {eta_bins[eta_bin + 1]})')
    conditions.append(f'All(abs({muon_type}.phi) >= {phi_bins[phi_bin]} && abs({muon_type}.phi) < {phi_bins[phi_bin + 1]})')
    
    # 조건을 '&&'로 결합
    filter_condition = " && ".join(conditions)
    print(f"Filter condition: {filter_condition}")
    filtered_df = df.Filter(filter_condition)  # 수정된 부분: RVec 대신 bool 반환
    avg_inv_pt = filtered_df.Define("inv_pt", f"1.0 / {muon_type}.pt").Mean("inv_pt")
    return avg_inv_pt.GetValue()

def calculate_inverse_pt_average_parallel(df, charge, eta_bin, phi_bin, is_gen=False):
    # 기존의 calculate_inverse_pt_average 함수를 호출
    return calculate_inverse_pt_average(df, charge, eta_bin, phi_bin, is_gen)

# Correction factors 계산
def calculate_correction_factors_parallel():
    correction_factor_Cs = {}
    correction_factor_D_ms = {}
    correction_factor_D_as = {}
    correction_factor_As = {}
    correction_factor_Ms = {}
    data_inv_pt_avg = {}

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = {}
        for eta_bin in range(len(eta_bins) - 1):
            for phi_bin in range(len(phi_bins) - 1):
                for charge in charge_bins:
                    # 병렬로 calculate_inverse_pt_average 호출
                    futures[(charge, eta_bin, phi_bin)] = executor.submit(
                        calculate_inverse_pt_average_parallel, mc_df, charge, eta_bin, phi_bin, True
                    )

        for key, future in futures.items():
            charge, eta_bin, phi_bin = key
            mc_inv_pt_avg = future.result()
            data_inv_pt_avg[key] = fine_tune_factor[key] * calculate_inverse_pt_average(data_df, charge, eta_bin, phi_bin)
            correction_factor_Cs[key] = mc_inv_pt_avg - data_inv_pt_avg[key]

        for eta_bin in range(len(eta_bins) - 1):
            for phi_bin in range(len(phi_bins) - 1):
                correction_factor_D_ms[(eta_bin, phi_bin)] = (
                    correction_factor_Cs[(1, eta_bin, phi_bin)] + correction_factor_Cs[(-1, eta_bin, phi_bin)]) / 2
                
                correction_factor_D_as[(eta_bin, phi_bin)] = (
                    correction_factor_Cs[(1, eta_bin, phi_bin)] - correction_factor_Cs[(-1, eta_bin, phi_bin)]) / 2
                
                correction_factor_As[(eta_bin, phi_bin)] = (
                    correction_factor_D_as[(eta_bin, phi_bin)] -
                    ((correction_factor_D_ms[(eta_bin, phi_bin)] * (data_inv_pt_avg[(1, eta_bin, phi_bin)] - data_inv_pt_avg[(-1, eta_bin, phi_bin)])) /
                     (data_inv_pt_avg[(1, eta_bin, phi_bin)] + data_inv_pt_avg[(-1, eta_bin, phi_bin)]))
                )
                
                correction_factor_Ms[(eta_bin, phi_bin)] = 1 + (
                    2 * correction_factor_D_ms[(eta_bin, phi_bin)] /
                    (data_inv_pt_avg[(1, eta_bin, phi_bin)] + data_inv_pt_avg[(-1, eta_bin, phi_bin)])
                )

    return correction_factor_As, correction_factor_Ms


# Fine-tuning factor 업데이트
def update_fine_tune_factors(df):
    new_factors = {}
    for eta_bin in range(len(eta_bins) - 1):
        for phi_bin in range(len(phi_bins) - 1):
            for charge in charge_bins:
                zmass, dzmass = calculate_z_mass_stats(df, charge, eta_bin, phi_bin)
                if zmass != 0:  # 0으로 나누기 방지
                    new_factors[(charge, eta_bin, phi_bin)] = 1 + (2 * dzmass / zmass)
                else:
                    new_factors[(charge, eta_bin, phi_bin)] = 1.0
    return new_factors


# Fine-tuning 루프
max_iterations = 10
convergence_threshold = 0.001
iteration = 0
convergence = False

try :
    while not convergence and iteration < max_iterations:
        print(f"Iteration {iteration + 1}")
        
        # 현재 correction factors 계산
        correction_factors = calculate_correction_factors_parallel()
        
        # PT 보정 적용 및 p4 업데이트
        data_df = data_df.Define("corrected_pt",
            """
            ROOT::RVec<double> result;
            for(size_t i = 0; i < muons.pt.size(); ++i) {
                result.push_back(apply_correction(muons.pt[i], muons.charge[i], muons.eta[i], muons.phi[i], correction_factors));
            }
            return result;
            """)
        
        # 보정된 pt로 새로운 XYZTLorentzVector 좌표 계산
        data_df = data_df.Define("new_coordinates",
            """
            struct Coordinates {
                ROOT::RVec<double> x, y, z, t;
            };
            Coordinates result;
            for(size_t i = 0; i < muons.pt.size(); ++i) {
                double pt = corrected_pt[i];
                double eta = muons.eta[i];
                double phi = muons.phi[i];
                double p = pt * std::cosh(eta);
                result.x.push_back(pt * std::cos(phi));
                result.y.push_back(pt * std::sin(phi));
                result.z.push_back(pt * std::sinh(eta));
                result.t.push_back(std::sqrt(p*p + 0.105658*0.105658));
            }
            return result;
            """)
        
        # 기존 브랜치 덮어쓰기
        data_df = data_df.Redefine("muons.pt", "corrected_pt")
        data_df = data_df.Redefine("muons.p4.fCoordinates.fX", "new_coordinates.x")
        data_df = data_df.Redefine("muons.p4.fCoordinates.fY", "new_coordinates.y")
        data_df = data_df.Redefine("muons.p4.fCoordinates.fZ", "new_coordinates.z")
        data_df = data_df.Redefine("muons.p4.fCoordinates.fT", "new_coordinates.t")
        
        # Z mass 재계산
        data_df = data_df.Define("new_z_mass",
            """
            ROOT::RVec<double> result(zBoson.mass.size());
            ROOT::Math::XYZTLorentzVector mu1(muons.p4.fCoordinates.fX[0], 
                                            muons.p4.fCoordinates.fY[0],
                                            muons.p4.fCoordinates.fZ[0],
                                            muons.p4.fCoordinates.fT[0]);
            ROOT::Math::XYZTLorentzVector mu2(muons.p4.fCoordinates.fX[1],
                                            muons.p4.fCoordinates.fY[1],
                                            muons.p4.fCoordinates.fZ[1],
                                            muons.p4.fCoordinates.fT[1]);
            result[1] = (mu1 + mu2).M();
            return result;
            """)
        data_df = data_df.Redefine("zBoson.mass", "new_z_mass")
        
        # Fine-tune factors 업데이트
        previous_factors = fine_tune_factor.copy()
        fine_tune_factor = update_fine_tune_factors(data_df)
        
        # 수렴 확인
        max_diff = 0.0
        for key in fine_tune_factor:
            diff = abs(fine_tune_factor[key] - previous_factors[key])
            max_diff = max(max_diff, diff)
        
        convergence = max_diff < convergence_threshold
        iteration += 1
        
        print(f"Maximum difference: {max_diff}")
except KeyboardInterrupt :
    print("교정을 멈추고 현재 상태를 저장합니다.")

# 최종 결과 저장
output_file_data = "data_corrected.root"
data_df.Snapshot("Analysis/Analysis", output_file_data)
print(f"Data 교정 파일 저장: {output_file_data}")