import ROOT
from ROOT import RDataFrame
import numpy as np
import gc
from tqdm import tqdm
# ROOT 설정
ROOT.EnableThreadSafety()
ROOT.EnableImplicitMT(48)

gc.collect()
ROOT.gSystem.Load("libAnalysisClasses.so")

# 파일 로드
mc_file = "mc.root"
data_file = "data.root"

# RDataFrame 초기화
mc_df = RDataFrame("Analysis/Analysis", mc_file)
data_df = RDataFrame("Analysis/Analysis", data_file)

# Binning 설정
eta_bins = [0, 0.9, 1.2, 2.1, 2.4]
phi_bins = np.linspace(0, np.pi, 11)
charge_bins = [-1, 1]
fine_tune_factor = {}

# Fine-tuning 루프 변수
max_iterations = 10
convergence_threshold = 0.001
iteration = 0
convergence = False

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

# Z mass 평균 계산 함수
def calculate_z_mass_stats(df, charge, eta_bin, phi_bin):
    filtered_df = df.Filter(
        f"muons.isTightReco == 1 && muons.isTightRecoZ == 1 &&"
        f"muons.charge == {charge} && "
        f"abs(muons.eta) >= {eta_bins[eta_bin]} && abs(muons.eta) < {eta_bins[eta_bin + 1]} && "
        f"abs(muons.phi) >= {phi_bins[phi_bin]} && abs(muons.phi) < {phi_bins[phi_bin + 1]}"
    )
    
    mean_mass = filtered_df.Mean("zBoson.mass").GetValue()
    mean_dmass = filtered_df.Define("dMass", "zBoson.mass - 91.1876").Mean("dMass").GetValue()
    
    return mean_mass, mean_dmass

# <1/pt> 평균값 계산
def calculate_inverse_pt_average(df, charge, eta_bin, phi_bin):
    filtered_df = df.Filter(
        f"muons.isTightReco == 1 && muons.isTightRecoZ == 1 &&"
        f"muons.charge == {charge} && "
        f"abs(muons.eta) >= {eta_bins[eta_bin]} && abs(muons.eta) < {eta_bins[eta_bin + 1]} && "
        f"abs(muons.phi) >= {phi_bins[phi_bin]} && abs(muons.phi) < {phi_bins[phi_bin + 1]}"
    )
    avg_inv_pt = filtered_df.Define("inv_pt", "1.0 / muons.pt").Mean("inv_pt")
    return avg_inv_pt.GetValue()

# Correction factors 계산
def calculate_correction_factors():
    correction_factor_Cs = {}
    correction_factor_D_ms = {}
    correction_factor_D_as = {}
    correction_factor_As = {}
    correction_factor_Ms = {}
    data_inv_pt_avg = {}

    for eta_bin in tqdm(range(len(eta_bins) - 1)):
        for phi_bin in range(len(phi_bins) - 1):
            for charge in charge_bins:
                # fine-tune factor 적용
                mc_inv_pt_avg = calculate_inverse_pt_average(mc_df, charge, eta_bin, phi_bin)
                data_inv_pt_avg[(charge, eta_bin, phi_bin)] = fine_tune_factor[(charge, eta_bin, phi_bin)] * calculate_inverse_pt_average(data_df, charge, eta_bin, phi_bin)
                correction_factor_Cs[(charge, eta_bin, phi_bin)] = mc_inv_pt_avg - data_inv_pt_avg[(charge, eta_bin, phi_bin)]

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


# Correction 적용 함수
def apply_correction(pt, charge, eta, phi, correction_factors):
    eta_bin = get_bin(abs(eta), eta_bins)
    phi_bin = get_bin(abs(phi), phi_bins)
    if eta_bin == -1 or phi_bin == -1:
        return pt
    
    As = correction_factors[0].get((eta_bin, phi_bin), 0.0)
    Ms = correction_factors[1].get((eta_bin, phi_bin), 1.0)
    
    return 1 / (((1 / pt) * Ms) + (charge * As))

# 4-벡터 업데이트
def update_p4(muons, corrected_pts):
    from ROOT import Math
    new_p4s = []
    for muon, pt in zip(muons, corrected_pts):
        p4 = muon.p4()
        new_p4 = Math.PtEtaPhiMVector(pt, p4.Eta(), p4.Phi(), p4.M())
        new_p4s.append(new_p4)
    return new_p4s

# Z 보존 재구성
def reconstruct_Z(df):
    reconstructed_Z_mass = (
        df.Filter("Sum(muons.isTightRecoZ == 1) == 2")
        .Define("Z_mass", "ROOT::Math::InvariantMass(muons[0].p4() + muons[1].p4())")
        .Mean("Z_mass")
    )
    return reconstructed_Z_mass.GetValue()

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
while not convergence and iteration < max_iterations:
    print(f"Iteration {iteration + 1}")
    
    # 현재 correction factors 계산
    correction_factors = calculate_correction_factors()
    
    # PT 보정 적용
    data_df = data_df.ReDefine("muons.pt",
        "apply_correction(muons.pt, muons.charge, muons.eta, muons.phi, correction_factors)")
    
    # Z boson 재구성
    data_df = data_df.Redefine("muons.p4",
        "update_p4(muons, corrected_pt)")
    
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

# 최종 결과 저장
output_file_data = "data_corrected.root"
data_df.Snapshot("Analysis/Analysis", output_file_data)
print(f"Data 교정 파일 저장: {output_file_data}")