import ROOT

# 원본 ROOT 파일 열기
input_file = ROOT.TFile.Open("/data1/users/dndus0107/AnalysisResults/processed_data/Analysis_Data_22EFG.root", "READ")
input_tree = input_file.Get("Analysis/Analysis")  # "TreeName"을 원본 트리 이름으로 바꾸세요

# 출력 ROOT 파일 생성
output_file = ROOT.TFile.Open("/data1/users/dndus0107/AnalysisResults/processed_data/mini_Analysis_Data_22EFG.root", "RECREATE")
output_tree = input_tree.CloneTree(0)  # 빈 트리 생성

# 엔트리 범위 설정
start_entry = 0
end_entry = 10000

# 엔트리 복사
for i in range(start_entry, end_entry):
    input_tree.GetEntry(i)
    output_tree.Fill()

# 출력 파일 저장 및 닫기
output_file.Write()
output_file.Close()

# 입력 파일 닫기
input_file.Close()

print(f"엔트리 {start_entry}부터 {end_entry}까지 추출 완료!")