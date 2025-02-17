import os, glob
import ROOT
from tqdm import tqdm

# 1. ROOT의 ImplicitMT를 활성화하여 다중 스레드를 사용 (기본 스레드 수는 CPU 코어 수에 따라 자동 결정)
ROOT.ROOT.EnableImplicitMT()

# 2. 입력: 폴더 리스트와 추출할 tree/branch 리스트
folder_list = ["/data1/users/dndus0107/raw_data_samples/Run2022C/SingleMuon/NANOAOD/16Dec2023-v1/50000/"]
branch_list = [
    "Events/Jet_nMuons", 
    "Events/Jet_muEF", 
    "Events/Jet_eta",
]
output_path = './extract_branch'

# 3. branch_list를 tree 이름을 기준으로 그룹화
tree_branches = {}
for tb in branch_list:
    try:
        tree, branch = tb.split("/")
    except ValueError:
        print(f"입력 형식 오류: '{tb}'는 'tree/branch' 형식이어야 합니다.")
        continue
    tree_branches.setdefault(tree, []).append(branch)

# 4. 모든 지정 폴더에서 *.root 파일 리스트 생성
root_files = []
for folder in folder_list:
    root_files.extend(glob.glob(os.path.join(folder, "*.root")))

if not root_files:
    print("입력 폴더에 ROOT 파일이 없습니다.")
    exit(1)

# 5. 각 tree 별로 RDataFrame을 생성하고, 선택한 브랜치만을 Snapshot으로 저장
for tree, branches in tree_branches.items():
    print(f"Processing tree '{tree}' with branches: {branches}")
    for file in tqdm(root_files, desc=f"Processing {tree} in ROOT files"):
         # 각 파일에 대해 RDataFrame 생성
         df = ROOT.RDataFrame(tree, file)
         base = os.path.basename(file)
         out_file = f"output_{tree}_{base}"
         # 지정된 브랜치만 포함하는 Snapshot 저장
         df.Snapshot(tree, out_file, branches)
         # 개별 파일 처리 완료 메시지 (선택 사항)
         # print(f"Processed {file} -> {out_file}")