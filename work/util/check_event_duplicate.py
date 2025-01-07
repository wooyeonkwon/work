import ROOT
from ROOT import RDataFrame
import sys
from tqdm import tqdm
import numpy as np
ROOT.EnableImplicitMT(56)

def check_duplicate_eventnumbers_rdf(file_list):
    # ROOT 파일 열기
    for file_path in file_list:
        if not file_path.endswith(".root"):
            print(f"Skipping non-ROOT file: {file_path}")
            continue

        root_file = ROOT.TFile.Open(file_path)

        if not root_file or root_file.IsZombie():
            print(f"Error: Unable to open file {file_path}")
            continue

        print(f"Checking file: {file_path}")
        corrupted = False

        rdf = RDataFrame("Analysis/Analysis", file_path )

        # eventnumber 필드 읽기
        data = rdf.AsNumpy(["eventNumber", "runNumber"])
        event_numbers = data["eventNumber"]
        run_numbers = data["runNumber"]

        # 중복 확인
        event_run_pairs = {}
        duplicate_event_run_pairs = []
        try :
            for event_number, run_number in tqdm(zip(event_numbers, run_numbers)):
                if (event_number, run_number) in event_run_pairs:
                    duplicate_event_run_pairs.append((event_number, run_number))
                else:
                    event_run_pairs[(event_number, run_number)] = True
        except KeyboardInterrupt :
            print("이벤트 체크를 멈추고 진행상황을 보고합니다.")

        if duplicate_event_run_pairs:
            print("중복된 eventnumber가 발견되었습니다:")
            with open("duplicate_event.txt", "w") as file:
                for event_number, run_number in duplicate_event_run_pairs:
                    file.write(f"eventNumber: {event_number}, runNumber: {run_number}\n")
        else:
            print("중복된 eventnumber가 없습니다.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_root_files.py <file1.root> <file2.root> ...")
        sys.exit(1)

    file_list = sys.argv[1:]
    check_duplicate_eventnumbers_rdf(file_list)