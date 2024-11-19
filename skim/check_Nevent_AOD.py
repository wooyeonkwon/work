# -*- coding: utf-8 -*-

import os
import ROOT
import argparse
from concurrent.futures import ProcessPoolExecutor

def count_events_in_file(file_path):
    """
    단일 ROOT 파일의 이벤트 수를 계산하는 함수.
    Args:
        file_path (str): ROOT 파일 경로.
    Returns:
        int: ROOT 파일의 이벤트 수.
    """
    try:
        root_file = ROOT.TFile.Open(file_path)
        if not root_file or root_file.IsZombie():
            print(f"파일을 열 수 없음: {file_path}")
            return 0

        tree = root_file.Get("Events")  # CMS AOD에서는 "Events" 트리 사용
        if not tree:
            print(f"'Events' 트리를 찾을 수 없음: {file_path}")
            return 0
        
        event_count = tree.GetEntries()
        root_file.Close()
        return event_count

    except Exception as e:
        print(f"파일 처리 중 오류 발생: {file_path} - {e}")
        return 0


def count_events_in_root_files_parallel(directory, max_workers=4):
    """
    다중 코어를 사용하여 디렉토리 내 ROOT 파일의 총 이벤트 수를 계산하는 함수.
    Args:
        directory (str): ROOT 파일이 저장된 디렉토리 경로.
        max_workers (int): 병렬로 작업할 프로세스 수.
    Returns:
        int: 디렉토리 내 모든 ROOT 파일의 총 이벤트 수.
    """
    root_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.root')]

    total_events = 0
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # 각 파일에 대해 이벤트 수 계산 작업 병렬 처리
        results = executor.map(count_events_in_file, root_files)
        total_events = sum(results)

    return total_events


if __name__ == "__main__":
    # 명령줄 인자를 처리하기 위한 argparse 설정
    parser = argparse.ArgumentParser(description="Count total events in ROOT files using multiple cores.")
    parser.add_argument("-d", "--directory", type=str, required=True, 
                        help="The directory containing ROOT files to process.")
    parser.add_argument("-w", "--workers", type=int, default=4, 
                        help="Number of parallel workers (default: 4).")

    args = parser.parse_args()

    # 입력받은 디렉토리 경로 처리
    directory_path = args.directory
    max_workers = args.workers

    if not os.path.isdir(directory_path):
        print(f"유효하지 않은 디렉토리: {directory_path}")
    else:
        total_events = count_events_in_root_files_parallel(directory_path, max_workers)
        print(f"총 이벤트 수: {total_events}")
