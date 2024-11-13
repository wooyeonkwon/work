#!/bin/bash

directory_path="/pnfs/knu.ac.kr/data/cms/store/user/wkwon/AnalysisResults/Muon/crab_MuonSkimming_Run2022D/241108_143535/0000/"

root_files=($(find "$directory_path" -type f -name "*.root"))

num_groups=16
files_per_group=$(( ${#root_files[@]} / num_groups ))

for i in $(seq 0 $((num_groups - 1)))
do
    
    start_index=$(( i * files_per_group ))
    end_index=$(( start_index + files_per_group ))

    sublist=("${root_files[@]:$start_index:$files_per_group}")
    sublist_file="list${i}.txt"
    output_file="output_data_${i}.root"

   
    printf "%s\n" "${sublist[@]}" > "$sublist_file"

    echo "Starting cmsRun for list${i}.txt with output file output_data_${i}.root"
    
    (cmsRun conf_FMT.py "$sublist_file" "$output_file" > "output_${i}.log" 2>&1) &
done


wait

echo "All cmsRun conf_FMT.py works have finished"
