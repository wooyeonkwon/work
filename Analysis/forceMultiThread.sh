#!/bin/bash

#directory setting
directory_path="/pnfs/knu.ac.kr/data/cms/store/user/wkwon/AnalysisResults/Muon/crab_MuonSkimming_Run2022D/241108_143535/0000/"

root_files=($(find "$directory_path" -type f -name "*.root"))

#multithread setting
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


# Wait for all background processes to finish
wait

# Check each log file for errors
for i in $(seq 0 $((num_groups - 1)))
do
    if grep -q "Fatal Exception" "output_$i.log" || grep -q "Traceback" "output_$i.log"; then
        echo "Warning: Errors detected in output_$i.log. Please check this log for details."
    else
        echo "cmsRun for list${i}.txt completed successfully."
    fi
done

echo "All cmsRun conf_FMT.py processes have finished."
