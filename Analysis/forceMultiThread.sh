#!/bin/bash

# Initialize variables
use_custom_lists=false
custom_lists=()
directory_path="/pnfs/knu.ac.kr/data/cms/store/user/wkwon/AnalysisResults/Muon/crab_MuonSkimming_Run2022D/241108_143535/0000/"

# Parse command line options
while getopts "d:l:" opt; do
    case $opt in
        d)
            directory_path="$OPTARG"
            ;;
        l)
            use_custom_lists=true
            IFS=',' read -r -a custom_lists <<< "$OPTARG"
            ;;
        *)
            echo "Usage: $0 [-d directory_path] [-l list1.txt,list2.txt,...]"
            exit 1
            ;;
    esac
done

# Find all root files only if custom lists are not used
if [ "$use_custom_lists" = false ]; then
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
        echo "Starting cmsRun for $sublist_file with output file $output_file"
        (cmsRun conf_FMT.py "$sublist_file" "$output_file" > "output_${i}.log" 2>&1) &
    done

else
    num_groups=${#custom_lists[@]}

    for i in $(seq 0 $((num_groups - 1)))
    do
        sublist_file="${custom_lists[$i]}"
        output_file="output_data_${i}.root"

        if [ ! -f "$sublist_file" ]; then
            echo "Error: File $sublist_file does not exist."
            continue
        fi

        echo "Starting cmsRun for $sublist_file with output file $output_file"
        (cmsRun conf_FMT.py "$sublist_file" "$output_file" > "output_${i}.log" 2>&1) &
    done
fi

wait

for i in $(seq 0 $((num_groups - 1)))
do
    if grep -q "Fatal Exception" "output_$i.log" || grep -q "Traceback" "output_$i.log"; then
        echo "Warning: Errors detected in output_$i.log. Please check this log for details."
    else
        echo "cmsRun for sublist ${custom_lists[$i]:-list$i.txt} completed successfully."
    fi
done

echo "All cmsRun conf_FMT.py processes have finished."
