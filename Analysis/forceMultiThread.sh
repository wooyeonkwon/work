#!/bin/bash

# Default directory paths
directory_path="/pnfs/knu.ac.kr/data/cms/store/user/wkwon/AnalysisResults/Muon/crab_MuonSkimming_Run2022D/241108_143535/0000/"
output_directory="/pnfs/knu.ac.kr/data/cms/store/user/wkwon/AnalysisResults/processed_data/"

# Parse command-line arguments
while getopts "d:l:o:" opt; do
    case $opt in
        d) directory_path="$OPTARG" ;;
        l) IFS=',' read -r -a custom_lists <<< "$OPTARG" ;;
        o) output_directory="$OPTARG" ;;
        *) echo "Usage: $0 [-d <directory_path>] [-l <list1.txt,list2.txt,...>] [-o <output_directory>]"; exit 1 ;;
    esac
done

# Ensure output directory exists
mkdir -p "$output_directory"

# Expand wildcard paths if directory_path contains a wildcard
if [[ "$directory_path" == *"*"* ]]; then
    expanded_paths=$(echo $directory_path)
    if [ -z "$expanded_paths" ]; then
        echo "Error: No directories match the wildcard pattern: $directory_path"
        exit 1
    fi
    root_files=($(find $expanded_paths -type f -name "*.root"))
else
    root_files=($(find "$directory_path" -type f -name "*.root"))
fi

# Find all root files if custom lists are not provided
if [ -z "${custom_lists+x}" ]; then
    # Number of groups to split the files into
    num_groups=16
    files_per_group=$(( ${#root_files[@]} / num_groups ))

    for i in $(seq 0 $((num_groups - 1)))
    do
        start_index=$(( i * files_per_group ))
        end_index=$(( start_index + files_per_group ))

        sublist=(${root_files[@]:$start_index:$files_per_group})
        sublist_file="$output_directory/list${i}.txt"
        output_file="$output_directory/output_data_${i}.root"

        # Write the sublist to a file
        printf "%s\n" "${sublist[@]}" > "$sublist_file"

        # Run cmsRun in the background
        echo "Starting cmsRun for $sublist_file with output file $output_file"
        (cmsRun conf_FMT.py "$sublist_file" "$output_file" > "$output_directory/output_${i}.log" 2>&1) &
    done
else
    # Use custom lists provided via -l
    for i in "${!custom_lists[@]}"
    do
        sublist_file="${custom_lists[$i]}"
        output_file="$output_directory/output_data_${i}.root"

        # Run cmsRun in the background
        echo "Starting cmsRun for $sublist_file with output file $output_file"
        (cmsRun conf_FMT.py "$sublist_file" "$output_file" > "$output_directory/output_${i}.log" 2>&1) &
    done
fi

# Wait for all background processes to finish
wait

# Check each log file for errors
if [ -z "${custom_lists+x}" ]; then
    for i in $(seq 0 $((num_groups - 1)))
    do
        if grep -q "Fatal Exception" "$output_directory/output_${i}.log" || grep -q "Traceback" "$output_directory/output_${i}.log"; then
            echo "Warning: Errors detected in $output_directory/output_${i}.log. Please check this log for details."
        else
            echo "cmsRun for list${i}.txt completed successfully."
        fi
    done
else
    for i in "${!custom_lists[@]}"
    do
        if grep -q "Fatal Exception" "$output_directory/output_${i}.log" || grep -q "Traceback" "$output_directory/output_${i}.log"; then
            echo "Warning: Errors detected in $output_directory/output_${i}.log. Please check this log for details."
        else
            echo "cmsRun for ${custom_lists[$i]} completed successfully."
        fi
    done
fi

echo "All cmsRun conf_FMT.py processes have finished."
