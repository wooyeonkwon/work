#!/bin/bash

# Default directory paths
directory_path="/data1/users/dndus0107/AnalysisResults/Muon/crab_MuonSkimming_Run2022C/241108_143513/0000/"
output_directory="/data1/users/dndus0107/AnalysisResults/processed_data_2022C/"

#thread
num_groups=54

# Parse command-line arguments
while getopts "d:l:o:" opt; do
    case $opt in
        d) directory_path="$OPTARG" ;;  # Set directory_path from -d option
        l) IFS=',' read -r -a custom_lists <<< "$OPTARG" ;;  # Parse -l option as a comma-separated list
        o) output_directory="$OPTARG" ;;  # Set output directory from -o option
        *) echo "Usage: $0 [-d <directory_path>] [-l <list1.txt,list2.txt,...>] [-o <output_directory>]"; exit 1 ;;
    esac
done

# Ensure output directory exists
mkdir -p "$output_directory"

# Expand wildcard paths if directory_path contains a wildcard
if [[ "$directory_path" == *"*"* ]]; then
    expanded_paths=($(ls -d $directory_path 2>/dev/null))  # Expand wildcard paths
    if [ ${#expanded_paths[@]} -eq 0 ]; then
        echo "Error: No directories match the wildcard pattern: $directory_path"
        exit 1
    fi
    root_files=()
    for path in "${expanded_paths[@]}"; do
        root_files+=( $(find "$path" -type f -name "*.root") )  # Collect .root files from all expanded paths
    done
else
    root_files=($(find "$directory_path" -type f -name "*.root"))  # Collect .root files from the specified directory
fi

# Process custom lists if provided
if [ ! -z "${custom_lists+x}" ]; then
    # Combine all custom lists into one array
    combined_root_files=()
    for custom_list in "${custom_lists[@]}"; do
        while IFS= read -r line; do
            combined_root_files+=("$line")
        done < "$custom_list"
    done

    # Split combined_root_files into num_groups
    files_per_group=$(( (${#combined_root_files[@]} + num_groups - 1) / num_groups ))  # Ensure proper division for uneven splits

    for i in $(seq 0 $((num_groups - 1)))
    do
        start_index=$(( i * files_per_group ))
        end_index=$(( start_index + files_per_group ))

        sublist=(${combined_root_files[@]:$start_index:$files_per_group})  # Extract a subset of files for this group
        if [ ${#sublist[@]} -eq 0 ]; then
            continue  # Skip empty sublists
        fi

        sublist_file="$output_directory/list${i}.txt"  # Sublist file path in the output directory
        output_file="$output_directory/output_data_${i}.root"  # Output ROOT file path

        # Write the sublist to a file
        printf "%s\n" "${sublist[@]}" > "$sublist_file"

        # Run cmsRun in the background
        echo "Starting cmsRun for $sublist_file with output file $output_file"
        (cmsRun conf_FMT.py "$sublist_file" "$output_file" > "$output_directory/output_${i}.log" 2>&1) &
    done
else
    # Use root_files directly if no custom lists are provided
    files_per_group=$(( (${#root_files[@]} + num_groups - 1) / num_groups ))  # Ensure proper division for uneven splits

    for i in $(seq 0 $((num_groups - 1)))
    do
        start_index=$(( i * files_per_group ))
        end_index=$(( start_index + files_per_group ))

        sublist=(${root_files[@]:$start_index:$files_per_group})  # Extract a subset of files for this group
        if [ ${#sublist[@]} -eq 0 ]; then
            continue  # Skip empty sublists
        fi

        sublist_file="$output_directory/list${i}.txt"  # Sublist file path in the output directory
        output_file="$output_directory/output_data_${i}.root"  # Output ROOT file path

        # Write the sublist to a file
        printf "%s\n" "${sublist[@]}" > "$sublist_file"

        # Run cmsRun in the background
        echo "Starting cmsRun for $sublist_file with output file $output_file"
        (cmsRun conf_FMT.py "$sublist_file" "$output_file" > "$output_directory/output_${i}.log" 2>&1) &
    done
fi

# Wait for all background processes to finish
wait

# Check each log file for errors
for i in $(seq 0 $((num_groups - 1)))
do
    if grep -q "Fatal Exception" "$output_directory/output_${i}.log" || grep -q "Traceback" "$output_directory/output_${i}.log"; then
        echo "Warning: Errors detected in $output_directory/output_${i}.log. Please check this log for details."
    else
        echo "cmsRun for list${i}.txt completed successfully."
    fi
done

echo "All cmsRun conf_FMT.py processes have finished."
