#!/bin/bash

# Default values
input_paths=()
output_directory=""
num_groups=56

# Function to print usage
print_usage() {
    echo "Usage: $0 -i <input_path1,input_path2,...> -o <output_directory> [-l <list1.txt,list2.txt,...>] [-w <num_workers>]"
    echo "Required options:"
    echo "  -i: Input directory paths (comma-separated)"
    echo "  -o: Output directory path"
    echo "Optional options:"
    echo "  -l: Comma-separated list of input files"
    echo "  -w: Number of worker groups (default: 56)"
    exit 1
}

# Parse command-line arguments
while getopts "i:l:o:w:" opt; do
    case $opt in
        i) IFS=',' read -r -a input_paths <<< "$OPTARG" ;;  # Parse -i option as comma-separated paths
        l) IFS=',' read -r -a custom_lists <<< "$OPTARG" ;;  # Parse -l option as a comma-separated list
        o) output_directory="$OPTARG" ;; # Set output directory from -o option
        w) num_groups="$OPTARG" ;;      # Set number of worker groups from -w option
        *) print_usage ;;
    esac
done

# Check for required options
if [ ${#input_paths[@]} -eq 0 ] || [ -z "$output_directory" ]; then
    echo "Error: -i (input paths) and -o (output directory) are required options"
    print_usage
fi

# Validate num_groups is a positive integer
if ! [[ "$num_groups" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: Number of workers (-w) must be a positive integer"
    exit 1
fi

# Ensure output directory exists
mkdir -p "$output_directory"

# Initialize empty array for all root files
root_files=()

# Process each input path
for input_path in "${input_paths[@]}"; do
    # Expand wildcard paths if input_path contains a wildcard
    if [[ "$input_path" == *"*"* ]]; then
        expanded_paths=($(ls -d $input_path 2>/dev/null))  # Expand wildcard paths
        if [ ${#expanded_paths[@]} -eq 0 ]; then
            echo "Warning: No directories match the wildcard pattern: $input_path"
            continue
        fi
        for path in "${expanded_paths[@]}"; do
            found_files=($(find "$path" -type f -name "*.root"))
            root_files+=("${found_files[@]}")
        done
    else
        found_files=($(find "$input_path" -type f -name "*.root"))
        root_files+=("${found_files[@]}")
    fi
done

# Check if any root files were found
if [ ${#root_files[@]} -eq 0 ]; then
    echo "Error: No .root files found in any of the input paths:"
    printf '%s\n' "${input_paths[@]}"
    exit 1
fi

echo "Found ${#root_files[@]} root files in total"

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
    files_per_group=$(( (${#root_files[@]} + num_groups - 1) / num_groups ))  # Calculate files per group with ceiling division

    for i in $(seq 0 $((num_groups - 1)))
    do
        start_index=$(( i * files_per_group ))
        current_group_size=$files_per_group
        
        # Adjust group size for the last group if necessary
        if [ $(( start_index + current_group_size )) -gt ${#root_files[@]} ]; then
            current_group_size=$(( ${#root_files[@]} - start_index ))
        fi

        # Skip if we've processed all files
        if [ $start_index -ge ${#root_files[@]} ]; then
            continue
        fi

        # Extract sublist for this group
        sublist=("${root_files[@]:$start_index:$current_group_size}")

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