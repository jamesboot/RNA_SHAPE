#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <directory> <output.csv>"
    exit 1
fi

directory=$1
output_file=$2

# Ensure the output file is empty before writing
> "$output_file"

declare -A file1_map

declare -A file2_map

#echo "CommonName,File1,File2" > "$output_file"

# Find files that match the expected patterns
while IFS= read -r file; do
    filename=$(basename "$file")
    path=$(realpath "$file")
    
    if [[ "$filename" =~ ^(.+?)_R1_001\.fastq\.gz$ ]]; then
        common_name="${BASH_REMATCH[1]}"
        file1_map["$common_name"]="$path"
    elif [[ "$filename" =~ ^(.+?)_R2_001\.fastq\.gz$ ]]; then
        common_name="${BASH_REMATCH[1]}"
        file2_map["$common_name"]="$path"
    fi
done < <(find -L "$directory" -type f -name "*_R1_001.fastq.gz" -o -name "*_R2_001.fastq.gz")

# Write paired files to CSV
for common_name in "${!file1_map[@]}"; do
    if [[ -n "${file1_map[$common_name]}" && -n "${file2_map[$common_name]}" ]]; then
        echo "$common_name,${file1_map[$common_name]},${file2_map[$common_name]}" >> "$output_file"
    fi
done

echo "Output saved to $output_file"