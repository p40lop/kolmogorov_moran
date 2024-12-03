#!/bin/bash

# Check for correct number of arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <folder_name> <file_name>"
    exit 1
fi

folder_name="$1"
file_name="$2"

# Find all matching files recursively
find "$folder_name" -type f -name "$file_name" | while read file; do
    awk '{a[$1] = $0} END {PROCINFO["sorted_in"] = "@ind_num_asc"; for (i in a) print a[i]}' "$file" > "${file}.sorted"
    rm "${file}"
    mv "${file}.sorted" "${file}"
done