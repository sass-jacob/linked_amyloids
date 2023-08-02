#!/bin/bash
#SBATCH -J make_directory
#SBATCH -N 1
#SBATCH -n 1

directory="~/md_proteins"
# Find the largest beginning number in existing directories
largest_number=0

for dir in "${directory}"/*; do
  echo "$dir"
  if [[ -d "$dir" ]]; then
    dir_name=$(basename "$dir")
    number="${dir_name%%[!0-9]*}"
    number=$(echo "$number" | sed 's/^0*//') # Remove leading zeros
    if (( number > largest_number )); then
      largest_number=$number
    fi
  fi
done

echo "$largest_number"

directory="~/md_proteins/pdb_files"
beginning_number=$((largest_number + 1))

# Iterate through the files in the directory
for file in "${directory}"/*; do
  # Extract the filename without the extension
  filename=$(basename "$file")
  filename_without_extension="${filename%.*}"
  
  # Extract the first 6 letters from the filename
  prefix="${filename_without_extension:0:6}"

  # Create the directory with a unique number
  dir_name=$(printf "%03d_%s" "$beginning_number" "${prefix}")
  mkdir -p "../${dir_name}"
  
  # Move the file to the newly created directory
  mv "$file" "../${dir_name}/${prefix}_afold.pdb"
  
  ((beginning_number++))
done