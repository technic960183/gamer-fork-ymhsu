#!/bin/bash

# Script to convert symlinks to actual files
# Usage: sh convert_symlinks_to_files.sh

echo "Converting symlinks to actual files..."

# Change to src directory
cd "$(dirname "$0")"

# Files to convert
FILES=(
    "Model_Hydro/GPU_Hydro/CUFLU_FluidSolver_MHM.cu"
    "Model_Hydro/GPU_Hydro/CUFLU_CR_TwoMoment.cu"
)

for file in "${FILES[@]}"; do
    if [ -L "$file" ]; then
        echo "Processing: $file"
        
        # Get the target of the symlink
        target=$(readlink -f "$file")
        
        if [ ! -f "$target" ]; then
            echo "Error: Target file $target does not exist"
            exit 1
        fi
        
        echo "  Symlink points to: $target"
        
        # Copy the target to a temporary file
        temp_file="${file}.tmp"
        cp "$target" "$temp_file"
        
        if [ $? -ne 0 ]; then
            echo "Error: Failed to copy $target to $temp_file"
            exit 1
        fi
        
        # Remove the symlink
        rm "$file"
        
        if [ $? -ne 0 ]; then
            echo "Error: Failed to remove symlink $file"
            # Cleanup temp file
            rm "$temp_file"
            exit 1
        fi
        
        # Move the temporary file to the original symlink location
        mv "$temp_file" "$file"
        
        if [ $? -ne 0 ]; then
            echo "Error: Failed to move $temp_file to $file"
            exit 1
        fi
        
        echo "  Converted: $file is now a regular file"
    else
        echo "Skipping: $file (not a symlink)"
    fi
done

echo "Done! All specified symlinks have been converted to regular files."
