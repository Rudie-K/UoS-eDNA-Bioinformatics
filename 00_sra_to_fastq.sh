#!/bin/bash
# ==============================================================================
# Stage 0: Data Acquisition & Standardization (Fully Automated)
# Purpose: Auto-locate benchmark directory and convert SRA data to 'fqconvert'.
# Usage: ./00_sra_to_fastq.sh
# ==============================================================================

# 1. Dependency Check
if ! command -v fasterq-dump &> /dev/null; then
    echo "[Error] SRA Toolkit not found. Please install (sudo apt install sra-toolkit)."
    exit 1
fi

# 2. Automated Target Directory Discovery
# This checks for the exact dataset folder name within your current directory
TARGET_NAME="SR normalised Clark et al 2024"

if [ -d "$TARGET_NAME" ]; then
    INPUT_DIR="$TARGET_NAME"
    echo "[Found] Target dataset identified: $INPUT_DIR"
else
    # Fallback/Safety Check: Create an empty template if not present
    echo "[Notice] '$TARGET_NAME' folder not detected in current working directory."
    mkdir -p "$TARGET_NAME"
    echo "         Created empty folder. Please place your raw .sra files inside '$TARGET_NAME' and rerun."
    exit 1
fi

# 3. Stable Target Directory Logic for Stage 1
OUTPUT_DIR="fqconvert"

if [ -d "$OUTPUT_DIR" ]; then
    echo "[Notice] Found existing '$OUTPUT_DIR' folder. Clearing old contents for fresh run..."
    rm -rf "$OUTPUT_DIR"/*
else
    mkdir -p "$OUTPUT_DIR"
fi

echo "-------------------------------------------------------"
echo "Starting Conversion Run"
echo "Source: $INPUT_DIR"
echo "Output Folder: $OUTPUT_DIR"
echo "-------------------------------------------------------"

# 4. Batch Conversion Loop (Handles files with potential spaces safely)
for sra_file in "$INPUT_DIR"/*; do
    if [ -f "$sra_file" ]; then
        filename=$(basename "$sra_file")
        echo "Standardizing: $filename"
        
        # --split-files is MANDATORY for paired-end synchronization
        fasterq-dump --split-files --outdir "$OUTPUT_DIR" "$sra_file"
        
        # Check return code for analytical transparency
        if [ $? -eq 0 ]; then
            echo " [Success] $filename processed."
        else
            echo " [Error] Failed to convert $filename."
        fi
    fi
done

echo "-------------------------------------------------------"
echo "Standardization complete. Files located in: $OUTPUT_DIR"
