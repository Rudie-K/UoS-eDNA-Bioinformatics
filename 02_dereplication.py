#!/usr/bin/env python3
import os
import subprocess
import shlex
import csv
import re
import sys
from datetime import datetime

# ==========================================
# CONFIGURATION: AUTOMATED RUN RESOLUTION
# ==========================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RUNS_DIR = os.path.join(BASE_DIR, "runs")

# 1. Dynamically locate the latest run folder from upstream processing
if not os.path.exists(RUNS_DIR) or not os.listdir(RUNS_DIR):
    print("[Error] No active run folders found. Please run Stage 1 (trim_diagnostic.py) first.")
    sys.exit(1)

all_runs = sorted([d for d in os.listdir(RUNS_DIR) if os.path.isdir(os.path.join(RUNS_DIR, d))])
TARGET_RUN = all_runs[-1]
RUN_PATH = os.path.join(RUNS_DIR, TARGET_RUN)

TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
OUTPUT_BASE = os.path.join(RUN_PATH, f"dereplicated_data_{TIMESTAMP}")
LOG_DIR = os.path.join(RUN_PATH, "logs", f"derep_{TIMESTAMP}")

# Report Files for the Institutional Audit Trail
SUMMARY_CSV = os.path.join(OUTPUT_BASE, f"derep_master_summary_{TIMESTAMP}.csv")
AUDIT_REPORT = os.path.join(OUTPUT_BASE, f"derep_technical_audit_{TIMESTAMP}.txt")

INPUT_BASE = os.path.join(RUN_PATH, "processed_data")
LOCI = ["MiFish_12S", "MarVer3_16S"]
MIN_FILE_SIZE = 1 * 1024  # Standardized 1KB guard limit to bypass completely blank files safely

# ==========================================
# HELPER: AUDIT LOGIC
# ==========================================
def count_fasta_stats(fasta_path):
    """Parses FASTA to count unique variants and total biological reads."""
    uniques = 0
    total_size = 0
    if os.path.exists(fasta_path):
        with open(fasta_path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    uniques += 1
                    # Extracts the required ;size=X; abundance metadata tag
                    match = re.search(r";size=(\d+);", line)
                    if match:
                        total_size += int(match.group(1))
    return uniques, total_size

def run_dereplication():
    print(f"\n>>> Stage 2A: Starting Automated Audited Dereplication")
    print(f">>> Target Run Folder: {TARGET_RUN}")
    print(f">>> Input Source:      {INPUT_BASE}")
    os.makedirs(OUTPUT_BASE, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)

    summary_data = []

    for locus in LOCI:
        locus_in = os.path.join(INPUT_BASE, locus)
        locus_out = os.path.join(OUTPUT_BASE, locus)
        if not os.path.exists(locus_in): 
            continue
        os.makedirs(locus_out, exist_ok=True)

        files = [f for f in os.listdir(locus_in) if f.endswith("_R1.fq")]

        for f in files:
            sample_id = f.replace("_R1.fq", "")
            r1 = os.path.join(locus_in, f)
            fa_out = os.path.join(locus_out, f"{sample_id}_unique.fa")
            log_file = os.path.join(LOG_DIR, f"{sample_id}_{locus}_derep.log")
            
            # 1. Quality Guard: Size Validation Check
            file_size = os.path.getsize(r1)
            if file_size < MIN_FILE_SIZE:
                summary_data.append([sample_id, locus, "Skipped", f"File size too minimal ({file_size} bytes)", 0, 0])
                continue

            # 2. Safe Cross-Platform Fastq Sanitization Block
            sanitized_fq = r1 + ".tmp"
            print(f"  Sanitizing & Cleaning: {sample_id} ({locus})...")
            
            # Formatted clean-read selection pattern filter
            awk_script = r'NR%4==1{a=$0} NR%4==2{b=$0} NR%4==3{c=$0} NR%4==0{if(length(b)>0) print a"\n"b"\n"c"\n"$0}'
            try:
                with open(sanitized_fq, "w") as out_f:
                    subprocess.run(["awk", awk_script, r1], stdout=out_f, check=True)
            except (subprocess.CalledProcessError, FileNotFoundError):
                # Fallback to copy if system awk environment is missing
                import shutil
                shutil.copy2(r1, sanitized_fq)

            # 3. USEARCH Fastx Uniques Abundance Execution Matrix
            print(f"  Dereplicating variants via USEARCH...")
            cmd = [
                "usearch", "-fastx_uniques", sanitized_fq, "-fastaout", fa_out,
                "-sizeout", "-relabel", f"{sample_id}.", "-strand", "both"
            ]

            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            
            # Clean up the large text temporary file to preserve drive space
            if os.path.exists(sanitized_fq): 
                os.remove(sanitized_fq)
                
            with open(log_file, "w") as lh:
                lh.write(result.stdout)

            if result.returncode == 0:
                uniques, total_reads = count_fasta_stats(fa_out)
                print(f"    [Success] Formatted {uniques} variants ({total_reads} total biological reads).")
                summary_data.append([sample_id, locus, "Success", "N/A", total_reads, uniques])
            else:
                error_context = result.stdout.strip()[:50].replace("\n", " ")
                print(f"    [Error Failed] Check log trail file: {os.path.basename(log_file)}")
                summary_data.append([sample_id, locus, "Failed", f"USEARCH: {error_context}", 0, 0])

    # ==========================================
    # FINAL REPORT AUDIT LOG TRAIL COMPILATION
    # ==========================================
    # Write Master Summary Data Matrix CSV File
    with open(SUMMARY_CSV, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Sample", "Locus", "Status", "Details", "Total_Biological_Reads", "Unique_Variants"])
        writer.writerows(summary_data)

    # Write Human-Readable Technical Evaluation Audit Text File
    with open(AUDIT_REPORT, "w") as report:
        report.write(f"Sussex Stage 2 Audit Trail Baseline Matrix: {TARGET_RUN}\n" + "="*60 + "\n")
        report.write(f"Execution Completed On: {datetime.now()}\n")
        report.write(f"Target Output Folder:    {OUTPUT_BASE}\n\n")
        for row in summary_data:
            # FIXED: Aligned indices to display data column items accurately
            report.write(f"Sample: {row[0]:<15} | Locus: {row[1]:<12} | Status: {row[2]:<8} | Reads: {row[4]:<8} | Variants: {row[5]:<6}\n")

    print(f"\n[Complete] Stage 2A Dereplication Concluded Successfully.")
    print(f"[Report] Master summary audit catalog written to: {SUMMARY_CSV}")

if __name__ == "__main__":
    run_dereplication()
