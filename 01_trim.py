#!/usr/bin/env python3
import csv
import gzip
import os
import re
import sys
import shutil
import subprocess
from datetime import datetime

# ==========================================================
# 1. STREAMLINED CONFIGURATION: Targets Stable 'fqconvert'
# ==========================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RUN_ROOT = os.path.join(SCRIPT_DIR, "runs")

# Stably points to the automated folder created in Stage 0
INPUT_DIR = os.path.join(SCRIPT_DIR, "fqconvert")

# Fail-safe check if Stage 0 has not been executed yet
if not os.path.exists(INPUT_DIR):
    os.makedirs(INPUT_DIR, exist_ok=True)
    print(f"[Notice] Created empty source folder at: {INPUT_DIR}")
    print("Please run 'sra_to_fastq.sh' first to populate data.")

def create_run_directory(run_root):
    run_date = datetime.now().strftime("%Y-%m-%d")
    run_index = 1
    while True:
        run_dir = os.path.join(run_root, f"{run_date}_{run_index}")
        if not os.path.exists(run_dir):
            os.makedirs(run_dir, exist_ok=True)
            return run_dir
        run_index += 1

RUN_DIR = create_run_directory(RUN_ROOT)
OUTPUT_BASE = os.path.join(RUN_DIR, "processed_data")
LOG_DIR = os.path.join(RUN_DIR, "logs")
REPORT_FILE = os.path.join(RUN_DIR, "final_processing_report.txt")
SUMMARY_CSV = os.path.join(RUN_DIR, "master_processing_summary.csv")

SUMMARY_FIELDS = ["Sample", "Locus", "Mode", "Variant", "Status", "Details", "Total", "Recovered", "Rate"]

# Ensure Cutadapt is available before kicking off the run matrix
if shutil.which("cutadapt") is None:
    raise RuntimeError("cutadapt tool was not detected on system PATH.")

primers = {
    "MiFish_12S": {
        "min_len": 100, "max_len": 300, 
        "q_score": 15, "error_rate": 0.20,
        "variants": {
            "Standard_Universal": {
                "fwd": "GTCGGTAAAACTCGTGCCAGC", "rev": "CATAGTGGGGTATCTAATCCCAGTTTG",
                "fwd_rc": "GCTGGCACGAGTTTTACCGAC", "rev_rc": "CAAACTGGGATTAGATACCCCACTATG"
            }
        }
    },
    "MarVer3_16S": {
        "min_len": 200, "max_len": 350,
        "q_score": 15, # Relaxed to catch primers at noisy 250bp ends
        "error_rate": 0.20, # Increased tolerance for short terminal matches
        "variants": {
            "Standard_Valsecchi": {
                "fwd": "AGACGAGAAGACCCTRTG", "rev": "GGATTGCGCTGTTATCCC",
                "fwd_rc": "CAYAGGGTCTTCTCGTCT", "rev_rc": "GGGATAACAGCGCAATCC"
            }
        }
    }
}

os.makedirs(LOG_DIR, exist_ok=True)
for locus in primers.keys(): os.makedirs(f"{OUTPUT_BASE}/{locus}", exist_ok=True)

if not os.path.isdir(INPUT_DIR):
    raise FileNotFoundError(f"Input directory not found: {INPUT_DIR}")

if shutil.which("cutadapt") is None:
    raise RuntimeError("cutadapt is not available on PATH")

# ==========================================
# 2. HELPER FUNCTIONS
# ==========================================
def scout_primers(r1, locus, locus_cfg):
    """Broadened Scout (50bp) to identify orientation."""
    for v_name, s in locus_cfg['variants'].items():
        # Use short signatures rather than the full primer so we can detect noisy ends quickly.
        fwd_sig, rev_rc_sig = s['fwd'][:12], s['rev_rc'][-12:]
        fwd_hits, rev_3p_hits = 0, 0
        try:
            # Support both plain FASTQ and gzipped FASTQ without changing the rest of the pipeline.
            opener = gzip.open if r1.endswith('.gz') else open
            with opener(r1, 'rt') as f:
                for i, line in enumerate(f):
                    if i >= 4000: break
                    if i % 4 == 1:
                        if fwd_sig in line: fwd_hits += 1
                        if rev_rc_sig in line: rev_3p_hits += 1
            if fwd_hits > 20: return "Strict-5p", v_name, s
            if rev_3p_hits > 20: return "3p-ReadThrough", v_name, s
        except OSError:
            continue
    return "Length-Only", "None", None

def extract_stats(log_file):
    """Robust regex to capture recovery stats for the CSV."""
    stats = {"Total": "0", "Written": "0", "Percent": "0.0%"}
    if not os.path.exists(log_file): return stats
    try:
        with open(log_file, 'r') as f:
            text = f.read()
            # Cutadapt wording differs a little across versions, so accept more than one pattern.
            total_patterns = [
                r"Total read pairs processed:\s+([\d,]+)",
                r"Total reads processed:\s+([\d,]+)",
            ]
            written_patterns = [
                r"Pairs written \(passing filters\):\s+([\d,]+)\s+\(([\d.]+%)\)",
                r"Reads written \(passing filters\):\s+([\d,]+)\s+\(([\d.]+%)\)",
            ]
            for pattern in total_patterns:
                total = re.search(pattern, text)
                if total:
                    stats["Total"] = total.group(1).replace(",", "")
                    break
            for pattern in written_patterns:
                written = re.search(pattern, text)
                if written:
                    stats["Written"] = written.group(1).replace(",", "")
                    stats["Percent"] = written.group(2)
                    break
    except OSError:
        pass
    return stats

def build_cutadapt_command(cmd_args, cfg, out1, out2, r1, r2):
    # Build an argument list instead of a shell string so quoting is predictable and safer.
    return [
        "cutadapt",
        *cmd_args,
        "--revcomp",
        "-e", str(cfg["error_rate"]),
        "-O", "3",
        "-q", str(cfg["q_score"]),
        "-o", out1,
        "-p", out2,
        r1,
        r2,
    ]

def safe_label(text):
    # Normalize labels so they can safely appear in filenames.
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text)

def expected_mate_name(filename):
    # Style 1: Standard Illumina style (_R1_001 / _R2_001)
    if "_R1_001.fastq" in filename:
        return filename.replace("_R1_001.fastq", "_R2_001.fastq")
        
    # Style 2: SRA / Simple style (_1 / _2)
    if filename.endswith("_1.fastq.gz"):
        return filename[:-11] + "_2.fastq.gz"
    if filename.endswith("_1.fastq"):
        return filename[:-8] + "_2.fastq"
        
    return None

def find_mate_file(input_dir, filename):
    # Allow either compression style on the mate so mixed input sets do not get skipped unnecessarily.
    mate_name = expected_mate_name(filename)
    if mate_name is None:
        return None

    mate_path = os.path.join(input_dir, mate_name)
    if os.path.exists(mate_path):
        return mate_path

    if mate_name.endswith(".gz"):
        alternate_mate = mate_name[:-3]
    else:
        alternate_mate = mate_name + ".gz"

    alternate_path = os.path.join(input_dir, alternate_mate)
    if os.path.exists(alternate_path):
        return alternate_path

    return None

# ==========================================
# 3. MAIN PROCESSING LOOP
# ==========================================
summary_data = []

with open(REPORT_FILE, "w") as report:
    report.write("Sussex eDNA Pipeline Audit: Multi-Variant Processing\n" + "="*55 + "\n\n")
    
    # Captures BOTH types of forward files seamlessly
    files = sorted(
        f for f in os.listdir(INPUT_DIR)
        if f.endswith("_1.fastq") or f.endswith("_1.fastq.gz") or 
           f.endswith("_R1_001.fastq") or f.endswith("_R1_001.fastq.gz")
    )
    
    for filename in files:
        # Dynamically extract clean sample_id based on naming style
        if "_R1_001.fastq" in filename:
            suffix_len = 16 if filename.endswith(".gz") else 13
            sample_id = filename[:-suffix_len]
        else:
            suffix_len = 11 if filename.endswith(".gz") else 8
            sample_id = filename[:-suffix_len]

        r1 = os.path.join(INPUT_DIR, filename)
        r2 = find_mate_file(INPUT_DIR, filename)
        print(f"\n>>> Processing Sample: {sample_id}")

        if r2 is None:
            message = f"Sample: {sample_id} | Skipped: missing mate file"
            report.write(message + "\n")
            print("  [Skipped] mate file not found")
            # Record skipped samples explicitly so the CSV still tells the full story.
            for locus in primers:
                summary_data.append({
                    "Sample": sample_id,
                    "Locus": locus,
                    "Mode": "Skipped",
                    "Variant": "None",
                    "Status": "SKIPPED",
                    "Details": "Missing mate file",
                    "Total": "0",
                    "Recovered": "0",
                    "Rate": "0.0%",
                })
            continue
        
        for locus, cfg in primers.items():
            mode, variant, s = scout_primers(r1, locus, cfg)
            
            # Use Locus-specific tuning for Cutadapt
            if mode == "Strict-5p":
                # Temporarily removed length constraint for diagnostic logging
                cmd_args = ["-g", f"{s['fwd']}", "-G", f"{s['rev']}", "--discard-untrimmed"]
            elif mode == "3p-ReadThrough":
                cmd_args = ["-a", s['rev_rc'], "-A", s['fwd_rc'], "--discard-untrimmed"]
            else:
                cmd_args = ["-m", str(cfg['min_len']), "-M", str(cfg['max_len'])]


            mode_tag = safe_label(mode)
            variant_tag = safe_label(variant)
            # Embed the trimming decision in the filename so downstream files keep their provenance.
            out1, out2 = (
                f"{OUTPUT_BASE}/{locus}/{sample_id}_{mode_tag}_{variant_tag}_R1.fq",
                f"{OUTPUT_BASE}/{locus}/{sample_id}_{mode_tag}_{variant_tag}_R2.fq",
            )
            log_file = f"{LOG_DIR}/{sample_id}_{locus}_{mode_tag}_{variant_tag}.log"

            cmd = build_cutadapt_command(cmd_args, cfg, out1, out2, r1, r2)

            with open(log_file, "w") as log_handle:
                completed = subprocess.run(cmd, stdout=log_handle, stderr=subprocess.STDOUT, check=False)

            if completed.returncode != 0:
                report.write(f"Sample: {sample_id} | Locus: {locus} | Mode: {mode} | Variant: {variant} | ERROR: Cutadapt failed with exit code {completed.returncode}\n")
                print(f"  [Error] {locus} ({variant}) via {mode} failed with exit code {completed.returncode}")
                # Keep failures in the summary so they can be audited later instead of disappearing.
                summary_data.append({
                    "Sample": sample_id,
                    "Locus": locus,
                    "Mode": mode,
                    "Variant": variant,
                    "Status": "ERROR",
                    "Details": f"Cutadapt exited with code {completed.returncode}",
                    "Total": "0",
                    "Recovered": "0",
                    "Rate": "0.0%",
                })
                continue
            
            # Collect and log data
            stats = extract_stats(log_file)
            recovered_count = int(stats["Written"])
            # so they don't break dereplication/denoising steps later.
            if recovered_count == 0:
                if os.path.exists(out1): os.remove(out1)
                if os.path.exists(out2): os.remove(out2)
                print(f"  [Dropped] {locus} ({variant}) - 0 reads recovered.")
                
                summary_data.append({
                    "Sample": sample_id, "Locus": locus, "Mode": mode, "Variant": variant, 
                    "Status": "DROPPED", "Details": "Zero reads matching primers", 
                    "Total": stats["Total"], "Recovered": "0", "Rate": "0.0%"
                })
                continue
            # Successful runs get a standard OK record plus the parsed recovery stats.
            summary_row = {"Sample": sample_id, "Locus": locus, "Mode": mode, "Variant": variant, 
                           "Status": "OK", "Details": "", "Total": stats["Total"], "Recovered": stats["Written"], "Rate": stats["Percent"]}
            summary_data.append(summary_row)

            report.write(f"Sample: {sample_id} | Locus: {locus} | Mode: {mode} | Variant: {variant} | Recovery: {stats['Percent']}\n")
            print(f"  [Complete] {locus} ({variant}) via {mode}. Recovery: {stats['Percent']}")

# Final Step: Generate the Master CSV
with open(SUMMARY_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS)
    writer.writeheader()
    writer.writerows(summary_data)

print(f"\nDone! Master Summary: {SUMMARY_CSV}")