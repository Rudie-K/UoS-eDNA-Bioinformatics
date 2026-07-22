#!/usr/bin/env python3
import os
import subprocess
import csv
import time
import re
import sys
from datetime import datetime
from pathlib import Path

# ==============================================================================
# 1. GOVERNANCE CONFIGURATION & TARGET METRIC BRACKETS
# ==============================================================================
THRESHOLDS = {"Species": 99.0, "Genus": 97.0, "Family": 95.0}
ABUNDANCE_FILTERS = {"MiFish_12S": 0.0002, "MarVer3_16S": 0.00025}

LENGTH_WINDOWS = {
    "MiFish_12S": {"min": 163, "max": 185},
    "MarVer3_16S": {"min": 232, "max": 274}
}

ABSOLUTE_MIN_SIZE = 5

# ==============================================================================
# 2. RUN PATHWAY ENVIRONMENT RESOLUTION
# ==============================================================================
BASE_DIR = Path(__file__).resolve().parent
RUNS_DIR = BASE_DIR / "runs"

if not RUNS_DIR.exists() or not any(RUNS_DIR.iterdir()):
    print("[Error] No run directories discovered. Execute Stage 1 and Stage 2 first.")
    sys.exit(1)

all_runs = sorted([d for d in RUNS_DIR.iterdir() if d.is_dir()])
RUN_PATH = all_runs[-1]
TARGET_RUN = RUN_PATH.name

def get_latest_folder(prefix):
    folders = sorted([d for d in RUN_PATH.iterdir() if d.is_dir() and d.name.startswith(prefix)])
    return folders[-1] if folders else None

def load_stage1_counts():
    summary_csv_path = RUN_PATH / "master_processing_summary.csv"
    read_counts = {}
    if not summary_csv_path.exists():
        return read_counts
    with summary_csv_path.open('r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                sample_id = row['Sample']
                read_counts[sample_id] = read_counts.get(sample_id, 0) + int(row['Total'])
            except (KeyError, ValueError):
                continue
    return read_counts

def parse_fasta_sequences(fasta_path):
    """Extract headers, size annotations, and raw sequences from a FASTA file.
    
    Translates USEARCH amp designations inside the report file into Zotu IDs
    to recover and assign the true sequence abundances.
    """
    sequences = []
    current_header = None
    current_seq = []
    
    if not fasta_path.exists():
        return sequences

    # Reconstruct the correct matching path to the adjacent report text file
    size_lookup = {}
    report_path = Path(str(fasta_path).replace("_zotus.fa", "_unoise_report.txt"))
    
    if report_path.exists():
        try:
            with report_path.open('r') as rf:
                for line in rf:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        # Extract the unique unique ID string (e.g., "SRR..._3p..._Valsecchi.1;size=15297;")
                        raw_id_field = parts[0]
                        # Extract the final clustering mapping token (e.g., "amp1" or "bad")
                        mapping_decision = parts[2]
                        
                        # Process only rows that were successfully retained as active operational units
                        if mapping_decision.startswith("amp"):
                            # Translate "amp1" -> "1"
                            amp_number = mapping_decision.replace("amp", "").strip()
                            zotu_key = f"Zotu{amp_number}"
                            
                            # Parse out the actual sequence read abundance
                            size_match = re.search(r";size=(\d+);", raw_id_field)
                            if size_match:
                                size_lookup[zotu_key] = int(size_match.group(1))
        except Exception:
            pass

    # Read the FASTA file arrays sequentially
    with fasta_path.open('r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    seq_str = "".join(current_seq)
                    # Isolate the clean header key (e.g., "Zotu1")
                    zotu_id = current_header.split(';')[0].strip()
                    
                    # Pull the true recovered size, default to 1 if missing
                    size = size_lookup.get(zotu_id, 1)
                    sequences.append({"header": current_header, "seq": seq_str, "size": size})
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            seq_str = "".join(current_seq)
            zotu_id = current_header.split(';')[0].strip()
            size = size_lookup.get(zotu_id, 1)
            sequences.append({"header": current_header, "seq": seq_str, "size": size})
            
    return sequences

# ==============================================================================
# 3. OPTIMIZED EXECUTION PROCESS MATRIX
# ==============================================================================
def main():
    print(f"\n>>> Stage 3: Starting Optimized Combined Taxonomic Pipeline")
    print(f">>> Target Run Folder: {TARGET_RUN}")

    denoised_dir = get_latest_folder("denoised_data_")
    if not denoised_dir:
        print("[Fatal Error] Missing required Stage 2 denoised data subfolders.")
        sys.exit(1)

    read_counts = load_stage1_counts()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = RUN_PATH / f"final_species_list_{timestamp}"
    log_dir = RUN_PATH / f"logs" / f"taxonomy_{timestamp}"
    output_dir.mkdir(exist_ok=True)
    log_dir.mkdir(exist_ok=True, parents=True)

    loc_totals = {"MiFish_12S": {"total": 0, "kept": 0}, "MarVer3_16S": {"total": 0, "kept": 0}}
    
    global_uniques = {}  
    sample_mappings = [] 
    mode_lookup = {}
    
    summary_csv_path = RUN_PATH / "master_processing_summary.csv"
    if summary_csv_path.exists():
        with summary_csv_path.open('r', newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                mode_lookup[(row['Sample'], row['Locus'])] = row['Mode']

    loci = ["MiFish_12S", "MarVer3_16S"]
    for locus in loci:
        locus_in = denoised_dir / locus
        if not locus_in.exists():
            continue
        
        zotu_files = list(locus_in.glob("*_zotus.fa"))
        for zf in zotu_files:
            sample_id = zf.name.replace("_zotus.fa", "")
            records = parse_fasta_sequences(zf)
            processing_mode = mode_lookup.get((sample_id, locus), "Length-Only")
            
            print(f"\nProcessing sample: {sample_id} ({locus})")
            print(f"Found {len(records)} sequences")

            for r in records:
                loc_totals[locus]["total"] += 1
                seq_len = len(r["seq"])

                print(f"  {r['header']}")
                print(f"     Length = {seq_len}")
                print(f"     Size   = {r['size']}")

                if processing_mode == "Length-Only":
                    bounds = LENGTH_WINDOWS[locus]
                    if not (bounds["min"] <= seq_len <= bounds["max"]):
                        print("     FAILED length filter")
                        continue

                if r["size"] < ABSOLUTE_MIN_SIZE:
                    print("     FAILED abundance filter")
                    continue

                print("     PASSED local filters")

                loc_totals[locus]["kept"] += 1
                seq_text = r["seq"].upper()
                
                if seq_text not in global_uniques:
                    global_id = f"GlobalSeq_{len(global_uniques) + 1:06d}"
                    global_uniques[seq_text] = {"id": global_id, "locus": locus}
                else:
                    global_id = global_uniques[seq_text]["id"]
                    
                sample_mappings.append({
                    "sample_id": sample_id, "locus": locus, "zotu_id": r["header"],
                    "seq": seq_text, "global_id": global_id, "size": r["size"]
                })


    if not global_uniques:
        print("[Notice] Zero target sequences passed local pre-filters. Aborting remote call.")
        sys.exit(0)

    print(f"  [12S MiFish] Sequences processed: {loc_totals['MiFish_12S']['total']} -> {loc_totals['MiFish_12S']['kept']} passed.")
    print(f"  [16S MarVer] Sequences processed: {loc_totals['MarVer3_16S']['total']} -> {loc_totals['MarVer3_16S']['kept']} passed.")
    print(f"  🚀 Global Run Pooling compressed queries down to {len(global_uniques)} unique NCBI uploads.")

    global_query_fa = output_dir / "global_query_matrix.fa"
    with global_query_fa.open('w') as out_fa:
        for seq_text, info in global_uniques.items():
            out_fa.write(f">{info['id']}\n{seq_text}\n")

    global_blast_txt = output_dir / "global_blast_results.txt"
    blast_log = log_dir / "global_ncbi_connection.log"
    
    cmd = [
        "blastn", "-query", str(global_query_fa), "-db", "nt", "-remote",
        "-max_target_seqs", "1",
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore sscinames sskingdoms",
        "-out", str(global_blast_txt)
    ]

    success = False
    for attempt in range(1, 4):
        print(f"  Contacting NCBI Endpoint Cloud servers... (Attempt {attempt}/3)")
        completed = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        
        if completed.returncode == 0:
            print("    [Success] Global alignment table received safely from NCBI.")
            success = True
            break
        else:
            with blast_log.open('a') as lh:
                lh.write(f"--- Attempt {attempt} Error Output ---\n{completed.stdout}\n")
            if attempt < 3:
                print("    [Warning] Cloud endpoint rate-limit or drop detected. Retrying in 30s...")
                time.sleep(30)

    if not success:
        print("[Fatal Error] Remote NCBI connection dropped permanently. Run diagnostic check.")
        sys.exit(1)

    print("\n--- Step 3: Re-Mapping Alignments to Institutional Framework ---")
    blast_map = {}
    if global_blast_txt.exists():
        with global_blast_txt.open('r') as bf:
            reader = csv.reader(bf, delimiter='\t')
            for row in reader:
                if len(row) >= 3:
                    # FIX: Corrected row index assignments to solve previous syntax/parsing bugs
                    raw_id = row[0].replace("lcl|", "").strip()
                    blast_map[raw_id] = {
                        "accession": row[1], 
                        "identity": float(row[2]),
                        "sci_name": row[6] if len(row) > 6 and row[6].strip() != "" and row[6] != "N/A" else f"Unknown ({row[1]})"
                    }

    master_composition = []

    global_seq_lookup = {info["id"]: seq_text for seq_text, info in global_uniques.items()}

    for locus in loci:
        locus_maps = [m for m in sample_mappings if m["locus"] == locus]
        distinct_samples = sorted(list(set(m["sample_id"] for m in locus_maps)))
        
        for sample_id in distinct_samples:
            s_records = [m for m in locus_maps if m["sample_id"] == sample_id]
            total_reads = read_counts.get(sample_id, 0)
            if total_reads == 0:
                total_reads = sum(r["size"] for r in s_records)
                
            min_required = total_reads * ABUNDANCE_FILTERS[locus]
            audit_log = output_dir / f"{sample_id}_{locus}_filtering_audit.txt"
            
            with audit_log.open('w') as log:
                log.write(f"Sussex Stage 3 Consolidated Audit: {sample_id} ({locus})\n")
                log.write(f"Calculated Sequence Total: {total_reads} | Abundance Gate Filter: {min_required:.2f} reads\n\n")
                
                for r in s_records:
                    size = r["size"]
                    zotu_id = r["zotu_id"]
                    global_id = r["global_id"]
                    
                    # FILTER A: Proportional Abundance Screen
                    # Retrieve BLAST result for this sequence
                    blast_hit = blast_map.get(global_id)

                    if blast_hit is None:
                        acc_id = "N/A"
                        pident = 0.0
                        sci_name = "No Match"
                        assigned_rank = "Unassigned"

                        log.write(f"[Discard] {zotu_id}: No BLAST match found.\n")
                        continue

                    acc_id = blast_hit["accession"]
                    pident = blast_hit["identity"]
                    sci_name = blast_hit["sci_name"]

                    # FILTER A: Proportional Abundance Screen
                    if size < min_required:
                        log.write(
                            f"[Discard] {zotu_id}: Abundance {size} below minimum "
                            f"threshold ({min_required:.2f} reads)\n"
                        )
                        continue

                    # FILTER B: Identity Threshold Screen
                    if pident >= THRESHOLDS["Species"]:
                        assigned_rank = "Species"
                    elif pident >= THRESHOLDS["Genus"]:
                        assigned_rank = "Genus"
                    elif pident >= THRESHOLDS["Family"]:
                        assigned_rank = "Family"
                    else:
                        assigned_rank = "Below Threshold"

                    if assigned_rank == "Below Threshold":
                        log.write(
                            f"[Discard] {zotu_id}: Identity falls below "
                            f"institutional 95% threshold ({pident:.2f}%)\n"
                        )
                        continue

                    raw_sequence_string = global_seq_lookup.get(global_id, "N/A")


                    master_composition.append([
                        sample_id,
                        locus,
                        raw_sequence_string,
                        sci_name,
                        acc_id,
                        assigned_rank,
                        pident,
                        size
                    ])

                    log.write(
                        f"[KEEP] {sci_name} ({assigned_rank}): "
                        f"{pident:.2f}% identity, {size} total reads\n"
                    )

    master_csv = output_dir / "master_species_composition.csv"
    with master_csv.open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Sample", "Locus", "Sequence", "Scientific_Name", "Accession", "Rank", "Identity_Pct", "Reads"])
        writer.writerows(master_composition)

    if global_query_fa.exists(): 
        os.remove(global_query_fa)

    print(f"\n[Complete] Combined Stage 3 Processing Finished Successfully.")
    print(f"[Report] Master composition table catalog generated here: {master_csv}")

if __name__ == "__main__":
    main()
