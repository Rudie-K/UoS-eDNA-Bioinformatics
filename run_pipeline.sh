#!/usr/bin/env bash
# ==============================================================================
# Sussex eDNA Pipeline
#
# Complete automated workflow
#
# Stage 0  - SRA → FASTQ
# Stage 1  - Primer trimming
# Stage 2  - Dereplication
# Stage 3  - Denoising
# Stage 4  - BLAST annotation
# Stage 5  - Taxonomic annotation
#
# Usage:
#     chmod +x run_pipeline.sh
#     ./run_pipeline.sh
# ==============================================================================

python3 - <<'EOF'
import importlib
packages = [
    "pandas",
    "requests",
    "Bio",
]

missing = []

for p in packages:
    try:
        importlib.import_module(p)
    except ImportError:
        missing.append(p)

if missing:
    print("\nMissing Python packages:")
    for m in missing:
        print("  -", m)
    raise SystemExit(1)
EOF

set -euo pipefail

###############################################################################
# Pretty printing
###############################################################################

line() {
    printf '%*s\n' 80 '' | tr ' ' '='
}

stage() {
    echo
    line
    echo "$1"
    line
}

success() {
    echo
    echo "[SUCCESS] $1"
}

fail() {
    echo
    echo "[FAILED] $1"
    echo
    echo "Pipeline terminated."
    exit 1
}

###############################################################################
# Dependency checks
###############################################################################

stage "Checking required software..."

dependencies=(
    python3
    fasterq-dump
    cutadapt
    usearch
)

for program in "${dependencies[@]}"; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "[Missing] $program"
        exit 1
    fi
done

echo "[OK] All required command line software detected."

###############################################################################
# Stage 0
###############################################################################

stage "Stage 0 / 5 : SRA → FASTQ"

bash 00_sra_to_fastq.sh || fail "Stage 0 (SRA conversion)"

success "Stage 0 complete."

###############################################################################
# Stage 1
###############################################################################

stage "Stage 1 / 5 : Primer trimming"

python3 01_trim.py || fail "Stage 1 (Primer trimming)"

success "Stage 1 complete."

###############################################################################
# Stage 2
###############################################################################

stage "Stage 2 / 5 : Dereplication"

python3 02_dereplicate.py || fail "Stage 2 (Dereplication)"

success "Stage 2 complete."

###############################################################################
# Stage 3
###############################################################################

stage "Stage 3 / 5 : Denoising"

python3 03_denoise.py || fail "Stage 3 (UNOISE3)"

success "Stage 3 complete."

###############################################################################
# Stage 4
###############################################################################

stage "Stage 4 / 5 : BLAST annotation"

python3 04_blast.py || fail "Stage 4 (BLAST)"

success "Stage 4 complete."

###############################################################################
# Stage 5
###############################################################################

stage "Stage 5 / 5 : Taxonomy annotation"

python3 05_taxonomy.py || fail "Stage 5 (Taxonomy)"

success "Stage 5 complete."

###############################################################################
# Finished
###############################################################################

line
echo "PIPELINE COMPLETE"
line

echo
echo "All stages finished successfully."
echo
echo "Results are available in the newest folder under:"
echo
echo "    runs/"
echo