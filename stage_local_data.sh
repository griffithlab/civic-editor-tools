#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

set -euo pipefail

run_step() {
    local description="$1"
    local script="$2"

    echo ""
    echo ">>> $description"
    "$script"
}

# Run each data download script one at a time. If files already exist, nothing will be downloaded

run_step "Ensembl protein fastas" ./data/ensembl/get_ensembl_fastas.sh
run_step "Ensembl build37 transcript annotation info" ./data/ensembl/get_ensembl_grch37_info.sh
run_step "Ensembl v75 version numbers" ./data/ensembl/get_ensembl_v75_version_numbers.py
run_step "Entrez gene to RefSeq mappings" ./data/entrez/get_refseq_mappings.sh
run_step "MANE select transcript info" ./data/refseq/get_mane_summary.sh
run_step "RefSeq protein fasta" ./data/refseq/get_refseq_protein.sh

echo ""
echo "All downloads complete."



