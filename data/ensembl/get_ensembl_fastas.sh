#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

versions=()

while IFS= read -r v; do
    versions+=("$v")
done < ensembl_versions.txt

mkdir -p version_data/indexed

for v in "${versions[@]}"; do
    echo "Processing $v"

    outfile_cdna="version_data/Homo_sapiens.GRCh38.cdna.all.v${v}.fa.gz"
    if [[ ! -f "$outfile_cdna" ]]; then
        wget https://ftp.ensembl.org/pub/release-${v}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O $outfile_cdna
    else
	    echo "  $outfile_cdna already exists"
    fi 

    outfile_ncrna="version_data/Homo_sapiens.GRCh38.ncrna.v${v}.fa.gz"
    if [[ ! -f "$outfile_ncrna" ]]; then
        wget https://ftp.ensembl.org/pub/release-${v}/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz -O $outfile_ncrna
    else
        echo "  $outfile_ncrna already exists"
    fi

    outfile_pep="Homo_sapiens.GRCh38.pep.all.v${v}.fa.gz"
    if [[ ! -f "version_data/${outfile_pep}" ]]; then
        wget https://ftp.ensembl.org/pub/release-${v}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz -O version_data/$outfile_pep
    else
        echo "  version_data/$outfile_pep already exists"
        outfile_pep2="Homo_sapiens.GRCh38.pep.all.v${v}.fa"
        if [[ ! -f "version_data/indexed/${outfile_pep2}" ]]; then
            cp version_data/$outfile_pep version_data/indexed/
	        gunzip version_data/indexed/$outfile_pep
        else
            echo "  version_data/indexed/$outfile_pep2 already exists"
        fi
    fi

done

