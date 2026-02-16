#!/usr/bin/env bash

versions=(109 110 111 112 113 114 115)

for v in "${versions[@]}"; do
    echo "Processing $v"

    outfile_cdna="Homo_sapiens.GRCh38.cdna.all.v${v}.fa.gz"
    if [[ ! -f "$outfile_cdna" ]]; then
        wget https://ftp.ensembl.org/pub/release-${v}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
        mv Homo_sapiens.GRCh38.cdna.all.fa.gz $outfile_cdna
    else
	echo "  $outfile_cdna already exists"
    fi 

    outfile_ncrna="Homo_sapiens.GRCh38.ncrna.v${v}.fa.gz"
    if [[ ! -f "$outfile_ncrna" ]]; then
        wget https://ftp.ensembl.org/pub/release-${v}/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
        mv Homo_sapiens.GRCh38.ncrna.fa.gz $outfile_ncrna
    else
        echo "  $outfile_ncrna already exists"
    fi

    outfile_pep="Homo_sapiens.GRCh38.pep.all.v${v}.fa.gz"
    if [[ ! -f "$outfile_pep" ]]; then
        wget https://ftp.ensembl.org/pub/release-${v}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
        mv Homo_sapiens.GRCh38.pep.all.fa.gz $outfile_pep
    else
        echo "  $outfile_pep already exists"
    fi

done

