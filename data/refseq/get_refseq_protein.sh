#!/usr/bin/env bash

fasta_file_url="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.40_GRCh38.p14"
fasta_file_name_compressed="GCF_000001405.40_GRCh38.p14_protein.faa.gz"
fasta_file_name_uncompressed="GCF_000001405.40_GRCh38.p14_protein.faa"

if [[ ! -f "${fasta_file_name_compressed}" ]]; then
    wget ${fasta_file_url}/${fasta_file_name_compressed}
else
	echo "  ${fasta_file_name_compressed} already exists"
fi

mkdir -p indexed
if [[ ! -f "indexed/${fasta_file_name_uncompressed}" ]]; then
	cp GCF_000001405.40_GRCh38.p14_protein.faa.gz indexed/
	gunzip indexed/GCF_000001405.40_GRCh38.p14_protein.faa.gz
else
	echo "  ${fasta_file_name_uncompressed} already exists"
fi


