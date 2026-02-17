#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

entrez_file_url="ftp://ftp.ncbi.nih.gov/gene/DATA"
entrez_file_name_compressed="gene2refseq.gz"
entrez_file_name_final="gene2refseq_human.tsv.gz"

if [[ ! -f "${entrez_file_name_final}" ]]; then
	wget ${entrez_file_url}/${entrez_file_name_compressed}
	gzcat $entrez_file_name_compressed | head -n 1 > header.tsv
	gzcat $entrez_file_name_compressed | grep '^9606\t' > gene2refseq_human.tmp.tsv
	cat header.tsv gene2refseq_human.tmp.tsv > gene2refseq_human.tsv
	rm -f header.tsv gene2refseq_human.tmp.tsv $entrez_file_name_compressed
	gzip gene2refseq_human.tsv
else
	echo "  ${entrez_file_name_final} already exists"
fi

