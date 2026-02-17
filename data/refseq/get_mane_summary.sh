#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

mane_file_url="https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5"
mane_file_name_compressed="MANE.GRCh38.v1.5.summary.txt.gz"
mane_file_name_uncompressed="MANE.GRCh38.v1.5.summary.txt"

if [[ ! -f "${mane_file_name_uncompressed}" ]]; then
    wget ${mane_file_url}/${mane_file_name_compressed}
    gunzip ${mane_file_name_uncompressed}
else
	echo "  ${mane_file_name_uncompressed} already exists"
fi

