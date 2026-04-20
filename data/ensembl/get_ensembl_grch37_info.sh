#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p build37

#Get the last version of Ensembl annotations that were actually predicted for build37 (v75)
GTF_URL="https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz"
GTF_FILE=$(basename "$GTF_URL")
GTF_FILE_OUT="build37/${GTF_FILE}"

if [[ ! -f "$GTF_FILE_OUT" ]]; then
    wget $GTF_URL -O $GTF_FILE_OUT
else
    echo "$GTF_FILE_OUT already exists"
fi

#Get the most updated transfer of ensembl annotations onto build37
#Even though this says v115 of ensembl, the "87" means that was the version of annotations imported to build37
#https://ftp.ensembl.org/pub/grch37/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz

GTF_URL="https://ftp.ensembl.org/pub/grch37/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
GTF_FILE=$(basename "$GTF_URL")
GTF_FILE_OUT="build37/${GTF_FILE}"

if [[ ! -f "$GTF_FILE_OUT" ]]; then
    wget $GTF_URL -O $GTF_FILE_OUT
else
    echo "$GTF_FILE_OUT already exists"
fi


