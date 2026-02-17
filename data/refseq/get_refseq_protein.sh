#!/usr/bin/env bash

wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_protein.faa.gz

mkdir indexed
cp GCF_000001405.40_GRCh38.p14_protein.faa.gz indexed/
gunzip indexed/GCF_000001405.40_GRCh38.p14_protein.faa.gz

