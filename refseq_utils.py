#!/usr/bin/env python3

import os
import gzip
import pickle
from Bio import SeqIO

def build_refseq_fasta_index(refseq_fasta_path, index_path):
    """Creates a persistent index file (.idx)"""

    SeqIO.index_db(index_path, refseq_fasta_path, "fasta")


def get_refseq_protein_indexed(refseq_protein_id, index_path):
    """Given a refseq protein id and a prebuilt SeqIO index, retrieve the protein sequence"""
    index = SeqIO.index_db(index_path)
    
    if refseq_protein_id not in index:
        raise ValueError(f"{refseq_protein_id} not found in index")

    return str(index[refseq_protein_id].seq)


def main():

    test_refseq_protein_id = "NP_006222.2"

    # build the index
    refseq_fasta_path = "data/refseq/indexed/GCF_000001405.40_GRCh38.p14_protein.faa"
    refseq_fasta_index_path = "data/refseq/indexed/GCF_000001405.40_GRCh38.p14_protein.faa.idx"
    build_refseq_fasta_index(refseq_fasta_path, refseq_fasta_index_path)

    # Later
    seq = get_refseq_protein_indexed(test_refseq_protein_id, refseq_fasta_index_path)
    
    print(f"Found sequence of {test_refseq_protein_id} with length: {len(seq)}:")
    print(seq)


if __name__ == "__main__":
    main()




