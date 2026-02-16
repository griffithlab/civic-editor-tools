#!/usr/bin/env python3

import os
import gzip
import pickle

def get_refseq_protein(refseq_id, refseq_fasta_path):
    """
    Parse a gzipped RefSeq FASTA file and return the protein sequence
    for an exact RefSeq ID match (e.g. 'NM_006231.4').
    """

    if not refseq_id:
        raise ValueError("refseq_id must be provided")

    found = False
    sequence_lines = []

    with gzip.open(refseq_fasta_path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                header_id = line[1:].split()[0]

                # If we were already collecting and hit a new header, stop
                if found:
                    break

                if header_id == refseq_id:
                    found = True
                    continue
            else:
                if found:
                    sequence_lines.append(line.strip())

    if not found:
        raise ValueError(f"RefSeq ID '{refseq_id}' not found in {refseq_fasta_path}")

    return "".join(sequence_lines)


def main():

    refseq_fasta_path = "data/refseq/GCF_000001405.40_GRCh38.p14_protein.faa.gz"
    test_refseq_protein_id = "NP_006222.2"

    seq = get_refseq_protein(test_refseq_protein_id, refseq_fasta_path)
    
    print(f"Found sequence of {test_refseq_protein_id} with length: {len(seq)}:")
    print(seq)


if __name__ == "__main__":
    main()




