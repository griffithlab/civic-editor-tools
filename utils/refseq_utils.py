#!/usr/bin/env python3

import sys
import os
import gzip
import pickle
from Bio import SeqIO
from pathlib import Path

base_dir = Path(__file__).resolve().parent

def build_refseq_fasta_index(refseq_fasta_path, index_path):
    """Creates a persistent index file (.idx)"""

    for path in [index_path, 
                 index_path.with_suffix(".idx-wal"), 
                 index_path.with_suffix(".idx-shm")]:
        if path.exists():
            path.unlink()

    print(f"Building index for refseq fasta...\n{index_path}")
    SeqIO.index_db(str(index_path), str(refseq_fasta_path), "fasta")


def get_refseq_protein_indexed(refseq_protein_id, index_path):
    """Given a refseq protein id and a prebuilt SeqIO index, retrieve the protein sequence from a fasta"""
    
    index = SeqIO.index_db(str(index_path))
    
    if refseq_protein_id not in index:
        error_message = (
            f"{refseq_protein_id} not found in protein sequence index. Do the following two commands to get the missing fasta sequence and rebuild fasta the index:\n"
            f"  ./data/refseq/get_missing_refseq_protein_fasta.py {refseq_protein_id}\n"
            f"  ./utils/refseq_utils.py"
        )
        sys.exit(error_message)

    return str(index[refseq_protein_id].seq)


def merge_refseq_fastas(refseq_fasta_path: Path, missing_refseq_fasta_path: Path, merged_refseq_fasta_path: Path) -> None:
    """
    Merge two RefSeq FASTA files into one, failing if any accession appears in both inputs.
    """

    # Load both files into dicts keyed by accession
    print(f"Reading {refseq_fasta_path} ...")
    primary = SeqIO.to_dict(SeqIO.parse(refseq_fasta_path, "fasta"))
    print(f"  {len(primary):,} records")

    print(f"Reading {missing_refseq_fasta_path} ...")
    supplemental = {}
    internal_duplicates = []
    for record in SeqIO.parse(missing_refseq_fasta_path, "fasta"):
        if record.id in supplemental:
            internal_duplicates.append(record.id)
        else:
            supplemental[record.id] = record

    if internal_duplicates:
        dup_list = "\n  ".join(sorted(internal_duplicates))
        raise ValueError(
            f"{len(internal_duplicates)} duplicate accession(s) found within "
            f"{missing_refseq_fasta_path.name} — merge aborted:\n  {dup_list}"
        )
    print(f"  {len(supplemental):,} records (no internal duplicates)")

    # Check for duplicates before touching the output file
    duplicates = primary.keys() & supplemental.keys()
    if duplicates:
        dup_list = "\n  ".join(sorted(duplicates))
        raise ValueError(
            f"{len(duplicates)} duplicate accession(s) found across both FASTA files "
            f"— merge aborted:\n  {dup_list}"
        )

    # Write merged output
    merged_refseq_fasta_path.parent.mkdir(parents=True, exist_ok=True)
    merged = {**primary, **supplemental}
    print(f"Writing {len(merged):,} records → {merged_refseq_fasta_path} ...")
    written = SeqIO.write(merged.values(), merged_refseq_fasta_path, "fasta")
    print(f"  {written:,} records written")


def main():

    test_refseq_protein_id = "NP_006222.2"

    # build the index
    refseq_fasta_path = base_dir / f"../data/refseq/indexed/GCF_000001405.40_GRCh38.p14_protein.faa"
    missing_refseq_fasta_path = base_dir / f"../data/refseq/indexed/missing_refseqs_protein.faa"
    merged_refseq_fasta_path = base_dir / f"../data/refseq/indexed/merged.faa"
    merged_refseq_fasta_index_path = merged_refseq_fasta_path.with_suffix(".faa.idx")

    # combine the baseline and missing refseq fasta files together
    merge_refseq_fastas(refseq_fasta_path, missing_refseq_fasta_path, merged_refseq_fasta_path)   

    # build a fasta index for the combined fasta with both baseline and added fasta records
    build_refseq_fasta_index(merged_refseq_fasta_path, merged_refseq_fasta_index_path)

    seq = get_refseq_protein_indexed(test_refseq_protein_id, merged_refseq_fasta_index_path)
    
    print(f"Found sequence of {test_refseq_protein_id} with length: {len(seq)}:")
    print(seq)


if __name__ == "__main__":
    main()

