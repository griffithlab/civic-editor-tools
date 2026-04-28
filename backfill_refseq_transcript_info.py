#!/usr/bin/env python3

""" 
backfill_refseq_transcript_info.py

Helper script to gather RefSeq transcript information for old transcripts in ClinGen for genes in CIViC.

Other tools may need this transcript info (refseq transcript to protein ids mappings)
   
"""     
        
import os
import sys
import argparse
from civicpy import civic
from pathlib import Path
import subprocess
import time
    
# Local utility imports
from utils import generic_utils
from utils import civicpy_utils
from utils import clingen_ar_utils
from utils import entrez_utils
from utils import refseq_utils

NCBI_REQUEST_DELAY = 0.5  # seconds between efetch calls

base_dir = Path(__file__).resolve().parent


def load_processed_genes(checkpoint_file: Path) -> set:
    """Load the set of already-processed gene names from the checkpoint file."""
    if not checkpoint_file.exists():
        return set()
    with checkpoint_file.open() as f:
        return {line.strip() for line in f if line.strip()}


def mark_gene_processed(checkpoint_file: Path, gene_name: str) -> None:
    """Append a gene name to the checkpoint file."""
    with checkpoint_file.open("a") as f:
        f.write(gene_name + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Backfill RefSeq transcript info for CIViC genes via ClinGen."
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=None,
        metavar="FILE",
        help=(
            "Optional path to a checkpoint file. Genes listed in this file will be "
            "skipped. Each newly completed gene is appended to the file so the run "
            "can be resumed if interrupted."
        ),
    )
    args = parser.parse_args()

    include_list = ['accepted', 'submitted']

    genes = civic.get_all_genes(include_status=include_list, allow_cached=True)
    print(f"Total genes in CIViC: {len(genes)}\n")

    # Load checkpoint if provided
    processed_genes = set()
    if args.checkpoint:
        processed_genes = load_processed_genes(args.checkpoint)
        if processed_genes:
            print(f"Checkpoint loaded from '{args.checkpoint}': "
                  f"skipping {len(processed_genes)} already-processed gene(s).\n")

    refseq_to_protein_file = base_dir / "data/entrez/gene2refseq_human.tsv.gz"
    refseq_to_protein_missing_file = base_dir / "data/entrez/gene2refseq_human_missing.tsv"
    refseq_fasta_index_path = base_dir / "data/refseq/indexed/merged.faa.idx"

    refseq_transcript_to_protein_map = entrez_utils.load_refseq_transcript_to_protein_map(
        refseq_to_protein_file, refseq_to_protein_missing_file
    )

    for g in genes:
        gene_name = g.name

        if gene_name in processed_genes:
            print(f"Skipping already-processed gene: {gene_name}")
            continue

        print(f"Working on gene: {gene_name}")
        
        clingen_transcript_ids = {}
        
        clingen_gene_transcripts_json = clingen_ar_utils.get_reference_sequences_by_gene(gene_name)

        if not clingen_gene_transcripts_json:
            continue

        clingen_reference_sequence_ids = clingen_ar_utils.extract_reference_sequences(clingen_gene_transcripts_json)

        for clingen_reference_sequence_id in clingen_reference_sequence_ids:
            if clingen_reference_sequence_id.startswith("ENST"):
                continue
    
            if clingen_reference_sequence_id.startswith(("NR_", "XM_", "XR_")):
                continue

            print(f"  ClinGen Reference Sequence ID: {clingen_reference_sequence_id}")

            protein_id = None
            protein_seq = None
            if clingen_reference_sequence_id in refseq_transcript_to_protein_map:
                protein_id = refseq_transcript_to_protein_map[clingen_reference_sequence_id]
                print(f"    Protein ID: {protein_id}")
                protein_seq = refseq_utils.get_refseq_protein_indexed(protein_id, refseq_fasta_index_path)
            else:
                missing_refseq_mappings_script = base_dir / "data/entrez/get_missing_refseq_mappings.py"
                print(f"Transcript {clingen_reference_sequence_id} not in map — fetching via NCBI...")

                try:
                    result = subprocess.run(
                        [sys.executable, str(missing_refseq_mappings_script), clingen_reference_sequence_id],
                        capture_output=True,
                        text=True,
                        check=True,
                    )
                    print(result.stdout)

                    refseq_transcript_to_protein_map = entrez_utils.load_refseq_transcript_to_protein_map(
                        refseq_to_protein_file, refseq_to_protein_missing_file
                    )

                    if clingen_reference_sequence_id in refseq_transcript_to_protein_map:
                        protein_id = refseq_transcript_to_protein_map[clingen_reference_sequence_id]
                        protein_seq = refseq_utils.get_refseq_protein_indexed(protein_id, refseq_fasta_index_path)
                    else:
                        print(
                            f"Warning: {clingen_reference_sequence_id} still missing after fetch — skipping.",
                            file=sys.stderr,
                        )

                except subprocess.CalledProcessError as e:
                    print(
                        f"Error fetching {clingen_reference_sequence_id}:\n{e.stderr}",
                        file=sys.stderr,
                    )
                finally:
                    time.sleep(NCBI_REQUEST_DELAY)

        # Gene fully processed — record it in the checkpoint file if one was specified
        if args.checkpoint:
            mark_gene_processed(args.checkpoint, gene_name)
            processed_genes.add(gene_name)


if __name__ == "__main__":
    main()

