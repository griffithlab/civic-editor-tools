#!/usr/bin/env python3

""" 
backfill_refseq_transcript_info.py

Helper script to gather RefSeq transcript information for old transcripts in ClinGen for genes in CIViC.

Other tools may need this transcript info (refseq transcript to protein ids mappings)
   
"""     
        
import os
import sys
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

def main():

    #include_list=None #to use every single variant regardless of status
    include_list = ['accepted', 'submitted']

    #use civicpy to get all genes in civic
    genes = civic.get_all_genes(include_status=include_list, allow_cached=True)
    print(f"Total genes in CIViC: {len(genes)}\n")

    #get mappings of transcript to protein identifiers for refseq transcripts
    refseq_to_protein_file = base_dir / f"data/entrez/gene2refseq_human.tsv.gz"
    refseq_to_protein_missing_file = base_dir / f"data/entrez/gene2refseq_human_missing.tsv"
    refseq_fasta_index_path = base_dir / f"data/refseq/indexed/GCF_000001405.40_GRCh38.p14_protein.faa.idx"

    refseq_transcript_to_protein_map = entrez_utils.load_refseq_transcript_to_protein_map(refseq_to_protein_file, refseq_to_protein_missing_file)

    for g in genes:
        gene_name = g.name
        print(f"Working on gene: {gene_name}")
        
        #create a data structure that will store all clingen supported transcripts ids per gene
        clingen_transcript_ids = {}
        
        #for each gene, get transcripts from clingen allele registry
        clingen_gene_transcripts_json = clingen_ar_utils.get_reference_sequences_by_gene(gene_name)

        #get the reference sequence ids associated with this gene
        clingen_reference_sequence_ids = clingen_ar_utils.extract_reference_sequences(clingen_gene_transcripts_json)


        for clingen_reference_sequence_id in clingen_reference_sequence_ids:
            #skip ensembl transcripts
            if clingen_reference_sequence_id.startswith("ENST"):
                continue
    
            #skip invalid refseq transcripts
            if clingen_reference_sequence_id.startswith(("NR_", "XM_", "XR_")):
                continue
        
            #get the protein ID for the current transcript id
            protein_id = None
            protein_seq = None
            if clingen_reference_sequence_id in refseq_transcript_to_protein_map:
                protein_id = refseq_transcript_to_protein_map[clingen_reference_sequence_id]
                #make sure the protein sequence is available
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

                    # Reload the map so the freshly fetched entry is available
                    refseq_transcript_to_protein_map = entrez_utils.load_refseq_transcript_to_protein_map(
                        refseq_to_protein_file, refseq_to_protein_missing_file
                    )

                    # Retry the lookup now that the map has been refreshed
                    if clingen_reference_sequence_id in refseq_transcript_to_protein_map:
                        protein_id = refseq_transcript_to_protein_map[clingen_reference_sequence_id]
                        protein_seq = refseq_utils.get_refseq_protein_indexed(
                            protein_id, refseq_fasta_index_path
                        )
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
                    # Always pause after any efetch call, success or failure
                    time.sleep(NCBI_REQUEST_DELAY)


if __name__ == "__main__":
    main()





