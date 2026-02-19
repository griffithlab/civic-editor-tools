#!/usr/bin/env python3

"""
review_variant_coordinates.py

Command-line tool to review and validate genomic coordinates for a CIViC variant.

This script retrieves variant data from CIViC, cross-references external
annotation sources (e.g., ClinGen Allele Registry, NCBI Entrez), and supports
review workflows tied to a specific CIViC contributor.
"""

import os
import argparse
import sys
from civicpy import civic
from pathlib import Path

# Local utility imports
from utils import generic_utils
from utils import civic_graphql_utils
from utils import civicpy_utils
from utils import clingen_ar_utils
from utils import entrez_utils
from utils import ensembl_utils
from utils import refseq_utils

base_dir = Path(__file__).resolve().parent

def parse_args():
    """Obtain command line arguments from the user"""
    parser = argparse.ArgumentParser(
        description=(
            "Review pending genomic coordinate revisions for one or more CIViC variant using "
            "multiple external annotation sources."
        )
    )

    parser.add_argument(
        "--contributor-id",
        dest="contributor_id",
        type=int,
        required=True,
        help="CIViC contributor ID performing the review (integer, e.g. 15)"
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--variant-id",
        dest="variant_id",
        type=int,
        help="CIViC variant ID to review (integer, e.g. 1832)"
    )
    group.add_argument(
        "--all-variants",
        dest="all_variants",
        action="store_true",
        help="Review all CIViC variants"
    )

    return parser.parse_args()

def get_variant_ids_to_process(variant_id, all_variants):
    """Determine which variant IDs to work with based on user supplied choices """
    variant_ids_to_process = []
    if variant_id:
        variant_ids_to_process.append(variant_id)
    
    if all_variants:
        include_list = ['accepted', 'submitted']
        variants = civic.get_all_gene_variants(include_status=include_list, allow_cached=True)
        variant_ids_to_process = civicpy_utils.extract_variant_id_list(variants)
        print(f"Total variant ids obtained from CIViCpy: {len(variant_ids_to_process)}\n")

    return variant_ids_to_process

def variant_is_black_listed(vid, black_listed_variant_ids, black_list_path, contributor_id):
    """Skip variants that are in a manually maintainted black list file"""
    if vid in black_listed_variant_ids:
        print(f"\nSkipping CIViC variant {vid} because it was found in the blacklist: {black_list_path}")
        return True

    return False


def variant_is_deprecated(vid, variant_data_basic):
    """Skip variants that have a deprecated status"""
    if (variant_data_basic['deprecated'] == "True"):
        print(f"Skipping CIViC variant {vid} because it has a deprecated status")
        return True
    
    return False


def variant_type_is_unsupported(civic_variant_name, guessed_gene_variant_type, target_variant_type):
    """Skip variants that do not have the expected variant type guessed from the name"""
    if (guessed_gene_variant_type != target_variant_type):
        print(f"Guessed variant type for: {civic_variant_name!r} -> {guessed_gene_variant_type} is not supported - skipping")
        return True
    else:
        print(f"Guessed variant type for: {civic_variant_name!r} -> {guessed_gene_variant_type} is supported")

    return False


def main(variant_id: int, contributor_id: int, all_variants: bool):

    #define input data files
    black_list_path = base_dir / f"data/civic_variant_blacklist.tsv"
    refseq_fasta_index_path = base_dir / f"data/refseq/indexed/GCF_000001405.40_GRCh38.p14_protein.faa.idx"
    ensembl_versions_file = base_dir / f"data/ensembl/ensembl_versions.txt"
    refseq_to_protein_file = base_dir / f"data/entrez/gene2refseq_human.tsv.gz"

    #load black listed variants file
    black_listed_variant_ids = civic_graphql_utils.load_blacklisted_variant_ids(black_list_path)
    
    #get mappings of transcript to protein identifiers for refseq transcripts
    refseq_transcript_to_protein_map = entrez_utils.load_refseq_transcript_to_protein_map(refseq_to_protein_file)

    #get mappings of transcript to protein identifier for ensembl transcripts
    ensembl_transcript_to_protein_map = ensembl_utils.load_ensembl_transcript_to_protein_map(ensembl_versions_file)

    #get ensembl transcript biotype to transcript identifiers
    ensembl_transcript_to_biotype_map = ensembl_utils.load_ensembl_transcript_to_biotype_map(ensembl_versions_file)

    #create a data structure that will store all clingen supported transcripts ids per gene
    clingen_transcript_ids = {}

    #get civic variant IDs to evaluate, either from the user, or by querying CIViCpy
    variant_ids_to_process = get_variant_ids_to_process(variant_id, all_variants)

    #iterate over each variant and examine revisions associated with it
    for vid in variant_ids_to_process:

        print(f"\nReviewing CIViC variant {vid} for revisions that could be reviewed by contributor: {contributor_id}")

        #skip if this variant is black listed
        if variant_is_black_listed(vid, black_listed_variant_ids, black_list_path, contributor_id): continue

        #query the graphql api for basic variant info
        variant_data_basic = civic_graphql_utils.gather_variant_details(vid)

        #get the civic variant name - much will be assumed based on this name
        civic_variant_name = variant_data_basic['variant_name']

        #skip variants with deprecated status
        if variant_is_deprecated(vid, variant_data_basic): continue

        #guess the variant type based on the CIViC variant name - skip unless it is a coding snv
        guessed_gene_variant_type = generic_utils.guess_variant_type(civic_variant_name)

        if variant_type_is_unsupported(civic_variant_name, guessed_gene_variant_type, "snv_coding"): continue

        #get the three components of a simple snv_coding variant: ref_aa_1, pos, var_aa_1
        ref_aa_1, pos, var_aa_1 = generic_utils.parse_snv_coding_name_components(civic_variant_name)
        print(f"{civic_variant_name} -> ref_aa_1: {ref_aa_1} ->  pos: {pos} -> var_aa_1: {var_aa_1}")

        #query the graphql api for more detailed variant revision info
        variant_data = civic_graphql_utils.gather_variant_revisions(vid, contributor_id)

        print(
            f"Variant revision info from gather_variant_revisions()\n"
            f"Variant ID used for graphql query: {variant_data['variant_id']}\n"
            f"  Variant name: {variant_data['variant_name']}\n"
            f"  Feature name: {variant_data['feature_name']}\n"
            f"  Open gene-variant revisions (total): {variant_data['open_revision_count_variant']}\n"
            f"  Open gene-variant revisions from specified contributor: {variant_data['contributor_revisions']}\n"
            f"  Open gene-variant revisions from all others users: {variant_data['open_revisions_non_contributor']}\n"
            f"  Variant coordinates id: {variant_data['variant_coordinates_id']}"
        )

        variant_revisions = variant_data['variant_revisions']

        #if there is an outstanding revision to the variant name itself, warn the user
        if variant_data['name_change']:
            print(f"WARNING. The variant name itself has an outstanding revision!")
            print(f"  Since this entire exercise derives from that name, this must be resolved first\n")

        #create the p. notation for the variant name (e.g. 'S459F' -> 'p.Ser459Phe')
        civic_variant_name_p_3letter = generic_utils.snv_coding_to_p_3letter(civic_variant_name)
        print(f"\nVariant name in p. notation: {civic_variant_name_p_3letter}")

        #skip a variant if it has 0 pending revision from other users
        if variant_data['open_revisions_non_contributor'] == 0:
            print(f"No open revision for this variant - skipping")
            continue

        #get all clingen allele registry transcripts supported for the gene of this variant
        #only query the clingen API if we don't already have transcripts for this gene
        gene_name = variant_data['feature_name']
        clingen_gene_transcript_ids = clingen_transcript_ids.get(gene_name)

        if clingen_gene_transcript_ids is None:
            clingen_gene_transcript_ids = clingen_ar_utils.get_reference_sequences_by_gene(gene_name)
            clingen_transcript_ids[gene_name] = clingen_gene_transcript_ids
        
        clingen_reference_sequence_ids = clingen_ar_utils.extract_reference_sequences(clingen_gene_transcript_ids)

        #when there are multiple versions of the same transcript, keep only the latest one
        clingen_reference_sequence_ids_latest = clingen_ar_utils.keep_latest_transcript_versions(clingen_reference_sequence_ids)

        #go through each transcript ID from clingen allele registry and see if it is worth checking for the current CIViC variant name
        for clingen_reference_sequence_id in clingen_reference_sequence_ids_latest:
            #skip invalid ensembl transcripts
            if clingen_reference_sequence_id.startswith("ENST"):
                #if the transcript ID is NOT in the biotype map, assume it is old, skip it
                if clingen_reference_sequence_id not in ensembl_transcript_to_biotype_map:
                    continue
                #if the transcript ID is an ensembl ID and its biotype is NOT protein_coding, skip it
                if ensembl_transcript_to_biotype_map[clingen_reference_sequence_id] != "protein_coding":
                    continue

            #get the protein ID for the current transcript id
            protein_id = None
            protein_seq = None
            if clingen_reference_sequence_id in refseq_transcript_to_protein_map:
                protein_id = refseq_transcript_to_protein_map[clingen_reference_sequence_id]
                protein_seq = refseq_utils.get_refseq_protein_indexed(protein_id, refseq_fasta_index_path)

            elif clingen_reference_sequence_id in ensembl_transcript_to_protein_map:
                protein_id = ensembl_transcript_to_protein_map[clingen_reference_sequence_id]
                protein_seq = ensembl_utils.get_ensembl_protein_indexed(protein_id, ensembl_versions_file)

            else:
                raise ValueError(
                    f"Transcript ID {clingen_reference_sequence_id} not found in "
                    "RefSeq or Ensembl transcript-to-protein maps"
                )

            #unless the protein sequence has the expected reference amino acid at the expected position, skip it
            if not generic_utils.reference_aa_positions_matches(ref_aa_1, pos, protein_seq):
                continue

            print(f"  transcript id: {clingen_reference_sequence_id} -> protein_id: {protein_id} -> {protein_id}:{civic_variant_name_p_3letter}")

        #- Variant ambiguity check (consider an example variant "BRAF V600E"
        #  - For a given gene get all transcripts (RefSeq and Ensembl) in ClinGen Allele Registry (CAR)
        #  - Check which of these transcripts have the expected ref AA at the expected position
        #  - Starting from variant name, contruct possible p. hgvs expressions for all RefSeq and Ensembl transcript in CAR
        #  - Check each of these p. hgvs expressions and get PAIDs. Get all CAIDs associated with these
        #  - Skip CAIDs that are not a simple SNV?
        #  - Get the g. HGVS expression associated with all remaining CAIDs (make not of the MANE select)
        #  - Are there multiple distinct g. HGVS values that the variant name could refer to? If so, warn the user

if __name__ == "__main__":
    args = parse_args()

    main(
        contributor_id=args.contributor_id,
        variant_id=args.variant_id,
        all_variants=args.all_variants
    )

