#!/usr/bin/env python3

"""
review_variant_coordinates.py

Command-line tool to review and validate genomic coordinates for a CIViC variant.

This script retrieves variant data from CIViC, cross-references external
annotation sources (e.g., ClinGen Allele Registry, NCBI Entrez), and supports
review workflows tied to a specific CIViC contributor.
"""

import argparse
import sys
from civicpy import civic

# Local utility imports
import generic_utils
import civic_graphql_utils
import civicpy_utils
import clingen_ar_utils
import entrez_utils

def parse_args():
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


def main(variant_id: int, contributor_id: int, all_variants: bool):

    #load black listed variants file
    black_list_path = "data/civic_variant_blacklist.tsv"
    black_listed_variant_ids = civic_graphql_utils.load_blacklisted_variant_ids(black_list_path)

    clingen_transcript_ids = {}

    variant_ids_to_process = []
    if variant_id:
        variant_ids_to_process.append(variant_id)
    
    if all_variants:
        include_list = ['accepted', 'submitted']
        variants = civic.get_all_gene_variants(include_status=include_list, allow_cached=True)
        variant_ids_to_process = civicpy_utils.extract_variant_id_list(variants)
        print(f"Total variant ids obtained from CIViCpy: {len(variant_ids_to_process)}\n")

    for vid in variant_ids_to_process:
        #skip if this variant is black listed
        if vid in black_listed_variant_ids:
            print(f"\nSkipping CIViC variant {vid} because it was found in the blacklist: {black_list_path}")
            continue

        print(f"\nReviewing CIViC variant {vid} for revisions that could be reviewed by contributor: {contributor_id}")

        #query the graphql api for basic variant info
        variant_data_basic = civic_graphql_utils.gather_variant_details(vid)

        civic_variant_name = variant_data_basic['variant_name']

        #skip variants with deprecated status
        if (variant_data_basic['deprecated'] == "True"):
            print(f"Skipping CIViC variant {vid} because it has a deprecated status")
            continue

        #guess the variant type based on the CIViC variant name - skip unless it is a coding snv
        guessed_gene_variant_type = generic_utils.guess_variant_type(civic_variant_name)
    
        if (guessed_gene_variant_type == "snv_coding"):
            print(f"Guessed variant type for: {civic_variant_name!r} -> {guessed_gene_variant_type}")
        else:
            print(f"Guessed variant type for: {civic_variant_name!r} -> {guessed_gene_variant_type} is not supported here - skipping")
            continue

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
        print(f"Variant name in p. notation: {civic_variant_name_p_3letter}")

        #skip a variant if it has 0 pending revision from other users
        if variant_data['open_revisions_non_contributor'] == 0:
            print(f"No open revision for this variant - skipping")
            continue

        #get all clingen allele registry transcripts supported for the gene of this variant
        gene_name = variant_data['feature_name']
        clingen_gene_transcript_ids = clingen_transcript_ids.get(gene_name)

        if clingen_gene_transcript_ids is None:
            clingen_gene_transcript_ids = clingen_ar_utils.get_reference_sequences_by_gene(gene_name)
            clingen_transcript_ids[gene_name] = clingen_gene_transcript_ids

    
        #- Variant ambiguity check (consider an example variant "BRAF V600E"
        #  - For a given gene get all transcripts (RefSeq and Ensembl) in ClinGen Allele Registry (CAR)
        #  - Check which of these transcripts have the expected ref AA at the expected position
        #  - Starting from the name of a variant, contruct possible p. hgvs expressions for all RefSeq and Ensembl transcript is CAR
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

