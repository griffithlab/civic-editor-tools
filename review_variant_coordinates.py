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

# Local utility imports
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
        "--variant-id",
        dest="variant_id",
        type=int,
        required=True,
        help="CIViC variant ID to review (integer)"
    )

    parser.add_argument(
        "--contributor-id",
        dest="contributor_id",
        type=int,
        required=True,
        help="CIViC contributor ID performing the review (integer)"
    )

    return parser.parse_args()


def main(variant_id: int, contributor_id: int):

    print(f"Reviewing CIViC variant {variant_id} for revisions that could be reviewed by contributor: {contributor_id}")

    variant_data = civic_graphql_utils.gather_variant_revisions(variant_id, contributor_id)

    print(
        f"\nVariant revision info from gather_variant_revisions()\n"
        f"Variant ID used for graphql query: {variant_data['variant_id']}\n"
        f"  Variant name: {variant_data['variant_name']}\n"
        f"  Feature name: {variant_data['feature_name']}\n"
        f"  Open gene-variant revisions (total): {variant_data['open_revision_count_variant']}\n"
        f"  Open gene-variant revisions from specified contributor: {variant_data['contributor_revisions']}\n"
        f"  Open gene-variant revisions from all others users: {variant_data['open_revisions_non_contributor']}\n"
        f"  Variant coordinates id: {variant_data['variant_coordinates_id']}"
    )



if __name__ == "__main__":
    args = parse_args()

    main(
        variant_id=args.variant_id,
        contributor_id=args.contributor_id
    )

