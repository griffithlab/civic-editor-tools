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
from utils import compare_utils

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
    group.add_argument(
        "--allow-variants-without-revisions",
        dest="allow_variants_without_revisions",
        action="store_true",
        help="Even if a variant has no outstanding revisions, process it anyway"
    )

    return parser.parse_args()

def verify_connectivity():
    """Tests whether internet access is working and then if we can actually access the APIs needed before attempting anything """
    if not generic_utils.check_connection():
        print("No internet access. Aborting.")
        sys.exit(1)

    if not generic_utils.check_apis():
        print("Required APIs are unavailable. Aborting.")
        sys.exit(1)

    print("Internet and API connectivity verified.")


def prompt_to_proceed(message: str = None) -> None:
    print("\n" + "=" * 80)
    if message:
        print(message)
    print("Press Enter to proceed or Ctrl+C to cancel...")
    print("=" * 80)
    try:
        input()
    except KeyboardInterrupt:
        print("\nAborted.")
        sys.exit(0)


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
    """Test if variant is in a manually maintainted black list file"""
    if vid in black_listed_variant_ids:
        print(f"\nSkipping CIViC variant {vid} because it was found in the blacklist: {black_list_path}")
        return True

    return False


def variant_is_deprecated(vid, variant_data_basic):
    """Test if variant has a deprecated status"""
    if (variant_data_basic['deprecated'] == "True"):
        print(f"Skipping CIViC variant {vid} because it has a deprecated status")
        return True
    
    return False


def variant_type_is_unsupported(civic_variant_name, guessed_gene_variant_type, target_variant_type):
    """Test if variant has the expected variant type guessed from the name"""
    if (guessed_gene_variant_type != target_variant_type):
        print(f"Guessed variant type for variant name {civic_variant_name}: {guessed_gene_variant_type} (not supported - skipping)")
        return True
    else:
        print(f"Guessed variant type for variant name {civic_variant_name}: {guessed_gene_variant_type} (is supported)")

    return False


def variant_name_has_revision(variant_data):
    """Test is variant has a pending revision on the variant name itself"""
    if variant_data['name_change']:
        print(f"WARNING. The variant name itself has an outstanding revision!")
        print(f"  Since this entire exercise derives from that name, this should be resolved first. Skipping this variant\n")
        return True

    return False


def variant_has_no_open_revisions(variant_data, allow_variants_without_revisions):
    """Test is variants has no open revisions """
    if variant_data['open_revisions_non_contributor'] == 0 and not allow_variants_without_revisions:
        print(f"No open revision(s) by other contributors for this variant - skipping")
        return True
    
    return False


def get_clingen_gene_transcripts_json (gene_name, clingen_transcript_ids):
    """Get a return a json object of clingen transcripts for a gene, store them so that the gene API will only be queried once"""
    clingen_gene_transcripts_json = clingen_transcript_ids.get(gene_name)

    if clingen_gene_transcripts_json is None:
        clingen_gene_transcripts_json = clingen_ar_utils.get_reference_sequences_by_gene(gene_name)
        clingen_transcript_ids[gene_name] = clingen_gene_transcripts_json

    return clingen_gene_transcripts_json


def get_compatible_clingen_transcripts(clingen_gene_transcripts_json, refseq_transcript_to_protein_map, ensembl_transcript_to_protein_map, ensembl_transcript_to_biotype_map, refseq_fasta_index_path, ensembl_versions_file, ref_aa_1, pos, var_aa_1, civic_variant_name_p_3letter):
    """"Filter clingen allele registry transcripts to those that are useful/compatible with this variant"""
    
    print(f"\nIdentifying compatible ClinGen Allele Registry transcript IDs:")
    clingen_transcript_sequence_ids_final = []
    clingen_protein_sequence_ids_final = []

    #from the json of transcript info, get a list of the transcript IDs supported by clingen allele regsitry
    clingen_reference_sequence_ids = clingen_ar_utils.extract_reference_sequences(clingen_gene_transcripts_json)

    #when there are multiple versions of the same transcript, keep only the latest one
    clingen_reference_sequence_ids_latest = clingen_ar_utils.keep_latest_transcript_versions(clingen_reference_sequence_ids)

    #go through each transcript ID from CAR and see if it is worth checking for the current CIViC variant name
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
        #also check for methionine counting ambiguity (is there a matching ref AA one position to the right of the named position?)
        if not generic_utils.reference_aa_positions_matches(ref_aa_1, pos, protein_seq, protein_id):
            continue

        print(f"  Transcript ID: {clingen_reference_sequence_id} -> Protein ID: {protein_id} -> HGVS: {protein_id}:{civic_variant_name_p_3letter}")

        clingen_transcript_sequence_ids_final.append(clingen_reference_sequence_id)
        clingen_protein_sequence_ids_final.append(protein_id)

    return clingen_transcript_sequence_ids_final, clingen_protein_sequence_ids_final

def get_clingen_alleles_from_p_hgvs(clingen_protein_sequence_ids_final, civic_variant_name_p_3letter):
    """Get clingen allele info starting from a list of protein HGVS expressions"""
    clingen_alleles = []
    for protein_id in clingen_protein_sequence_ids_final:

        #create the protein HGVS expression
        protein_hgvs = f"{protein_id}:{civic_variant_name_p_3letter}"

        #query clingen API with protein_hgvs and get a json object of protein allele info
        pa_json = clingen_ar_utils.get_allele_by_hgvs(protein_hgvs)

        #extract transcript CAID and transcript HGVS values associated with the protein allele
        transcript_cas = clingen_ar_utils.extract_transcript_cas(pa_json)

        for tx in transcript_cas:
            clingen_alleles.append(tx['caid'])

    return list(set(clingen_alleles))

def get_clingen_allele_jsons(clingen_allele_ids):
    """ 
    For a list of CAIDs get the json object for each from clingen API
    Returns a dict of {caid: json_object}
    """     
    allele_jsons = {}
    
    for caid in clingen_allele_ids:
        ca_json = clingen_ar_utils.get_allele_by_id(caid)
        allele_jsons[caid] = ca_json

    return allele_jsons

def variant_is_ambiguous_in_genome(clingen_allele_info):
    """
    Check if a set of clingen alleles resolve to multiple distinct genomic variants
    The goal here is to find cases where a variant like BRAF V600E could be two distinct things
    This can happen when distinct transcripts both have a V at their position 600 ...
    but these correspond to distinct genomic positions
    """

    build_37_positions = []
        
    for caid, ca_json in clingen_allele_info.items():
        ca_json = clingen_ar_utils.get_allele_by_id(caid)

        clingen_coords = clingen_ar_utils.extract_genomic_coords(ca_json)
        for coord in clingen_coords:
            assembly = coord['assembly']
            if assembly != "GRCh37":
                continue
            genomic_variant = f"chr{coord['chr']}:{coord['start']}-{coord['end']}{coord['ref']}>{coord['alt']}"
            build_37_positions.append(genomic_variant)

    unique_positions = list(set(build_37_positions))

    if len(unique_positions) > 1:
        print(f"\nWARNING! Ambiguous genomic positions found: {unique_positions}")
        return True

    if len(clingen_allele_info) > 1:
        print(f"\nChecked {len(clingen_allele_info)} CAID(s) - no ambiguous genomic variants found.")

    return False


def get_build37_ensembl_transcripts_for_variant(clingen_transcript_sequence_ids_final, build37_ensembl_transcripts):
    """
    For each ENST transcript in the input list, find matching entries in v75 and v87
    annotations, allowing version-number mismatches.
    Returns a list of (clingen_id, v75_match, v87_match) tuples.
    """

    variant_build37_ensembl_transcripts = []

    for clingen_id in clingen_transcript_sequence_ids_final:
        if not clingen_id.startswith("ENST"):
            continue

        base_id = clingen_id.split(".")[0]  # strip version for comparison

        #search for matches, ignoring version numbers, and stop seaching once one is found
        v75_match = next(
            (k for k in build37_ensembl_transcripts["v75"] if k.split(".")[0] == base_id), None
        )
        v87_match = next(
            (k for k in build37_ensembl_transcripts["v87"] if k.split(".")[0] == base_id), None
        )

        #store each result as a tuple of (clingen_id, v75_match, v87_match) where the match values are the versioned keys as they appear in the annotation dicts
        #this makes it easy to then do a follow-up lookup like: build37_ensembl_transcripts["v75"][v75_match] if you need the full annotation entry
        variant_build37_ensembl_transcripts.append((clingen_id, v75_match, v87_match))

    return variant_build37_ensembl_transcripts


def display_accepted_variant_info(variant_id, accepted_variant_data):
    """Create a human readable summary of variant info already accepted in CIViC """
    print(
        f"\nAlready accepted variant details for variant id: {variant_id}\n"
        f"  Allele Registry ID: {accepted_variant_data['allele_registry_id']}\n"
        f"  Name: {accepted_variant_data['name']}\n"
        f"  Variant Aliases: {accepted_variant_data['variant_aliases']}\n"
        f"  HGVS Descriptions: {accepted_variant_data['hgvs_descriptions']}\n"
        f"  ClinVar IDs: {accepted_variant_data['clinvar_ids']}\n"
        f"  Reference Build: {accepted_variant_data['reference_build']}\n"
        f"  Chromosome: {accepted_variant_data['chromosome']}\n"
        f"  Start: {accepted_variant_data['start']}\n"
        f"  Stop: {accepted_variant_data['stop']}\n"
        f"  Reference Bases {accepted_variant_data['reference_bases']}\n"
        f"  Variant Bases: {accepted_variant_data['variant_bases']}\n"
        f"  Representative Transcript {accepted_variant_data['representative_transcript']}\n"
        f"  Ensembl Version: {accepted_variant_data['ensembl_version']}\n"
    )
    return

def main(variant_id: int, contributor_id: int, all_variants: bool, allow_variants_without_revisions: bool):

    #define input data files
    black_list_path = base_dir / f"data/civic_variant_blacklist.tsv"
    refseq_fasta_index_path = base_dir / f"data/refseq/indexed/GCF_000001405.40_GRCh38.p14_protein.faa.idx"
    ensembl_versions_file = base_dir / f"data/ensembl/ensembl_versions.txt"
    refseq_to_protein_file = base_dir / f"data/entrez/gene2refseq_human.tsv.gz"

    revision_value_key = {
        "variant": "revision_values_list",
        "coordinate": "suggested_value",
    }

    #make sure internet and API access is working before attempting anything
    verify_connectivity()

    #load black listed variants file
    black_listed_variant_ids = civic_graphql_utils.load_blacklisted_variant_ids(black_list_path)
 
    #get mappings of transcript to protein identifiers for refseq transcripts
    refseq_transcript_to_protein_map = entrez_utils.load_refseq_transcript_to_protein_map(refseq_to_protein_file)

    #get mappings of transcript to protein identifier for ensembl transcripts
    ensembl_transcript_to_protein_map = ensembl_utils.load_ensembl_transcript_to_protein_map(ensembl_versions_file)

    #get ensembl transcript biotype to transcript identifiers
    ensembl_transcript_to_biotype_map = ensembl_utils.load_ensembl_transcript_to_biotype_map(ensembl_versions_file)

    #get build37 ensembl transcript IDs with version numbers for (ensembl v75 and build37 imported ensembl v87)
    build37_ensembl_transcripts = ensembl_utils.load_build37_ensembl_transcripts()

    #summarize user info based on contributor id
    user_details = civic_graphql_utils.gather_user_details(contributor_id)
    print(f"\nContributor (id: {contributor_id}) is {user_details['user_name']} aka {user_details['user_display_name']} ({user_details['user_role']})")

    prompt_to_proceed("Verify your user info above. Revisions by this user will be ignored/skipped. \nYou can't moderate your own submissions.")

    #create a data structure that will store all clingen supported transcripts ids per gene
    clingen_transcript_ids = {}

    #get civic variant IDs to evaluate, either from the user, or by querying CIViCpy
    variant_ids_to_process = get_variant_ids_to_process(variant_id, all_variants)

    #iterate over each variant and examine revisions associated with it
    for vid in variant_ids_to_process:

        print(f"\nReviewing CIViC variant ID {vid} for revisions that could be reviewed by contributor ID: {contributor_id}")

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

        if variant_type_is_unsupported(civic_variant_name, guessed_gene_variant_type, "Missense Variant"): continue

        #get the three components of a simple "Missense Variant" variant: ref_aa_1, pos, var_aa_1
        ref_aa_1, pos, var_aa_1 = generic_utils.parse_snv_coding_name_components(civic_variant_name)

        #query the graphql api for more detailed variant revision info
        variant_data = civic_graphql_utils.gather_variant_revisions(vid, contributor_id)
        variant_data = civic_graphql_utils.merge_revision_data(variant_data)

        #if there is an outstanding revision to the variant name itself, warn the user
        if variant_name_has_revision(variant_data): continue

        #skip a variant if it has 0 pending revisions from other users - unless the user wishes to bypass this
        if variant_has_no_open_revisions(variant_data, allow_variants_without_revisions): continue

        #create the p. notation for the variant name (e.g. 'S459F' -> 'p.Ser459Phe')
        civic_variant_name_p_3letter = generic_utils.snv_coding_to_p_3letter(civic_variant_name)

        #get the gene name for the current variant
        gene_name = variant_data['feature_name']

        #provide a basic summary of variant info from civic
        print(
            f"\nVariant revision info:\n"
            f"  Feature name: {variant_data['feature_name']}\tVariant name in p. notation: {civic_variant_name_p_3letter}\n"
            f"  Open gene-variant revisions: {variant_data['open_revision_count_variant']} (total);"
            f"  {variant_data['contributor_revisions']} by you; {variant_data['open_revisions_non_contributor']} by others"
        )

        #diplay sources associated with the evidence for the molecular profiles involving this variant
        sources = civicpy_utils.get_sources_for_variant(vid)
        print(f"\nSources for all evidence associated with this variant (check for validity of variant name/identity)")
        for url, source in sources.items():
            print(f"  {url} ({source['citation']}. {source['source_type']})")

        #get currently *accepted* variant and variant coordinate info in CIViC for this variant
        accepted_variant_data = civic_graphql_utils.gather_accepted_variant_data(vid)

        #print out a summary of accepted variant info
        display_accepted_variant_info(variant_id, accepted_variant_data)

        #get all clingen allele registry transcripts supported for the gene of this variant
        #only query the clingen API if we don't already have transcripts for this gene
        clingen_gene_transcripts_json = get_clingen_gene_transcripts_json(gene_name, clingen_transcript_ids)

        #filter the transcripts to those that are useful/compatible with this variant, return protein IDs
        clingen_transcript_sequence_ids_final, clingen_protein_sequence_ids_final = get_compatible_clingen_transcripts(clingen_gene_transcripts_json, refseq_transcript_to_protein_map, ensembl_transcript_to_protein_map, ensembl_transcript_to_biotype_map, refseq_fasta_index_path, ensembl_versions_file, ref_aa_1, pos, var_aa_1, civic_variant_name_p_3letter)

        #for the list of valid transcript IDs from ClinGen, obtain old build37 compatible versions (from v75, and v87 import)
        #returns a list of tuples: clingen_transcript_id, ensembl_v75_match, ensembl_v87_match
        variant_build37_ensembl_transcripts = get_build37_ensembl_transcripts_for_variant(clingen_transcript_sequence_ids_final, build37_ensembl_transcripts)

        #get a unique list of useful/compatible CAIDs for the list of protein HGVS expressions
        clingen_allele_ids = get_clingen_alleles_from_p_hgvs(clingen_protein_sequence_ids_final, civic_variant_name_p_3letter)

        #using the list of CAIDs, get the json allele info object for each from the clingen API
        clingen_allele_info = get_clingen_allele_jsons(clingen_allele_ids)

        #look across the clingen alleles for ambiguous genomic variants - warn the user if found
        variant_is_ambiguous_in_genome(clingen_allele_info)
       
        #iterate through each useful/compatible CAID and display information that helps the user review outstanding edits
        for caid, ca_json in clingen_allele_info.items():
            print(f"\nCAID: {caid}")
            
            #query the clingen API with a CAID and get a json of relevant info
            ca_json = clingen_ar_utils.get_allele_by_id(caid)

            #for each clingen CAID get info that we would expect to be submited to CIViC:
            #variant aliases, clinvar ids, hgvs expressions, genomic coordinates (chr, start, stop, ref var)

            #get and disply genomic coordinate information for this CAID (all builds)
            #save coord info for build37 specifically
            clingen_assembly = clingen_chromosome = clingen_start = clingen_end = clingen_ref_bases = clingen_alt_bases = None
            clingen_coords = clingen_ar_utils.extract_genomic_coords(ca_json)
            clingen_genomic_hgvs_expressions = []
            for clingen_coord in clingen_coords:
                clingen_genomic_hgvs_expressions.append(clingen_coord['genomic_hgvs'])
                print(f"  {clingen_coord['assembly']} chr{clingen_coord['chr']}:{clingen_coord['start']}-{clingen_coord['end']} {clingen_coord['ref']}>{clingen_coord['alt']}")
                if clingen_coord['assembly'] == "GRCh37":
                    clingen_assembly = clingen_coord['assembly']
                    clingen_chromosome = clingen_coord['chr']
                    clingen_start = clingen_coord['start']
                    clingen_end = clingen_coord['end']
                    clingen_ref_bases = clingen_coord['ref']
                    clingen_alt_bases = clingen_coord['alt']

            #show the genomic HGVS expression for this CAID
            print(f"  Genomic HGVS expressions: {', '.join(clingen_genomic_hgvs_expressions)}")

            #get the list of MANE Select HGVS expressions for this CAID
            clingen_mane_select_hgvs_expressions = clingen_ar_utils.extract_mane_select_hgvs_expressions(ca_json)
            print(f"  MANE Select HGVS expressions: {', '.join(clingen_mane_select_hgvs_expressions)}")
 
            #compare the CIViC variant name to the MANE select variant name and warning if it doesn't match
            mane_select_names = clingen_ar_utils.extract_mane_select_names_and_compare(ca_json, civic_variant_name_p_3letter)

            #combine genomic and MANE select HGVS expressions into a single list of valid options
            clingen_combined_hgvs_expressions = clingen_genomic_hgvs_expressions + clingen_mane_select_hgvs_expressions

            #extract possible variant aliases across all transcripts for this CAID
            clingen_variant_aliases = clingen_ar_utils.extract_possible_variant_aliases(ca_json)
            print(f"  Possible variant aliases: {','.join(clingen_variant_aliases)}")

            #extract ClinVar IDs for this CAID
            clingen_clinvar_ids = clingen_ar_utils.extract_clinvar_ids(ca_json)
            print(f"  ClinVar IDs: {','.join(str(id) for id in clingen_clinvar_ids)}")

            #TODO: what to do about "ensembl_version". Ignore?
            clingen_data = {
                "variant_type": guessed_gene_variant_type,
                "variant_aliases": clingen_variant_aliases,
                "hgvs_expressions": clingen_combined_hgvs_expressions,
                "clinvar_ids": clingen_clinvar_ids,
                "assembly": clingen_assembly,
                "chromosome": clingen_chromosome, 
                "start": clingen_start,
                "end": clingen_end,
                "ref_bases": clingen_ref_bases,
                "alt_bases": clingen_alt_bases,
                "representative_transcript": variant_build37_ensembl_transcripts,
                "ensembl_version": None
            }

            #############################################################################################
            #Perform comparison between the ClinGen Allele Info for this CAID and CIViC Variant Revisions
            print(f"\n  Comparing existing CIViC revisions to the values of this CAID ({caid}):")
            comparator = compare_utils.RevisionComparator(clingen_data)
            all_revisions = variant_data['all_revisions']
            for revision in all_revisions:
                #get revision value from correct dictionary slot according to revision type
                revision_value = revision[revision_value_key[revision["revision_type"]]]
                revision_id = revision['revision_id']
                user_display_name = revision['user_display_name']
                user_id = revision['user_id']
                field_name = revision['field_name']

                is_consistent = comparator.compare(field_name, revision_value, revision_id)





        #Pause before moving on to the next CIViC variant
        prompt_to_proceed("Processing complete for variant ({vid}: civic_variant_name)")


if __name__ == "__main__":
    args = parse_args()

    main(
        contributor_id=args.contributor_id,
        variant_id=args.variant_id,
        all_variants=args.all_variants,
        allow_variants_without_revisions=args.allow_variants_without_revisions
    )

