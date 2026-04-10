#!/usr/bin/env python3

import requests
import urllib.parse
import sys
import re

# handle both package and standalone execution
try:
    from . import generic_utils
except ImportError:
    import generic_utils

BASE_URL = "https://reg.genome.network"

def get_allele_by_hgvs(hgvs):
    """query clingen API using a protein HGVS (e.g. NP_004324.2:p.Val600Glu)"""
    url = f"{BASE_URL}/allele?hgvs={urllib.parse.quote(hgvs, safe='')}"
    r = requests.get(url, headers={"Accept": "application/json"})
    r.raise_for_status()
    return r.json()

def get_allele_by_id(allele_id):
    """query clingen API using a transcript level CAID
       e.g. https://reg.genome.network/allele/CA387358756
    """
    url = f"{BASE_URL}/allele/{allele_id}"
    r = requests.get(url, headers={"Accept": "application/json"})
    r.raise_for_status()
    return r.json()

def get_reference_sequences_by_gene(gene_name):
    """query clinen API using a gene name
       e.g. https://reg.genome.network/refseqs?gene=POLE
    """
    url = f"{BASE_URL}/refseqs?gene={gene_name}"
    r = requests.get(url, headers={"Accept": "application/json"})
    r.raise_for_status()
    return r.json()

def extract_transcript_cas(pa_json):
    """extract CAIDs from a protein PAID json object"""
    cas = []
    for aa in pa_json.get("aminoAcidAlleles", []):
        for tx in aa.get("matchingRegisteredTranscripts", []):
            ca_id = tx["@id"].split("/")[-1] #extract CAID from url
            cas.append({
                "caid": ca_id,
                "hgvs": tx["hgvs"]
            })
    return cas

def extract_genomic_coords(ca_json):
    """extract genomic coordinate info from a transcript CAID json object
       return a list of coordinate objects, each with:
       assembly, chr, start, end, ref, alt, genomic hgvs 
    """
    coords = []
    allowed_assemblies = {"GRCh37", "GRCh38"}

    for g in ca_json.get("genomicAlleles", []):
        assembly = g.get("referenceGenome")

        # Skip anything not GRCh37/38
        if assembly not in allowed_assemblies:
            continue

        chrom = g.get("chromosome")
        genomic_hgvs_expressions = []

        for hgvs in g.get("hgvs", []):
            if not hgvs.startswith("NC_"):
                continue
            genomic_hgvs_expressions.append(hgvs)

        for c in g.get("coordinates", []):
            coords.append({
                "assembly": assembly,
                "chr": chrom,
                "start": c.get("start"),
                "end": c.get("end"),
                "ref": c.get("referenceAllele"),
                "alt": c.get("allele"),
                "genomic_hgvs": ",".join(genomic_hgvs_expressions),
            })

    return coords


def extract_possible_variant_aliases(ca_json):
    """extract possible variant aliases from a transcript CAID json object
       return a list of aliases to check CIViC proposed aliases against
    """
    possible_variant_aliases = []
    for transcript in ca_json.get("transcriptAlleles", []):
        protein_effect = transcript.get("proteinEffect")
        if protein_effect:
            hgvs = protein_effect.get("hgvs")
            possible_variant_aliases.append(hgvs)


    possible_variant_aliases_filtered = []

    for alias in possible_variant_aliases:
        seq_id, rest = alias.split(':', 1)
    
        # Skip unwanted sequence IDs
        if seq_id.startswith(('XP_', 'XR_')):
            continue
    
        # Skip nucleotide variants
        if rest.startswith('n.'):
            continue
    
        # Extract variant name
        variant_name = rest.split('.', 1)[1]
        possible_variant_aliases_filtered.append(variant_name)

        # Convert 3-letter AA codes to 1-letter and append as additional alias
        # Matches patterns like Ser432Phe, Ala123Val, etc.
        variant_name_1_aa = re.sub(
            r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})',
            lambda m: f"{generic_utils.aa_3_to_1(m.group(1))}{m.group(2)}{generic_utils.aa_3_to_1(m.group(3))}",
            variant_name
        )
        if variant_name_1_aa != variant_name:
            possible_variant_aliases_filtered.append(variant_name_1_aa)

    return list(dict.fromkeys(possible_variant_aliases_filtered))


def extract_mane_select_hgvs_expressions(ca_json):
    """extract mane select hgvs expression from a transcript CAID json object"""
    mane_select_hgvs_expressions = []

    for transcript in ca_json.get("transcriptAlleles", []):
        mane = transcript.get("MANE")
        if mane:
            nucleotide = mane.get("nucleotide", {})
            mane_select_hgvs_expressions.append(nucleotide.get("Ensembl", {}).get("hgvs"))
            mane_select_hgvs_expressions.append(nucleotide.get("RefSeq", {}).get("hgvs"))

            protein = mane.get("protein", {})
            mane_select_hgvs_expressions.append(protein.get("Ensembl", {}).get("hgvs"))
            mane_select_hgvs_expressions.append(protein.get("RefSeq", {}).get("hgvs"))
            
    return list(set(mane_select_hgvs_expressions))


def extract_clinvar_ids(ca_json):
    """"extract clinvar ids from a transcript CAID json object"""
    clinvar_ids = []
    external_records = ca_json.get("externalRecords")
    for clinvar_allele in external_records.get("ClinVarAlleles", []):
        clinvar_ids.append(clinvar_allele.get("alleleId"))

    return list(set(clinvar_ids))


def extract_reference_sequences(ref_seqs_json):
    """get all supported ensembl and refseq transcript ids for a gene from clingen allele registry"""
    tid_list = []
    for reference_sequence in ref_seqs_json:
        type = reference_sequence['type']
        if type != "transcript":
            continue
        external_records = reference_sequence['externalRecords']

        tid = (
            external_records.get('NCBI', {}).get('id')
            or
            external_records.get('Ensembl', {}).get('id')
        )
        # reject predicted/refseq model transcripts
        if tid.startswith(("XM_", "XR_")):
            continue

        if tid is None:
            sys.exit("Error: No transcript ID found in NCBI or Ensembl external records.")
        tid_list.append(tid)
    return tid_list


def keep_latest_transcript_versions(reference_sequence_ids):
    latest = {}

    for tid in reference_sequence_ids:
        if "." not in tid:
            raise ValueError(f"Transcript ID missing version: {tid}")

        base_id, version_str = tid.rsplit(".", 1)

        if not version_str.isdigit():
            raise ValueError(f"Invalid version in transcript ID: {tid}")

        version = int(version_str)

        if base_id not in latest or version > latest[base_id][1]:
            latest[base_id] = (tid, version)

    return [v[0] for v in latest.values()]


if __name__ == "__main__":
    """given a protein level civic variant (e.g. BRAF V600E) get possible coords from clingen"""

    gene_symbol = "POLE"
    protein_id = "NP_006222.2"
    p_dot_var = "p.Ser459Phe" #S459F
    hgvs_protein = f"{protein_id}:{p_dot_var}"

    print("\nHGVS protein query:", hgvs_protein)

    # query clingen api with a protein allele hgvs
    # Example: http://reg.genome.network/allele/PA094029
    pa = get_allele_by_hgvs(hgvs_protein)
    print("Protein allele:", pa["@id"])

    # extract transcript-level CA alleles that were found for the protein allele
    # example: http://reg.genome.network/allele/CA123643
    transcript_cas = extract_transcript_cas(pa)

    # for each transcript-level CA allele, get additional information 
    for tx in transcript_cas:
        print(f"\nTranscript {tx['hgvs']} ({tx['caid']})")

        # extract genomic coordinates and hgvs expressions for each transcript caid found
        ca_json = get_allele_by_id(tx["caid"])
        for coord in extract_genomic_coords(ca_json):
            print(
                f"  {coord['assembly']} chr{coord['chr']}:{coord['start']}-{coord['end']} "
                f"{coord['ref']}>{coord['alt']} {coord['genomic_hgvs']}"
            )
        
        # extract MANE select transcript hgvs expressions for each transcript caid found
        mane_select_hgvs_expressions = extract_mane_select_hgvs_expressions(ca_json)
        print(f"\nMANE Select HGVS expressions:\n  {','.join(mane_select_hgvs_expressions)}")

        # extract possible variant aliases across all transcripts
        variant_aliases = extract_possible_variant_aliases(ca_json)
        print(f"\nPossible variant aliases:\n {','.join(variant_aliases)}")

        # extract ClinVar IDs for each transcript caid found
        clinvar_ids = extract_clinvar_ids(ca_json)
        print(f"\nClinVar IDs:\n  {','.join(str(id) for id in clinvar_ids)}")

    # for a single gene, get all the CAR supported transcript identifiers
    rs = get_reference_sequences_by_gene(gene_symbol)
    reference_sequence_ids = extract_reference_sequences(rs)
    reference_sequence_ids_string = ",".join(reference_sequence_ids)
    print(f"\nClinGen supported transcript ids (excluding XR_ and NR_ transcripts):\n{reference_sequence_ids_string}")

    # produce a filtered list that keeps only the most recent version of each transcript support by CAR
    reference_sequence_ids_latest = keep_latest_transcript_versions(reference_sequence_ids)
    reference_sequence_ids_latest_string = ",".join(reference_sequence_ids_latest)
    print(f"\nClinGen supported transcript ids (limited to only the most recent versions):\n{reference_sequence_ids_latest_string}")


