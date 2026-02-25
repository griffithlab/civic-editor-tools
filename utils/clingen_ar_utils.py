#!/usr/bin/env python3

import requests
import urllib.parse
import sys

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
    """query clinen API using a gene name"""
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
    """extract genomic coordinate info from a transcript CAID json object"""
    coords = []
    allowed_assemblies = {"GRCh37", "GRCh38"}

    for g in ca_json.get("genomicAlleles", []):
        assembly = g.get("referenceGenome")

        # Skip anything not GRCh37/38
        if assembly not in allowed_assemblies:
            continue

        chrom = g.get("chromosome")

        for c in g.get("coordinates", []):
            coords.append({
                "assembly": assembly,
                "chr": chrom,
                "start": c.get("start"),
                "end": c.get("end"),
                "ref": c.get("referenceAllele"),
                "alt": c.get("allele"),
            })

    return coords

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

        # extract genomic coordinates for each transcript caid found
        ca = get_allele_by_id(tx["caid"])
        for c in extract_genomic_coords(ca):
            print(
                f"  {c['assembly']} chr{c['chr']}:{c['start']}-{c['end']} "
                f"{c['ref']}>{c['alt']}"
            )

    # for a single gene, get all the CAR supported transcript identifiers
    rs = get_reference_sequences_by_gene(gene_symbol)
    reference_sequence_ids = extract_reference_sequences(rs)
    print(f"\nClinGen supported transcript ids (excluding XR_ and NR_ transcripts):\n {reference_sequence_ids}")

    # produce a filtered list that keeps only the most recent version of each transcript support by CAR
    reference_sequence_ids_latest = keep_latest_transcript_versions(reference_sequence_ids)
    print(f"\nClinGen supported transcript ids (only most recent versions):\n {reference_sequence_ids_latest}")


