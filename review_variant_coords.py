#!/usr/bin/env python3

import requests
import urllib.parse

BASE_URL = "https://reg.genome.network"

# query clingen API using a protein HGVS (e.g. NP_004324.2:p.Val600Glu)
def get_allele_by_hgvs(hgvs):
    url = f"{BASE_URL}/allele?hgvs={urllib.parse.quote(hgvs, safe='')}"
    r = requests.get(url, headers={"Accept": "application/json"})
    r.raise_for_status()
    return r.json()

# query clingen API using a transcript level CAID
def get_allele_by_id(allele_id):
    url = f"{BASE_URL}/allele/{allele_id}"
    r = requests.get(url, headers={"Accept": "application/json"})
    r.raise_for_status()
    return r.json()

# extract CAIDs from a protein PAID json object
def extract_transcript_cas(pa_json):
    cas = []
    for aa in pa_json.get("aminoAcidAlleles", []):
        for tx in aa.get("matchingRegisteredTranscripts", []):
            ca_id = tx["@id"].split("/")[-1] #extract CAID from url
            cas.append({
                "caid": ca_id,
                "hgvs": tx["hgvs"]
            })
    return cas

# extract genomic coordinate info from a transcript CAID json object
def extract_genomic_coords(ca_json):
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


# given a protein level civic variant (e.g. BRAF V600E) get possible coords from clingen
if __name__ == "__main__":

    hgvs_protein = "NP_004324.2:p.Val600Glu"
    print("\nHGVS protein query:", hgvs_protein)

    # query clingen api with a protein allele hgvs
    # Example: http://reg.genome.network/allele/PA094029
    pa = get_allele_by_hgvs(hgvs_protein)
    print("Protein allele:", pa["@id"])

    # extract transcript-level CA alleles that were found for the protein allele
    # example: http://reg.genome.network/allele/CA123643
    transcript_cas = extract_transcript_cas(pa)

    for tx in transcript_cas:
        print(f"\nTranscript {tx['hgvs']} ({tx['caid']})")

        # extract genomic coordinates for each transcript caid found
        ca = get_allele_by_id(tx["caid"])
        for c in extract_genomic_coords(ca):
            print(
                f"  {c['assembly']} chr{c['chr']}:{c['start']}-{c['end']} "
                f"{c['ref']}>{c['alt']}"
            )

