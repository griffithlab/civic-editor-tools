#!/usr/bin/env python3

from Bio import Entrez

import ssl
import certifi
from Bio import Entrez

def _patched_https_context():
    return ssl.create_default_context(cafile=certifi.where())

ssl._create_default_https_context = _patched_https_context

# REQUIRED by NCBI
Entrez.email = "mgriffit@wustl.edu"


def get_gene_id(gene_symbol, organism="human"):
    """Resolve gene symbol to NCBI Gene ID"""
    handle = Entrez.esearch(
        db="gene",
        term=f"{gene_symbol}[Symbol] AND {organism}[Organism]"
    )
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        raise ValueError(f"No Gene ID found for {gene_symbol}")

    return record["IdList"][0]


def get_mane_select_transcript(gene_id):
    """Extract MANE Select transcript (NM_) from Gene XML"""
    handle = Entrez.efetch(
        db="gene",
        id=gene_id,
        rettype="xml"
    )
    root = ET.fromstring(handle.read())
    handle.close()

    for comment in root.iter("Gene-commentary"):
        label = comment.findtext("Gene-commentary_label")
        if label == "MANE Select":
            accession = comment.findtext(
                "Gene-commentary_accession"
            )
            return accession

    raise ValueError("MANE Select transcript not found")


def get_protein_from_transcript(transcript_id):
    """Map NM_ transcript to NP_ protein"""
    handle = Entrez.elink(
        dbfrom="nucleotide",
        db="protein",
        id=transcript_id
    )
    record = Entrez.read(handle)
    handle.close()

    links = record[0]["LinkSetDb"]
    if not links:
        raise ValueError(f"No protein linked to {transcript_id}")

    protein_id = links[0]["Link"][0]["Id"]

    handle = Entrez.efetch(
        db="protein",
        id=protein_id,
        rettype="acc"
    )
    protein_acc = handle.read().strip()
    handle.close()

    return protein_acc


def main(gene_symbol):
    gene_id = get_gene_id(gene_symbol)
    mane_nm = get_mane_select_transcript(gene_id)
    mane_np = get_protein_from_transcript(mane_nm)

    print(f"Gene symbol: {gene_symbol}")
    print(f"NCBI Gene ID: {gene_id}")
    print(f"MANE Select transcript: {mane_nm}")
    print(f"Corresponding protein: {mane_np}")


if __name__ == "__main__":
    main("BRAF")

