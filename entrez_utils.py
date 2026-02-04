#!/usr/bin/env python3

from Bio import Entrez
import ssl
import certifi
from Bio import Entrez
import xml.etree.ElementTree as ET

def _patched_https_context():
    return ssl.create_default_context(cafile=certifi.where())

ssl._create_default_https_context = _patched_https_context

#Email required to use Entrez
Entrez.email = "mgriffit@wustl.edu"

# given a gene symbol (e.g. "BRAF" get entrez gene id)
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

# create a dictionary of gene symbol to gene/protein ids (read from flat file)
def load_mane_summary(path):
    mane = {}
    with open(path) as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            fields = dict(zip(header, line.rstrip("\n").split("\t")))
            if fields["MANE_status"] == "MANE Select":
                mane[fields["symbol"]] = {
                    "transcript": fields["RefSeq_nuc"],
                    "protein": fields["RefSeq_prot"]
                }
    return mane

def get_protein_from_transcript(transcript_id):
    """Map NM_ transcript to NP_ protein by parsing RefSeq record"""

    # NCBI efetch is unreliable with versioned NM_ accessions
    base_id = transcript_id.split(".")[0]

    handle = Entrez.efetch(
        db="nucleotide",
        id=base_id,
        rettype="gbwithparts",
        retmode="text"
    )
    record = handle.read()
    handle.close()

    for line in record.splitlines():
        if line.strip().startswith("/protein_id="):
            return line.split('"')[1]

    raise ValueError(f"No protein_id found for {transcript_id}")


def main(gene_symbol):
    gene_id = get_gene_id(gene_symbol)
    mane = load_mane_summary("data/MANE.GRCh38.v1.5.summary.txt")
    mane_nm = mane[gene_symbol]["transcript"]
    mane_np = mane[gene_symbol]["protein"]
    #mane_np = get_protein_from_transcript(mane_nm)

    print(f"Gene symbol: {gene_symbol}")
    print(f"NCBI Gene ID: {gene_id}")
    print(f"MANE Select transcript: {mane_nm}")
    #print(f"Corresponding protein: {mane_np}")


if __name__ == "__main__":
    main("BRAF")

