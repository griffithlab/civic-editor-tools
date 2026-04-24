#!/usr/bin/env python3

from Bio import Entrez
import ssl
import certifi
from Bio import Entrez
import xml.etree.ElementTree as ET
import gzip
from pathlib import Path

base_dir = Path(__file__).resolve().parent

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

def load_mane_summary(path):
    """create a dictionary of gene symbol to gene/protein ids for mane select transcripts only"""
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


def load_refseq_transcript_to_protein_map(gene2refseq_path: Path, gene2refseq_missing_path: Path) -> dict:
    """load refseq transcript to protein id mappings from a file"""
    tx_to_protein = {}

    #first parse the main source of refseq transcript to protein ids
    open_func = gzip.open if gene2refseq_path.suffix == ".gz" else open
    with open_func(gene2refseq_path, "rt") as f:
        for line in f:
            line = line.strip()

            # Skip header or empty lines
            if not line or line.startswith("#"):
                continue

            fields = line.split()

            # Defensive check (make sure enough columns exist)
            if len(fields) < 6:
                continue

            transcript_id = fields[3]
            protein_id = fields[5]

            if not transcript_id.startswith("NM_"):
                continue

            tx_to_protein[transcript_id] = protein_id

    #next add any missing values that have been gathered manually
    open_func = gzip.open if gene2refseq_missing_path.suffix == ".gz" else open
    with open_func(gene2refseq_missing_path, "rt") as f:
        for line in f:
            line = line.strip()

            # Skip header or empty lines
            if not line or line.startswith("#"):
                continue

            fields = line.split()

            # Defensive check (make sure enough columns exist)
            if len(fields) < 2:
                continue

            transcript_id = fields[0]
            protein_id = fields[1]

            if not transcript_id.startswith("NM_"):
                continue

            tx_to_protein[transcript_id] = protein_id

    return tx_to_protein


def main(gene_symbol):
    gene_id = get_gene_id(gene_symbol)
    mane = load_mane_summary(base_dir / f"../data/refseq/MANE.GRCh38.v1.5.summary.txt")
    mane_nm = mane[gene_symbol]["transcript"]
    mane_np = mane[gene_symbol]["protein"]
    #mane_np = get_protein_from_transcript(mane_nm)

    print(f"Gene symbol: {gene_symbol}")
    print(f"NCBI Gene ID: {gene_id}")
    print(f"MANE Select transcript: {mane_nm}")
    print(f"Corresponding protein: {mane_np}")

    gene2refseq_path = base_dir / f"../data/entrez/gene2refseq_human.tsv.gz"
    gene2refseq_missing_path = base_dir / f"../data/entrez/gene2refseq_human_missing.tsv"
    refseq_transcript_to_protein_map = load_refseq_transcript_to_protein_map(gene2refseq_path, gene2refseq_missing_path)
    mapped_np = refseq_transcript_to_protein_map.get(mane_nm)
    print(f"Mapped protein id: {mapped_np}")



if __name__ == "__main__":
    main("POLE")

