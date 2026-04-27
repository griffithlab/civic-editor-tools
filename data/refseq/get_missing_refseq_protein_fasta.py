#!/usr/bin/env python3

"""
get_missing_refseq_protein.py

Fetches a RefSeq protein sequence using the NCBI efetch CLI tool and appends
it to a local FASTA file. Creates the file if it does not exist.
After any write, rebuilds the SeqIO index alongside the FASTA.

Usage:
    python get_missing_refseq_protein.py NP_059989.2

Requires:
    NCBI Entrez Direct (efetch) and Biopython (for indexing only)
"""

import os
import sys
import time
import shutil
import argparse
import subprocess
from pathlib import Path

from Bio import SeqIO

REQUEST_DELAY = 0.5

BASE_DIR = Path(__file__).resolve().parent.parent.parent
FASTA_PATH = BASE_DIR / "data/refseq/indexed/missing_refseqs_protein.faa"
INDEX_PATH = FASTA_PATH.with_suffix(".faa.idx")

def check_efetch_available():
    """Exit early with a helpful message if efetch is not on PATH."""
    if shutil.which("efetch") is None:
        sys.exit(
            "Error: 'efetch' not found on PATH.\n"
            "Install NCBI Entrez Direct: https://www.ncbi.nlm.nih.gov/books/NBK179288/\n"
            "  sh -c \"$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)\""
        )


def fetch_protein_fasta(accession: str) -> str:
    """Fetch a protein record from NCBI in FASTA format using the efetch CLI."""
    print(f"Fetching {accession} from NCBI...")

    cmd = ["efetch", "-db", "protein", "-id", accession, "-format", "fasta"]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        sys.exit(f"efetch failed for '{accession}':\n{e.stderr.strip()}")
    finally:
        time.sleep(REQUEST_DELAY)

    fasta_text = result.stdout.strip()
    if not fasta_text:
        sys.exit(f"Error: NCBI returned an empty response for '{accession}'.")

    return fasta_text


def accession_in_fasta(accession: str, fasta_path: Path) -> bool:
    """Return True if any record in the FASTA has a matching accession."""
    for record in SeqIO.parse(fasta_path, "fasta"):
        if accession in record.id or accession in record.description:
            return True
    return False


def build_index(fasta_path: Path, index_path: Path) -> None:
    """Build (or rebuild) a SeqIO index for fast random access."""
    if index_path.exists():
        index_path.unlink()
    SeqIO.index_db(str(index_path), str(fasta_path), "fasta")
    print(f"Index written → {index_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Fetch a RefSeq protein sequence by accession (e.g. NP_059989.2) and append it to a local FASTA file."
    )
    parser.add_argument("accession", help="RefSeq protein accession, e.g. NP_059989.2")
    args = parser.parse_args()
    accession = args.accession.strip()

    if not (accession.startswith("NP_") or accession.startswith("XP_")):
        sys.exit(
            f"Error: '{accession}' does not look like a RefSeq protein accession "
            f"(expected NP_* or XP_*)."
        )

    check_efetch_available()
    FASTA_PATH.parent.mkdir(parents=True, exist_ok=True)

    if FASTA_PATH.exists():
        if accession_in_fasta(accession, FASTA_PATH):
            print(f"'{accession}' is already present in {FASTA_PATH} — nothing to do.")
            sys.exit(0)
        print(f"'{accession}' confirmed missing from {FASTA_PATH}. Will append.")
    else:
        print(f"{FASTA_PATH} does not exist — will create it.")

    fasta_text = fetch_protein_fasta(accession)

    if accession not in fasta_text:
        sys.exit(
            f"Error: fetched FASTA does not contain '{accession}'.\n"
            f"Header: {fasta_text.splitlines()[0]!r}"
        )

    mode = "a" if FASTA_PATH.exists() else "w"
    with FASTA_PATH.open(mode) as fh:
        if mode == "a":
            fh.write("\n")
        fh.write(fasta_text + "\n")

    action = "Appended to" if mode == "a" else "Created"
    print(f"{action} {FASTA_PATH}")

    build_index(FASTA_PATH, INDEX_PATH)


if __name__ == "__main__":
    main()


