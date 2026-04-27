#!/usr/bin/env python3
"""
Fetch RefSeq protein IDs from RefSeq transcript IDs using NCBI's efetch utility.
Mappings are cached in gene2refseq_human_missing.tsv to avoid redundant lookups.
"""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

#set data path location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_FILE = Path(os.path.join(SCRIPT_DIR, "gene2refseq_human_missing.tsv"))
CACHE_HEADER = "#transcript_id\tprotein_id\n"

def check_efetch_installed():
    """Verify that efetch (NCBI E-utilities) is available on PATH."""
    if shutil.which("efetch") is None:
        sys.exit(
            "Error: 'efetch' is not installed or not on your PATH.\n"
            "Install the NCBI E-utilities: "
            "sh -c \"$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)\""
        )


def load_cache() -> dict[str, str]:
    """Load existing transcript→protein mappings from the cache file."""
    cache: dict[str, str] = {}
    if not CACHE_FILE.exists():
        return cache

    with CACHE_FILE.open() as fh:
        for line in fh:
            line = line.strip()
            # Skip header and blank lines
            if not line or line.startswith("transcript_id"):
                continue
            parts = line.split("\t")
            if len(parts) == 2:
                cache[parts[0]] = parts[1]

    return cache


def save_mapping(transcript_id: str, protein_id: str):
    """Append a new transcript→protein mapping to the cache file."""
    write_header = not CACHE_FILE.exists() or CACHE_FILE.stat().st_size == 0

    with CACHE_FILE.open("a") as fh:
        if write_header:
            fh.write(CACHE_HEADER)
        fh.write(f"{transcript_id}\t{protein_id}\n")


def fetch_protein_id(transcript_id: str) -> str:
    """
    Run efetch for the given transcript ID and parse out the protein_id field.

    Returns the protein ID string (e.g. 'NP_065681.1'), or raises RuntimeError
    if the protein ID cannot be found in the returned record.
    """
    cmd = ["efetch", "-db", "nucleotide", "-id", transcript_id, "-format", "gb"]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"efetch failed for '{transcript_id}' (exit {exc.returncode}):\n"
            f"{exc.stderr.strip()}"
        ) from exc

    # Match lines like:   /protein_id="NP_065681.1"
    match = re.search(r'/protein_id="([^"]+)"', result.stdout)
    if not match:
        raise RuntimeError(
            f"No protein_id found in GenBank record for '{transcript_id}'.\n"
            "The transcript may be non-coding or the accession may be invalid."
        )

    return match.group(1)


def get_protein_id(transcript_id: str) -> str:
    """
    Return the protein ID for a transcript, using the cache when possible.
    Updates the cache file with any newly fetched mapping.
    """
    cache = load_cache()

    if transcript_id in cache:
        print(f"[cache] {transcript_id} → {cache[transcript_id]}")
        return cache[transcript_id]

    print(f"[fetch] Querying NCBI for {transcript_id} ...")
    protein_id = fetch_protein_id(transcript_id)
    save_mapping(transcript_id, protein_id)
    print(f"[fetch] {transcript_id} → {protein_id}  (saved to {CACHE_FILE})")
    return protein_id


def main():
    if len(sys.argv) < 2:
        sys.exit(f"Usage: {sys.argv[0]} <transcript_id> [transcript_id ...]")

    check_efetch_installed()

    for transcript_id in sys.argv[1:]:
        try:
            protein_id = get_protein_id(transcript_id)

        except RuntimeError as exc:
            print(f"Error: {exc}", file=sys.stderr)


if __name__ == "__main__":
    main()
