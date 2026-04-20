#!/usr/bin/env python3

import re
import urllib.request
import urllib.error
import socket
import sys
import ssl
import requests

def aa_1_to_3(aa):
    """Convert single letter AA abbreviations into three letter AA abbreviations"""
    aa_map = {
        "A": "Ala",
        "R": "Arg",
        "N": "Asn",
        "D": "Asp",
        "C": "Cys",
        "E": "Glu",
        "Q": "Gln",
        "G": "Gly",
        "H": "His",
        "I": "Ile",
        "L": "Leu",
        "K": "Lys",
        "M": "Met",
        "F": "Phe",
        "P": "Pro",
        "S": "Ser",
        "T": "Thr",
        "W": "Trp",
        "Y": "Tyr",
        "V": "Val",
        "*": "Ter"   # stop codon (optional)
    }
    try:
        return aa_map[aa.upper()]
    except KeyError:
        raise ValueError(f"Invalid amino acid code: {aa}")


def aa_3_to_1(aa_3):
    """Convert three letter AA abbreviations into single letter AA abbreviations"""
    # Mapping in reverse: 3-letter -> 1-letter
    aa_map_rev = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLU": "E",
        "GLN": "Q",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
        "TER": "*",  # Stop codon
        "STOP": "*"  # Common alternative for stop
    }
    
    # Normalize input to uppercase
    aa_3_upper = aa_3.upper()
    
    try:
        return aa_map_rev[aa_3_upper]
    except KeyError:
        raise ValueError(f"Invalid three-letter amino acid code: {aa_3}")


def is_valid_aa(aa):
    """Test is a single AA value is valid"""
    valid_aas = {
        "A", "R", "N", "D", "C",
        "E", "Q", "G", "H", "I",
        "L", "K", "M", "F", "P",
        "S", "T", "W", "Y", "V"
    }
    return aa.upper() in valid_aas


def guess_variant_type(variant):
    """Based only on a variant name from civic attempt to guess the likely variant type"""
    if not variant:
        return None

    variant = variant.strip()

    # 1) frameshift: e.g. A321fs
    match = re.match(r'^([A-Za-z])(\d+)(fs)', variant, re.IGNORECASE)
    if match:
        aa = match.group(1)
        if is_valid_aa(aa):
            return "frameshift"

    # 2) coding SNV: e.g. S459F
    snv_pattern = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])(?=\s|$)")

    match = snv_pattern.match(variant)
    if match:
        ref_aa, pos, alt_aa = match.groups()

        if is_valid_aa(ref_aa) and is_valid_aa(alt_aa):
            return "Missense Variant"
    
    # 3) splice site
    if variant.startswith("Splice Site"):
        return "splice_site"

    # None of the categories above were guessed from the variant name
    return None


def snv_coding_to_p_3letter(variant):
    """Convert a coding SNV like 'S459F' to 'p.Ser459Phe'"""
    
    if guess_variant_type(variant) != "Missense Variant":
        raise ValueError(f"Expected 'Missense Variant', got '{variant_type}' for variant '{variant}'")

    match = re.match(r"^([A-Za-z])(\d+)([A-Za-z])", variant)
    if not match:
        raise ValueError(f"Variant '{variant}' does not match expected coding SNV format (e.g. 'S459F')")

    ref_aa, pos, alt_aa = match.groups()

    return f"p.{aa_1_to_3(ref_aa)}{pos}{aa_1_to_3(alt_aa)}"


def parse_snv_coding_name_components(snv_coding_name):
    """
    Parse a CIViC-style protein variant like 'S459F'
    into (ref_aa_1, pos, var_aa_1).

    Raises ValueError if the format is invalid.
    """
    pattern = re.compile(r"^([A-Z])(\d+)([A-Z])$")

    match = pattern.match(snv_coding_name)
    if not match:
        raise ValueError(f"Invalid variant format: {snv_coding_name}")

    ref_aa_1, pos, var_aa_1 = match.groups()

    return ref_aa_1, int(pos), var_aa_1


def reference_aa_positions_matches(ref_aa_1, pos, protein_seq, protein_id):
    """
    Check whether the amino acid at position `pos` in `protein_seq`
    matches `ref_aa_1`.

    Also checks for methionine counting ambiguity: if the amino acid at
    pos+1 matches ref_aa_1, a warning is issued. This can occur when the 
    variant was named without counting the initiator methionine, shifting 
    all positions by one.
    """
    # Basic validation
    if not isinstance(ref_aa_1, str) or len(ref_aa_1) != 1:
        raise ValueError("ref_aa_1 must be a single-character string")
    if not isinstance(pos, int) or pos < 1:
        raise ValueError("pos must be a positive integer (1-based)")
    if not is_valid_aa(ref_aa_1):
        raise ValueError("ref_aa_1 must be a valid amino acid")
    if pos > len(protein_seq):
        return False  # position out of range

    # Check for AA match at the position used in the name
    primary_match = protein_seq[pos - 1] == ref_aa_1

    # Check for methionine counting ambiguity at pos+1
    secondary_pos = pos + 1
    if secondary_pos <= len(protein_seq):
        secondary_match = protein_seq[secondary_pos - 1] == ref_aa_1
        if secondary_match:
            print(f"\nWARNING! Possible methionine counting ambiguity detected for AA:{ref_aa_1} in {protein_id}:")
            print(f"\nPrimary position {pos} match: {primary_match}\tSecondary position {secondary_pos} match: {secondary_match}")

    return primary_match


def check_connection(timeout: int = 5) -> bool:
    """
    Returns True if internet access is available, False otherwise.
    Tests by making a lightweight HEAD request to relevant hosts.
    """
    test_urls = [
        "https://www.google.com",
        "https://www.civicdb.org",
        "https://reg.genome.network",
    ]
    for url in test_urls:
        try:
            context = ssl._create_unverified_context()
            req = urllib.request.Request(url, method="HEAD")
            with urllib.request.urlopen(req, timeout=timeout, context=context):
                return True
        except (urllib.error.URLError, socket.timeout, OSError):
            continue
    return False


def check_apis(timeout: int = 5) -> bool:
    """
    Returns True if both required APIs are reachable, False otherwise.
    """
    api_urls = [
        "https://www.civicdb.org",
        "https://reg.genome.network",
    ]
    for url in api_urls:
        try:
            response = requests.head(url, timeout=timeout)
            response.raise_for_status()
        except requests.RequestException as e:
            print(f"API unreachable: {url} — {e}")
            return False
    return True


def main():
    
    if not check_connection():
        print("No internet access. Aborting.")
        sys.exit(1)

    if not check_apis():
        print("Required APIs are unavailable. Aborting.")
        sys.exit(1)

    test_aa_1 = "W"
    test_aa_3 = aa_1_to_3(test_aa_1)
    print (f"test_aa_1: {test_aa_1}\ttest_aa_3: {test_aa_3}")
    
    test_aa_1 = aa_3_to_1(test_aa_3)
    print (f"test_aa_3: {test_aa_3}\ttest_aa_1: {test_aa_1}")

    test_aa = "K"
    if len(test_aa) == 1 and is_valid_aa(test_aa):
        print (f"test_aa: {test_aa} is valid")
    else:
        print (f"test_aa: {test_aa} is not valid")

    test_aa = "Z"
    if len(test_aa) == 1 and is_valid_aa(test_aa):
        print (f"test_aa: {test_aa} is valid")
    else:
        print (f"test_aa: {test_aa} is not valid")

    variants = [
        "S459F",
        "R161Q (c.482G>A)",
        "Splice Site (c.3028+1G>T)",
        "Mutation"
    ]

    for v in variants:
        vtype = guess_variant_type(v)
        vname_long = snv_coding_to_p_3letter(v)
        print(f"{v!r} -> {vtype}. Long form name: {vname_long}")

    test_snv_coding_variant_name = "S459F"
    ref_aa_1, pos, var_aa_1 = parse_snv_coding_name_components(test_snv_coding_variant_name)
    print(f"{test_snv_coding_variant_name} -> ref_aa_1: {ref_aa_1} / pos: {pos} / var_aa_1: {var_aa_1}")

    protein = "MSTNPKPQRKTKRNTNRRPQDVKFPGGG"
    print(f"\nTesting amino acid position matches for sequence: {protein}")
    print(f"  Does P at position 5 match? {reference_aa_positions_matches('P', 5, protein, 'test_id')}") #True
    print(f"  Does A at position 5 match? {reference_aa_positions_matches('A', 5, protein, 'test_id')}") #False
    print(f"  Does P at position 100 match? {reference_aa_positions_matches('P', 100, protein, 'test_id')}") #False


if __name__ == "__main__":
    
    main()




