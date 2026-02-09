#!/usr/bin/env python3

import re

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
    snv_pattern = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])")

    match = snv_pattern.match(variant)
    if match:
        ref_aa, pos, alt_aa = match.groups()

        if is_valid_aa(ref_aa) and is_valid_aa(alt_aa):
            return "snv_coding"
    
    # 3) splice site
    if variant.startswith("Splice Site"):
        return "splice_site"

    # None of the categories above were guessed from the variant name
    return None


def snv_coding_to_p_3letter(variant):
    """Convert a coding SNV like 'S459F' to 'p.Ser459Phe'"""
    
    if guess_variant_type(variant) != "snv_coding":
        return None

    match = re.match(r"^([A-Za-z])(\d+)([A-Za-z])", variant)
    if not match:
        return None

    ref_aa, pos, alt_aa = match.groups()

    return f"p.{aa_1_to_3(ref_aa)}{pos}{aa_1_to_3(alt_aa)}"


def main():
    
    test_aa_1 = "W"
    test_aa_3 = aa_1_to_3(test_aa_1)
    print (f"test_aa_1: {test_aa_1}\ttest_aa_3: {test_aa_3}")

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


if __name__ == "__main__":
    
    main()




