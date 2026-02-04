#!/usr/bin/env python3

from civicpy import civic

def extract_variant_info(v):
    return {
        "variant_id": v.id,
        "allele_registry_id": v.allele_registry_id,
        "variant_name": v.name,
        "variant_aliases": v.variant_aliases,
        "hgvs_expressions": v.hgvs_expressions,
        "gene": v.gene.name if v.gene else None,
        "variant_types": [vt.name for vt in v.variant_types],
        "subtype": v.subtype,
        "clinvar_entries": v.clinvar_entries,
        "ref_build": v.coordinates.reference_build,
        "chromosome": v.coordinates.chromosome,
        "start": v.coordinates.start,
        "stop": v.coordinates.stop,
        "reference_bases": v.coordinates.reference_bases,
        "variant_bases": v.coordinates.variant_bases,
        "representative_transcript": v.coordinates.representative_transcript,

    }

def main():
    variants = civic.get_all_gene_variants(include_status=['accepted', 'submitted'], allow_cached=True)
    print(f"Total variants in CIViC: {len(variants)}\n")

    for v in variants[:100]:
        info = extract_variant_info(v)
        print(info)
        print("-" * 80)

if __name__ == "__main__":
    main()

