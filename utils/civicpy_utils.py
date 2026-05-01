#!/usr/bin/env python3

from civicpy import civic

def extract_variant_info(v: list) -> dict:
    """Extract useful variant info from a civicpy variant object"""
    return {
        "variant_id": v.id,
        "feature_id": v.feature_id,
        "allele_registry_id": v.allele_registry_id,
        "variant_name": v.name,
        "variant_aliases": v.variant_aliases,
        "hgvs_expressions": v.hgvs_expressions,
        "gene": v.gene.name if v.gene else None,
        "entrez_name": v.entrez_name,
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

def extract_variant_id_list(variants: list) -> list[int]:
    """Produce a list of all variant IDs in a civicpy variants object, order by feature ID and return it"""
    variant_info = []

    for v in variants:
        info = extract_variant_info(v)
        variant_id = info['variant_id']
        feature_id = info['feature_id']
        variant_info.append((feature_id, variant_id))
    
    # Sort by feature_id (first element of tuple)
    variant_info.sort(key=lambda x: x[0])

    # Extract ordered variant_ids
    return [int(variant_id) for _, variant_id in variant_info]


def get_sources_for_variant(variant_id: int) -> dict[str, dict[str, str]]:
    """Starting with a variant id, get associated molecular profiles, and their associated evidence, and their associated sources """
    sources = {}
    
    variant = civic.get_variant_by_id(variant_id)
    molecular_profiles = variant.molecular_profiles
    for mp in molecular_profiles:
        evidence_items = mp.evidence_items
        for ei in evidence_items:
            source = ei.source
            source_url = source.source_url

            if source_url not in sources:
                sources[source_url] = {
                    "citation": source.citation,
                    "source_type": source.source_type,
                    "source_url": source_url,
                }

    return sources


def main():
    #include_list=None #to use every single variant regardless of status
    include_list = ['accepted', 'submitted']

    variants = civic.get_all_gene_variants(include_status=include_list, allow_cached=True)
    print(f"Total variants in CIViC: {len(variants)}\n")

    variant_ids = extract_variant_id_list(variants)

    for v in variants[:5]:
        info = extract_variant_info(v)
        print(info)
        print("-" * 80)
    print(len(variant_ids))

    #for a specific variant get all associated evidence sources
    variant_id = 1832
    sources = get_sources_for_variant(variant_id)



if __name__ == "__main__":
    main()

