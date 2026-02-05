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

def extract_variant_id_list(variants):
	variant_ids = []
	for v in variants:
		id = v.id
		variant_ids.append(id)
	return variant_ids

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

if __name__ == "__main__":
    main()

