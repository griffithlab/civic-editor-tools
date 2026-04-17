#!/usr/bin/env python3

class RevisionComparator:

    def __init__(self, clingen_data):
        self.clingen_data = clingen_data  # your externally gathered info

        self._dispatch = {
            "variant_type_ids":          self.compare_variant_types,
            "variant_alias_ids":         self.compare_variant_aliases,
            "hgvs_description_ids":      self.compare_hgvs_expressions,
            "clinvar_entry_ids":         self.compare_clinvar_ids,
            "reference_build":           self.compare_reference_build,
            "chromosome":                self.compare_chromosome,
            "start":                     self.compare_start,
            "stop":                      self.compare_stop,
            "reference_bases":           self.compare_reference_bases,
            "variant_bases":             self.compare_variant_bases,
            "representative_transcript": self.compare_representative_transcript,
            "ensembl_version":           self.compare_ensembl_version
            # add more field_names here as needed
        }

    def compare(self, field_name, revision_value):
        """method that matches a civic revision field to appropriate comparison logic method below"""
        handler = self._dispatch.get(field_name)
        if handler is None:
            raise NotImplementedError(f"No comparator defined for field: '{field_name}'")
        return handler(revision_value)

    def compare_variant_types(self, revision_value):
        # field specific logic here
        return revision_value == self.clingen_data["variant_type"]

    def compare_variant_aliases(self, revision_value):
        # field specific logic here
        return revision_value == self.clingen_data["variant_aliases"]

    def compare_hgvs_expressions(self, revision_value):
        # field specific logic here
        return revision_value == self.clingen_data["hgvs_expressions"]

    def compare_clinvar_ids(self, revision_value):
        # field specific logic here
        return revision_value == self.clingen_data["clinvar_ids"]

    def compare_reference_build(self, revision_value):
        # field specific logic here
        return revision_value == self.clingen_data["assembly"]

    def compare_chromosome(self, civic_chromosome):
        clingen_chromosome = self.clingen_data["chromosome"]
        if clingen_chromosome == civic_chromosome:
        	return True
        else:
            return False

    def compare_start(self, revision_value):
        # field specific logic here
        return revision_value == self.clingen_data["start"]

    def compare_stop(self, revision_value):
        # field specific logic here
        return revision_value == self.clingen_data["end"]

    def compare_reference_bases(self, revision_value):
        # field specific logic here
        return revision_value == self.clingen_data["ref_bases"]

    def compare_variant_bases(self, revision_value):
        # field specific logic here
        return revision_value == self.clingen_data["alt_bases"]

    def compare_representative_transcript(self, revision_value):
        # field specific logic here
        return revision_value == self.clingen_data["representative_transcript"]

    def compare_ensembl_version(self, revision_value):
        # field specific logic here
        return revision_value == self.clingen_data["ensembl_version"]


def main():

	#example info coming from ClinGen Allele Registry
    clingen_data = {
        "variant_aliases": ["Ser459Phe","S459F","Ser432Phe","S432F"],
        "chromosome": "chr12", 
        "hgvs_expressions": ["NP_006222.2:p.Ser459Phe", "ENST00000320574.10:c.1376C>T", "ENSP00000322570.5:p.Ser459Phe", "NM_006231.4:c.1376C>T"],
    }


    clingen_value = "chr1"

    #example revision info from CIViC
    civic_field_name = "chromosome"
    civic_revision_value = "1"

    comparator = RevisionComparator(clingen_data)

    is_consistent = comparator.compare(civic_field_name, civic_revision_value)

    if is_consistent:
        print(f"{civic_field_name}: clingen_value ({clingen_value}) matches civic_value ({civic_revision_value})")
    else:
        print(f"{civic_field_name}: clingen_value ({clingen_value}) does NOT match civic_value ({civic_revision_value})")

    return

if __name__ == "__main__":
    main()

