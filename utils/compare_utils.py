#!/usr/bin/env python3

from enum import Enum

from enum import Enum

class MatchLevel(Enum):
    """A class that defines colors to print out summary text based on the quality of matching information"""
    MATCH = "match"
    QUALIFIED_MATCH = "qualified_match"
    MISMATCH = "mismatch"

_MATCH_COLORS = {
    MatchLevel.MATCH: "\033[32m",           # green
    MatchLevel.QUALIFIED_MATCH: "\033[33m", # yellow
    MatchLevel.MISMATCH: "\033[31m",        # red
}

class RevisionComparator:
    """A class that facilitates use of an arbitrary set of methods that compare civic revision info to clingen info"""

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
    def _print_match(self, level: MatchLevel, message: str):
        """Helper method that prints out a message with color matched to the quality of the matching information"""
        color = _MATCH_COLORS[level]
        print(f"{color}{message}\033[0m")

    def compare(self, field_name, revision_value, revision_id):
        """Method that matches a civic revision field to appropriate comparison logic method below"""
        self.current_field_name = field_name
        self.current_revision_id = revision_id
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
        """Method that compares the chromosome values from a CIViC revision to one from ClinGen Allele Registry """
        clingen_chromosome = self.clingen_data["chromosome"]
        clingen_chromosome_normalized = clingen_chromosome.removeprefix("chr")
        field_name = self.current_field_name
        rid = self.current_revision_id

        if clingen_chromosome_normalized == civic_chromosome:
            self._print_match(MatchLevel.MATCH, f"  {self.current_field_name} (revision: {self.current_revision_id}). clingen_value: ({clingen_chromosome}) matches civic_value: ({civic_chromosome})")
            return True
        else:
            self._print_match(MatchLevel.MISMATCH, f"  {self.current_field_name} (revision: {self.current_revision_id}). clingen_value: ({clingen_chromosome}) mismatch civic_value: ({civic_chromosome})")
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


    #example revision info from CIViC
    civic_revision_id = 1
    civic_field_name = "chromosome"
    civic_revision_value = "12"

    comparator = RevisionComparator(clingen_data)

    is_consistent = comparator.compare(civic_field_name, civic_revision_value, civic_revision_id)

    return

if __name__ == "__main__":
    main()

