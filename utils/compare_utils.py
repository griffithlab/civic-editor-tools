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
    """A class that facilitates use of an arbitrary set of methods that compare civic values (accepted or revisions) to clingen info"""

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
        rid = self.current_revision_id
        user = self.current_user_display_name
        color = _MATCH_COLORS[level]

        if rid:
            revision_text = f" [revision: {rid}, user: {user}]"
            message = message + revision_text

        print(f"{color}{message}\033[0m")

    def compare(self, field_name, comparison_value, revision_id, user_display_name):
        """Method that matches a civic field to an appropriate comparison logic method below"""
        self.current_field_name = field_name
        self.current_revision_id = revision_id
        self.current_user_display_name = user_display_name
        handler = self._dispatch.get(field_name)
        if handler is None:
            raise NotImplementedError(f"No comparator defined for field: '{field_name}'")
        return handler(comparison_value)

    def compare_variant_types(self, civic_variant_type):
        guessed_gene_variant_type = self.clingen_data["variant_type"]
        field_name = self.current_field_name
        civic_variant_type_string = ','.join(civic_variant_type)

        if guessed_gene_variant_type in civic_variant_type:
            self._print_match(
                MatchLevel.MATCH, 
                f"    {field_name}. expected value: {guessed_gene_variant_type} "
                f"matches civic value: {civic_variant_type_string}."
            )
            return True
        else:
            self._print_match(
                 MatchLevel.MISMATCH, 
                 f"    {field_name}. expected value: {guessed_gene_variant_type} "
                 f"mismatches civic value: {civic_variant_type_string}."
            )
            return False

    def compare_variant_aliases(self, comparison_value):
        # field specific logic here
        return comparison_value == self.clingen_data["variant_aliases"]

    def compare_hgvs_expressions(self, comparison_value):
        # field specific logic here
        return comparison_value == self.clingen_data["hgvs_expressions"]

    def compare_clinvar_ids(self, comparison_value):
        # field specific logic here
        return comparison_value == self.clingen_data["clinvar_ids"]

    def compare_reference_build(self, comparison_value):
        # field specific logic here
        return comparison_value == self.clingen_data["assembly"]

    def compare_chromosome(self, civic_chromosome):
        """Method that compares the chromosome values from a CIViC submission (revision or accepted) to one from ClinGen Allele Registry """
        clingen_chromosome = self.clingen_data["chromosome"]
        clingen_chromosome_normalized = clingen_chromosome.removeprefix("chr")
        field_name = self.current_field_name

        if clingen_chromosome_normalized == civic_chromosome:
            self._print_match(
                MatchLevel.MATCH, 
                f"    {field_name}. clingen value: {clingen_chromosome} "
                f"matches civic value: {civic_chromosome}."
            )
            return True
        else:
            self._print_match(
                 MatchLevel.MISMATCH, 
                 f"    {field_name}. clingen value: {clingen_chromosome} "
                 f"mismatches civic value: {civic_chromosome}."
            )
            return False

    def compare_start(self, comparison_value):
        # field specific logic here
        return comparison_value == self.clingen_data["start"]

    def compare_stop(self, comparison_value):
        # field specific logic here
        return comparison_value == self.clingen_data["end"]

    def compare_reference_bases(self, comparison_value):
        # field specific logic here
        return comparison_value == self.clingen_data["ref_bases"]

    def compare_variant_bases(self, comparison_value):
        # field specific logic here
        return comparison_value == self.clingen_data["alt_bases"]

    def compare_representative_transcript(self, civic_representative_transcript):
        """
        Method that compares the representative transcript entry from CIViC (accepted or revision) to values from ClinGen Allele Registry.
        The ClinGen values come from "valid" transcripts that would make sense for the CIViC variant
        From these, matching build37 versioned ensembl transcripts are obtained from Ensembl v75 and build37 imported Ensembl v87
        """
        # field specific logic here

        variant_build37_ensembl_transcripts = self.clingen_data["representative_transcript"]
        field_name = self.current_field_name
        civic_base = civic_representative_transcript.split(".")[0]

        v75_partial_match_result = False
        v75_partial_match_transcript = None
        v87_partial_match_result = False
        v87_partial_match_transcript = None
        v75_match_result = False
        v75_match_transcript = None
        v87_match_result = False
        v87_match_transcript = None

        for clingen_id, v75_match, v87_match in variant_build37_ensembl_transcripts:
            #check for exact matches against either v75 or v87
            if civic_representative_transcript == v75_match:
                v75_match_result = True
                v75_match_transcript = v75_match
            if civic_representative_transcript == v87_match:
                v87_match_result = True
                v87_match_transcript = v87_match

            #check for partial matches (base ID only) against either v75 or v87
            if v75_match is not None and v75_match.split(".")[0] == civic_base:
                v75_partial_match_result = True
                v75_partial_match_transcript = v75_match
            if v87_match is not None and v87_match.split(".")[0] == civic_base:
                v87_partial_match_result = True
                v87_partial_match_transcript = v87_match

        #evaluate the match results
        if v75_match_result or v87_match_result:
            self._print_match(
                MatchLevel.MATCH, 
                f"    {field_name}. clingen value: {v75_match_transcript}(v75) or {v87_match_transcript}(v87) "
                f"matches civic value: {civic_representative_transcript}."
            )
            return True

        if v75_partial_match_result or v87_partial_match_result:
            self._print_match(
                MatchLevel.QUALIFIED_MATCH, 
                f"    {field_name}. clingen value: {v75_partial_match_transcript}(v75) or {v87_partial_match_transcript}(v87) "
                f"partially matches civic value: {civic_representative_transcript} (but no exact match)."
            )
            return True

        self._print_match(
            MatchLevel.MISMATCH, 
            f"    {field_name}. No clingen match found for civic value: {civic_representative_transcript}."
        )
        return False


    def compare_ensembl_version(self, civic_ensembl_version):
        """
        Method that compares CIViC revision value for Ensembl version against 
        a basic hard coded expectation of likely Ensembl version numbers
        """

        expected_ensembl_versions = self.clingen_data["ensembl_version"]
        expected_ensembl_versions_string = ','.join(expected_ensembl_versions)
        field_name = self.current_field_name

        if str(civic_ensembl_version) in expected_ensembl_versions:
            self._print_match(
                MatchLevel.MATCH, 
                f"    {field_name}. expected values: {expected_ensembl_versions_string} "
                f"matches civic value: {civic_ensembl_version}."
            )
            return True
        
        if str(civic_ensembl_version) not in expected_ensembl_versions:
            self._print_match(
                MatchLevel.MISMATCH, 
                f"    {field_name}. expected values: {expected_ensembl_versions_string} "
                f"mismatches civic value: {civic_ensembl_version}."
            )
 
        return False


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
    civic_user = "MalachiGriffith"

    comparator = RevisionComparator(clingen_data)

    is_consistent = comparator.compare(civic_field_name, civic_revision_value, civic_revision_id)

    return

if __name__ == "__main__":
    main()

