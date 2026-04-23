#!/usr/bin/env python3

import re
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

class ValueComparator:
    """
    A class that facilitates use of an arbitrary set of methods that compare civic values (accepted or revisions) to clingen info
    Each method will be called based on a field name that triggers use of a method with logic specific to that field
    """

    def __init__(self, clingen_data):
        self.clingen_data = clingen_data  # your externally gathered info

        self._dispatch = {
            "allele_registry_id":        self.compare_allele_registry_id,
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
    def _print_match(self, level: MatchLevel, message: str) -> None:
        """Helper method that prints out a message with color matched to the quality of the matching information"""
        rid = self.current_revision_id
        user = self.current_user_display_name
        color = _MATCH_COLORS[level]

        if rid:
            revision_text = f" [revision: {rid}, user: {user}]"
            print(f"{color}{message}\033[0m{revision_text}", end="")
        else:
            print(f"{color}{message}\033[0m", end="")

    def _normalize_hgvs(self, hgvs: str) -> str:
        """Helper method that strips sequence version number for loose comparison.
        e.g. 'ENST00000320574.10:c.1376C>T' -> 'ENST00000320574:c.1376C>T'
        """
        return re.sub(r'\.\d+(?=:)', '', hgvs)

    def compare(self, field_name: str, comparison_value: str | list, revision_id: int, user_display_name: str) -> bool:
        """Method that matches a civic field to an appropriate comparison logic method below"""
        self.current_field_name = field_name
        self.current_revision_id = revision_id
        self.current_user_display_name = user_display_name
        handler = self._dispatch.get(field_name)
        if handler is None:
            raise NotImplementedError(f"No comparator defined for field: '{field_name}'")
        return handler(comparison_value)

    def compare_allele_registry_id(self, civic_allele_registry_id: str) -> bool:
        """Method for comparison of accepted civic CAID to possible clingen CAID matched by variant name"""
        clingen_allele_registry_id = self.clingen_data["allele_registry_id"]
        field_name = self.current_field_name
        print(f"    {field_name}.",  end=" ")
        print(f"clingen value: {clingen_allele_registry_id}", end = " ")
        if clingen_allele_registry_id == civic_allele_registry_id:
            self._print_match(
                MatchLevel.MATCH, 
                f"matches civic value: {civic_allele_registry_id}.\n"
            )
            return True
        else:
            self._print_match(
                 MatchLevel.MISMATCH, 
                 f"mismatches civic value: {civic_allele_registry_id}.\n"
            )
            return False

    def compare_variant_types(self, civic_variant_type: list) -> bool:
        """Method for comparison of CIViC variant types (revision or accepted) to an expected variant type inferred from the variant name"""
        guessed_gene_variant_type = self.clingen_data["variant_type"]
        field_name = self.current_field_name
        civic_variant_type_string = ','.join(civic_variant_type)

        print(f"    {field_name}.", end=" ")
        print(f"expected value: {guessed_gene_variant_type}", end = " ")
        if guessed_gene_variant_type in civic_variant_type:
            self._print_match(
                MatchLevel.MATCH, 
                f"matches civic value: {civic_variant_type_string}.\n"
            )
            return True
        else:
            self._print_match(
                 MatchLevel.MISMATCH, 
                 f"mismatches civic value: {civic_variant_type_string}.\n"
            )
            return False

    def compare_variant_aliases(self, civic_variant_aliases: list) -> bool:
        """Method for comparison of CIViC variant aliases (revision or accepted) to a list of possible aliases obtained from ClinGen"""
        civic_variant_aliases_string = ','.join(civic_variant_aliases)
        clingen_variant_aliases = self.clingen_data["variant_aliases"]
        clingen_variant_aliases_string = ','.join(clingen_variant_aliases)
        clingen_variant_aliases_upper = {a.upper() for a in clingen_variant_aliases}
        field_name = self.current_field_name

        #find civic aliases that are in the clingen possible aliases list and those that are not
        matched_civic_aliases = []
        unmatched_civic_aliases = []

        for alias in civic_variant_aliases:
            if alias.upper() in clingen_variant_aliases_upper:
                matched_civic_aliases.append(alias)
            else:
                unmatched_civic_aliases.append(alias)

        matched_civic_aliases_string = ','.join(matched_civic_aliases)
        unmatched_civic_aliases_string = ','.join(unmatched_civic_aliases)

        print(f"    {field_name}.", end=" ")
        print(f"clingen values: {clingen_variant_aliases_string}", end=" ")
        if matched_civic_aliases:
            self._print_match(
                MatchLevel.MATCH, 
                f"\n      matches civic values: {matched_civic_aliases_string}."
            )
        if unmatched_civic_aliases:
            self._print_match(
                MatchLevel.MISMATCH, 
                f"\n      mismatches civic values: {unmatched_civic_aliases_string}."
            )

        return len(unmatched_civic_aliases) == 0

    def compare_hgvs_expressions(self, civic_hgvs_expressions: list) -> bool:
        """Method for comparison of CIViC hgvs expressions (revision or accepted) to a list of possible hgvs expressions obtained from ClinGen"""
        civic_hgvs_expressions_string = ','.join(civic_hgvs_expressions)
        clingen_hgvs_expressions = self.clingen_data["hgvs_expressions"]
        clingen_hgvs_expressions_string = ','.join(clingen_hgvs_expressions)
        field_name = self.current_field_name

        #find civic aliases that are in the clingen possible aliases list and those that are not
        matched_civic_hgvs_expressions = []
        unmatched_civic_hgvs_expressions = []

        for civic_hgvs in civic_hgvs_expressions:
            civic_match_found = False
            for clingen_hgvs in clingen_hgvs_expressions:
                if self._normalize_hgvs(civic_hgvs) == self._normalize_hgvs(clingen_hgvs):
                    civic_match_found = True
                    matched_civic_hgvs_expressions.append(civic_hgvs)
            if not civic_match_found:
                unmatched_civic_hgvs_expressions.append(civic_hgvs)

        matched_civic_hgvs_expressions_string = ','.join(matched_civic_hgvs_expressions)
        unmatched_civic_hgvs_expressions_string = ','.join(unmatched_civic_hgvs_expressions)

        print(f"    {field_name}.", end=" ")
        print(f"clingen values: {clingen_hgvs_expressions_string}", end=" ")
        if matched_civic_hgvs_expressions:
            self._print_match(
                MatchLevel.MATCH, 
                f"matches civic values: {matched_civic_hgvs_expressions_string}."
            )
        if unmatched_civic_hgvs_expressions:
            self._print_match(
                MatchLevel.MISMATCH, 
                f"mismatches civic values: {unmatched_civic_hgvs_expressions_string}."
            )

        return len(unmatched_civic_hgvs_expressions) == 0

    def compare_clinvar_ids(self, civic_clinvar_ids: list) -> bool:
        """Method that compares clinvar IDs from a CIViC submission (revision or accepted) to those from ClinGen Allele Registry"""
        #normalize clinvar ids to a list of strings
        clingen_clinvar_ids = [str(id) for id in self.clingen_data["clinvar_ids"]]
        civic_clinvar_ids = [str(id) for id in civic_clinvar_ids]
        clingen_clinvar_ids_string = ','.join(clingen_clinvar_ids)
        field_name = self.current_field_name

        #find civic clinvar IDs that are in the clingen clinvar list and those that are not
        matched_civic_clinvar_ids = []
        unmatched_civic_clinvar_ids = []

        for civic_clinvar_id in civic_clinvar_ids:
            if civic_clinvar_id in clingen_clinvar_ids:
                matched_civic_clinvar_ids.append(civic_clinvar_id)
            else:
                unmatched_civic_clinvar_ids.append(civic_clinvar_id)

        matched_civic_clinvar_ids_string = ','.join(matched_civic_clinvar_ids)
        unmatched_civic_clinvar_ids_string = ','.join(unmatched_civic_clinvar_ids)

        print(f"    {field_name}.", end=" ")
        print(f"clingen values: {clingen_clinvar_ids_string}", end = " ")
        if matched_civic_clinvar_ids:
            self._print_match(
                MatchLevel.MATCH, 
                f"matches civic values: {matched_civic_clinvar_ids_string}."
            )
        if unmatched_civic_clinvar_ids:
            self._print_match(
                MatchLevel.MISMATCH, 
                f"mismatches civic values: {unmatched_civic_clinvar_ids_string}."
            )

        return len(unmatched_civic_clinvar_ids) == 0

    def compare_reference_build(self, civic_reference_build: str) -> bool:
        """Method that compares the reference build in a CIViC submission (revision or accepted) to the expected value"""
        expected_reference_build = self.clingen_data["reference_build"]
        field_name = self.current_field_name

        print(f"    {field_name}.", end=" ")
        print(f"expected value: {expected_reference_build}", end=" ")
        if expected_reference_build.upper() == civic_reference_build.upper():
            self._print_match(
                MatchLevel.MATCH,
                f"matches civic value: {civic_reference_build}."
            )
            return True
        else:
            self._print_match(
                MatchLevel.MISMATCH,
                f"mismatches civic value: {civic_reference_build}."
            )
            return False

    def compare_chromosome(self, civic_chromosome: str) -> bool:
        """Method that compares the chromosome values from a CIViC submission (revision or accepted) to one from ClinGen Allele Registry"""
        clingen_chromosome = self.clingen_data["chromosome"]
        clingen_chromosome_normalized = clingen_chromosome.removeprefix("chr")
        field_name = self.current_field_name
 
        print(f"    {field_name}.", end=" ")
        print(f"clingen value: {clingen_chromosome}", end=" ")
        if clingen_chromosome_normalized == civic_chromosome:
            self._print_match(
                MatchLevel.MATCH, 
                f"matches civic value: {civic_chromosome}."
            )
            return True
        else:
            self._print_match(
                 MatchLevel.MISMATCH, 
                 f"mismatches civic value: {civic_chromosome}."
            )
            return False

    def compare_start(self, civic_start: int) -> bool:
        """
        Method that compares the genomic start position value from a CIViC submission (revision or accepted) to one from ClinGen Allele Registry
        Since clingen start positions are 0-based and civic start positions are 1 based, it adjusts for this
        """
        clingen_start = self.clingen_data["start"] + 1 
        field_name = self.current_field_name

        print(f"    {field_name}.", end=" ")
        print(f"clingen value (+1): {clingen_start}", end=" ")
        if clingen_start == civic_start:
            self._print_match(
                MatchLevel.MATCH,
                f"matches civic value: {civic_start}."
            )
            return True
        else:
            self._print_match(
                MatchLevel.MISMATCH,
                f"mismatches civic value: {civic_start}."
            )
            return False

    def compare_stop(self, civic_stop: int) -> bool:
        """Method that compares the genomic stop position value from a CIViC submission (revision or accepted) to one from ClinGen Allele Registry"""
        clingen_end = self.clingen_data["end"]
        field_name = self.current_field_name

        print(f"    {field_name}.", end=" ")
        print(f"clingen value: {clingen_end}", end=" ")
        if clingen_end == civic_stop:
            self._print_match(
                MatchLevel.MATCH,
                f"matches civic value: {civic_stop}."
            )
            return True
        else:
            self._print_match(
                MatchLevel.MISMATCH,
                f"mismatches civic value: {civic_stop}."
            )
            return False

    def compare_reference_bases(self, civic_reference_bases: str) -> bool:
        """Method that compares the reference bases value from a CIViC submission (revision or accepted) to one from ClinGen Allele Registry"""
        clingen_reference_bases = self.clingen_data["ref_bases"]
        field_name = self.current_field_name

        print(f"    {field_name}.", end=" ")
        print(f"clingen value: {clingen_reference_bases}", end=" ")
        if clingen_reference_bases == civic_reference_bases:
            self._print_match(
                MatchLevel.MATCH,
                f"matches civic value: {civic_reference_bases}."
            )
            return True
        else:
            self._print_match(
                MatchLevel.MISMATCH,
                f"mismatches civic value: {civic_reference_bases}."
            )
            return False

    def compare_variant_bases(self, civic_variant_bases: str) -> bool:
        """Method that compares the variant bases value from a CIViC submission (revision or accepted) to one from ClinGen Allele Registry"""
        clingen_variant_bases = self.clingen_data["alt_bases"]
        field_name = self.current_field_name

        print(f"    {field_name}.", end=" ")
        print(f"clingen value: {clingen_variant_bases}", end=" ")
        if clingen_variant_bases == civic_variant_bases:
            self._print_match(
                MatchLevel.MATCH,
                f"matches civic value: {civic_variant_bases}."
            )
            return True
        else:
            self._print_match(
                MatchLevel.MISMATCH,
                f"mismatches civic value: {civic_variant_bases}."
            )
            return False

    def compare_representative_transcript(self, civic_representative_transcript: str) -> bool:
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
        print(f"    {field_name}.", end=" ")
        print(f"clingen value: {v75_match_transcript}(v75) or {v87_match_transcript}(v87)", end=" ")
        if v75_match_result or v87_match_result:
            self._print_match(
                MatchLevel.MATCH, 
                f"matches civic value: {civic_representative_transcript}."
            )
            return True

        if v75_partial_match_result or v87_partial_match_result:
            self._print_match(
                MatchLevel.QUALIFIED_MATCH, 
                f"partially matches civic value: {civic_representative_transcript} (but no exact match)."
            )
            return True

        self._print_match(
            MatchLevel.MISMATCH, 
            f"No clingen match found for civic value: {civic_representative_transcript}."
        )
        return False

    def compare_ensembl_version(self, civic_ensembl_version: int) -> bool:
        """
        Method that compares CIViC revision value for Ensembl version against 
        a basic hard coded expectation of likely Ensembl version numbers
        """

        expected_ensembl_versions = self.clingen_data["ensembl_version"]
        expected_ensembl_versions_string = ','.join(expected_ensembl_versions)
        field_name = self.current_field_name

        print(f"    {field_name}.", end=" ")
        print(f"expected values: {expected_ensembl_versions_string}", end=" ")
        if str(civic_ensembl_version) in expected_ensembl_versions:
            self._print_match(
                MatchLevel.MATCH, 
                f"matches civic value: {civic_ensembl_version}."
            )
            return True
        
        if str(civic_ensembl_version) not in expected_ensembl_versions:
            self._print_match(
                MatchLevel.MISMATCH, 
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

    comparator = ValueComparator(clingen_data)

    is_consistent = comparator.compare(civic_field_name, civic_revision_value, civic_revision_id, civic_user)

    return

if __name__ == "__main__":
    main()

