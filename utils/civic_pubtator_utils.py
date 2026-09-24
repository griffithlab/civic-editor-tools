#!/usr/bin/env python3

import re
import csv
from pathlib import Path

base_dir = Path(__file__).resolve().parent

REPORT_FILENAME_PATTERN = re.compile(r"^report_(.+)\.tsv$")


def resolve_pub_reports_dir(civic_pubtator_data_path) -> Path:
    """
    Resolve the pub-reports directory from a supplied civic-pubtator-data path.
    Accepts either the pub-reports directory itself, or the repository root one
    level above it (in which case a 'pub-reports' subdirectory is expected).
    """
    path = Path(civic_pubtator_data_path)

    if path.name == "pub-reports" and path.is_dir():
        return path

    candidate = path / "pub-reports"
    if candidate.is_dir():
        return candidate

    raise FileNotFoundError(
        f"Could not find a 'pub-reports' directory at '{path}' or '{candidate}'"
    )


def find_report_tsv_files(pub_reports_dir) -> list:
    """
    Find all report TSV files in a pub-reports directory.
    Returns a list of (publication_id, file_path) tuples, sorted by publication id.
    """
    pub_reports_dir = Path(pub_reports_dir)
    report_files = []

    for tsv_path in pub_reports_dir.glob("report_*.tsv"):
        match = REPORT_FILENAME_PATTERN.match(tsv_path.name)
        if not match:
            continue
        publication_id = match.group(1)
        report_files.append((publication_id, tsv_path))

    report_files.sort(key=lambda entry: entry[0])

    return report_files


def variant_row_qualifies(entity_category: str, entity_type: str) -> bool:
    """
    Test whether a TSV row's entity_category/entity_type combination represents
    a candidate variant mention worth capturing.
    """
    if entity_category == "Variant":
        return True

    if entity_category == "NER_AIONER" and entity_type == "Mutation":
        return True

    return False


def parse_variant_entries(tsv_path) -> list:
    """
    Parse a single pub-reports TSV file and return the rows that represent candidate
    variant mentions: entity_category == "Variant", or entity_category == "NER_AIONER"
    with entity_type == "Mutation".

    Each returned entry captures the mention and hgvs columns (either may be blank,
    but not both), plus the source and count columns.
    """
    tsv_path = Path(tsv_path)
    variant_entries = []

    with tsv_path.open("r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        for row in reader:
            entity_category = row["entity_category"]
            entity_type = row["entity_type"]

            if not variant_row_qualifies(entity_category, entity_type):
                continue

            variant_entries.append({
                "entity_category": entity_category,
                "entity_type": entity_type,
                "mention": row["mention"] or None,
                "hgvs": row["hgvs"] or None,
                "source": row["source"] or None,
                "count": int(row["count"]) if row["count"] else None,
            })

    return variant_entries


def extract_variant_names(variant_entries: list) -> list:
    """
    Given a list of variant entries (as returned by parse_variant_entries), produce a
    deduplicated list of candidate variant name strings drawn from the mention and
    hgvs columns, preserving first-seen order.
    """
    seen = set()
    variant_names = []

    for entry in variant_entries:
        for name in (entry["mention"], entry["hgvs"]):
            if name and name not in seen:
                seen.add(name)
                variant_names.append(name)

    return variant_names


def load_all_publication_variants(civic_pubtator_data_path) -> dict:
    """
    Given a path to a civic-pubtator-data repository (or its pub-reports directory
    directly), parse every report TSV file and return a dict keyed by publication id,
    each value being the list of qualifying variant entries found in that publication.
    """
    pub_reports_dir = resolve_pub_reports_dir(civic_pubtator_data_path)
    report_files = find_report_tsv_files(pub_reports_dir)

    publication_variants = {}
    for publication_id, tsv_path in report_files:
        publication_variants[publication_id] = parse_variant_entries(tsv_path)

    return publication_variants


def build_report_file_index(civic_pubtator_data_path) -> dict:
    """
    Resolve the pub-reports directory once and return a dict of {publication_id: tsv_path}
    for every report TSV file found.

    Build this once per run (e.g. at the start of a batch review) and pass it into repeated
    calls to summarize_variant_name_matches(), rather than letting each call re-resolve the
    pub-reports directory and re-glob its files. The civic-pubtator-data repository this reads
    from is expected to grow over time, so re-listing it on every publication/variant checked
    would get more expensive as it does.
    """
    pub_reports_dir = resolve_pub_reports_dir(civic_pubtator_data_path)
    return dict(find_report_tsv_files(pub_reports_dir))


def summarize_variant_name_matches(report_file_index: dict, publication_ids: list, variant_names: list) -> dict:
    """
    For each requested publication, find exact matches (mention or hgvs) to each name in
    variant_names (typically the same underlying variant expressed in different naming
    styles, e.g. 'Pro81Ser', 'p.P81S', 'P81S'), and sum the count column across all
    qualifying entries that match each name. Each variant name's total is computed
    independently of the others.

    Each TSV row is a single named entity recognition (NER) record; its mention and hgvs
    values are both checked against a given variant name, but a row is credited to that
    variant name's total at most once even if both its mention and hgvs equal that same
    variant name.

    report_file_index should be built once via build_report_file_index() and reused across
    calls, rather than rebuilt for every publication/variant checked.

    Returns a dict of dicts: {publication_id: {variant_name: total_match_count}}.
    A publication_id with no report TSV file on disk is omitted from the returned dict
    (a warning is printed in that case).
    """
    match_summary = {}
    for publication_id in publication_ids:
        tsv_path = report_file_index.get(publication_id)

        if tsv_path is None:
            print(f"    WARNING: No report TSV file found for publication {publication_id}", end="")
            continue

        variant_name_counts = {variant_name: 0 for variant_name in variant_names}
        variant_entries = parse_variant_entries(tsv_path)
        for entry in variant_entries:
            for variant_name in variant_names:
                if entry["mention"] == variant_name or entry["hgvs"] == variant_name:
                    variant_name_counts[variant_name] += entry["count"] or 0

        match_summary[publication_id] = variant_name_counts

    return match_summary


def format_variant_name_match_summary(match_summary: dict, variant_names: list) -> str:
    """
    Format the output of summarize_variant_name_matches() into a human readable summary, e.g.:
    pmid1: varname1 (n=10 matches); varname2 (n=22 matches); varname3 (n=5 matches)
    pmid2: varname1 (n=0 matches); varname2 (n=1 matches); varname3 (n=0 matches)
    """
    lines = []
    for publication_id, variant_name_counts in match_summary.items():
        parts = [f"{variant_name} ({variant_name_counts[variant_name]})" for variant_name in variant_names]
        lines.append(f"PMID:{publication_id} matches: " + "; ".join(parts))

    return "\n".join(lines)


def main():

    #example path to a local checkout of civic-pubtator-data
    civic_pubtator_data_path = "/Users/mgriffit/git/civic-pubtator-data"

    pub_reports_dir = resolve_pub_reports_dir(civic_pubtator_data_path)
    print(f"Resolved pub-reports directory: {pub_reports_dir}")

    report_files = find_report_tsv_files(pub_reports_dir)
    print(f"Found {len(report_files)} report TSV file(s)")

    #find the first report file that actually has qualifying variant entries, for a useful demo
    test_publication_id, test_tsv_path, variant_entries = None, None, []
    for publication_id, tsv_path in report_files:
        entries = parse_variant_entries(tsv_path)
        if entries:
            test_publication_id, test_tsv_path, variant_entries = publication_id, tsv_path, entries
            break

    print(f"\nParsing variant entries for publication {test_publication_id} ({test_tsv_path.name})")
    print(f"Found {len(variant_entries)} qualifying variant entrie(s)")
    for entry in variant_entries[:10]:
        print(f"  {entry}")

    variant_names = extract_variant_names(variant_entries)
    print(f"\nUnique candidate variant names ({len(variant_names)}):")
    print(f"  {', '.join(variant_names[:20])}")

    print(f"\nLoading variant data for all publications (this may take a while)...")
    publication_variants = load_all_publication_variants(civic_pubtator_data_path)
    total_entries = sum(len(entries) for entries in publication_variants.values())
    print(f"Loaded {len(publication_variants)} publication(s), {total_entries} total qualifying variant entrie(s)")

    #summarize matches for the same underlying variant (VHL Pro81Ser) expressed in different naming styles
    #across a couple of publications known to mention it, plus one with no matches, to demonstrate n=0 output
    test_publication_ids = ["10340905", "11409863", "1"]
    test_variant_names = ["Pro81Ser", "p.P81S", "P81S"]

    report_file_index = build_report_file_index(civic_pubtator_data_path)

    print(f"\nSummarizing variant name matches for publications {test_publication_ids} and variant names {test_variant_names}:")
    match_summary = summarize_variant_name_matches(report_file_index, test_publication_ids, test_variant_names)
    #summarize_variant_name_matches() prints any "not found" warnings inline (no trailing newline), so start
    #the formatted results on a fresh line rather than running on from a warning left dangling above
    print(f"\n{format_variant_name_match_summary(match_summary, test_variant_names)}")


if __name__ == "__main__":
    main()
