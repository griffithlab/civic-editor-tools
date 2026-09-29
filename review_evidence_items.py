#!/usr/bin/env python3

"""
review_evidence_items.py

Command-line tool to review CIViC evidence items with outstanding (open) revisions.

This script retrieves evidence item data from CIViC (via CIViCpy and the CIViC
GraphQL API) and supports review workflows tied to a specific CIViC contributor.
"""

import os
import argparse
import sys
import webbrowser
from civicpy import civic
from pathlib import Path

# Local utility imports
from utils import generic_utils
from utils import civic_graphql_utils
from utils import civicpy_utils

base_dir = Path(__file__).resolve().parent

GREEN = "\033[32m"
YELLOW = "\033[33m"
RED = "\033[31m"
BLUE = "\033[34m"
RESET = "\033[0m"

def parse_args():
    """Obtain command line arguments from the user"""
    parser = argparse.ArgumentParser(
        description=(
            "Review pending revisions for one or more CIViC evidence items."
        )
    )

    parser.add_argument(
        "--contributor-id",
        dest="contributor_id",
        type=int,
        required=True,
        help="CIViC contributor ID performing the review (integer, e.g. 15)"
    )

    evidence_choice_group = parser.add_mutually_exclusive_group(required=True)
    evidence_choice_group.add_argument(
        "--evidence-id",
        dest="evidence_id",
        type=int,
        help="CIViC evidence item ID to review (integer, e.g. 1832)"
    )
    evidence_choice_group.add_argument(
        "--all-evidence",
        dest="all_evidence",
        action="store_true",
        help="Review all CIViC evidence items"
    )
    evidence_choice_group.add_argument(
        "--evidence-list-file",
        dest="evidence_list_file",
        type=str,
        help="Path to file with list of CIViC evidence item IDs to process, one per line (first column), rows with # will be ignored"
    )
    parser.add_argument(
        "--allow-evidence-without-revisions",
        dest="allow_evidence_without_revisions",
        action="store_true",
        help="Even if an evidence item has no outstanding revisions, process it anyway"
    )
    parser.add_argument(
        "--processed-evidence-file",
        dest="processed_evidence_file",
        type=str,
        help="Path to file with list of already processed evidence item IDs, one per line (first column), rows with # will be ignored, file will be updated as new evidence items are processed"
    )
    parser.add_argument(
        "--open-browser",
        dest="open_browser",
        action="store_true",
        help="Allow the user to control whether a browser view will be opened for each evidence item"
    )

    return parser.parse_args()

def verify_connectivity():
    """Tests whether internet access is working and then if we can actually access the APIs needed before attempting anything """
    if not generic_utils.check_connection():
        print("No internet access. Aborting.")
        sys.exit(1)

    if not generic_utils.check_apis(api_urls=["https://www.civicdb.org"]):
        print("Required APIs are unavailable. Aborting.")
        sys.exit(1)

    print("Internet and API connectivity verified.")


def open_evidence_revision(eid: int) -> None:
    url = f"https://civicdb.org/evidence/{eid}/revisions"
    webbrowser.open(url)


def prompt_to_proceed(message: str = None) -> None:
    print("\n" + "=" * 80)
    if message:
        print(message)
    print("Press Enter to proceed or Ctrl+C to cancel...")
    print("=" * 80)
    try:
        input()
    except KeyboardInterrupt:
        print("\nAborted.")
        sys.exit(0)


def load_processed_evidence(processed_evidence_file: str, already_processed_evidence: set) -> set:
    """Load a list of already processed evidence item IDs"""

    if os.path.exists(processed_evidence_file):
        with open(processed_evidence_file, "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    already_processed_evidence.add(int(line))

    print(f"Loaded {len(already_processed_evidence)} already-processed evidence items. These will be silently skipped")

    return already_processed_evidence


def mark_evidence_processed(processed_evidence_file, eid, out_f, already_processed_evidence):
    """Add an evidence item to the processed file"""

    if not processed_evidence_file:
        return

    if out_f:
        out_f.write(f"{eid}\n")
        out_f.flush()
    already_processed_evidence.add(eid)


def get_evidence_ids_to_process(evidence_id, all_evidence, evidence_list_file):
    """Determine which evidence item IDs to work with based on user supplied choices """
    evidence_ids_to_process = []
    include_list = ['accepted', 'submitted']

    #single id
    if evidence_id:
        evidence_ids_to_process.append(evidence_id)

    #all evidence items in civic
    if all_evidence:
        evidence_items = civic.get_all_evidence(include_status=include_list, allow_cached=True)
        evidence_ids_to_process = civicpy_utils.extract_evidence_id_list(evidence_items)
        print(f"Total evidence item ids obtained from CIViCpy: {len(evidence_ids_to_process)}\n")

    #a list of specific evidence item ids, optionally with a header row (detected, not assumed)
    if evidence_list_file:
        with open(evidence_list_file, "r") as f:
            seen_data_line = False
            for line_num, line in enumerate(f, start=1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                try:
                    evidence_id_value = int(parts[0])
                except ValueError:
                    #the first data line is allowed to be a non-numeric header (e.g. "id"); skip it silently
                    if not seen_data_line:
                        print(f"Skipping header row in {evidence_list_file}: {line}")
                        continue
                    sys.exit(
                        f"Fatal: {evidence_list_file} line {line_num} does not start with an integer evidence ID: {line!r}\n"
                        f"  Expected one evidence ID per line (first tab-delimited column). An optional single header "
                        f"row is allowed, and lines starting with # are ignored."
                    )
                evidence_ids_to_process.append(evidence_id_value)
                seen_data_line = True

    return evidence_ids_to_process


def evidence_has_no_open_revisions(open_revisions_non_contributor, allow_evidence_without_revisions):
    """Test if evidence item has no open revisions from other contributors"""
    if open_revisions_non_contributor == 0 and not allow_evidence_without_revisions:
        print(f"No open revision(s) by other contributors for this evidence item - skipping")
        return True

    return False


def display_evidence_summary(eid, evidence_data):
    """Create a human readable summary of the current live state of an evidence item """

    therapies = ', '.join(evidence_data['therapy_names']) if evidence_data['therapy_names'] else None
    phenotypes = ', '.join(evidence_data['phenotype_names']) if evidence_data['phenotype_names'] else None

    print(
        f"\nEvidence item details for evidence id: {eid}\n"
        f"  Name: {evidence_data['evidence_name']}"
        f" | Status: {evidence_data['status']}"
        f" | Open flags: {evidence_data['open_flag_count']}\n"
        f"  Molecular Profile: {evidence_data['molecular_profile_name']} (id: {evidence_data['molecular_profile_id']})"
        f" | Disease: {evidence_data['disease_name']}\n"
        f"  Evidence Type: {evidence_data['evidence_type']}"
        f" | Direction: {evidence_data['evidence_direction']}"
        f" | Level: {evidence_data['evidence_level']}"
        f" | Significance: {evidence_data['significance']}"
        f" | Rating: {evidence_data['evidence_rating']}\n"
        f"  Variant Origin: {evidence_data['variant_origin']}"
        f" | Therapy Interaction Type: {evidence_data['therapy_interaction_type']}"
        f" | Therapies: {therapies}\n"
        f"  Phenotypes: {phenotypes}\n"
        f"  Source: {evidence_data['source_citation']} ({evidence_data['source_type']}) {evidence_data['source_url']}\n"
        f"  Description: {evidence_data['description']}\n"
        f"  Open revision count: {evidence_data['open_revision_count']}"
    )


def display_evidence_revisions(revisions):
    """Display a human readable summary of open revisions for an evidence item"""

    print(f"\nOpen revisions to review:")
    if len(revisions) == 0:
        print(f"  No open revisions to be compared")
        return

    for revision in revisions:
        print(f"  [{revision['revision_id']}] {revision['field_name']} (by {revision['user_display_name']}): {revision['revision_values_string']}")


def main(evidence_id: int, contributor_id: int, all_evidence: bool, evidence_list_file: str, processed_evidence_file: str, allow_evidence_without_revisions: bool, open_browser: bool):

    #define input data files
    version_file = base_dir / f"RELEASE"

    #load the current release number for this software
    editor_tools_version = None
    with open(version_file) as f:
        editor_tools_version = f.read().strip()

    #make sure internet and API access is working before attempting anything
    verify_connectivity()

    #load a list of already processed evidence items supplied by the user
    already_processed_evidence = set()
    processed_evidence_out_f = None
    if processed_evidence_file:
        already_processed_evidence = load_processed_evidence(processed_evidence_file, already_processed_evidence)
        processed_evidence_out_f = open(processed_evidence_file, "a")

    #summarize user info based on contributor id
    user_details = civic_graphql_utils.gather_user_details(contributor_id)
    print(f"\nContributor (id: {contributor_id}) is {user_details['user_name']} aka {user_details['user_display_name']} ({user_details['user_role']})")

    prompt_to_proceed("Verify your user info above. Revisions by this user will be ignored/skipped. \nYou can't moderate your own submissions.")

    #get civic evidence item IDs to evaluate, either from the user, or by querying CIViCpy
    evidence_ids_to_process = get_evidence_ids_to_process(evidence_id, all_evidence, evidence_list_file)

    #iterate over each evidence item and examine revisions associated with it
    for eid in evidence_ids_to_process:

        #skip if this evidence item is stored as already processed
        if eid in already_processed_evidence:
            continue

        print(f"\n{BLUE}" + "=" * 80)
        print(f"Reviewing CIViC evidence item ID {eid} for revisions that could be reviewed by contributor ID: {contributor_id}{RESET}")

        #query the graphql api for the current live state of the evidence item
        evidence_data = civic_graphql_utils.gather_evidence_summary_data(eid)

        #provide a summary of the evidence item's current live state
        display_evidence_summary(eid, evidence_data)

        #query the graphql api for open revisions on this evidence item
        revision_data = civic_graphql_utils.gather_evidence_revisions(eid, contributor_id)
        open_revisions_non_contributor = evidence_data['open_revision_count'] - revision_data['contributor_revisions']

        #skip an evidence item if it has 0 pending revisions from other users - unless the user wishes to bypass this
        if evidence_has_no_open_revisions(open_revisions_non_contributor, allow_evidence_without_revisions):
            mark_evidence_processed(processed_evidence_file, eid, processed_evidence_out_f, already_processed_evidence)
            continue

        #display the open revisions for this evidence item
        display_evidence_revisions(revision_data['evidence_revisions'])

        #Open the evidence item revision view for the user
        if open_browser:
            open_evidence_revision(eid)

        #Display an example comment message in case the user is going to accept/reject something in CIViC
        print(f"\nTemplate comment for CIViC submission:\nEvidence item was reviewed with civic-editor-tools (v{editor_tools_version})")

        #If the user asked to keep track of already processed evidence items, add this one to the list
        mark_evidence_processed(processed_evidence_file, eid, processed_evidence_out_f, already_processed_evidence)

        #Pause before moving on to the next CIViC evidence item
        prompt_to_proceed(f"Processing complete for evidence item: {eid} ({evidence_data['evidence_name']})")


if __name__ == "__main__":
    args = parse_args()

    main(
        contributor_id=args.contributor_id,
        evidence_id=args.evidence_id,
        all_evidence=args.all_evidence,
        evidence_list_file=args.evidence_list_file,
        processed_evidence_file=args.processed_evidence_file,
        allow_evidence_without_revisions=args.allow_evidence_without_revisions,
        open_browser=args.open_browser
    )
