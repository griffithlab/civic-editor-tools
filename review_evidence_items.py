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
        "--civic-config",
        dest="civic_config",
        type=str,
        required=True,
        help="Path to a *.civic.conf file defining contributor_id, organization_id, api_key, and api_url (see example.civic.conf.template)"
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
    parser.add_argument(
        "--execute",
        dest="execute",
        action="store_true",
        help="Actually submit accept-revision mutations to CIViC. Without this flag, accepting revisions is dry-run only: the exact mutation payload is printed but nothing is sent"
    )
    parser.add_argument(
        "--llm-assist",
        dest="llm_assist",
        action="store_true",
        help="Show an LLM-generated advisory review of an evidence item's open revisions before asking which to accept. Purely advisory - never submits anything itself. Requires the llm/ package's dependencies (pip3 install -r install/requirements.txt) and a configured provider (see llm/profiles.toml)"
    )
    parser.add_argument(
        "--llm-profile",
        dest="llm_profile",
        type=str,
        default="claude",
        help="Provider profile (from llm/profiles.toml) to use with --llm-assist (default: claude)"
    )

    return parser.parse_args()

def verify_connectivity(civic_url: str):
    """Tests whether internet access is working and then if we can actually access the APIs needed before attempting anything """
    if not generic_utils.check_connection():
        print("No internet access. Aborting.")
        sys.exit(1)

    if not generic_utils.check_apis(api_urls=[civic_url]):
        print("Required APIs are unavailable. Aborting.")
        sys.exit(1)

    print("Internet and API connectivity verified.")


def open_evidence_revision(eid: int, civic_url: str) -> None:
    url = f"{civic_url}/evidence/{eid}/revisions"
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


def select_revision_ids_to_accept(revisions: list) -> list:
    """Ask which of the displayed open revisions (if any) to accept. Validates every typed
    ID against what was actually displayed, so a typo can't silently target the wrong revision."""
    if len(revisions) == 0:
        return []

    valid_ids = {revision['revision_id'] for revision in revisions}
    response = input("\nEnter revision ID(s) to accept (comma-separated), 'a'/'all', or Enter to skip: ").strip()

    if not response:
        return []

    if response.lower() in {"a", "all"}:
        return sorted(valid_ids)

    selected_revision_ids = []
    for token in response.split(","):
        token = token.strip()
        if not token:
            continue
        try:
            revision_id = int(token)
        except ValueError:
            sys.exit(f"Fatal: '{token}' is not a valid revision ID")
        if revision_id not in valid_ids:
            sys.exit(f"Fatal: revision ID {revision_id} is not one of the open revisions shown above: {sorted(valid_ids)}")
        selected_revision_ids.append(revision_id)

    return selected_revision_ids


def build_accept_comment(editor_tools_version: str) -> str:
    """Auto-generate the comment submitted alongside an AcceptRevisions mutation. Applies to
    every revision accepted in that call."""
    return f"Reviewed and accepted with assistance from civic-editor-tools (v{editor_tools_version})."


def build_llm_review_inputs(eid: int, evidence_data: dict, revisions: list) -> dict:
    """Build the `inputs` dict for the llm package's review_evidence_item task from data this
    script already has in hand (no extra API calls). source_full_text is always None today -
    this script has no mechanism to fetch the cited publication's full text, or even its
    abstract, from PubMed/PMC."""
    return {
        "evidence_item": {
            "id": eid,
            "statement": evidence_data['description'],
            "evidence_level": evidence_data['evidence_level'],
            #evidence_type/direction/significance come back from the API as raw GraphQL enums
            #(e.g. "NA", "SENSITIVITYRESPONSE") - humanized here so they match both the revision
            #current/proposed values below (already human-readable) and the terminology used in
            #knowledge/evidence_item_basics.md.
            "evidence_type": civic_graphql_utils.humanize_evidence_type(evidence_data['evidence_type']),
            "direction": civic_graphql_utils.humanize_evidence_direction(evidence_data['evidence_direction']),
            "significance": civic_graphql_utils.humanize_evidence_significance(evidence_data['significance']),
            "disease": evidence_data['disease_name'],
            "therapies": evidence_data['therapy_names'],
            "molecular_profile": evidence_data['molecular_profile_name'],
            "source": {
                "citation": evidence_data['source_citation'],
                "url": evidence_data['source_url'],
                "type": evidence_data['source_type'],
            },
            "rating": evidence_data['evidence_rating'],
        },
        "revisions": [
            {
                "revision_id": revision['revision_id'],
                "field_name": revision['field_name'],
                "current_value": revision['current_value'],
                "proposed_value": revision['proposed_value'],
                "comment": revision['comment'],
                "revisor_id": revision['user_id'],
            }
            for revision in revisions
        ],
        "source_full_text": None,
    }


def run_llm_review(inputs: dict, llm_profile: str):
    """
    Best-effort: returns an llm.providers.base.LLMResult, or None if the LLM-assisted review
    fails for any reason (llm package not installed, no configured provider/API key, network
    error, the model's output failing schema validation, a refusal). This step is purely
    advisory - it must never block the human review workflow that follows it, and it never
    itself calls accept_revisions() or any other mutation.
    """
    try:
        from llm import run_task
    except ImportError:
        print(
            f"{YELLOW}WARNING: --llm-assist requires the llm/ package's dependencies "
            f"(pip3 install -r install/requirements.txt). Continuing without it.{RESET}"
        )
        return None

    try:
        return run_task("review_evidence_item", inputs=inputs, profile=llm_profile)
    except Exception as exc:
        print(f"{YELLOW}WARNING: LLM-assisted review failed ({exc}); continuing without it.{RESET}")
        return None


def display_llm_review(result) -> None:
    """Print an llm.providers.base.LLMResult from the review_evidence_item task.
    result.parsed is an EvidenceItemReview - see
    llm/civic-editorial-review/tasks/review_evidence_item/schema.py."""
    review = result.parsed

    print(f"\n{BLUE}--- LLM-assisted review (advisory only - {result.provider}/{result.model}) ---{RESET}")
    print(f"Summary: {review.summary}")

    for assessment in review.revision_assessments:
        recommendation = assessment.recommendation.value
        color = GREEN if recommendation == "accept" else RED if recommendation == "reject" else YELLOW
        print(f"\n  [{assessment.revision_id}] {assessment.field_name}: {color}{recommendation}{RESET} (confidence: {assessment.confidence.value})")
        print(f"    {assessment.rationale}")

    if review.issues:
        print(f"\n{YELLOW}Cross-revision issues:{RESET}")
        for issue in review.issues:
            print(f"  [{issue.severity.value}] {issue.field}: {issue.description}")

    if review.suggested_revisions:
        print(f"\nSuggested alternative values:")
        for suggestion in review.suggested_revisions:
            print(f"  {suggestion.field}: '{suggestion.current_value}' -> '{suggestion.suggested_value}' ({suggestion.rationale})")

    print(f"{BLUE}--- end LLM-assisted review (nothing above was submitted to CIViC) ---{RESET}")


def main(civic_config: str, evidence_id: int, all_evidence: bool, evidence_list_file: str, processed_evidence_file: str, allow_evidence_without_revisions: bool, open_browser: bool, execute: bool, llm_assist: bool, llm_profile: str):

    #define input data files
    version_file = base_dir / f"RELEASE"

    #load the current release number for this software
    editor_tools_version = None
    with open(version_file) as f:
        editor_tools_version = f.read().strip()

    #load and validate the user's CIViC connection config (contributor/org id, api key, api url).
    #if it points at production, this forces a typed confirmation before any API call is made.
    config = civic_graphql_utils.load_civic_config(civic_config)
    civic_graphql_utils.confirm_production_usage(config['api_url'])
    civic_graphql_utils.apply_civic_config(config)
    civic_graphql_utils.summarize_civic_config(config)

    contributor_id = config['contributor_id']
    civic_url = config['api_url']

    #make sure internet and API access is working before attempting anything
    verify_connectivity(civic_url)

    #load a list of already processed evidence items supplied by the user
    already_processed_evidence = set()
    processed_evidence_out_f = None
    if processed_evidence_file:
        already_processed_evidence = load_processed_evidence(processed_evidence_file, already_processed_evidence)
        processed_evidence_out_f = open(processed_evidence_file, "a")

    #summarize user info based on contributor id
    user_details = civic_graphql_utils.gather_user_details(contributor_id)
    print(f"\nContributor (id: {contributor_id}) is {user_details['user_name']} aka {user_details['user_display_name']} ({user_details['user_role']})")

    #summarize organization info based on organization id
    organization_id = config['organization_id']
    organization_details = civic_graphql_utils.gather_organization_details(organization_id)
    print(f"Organization (id: {organization_id}) is {organization_details['organization_name']} ({organization_details['organization_member_count']} members)")

    if llm_assist:
        print(f"LLM-assisted review is ON (profile: {llm_profile}) - advisory only, never submits anything itself.")

    prompt_to_proceed("Verify your user and organization info above. Revisions by this user will be ignored/skipped. \nYou can't moderate your own submissions.")

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

        #Optionally show an LLM-generated advisory review of the open revisions, before asking
        #which (if any) to accept. Purely advisory: never blocks the human workflow below, and
        #never itself calls accept_revisions() or any other mutation.
        if llm_assist:
            llm_inputs = build_llm_review_inputs(eid, evidence_data, revision_data['evidence_revisions'])
            llm_result = run_llm_review(llm_inputs, llm_profile)
            if llm_result:
                display_llm_review(llm_result)

        #Open the evidence item revision view for the user
        if open_browser:
            open_evidence_revision(eid, civic_url)

        #ask which (if any) of the displayed revisions to accept, then accept them (dry-run
        #unless --execute was passed; accept_revisions() itself requires a typed confirmation
        #before actually sending, regardless of this flag)
        selected_revision_ids = select_revision_ids_to_accept(revision_data['evidence_revisions'])
        if selected_revision_ids:
            comment = build_accept_comment(editor_tools_version)
            civic_graphql_utils.accept_revisions(selected_revision_ids, organization_id, comment, dry_run=not execute)
        else:
            print(f"\nNo revisions selected to accept for this evidence item.")

        #If the user asked to keep track of already processed evidence items, add this one to the list
        mark_evidence_processed(processed_evidence_file, eid, processed_evidence_out_f, already_processed_evidence)

        #Pause before moving on to the next CIViC evidence item
        prompt_to_proceed(f"Processing complete for evidence item: {eid} ({evidence_data['evidence_name']})")


if __name__ == "__main__":
    args = parse_args()

    main(
        civic_config=args.civic_config,
        evidence_id=args.evidence_id,
        all_evidence=args.all_evidence,
        evidence_list_file=args.evidence_list_file,
        processed_evidence_file=args.processed_evidence_file,
        allow_evidence_without_revisions=args.allow_evidence_without_revisions,
        open_browser=args.open_browser,
        execute=args.execute,
        llm_assist=args.llm_assist,
        llm_profile=args.llm_profile
    )
