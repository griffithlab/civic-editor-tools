#!/usr/bin/env python3

from pathlib import Path
import requests
import pdb
import time
import sys
import stat
import json
import html
from requests.exceptions import ConnectionError, Timeout, RequestException, HTTPError

api_url = "https://civicdb.org/api/graphql"
api_key = None  # optional bearer token, set at runtime by apply_civic_config()
base_dir = Path(__file__).resolve().parent

RED = "\033[31m"
YELLOW = "\033[33m"
RESET = "\033[0m"

#HTTP statuses worth retrying: rate limiting and transient server errors
RETRYABLE_HTTP_STATUS_CODES = {429, 500, 502, 503, 504}

#default pacing applied before every GraphQL request, to stay under CIViC's rate limit
#regardless of caller or query pattern (single variant, --all-variants, --target-gene, etc.)
DEFAULT_GRAPHQL_REQUEST_DELAY = 0.25

#the only two CIViC deployments a *.civic.conf file is allowed to point at
PRODUCTION_CIVIC_URL = "https://civicdb.org"
STAGING_CIVIC_URL = "https://staging.civicdb.org"
VALID_CIVIC_URLS = {PRODUCTION_CIVIC_URL, STAGING_CIVIC_URL}

CIVIC_CONFIG_SUFFIX = ".civic.conf"
CIVIC_CONFIG_REQUIRED_KEYS = {"contributor_id", "organization_id", "api_key", "api_url"}

#the exact word a user must type to confirm they intend to query the live production database
PRODUCTION_CONFIRMATION_PHRASE = "production"

#the words (case-insensitive) a user may type to confirm a mutation should actually be sent (not dry-run)
ACCEPT_REVISIONS_CONFIRMATION_WORDS = {"y", "yes"}

def populate_variables_id(variables_template: str, graphql_id: int) -> str:
    """update a template graphql object string to inject the query ID to be used"""
    placeholder = "graphql_query_id1"

    count = variables_template.count(placeholder)
    if count != 1:
        raise ValueError(
            f"Expected exactly one '{placeholder}', found {count}"
    )

    if placeholder not in variables_template:
        raise ValueError(
            f"Placeholder '{placeholder}' not found in variables template"
        )

    populated = variables_template.replace(placeholder, str(graphql_id), 1)

    return populated


def populate_variables(variables_template: str, values: dict) -> str:
    """
    Update a template graphql variables string by substituting one or more named
    placeholders with their values. Each key in `values` must appear in the template
    exactly once. json.dumps() (not str()) is used for substitution so list/str/int
    values are all safely and correctly JSON-encoded - critical for a free-text field
    (e.g. a revision comment) that could otherwise contain characters (quotes, newlines)
    that would break the surrounding JSON structure if spliced in raw.
    """
    populated = variables_template

    for placeholder, value in values.items():
        count = populated.count(placeholder)
        if count != 1:
            raise ValueError(f"Expected exactly one '{placeholder}', found {count}")
        populated = populated.replace(placeholder, json.dumps(value), 1)

    return populated


def run_graphql_operation(api_url: str, operation_name: str, query_id: int, timeout: tuple = (20, 200),
                          retries: int = 4, backoff_factor: float = 2.0,
                          request_delay: float = DEFAULT_GRAPHQL_REQUEST_DELAY) -> requests.Response:
    """Load graphql query and variable json objects from file, update with a
    query id, and submit the query to the API. Paces every request by request_delay
    seconds (including retries) to reduce how often CIViC's rate limit is hit in the
    first place, and retries on transient errors."""
    query_path = base_dir / f"../graphql/{operation_name}_query.json"
    variables_path = base_dir / f"../graphql/{operation_name}_variables.json"

    if not query_path.exists():
        raise FileNotFoundError(f"Missing query file: {query_path}")
    if not variables_path.exists():
        raise FileNotFoundError(f"Missing variables file: {variables_path}")

    with query_path.open("r") as f:
        query = f.read()
    with variables_path.open("r") as f:
        variables = f.read()

    variables_updated = populate_variables_id(variables, query_id)

    last_exc: Exception | None = None

    #if a *.civic.conf file was loaded and applied (see apply_civic_config()), authenticate
    #requests with its api_key; otherwise requests are sent unauthenticated, as before
    headers = {"Content-Type": "application/json"}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"

    for attempt in range(1, retries + 1):
        time.sleep(request_delay)
        try:
            resp = requests.post(
                api_url,
                json={"query": query, "variables": variables_updated},
                timeout=timeout,
                headers=headers,
            )
            resp.raise_for_status()  # surface 4xx/5xx as exceptions
            return resp

        except (ConnectionError, Timeout) as exc:
            last_exc = exc
            if attempt < retries:
                wait = backoff_factor ** (attempt - 1)  # 1s, 2s, 4s …
                print(f"[Attempt {attempt}/{retries}] Network error: {exc}. Retrying in {wait:.1f}s…")
                time.sleep(wait)

        except HTTPError as exc:
            last_exc = exc
            status_code = exc.response.status_code if exc.response is not None else None

            # Non-retryable HTTP errors (auth failures, bad requests, not found, etc.)
            if status_code not in RETRYABLE_HTTP_STATUS_CODES:
                raise RuntimeError(f"GraphQL request failed: {exc}") from exc

            if attempt < retries:
                wait = backoff_factor ** (attempt - 1)  # 1s, 2s, 4s …
                retry_after = exc.response.headers.get("Retry-After") if exc.response is not None else None
                if retry_after:
                    try:
                        wait = max(wait, float(retry_after))
                    except ValueError:
                        pass
                print(f"[Attempt {attempt}/{retries}] HTTP {status_code} from CIViC API: {exc}. Retrying in {wait:.1f}s…")
                time.sleep(wait)

        except RequestException as exc:
            # Non-retryable errors (e.g. malformed request that requests itself rejects)
            raise RuntimeError(f"GraphQL request failed: {exc}") from exc

    raise RuntimeError(
        f"GraphQL operation '{operation_name}' failed after {retries} attempts. "
        f"Last error: {last_exc}"
    ) from last_exc


def run_graphql_mutation(api_url: str, operation_name: str, values: dict, timeout: tuple = (20, 200)) -> requests.Response:
    """
    Load a graphql mutation document and its variables template from file, substitute the
    given named values (see populate_variables()), and submit the mutation exactly once.

    Deliberately does NOT retry, unlike run_graphql_operation(). A failed or ambiguous
    (e.g. timed-out) mutation attempt must surface to the caller rather than being silently
    retried: unlike a read, the original attempt may have already succeeded server-side, so
    retrying risks double-submitting (e.g. accepting the same revisions twice).
    """
    query_path = base_dir / f"../graphql/{operation_name}_mutate.json"
    variables_path = base_dir / f"../graphql/{operation_name}_variables.json"

    if not query_path.exists():
        raise FileNotFoundError(f"Missing mutation file: {query_path}")
    if not variables_path.exists():
        raise FileNotFoundError(f"Missing variables file: {variables_path}")

    with query_path.open("r") as f:
        query = f.read()
    with variables_path.open("r") as f:
        variables = f.read()

    variables_updated = populate_variables(variables, values)

    headers = {"Content-Type": "application/json"}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"

    time.sleep(DEFAULT_GRAPHQL_REQUEST_DELAY)
    resp = requests.post(
        api_url,
        json={"query": query, "variables": variables_updated},
        timeout=timeout,
        headers=headers,
    )
    resp.raise_for_status()  # surface 4xx/5xx as exceptions - no retry, caller must handle

    return resp


def accept_revisions(revision_ids: list, organization_id: int, comment: str, dry_run: bool = True) -> dict:
    """
    Accept one or more open CIViC revisions via the AcceptRevisions mutation.

    dry_run (default True) prints the exact payload that would be sent and returns without
    sending it. When dry_run is False, prints the same payload and then requires an exact
    typed confirmation before the mutation is actually sent - this function is the single
    choke point every caller goes through, so this safety gate can't be skipped by a future
    caller forgetting to add it.
    """
    values = {
        "graphql_var_ids": revision_ids,
        "graphql_var_organization_id": organization_id,
        "graphql_var_comment": comment,
    }

    print(
        f"\nAccept revisions mutation payload:\n"
        f"  Revision IDs: {revision_ids}\n"
        f"  Organization ID: {organization_id}\n"
        f"  Comment: {comment}"
    )

    if dry_run:
        print(f"{YELLOW}Dry run - mutation NOT sent. Re-run with --execute to actually accept these revisions.{RESET}")
        return {"dry_run": True, "revision_ids": revision_ids}

    try:
        response = input(
            f"\nType 'y' or 'yes' to submit this mutation to {api_url}, "
            f"or anything else to abort: "
        )
    except KeyboardInterrupt:
        sys.exit("\nAborted.")

    if response.strip().lower() not in ACCEPT_REVISIONS_CONFIRMATION_WORDS:
        sys.exit("Aborted: accept confirmation did not match.")

    resp = run_graphql_mutation(api_url, "revision_AcceptRevisions", values)
    result = resp.json()

    if "errors" in result:
        sys.exit(f"Fatal: AcceptRevisions mutation returned errors: {result['errors']}")

    accepted_revision_ids = [r["id"] for r in result["data"]["acceptRevisions"]["revisions"]]

    return {"dry_run": False, "accepted_revision_ids": accepted_revision_ids}


def load_civic_config(config_path: str) -> dict:
    """
    Load, validate, and return CIViC connection settings from a user-supplied config file.

    The file must:
      - be named '*.civic.conf' (matched by this repo's .gitignore, so it can never be committed)
      - be readable/writable by its owner only, no group or other access (like an AWS .pem file)
      - contain simple 'key = value' lines (blank lines and lines starting with # are ignored)
        defining exactly: contributor_id, organization_id, api_key, api_url

    api_url must be exactly 'https://civicdb.org' or 'https://staging.civicdb.org'. Does not
    print or log the api_key value. Exits with a helpful message on any validation failure
    rather than raising, since this is always called directly from a CLI entry point.
    """
    path = Path(config_path)

    if not path.name.endswith(CIVIC_CONFIG_SUFFIX):
        sys.exit(
            f"Fatal: config file '{path}' must be named '*{CIVIC_CONFIG_SUFFIX}' "
            f"(this naming convention is what keeps it out of git via .gitignore)."
        )

    if not path.exists():
        sys.exit(f"Fatal: config file not found: {path}")

    #enforce AWS .pem-style permissions: owner access only, no group/other bits set.
    #this file contains an API key, so a loose mode is treated as fatal, not just a warning.
    mode = stat.S_IMODE(path.stat().st_mode)
    if mode & 0o077:
        sys.exit(
            f"Fatal: config file '{path}' is readable by group and/or other (mode {oct(mode)}).\n"
            f"  This file contains an API key and must be readable by its owner only. Fix it with:\n"
            f"    chmod 600 {path}"
        )

    values = {}
    with path.open("r") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                sys.exit(f"Fatal: {path} line {line_num} is not a 'key = value' line: {line!r}")
            key, _, value = line.partition("=")
            values[key.strip()] = value.strip()

    missing = CIVIC_CONFIG_REQUIRED_KEYS - values.keys()
    if missing:
        sys.exit(f"Fatal: {path} is missing required setting(s): {', '.join(sorted(missing))}")

    try:
        contributor_id = int(values["contributor_id"])
    except ValueError:
        sys.exit(f"Fatal: {path} 'contributor_id' must be an integer, got: {values['contributor_id']!r}")

    try:
        organization_id = int(values["organization_id"])
    except ValueError:
        sys.exit(f"Fatal: {path} 'organization_id' must be an integer, got: {values['organization_id']!r}")

    civic_api_key = values["api_key"]
    if not civic_api_key:
        sys.exit(f"Fatal: {path} 'api_key' must not be empty")

    civic_url = values["api_url"].rstrip("/")
    if civic_url not in VALID_CIVIC_URLS:
        valid_urls_string = " or ".join(sorted(VALID_CIVIC_URLS))
        sys.exit(f"Fatal: {path} 'api_url' must be exactly {valid_urls_string}, got: {values['api_url']!r}")

    return {
        "contributor_id": contributor_id,
        "organization_id": organization_id,
        "api_key": civic_api_key,
        "api_url": civic_url,
    }


def confirm_production_usage(civic_url: str) -> None:
    """If civic_url is the live production CIViC database, force an explicit typed confirmation
    before the caller is allowed to proceed to make any API calls."""
    if civic_url != PRODUCTION_CIVIC_URL:
        return

    print(
        f"\n{RED}{'!' * 80}\n"
        f"WARNING: This configuration points at the LIVE PRODUCTION CIViC database ({PRODUCTION_CIVIC_URL}).\n"
        f"Any revisions you review and comment on here are visible to real CIViC users.\n"
        f"{'!' * 80}{RESET}"
    )
    try:
        response = input(f"Type exactly: \"{PRODUCTION_CONFIRMATION_PHRASE}\" to proceed, or anything else to abort: ")
    except KeyboardInterrupt:
        sys.exit("\nAborted.")

    if response.strip() != PRODUCTION_CONFIRMATION_PHRASE:
        sys.exit("Aborted: production confirmation phrase did not match.")


def apply_civic_config(config: dict) -> None:
    """Point this module's GraphQL requests at the api_url/api_key from a loaded civic config"""
    global api_url, api_key
    api_url = f"{config['api_url']}/api/graphql"
    api_key = config['api_key']


def summarize_civic_config(config: dict) -> None:
    """Print the non-secret parts of a loaded civic config back to the user. Never prints the api_key value."""
    url_note = f"{YELLOW} (PRODUCTION){RESET}" if config['api_url'] == PRODUCTION_CIVIC_URL else " (staging)"
    print(
        f"\nCIViC configuration in use:\n"
        f"  API URL: {config['api_url']}{url_note}\n"
        f"  Contributor ID: {config['contributor_id']}\n"
        f"  Organization ID: {config['organization_id']}\n"
        f"  API Key: ******** (provided, not displayed)"
    )


def gather_user_details(user_id: int) -> dict:
    """execute graphql queries, parse json, build a simplied data structure with the user info needed"""
    #graphql template name: "user_Data"
    resp = run_graphql_operation(api_url, "user_Data", user_id)
    json = resp.json()
    
    user_name = json['data']['user']['name']
    user_display_name = json['data']['user']['displayName']
    user_role = json['data']['user']['role']
    user_data = {
        "user_name": user_name,
        "user_display_name": user_display_name,
        "user_role": user_role
    }
    #pdb.set_trace()

    return user_data


def gather_organization_details(organization_id: int) -> dict:
    """execute graphql queries, parse json, build a simplied data structure with the organization info needed"""
    #graphql template name: "organization_Data"
    resp = run_graphql_operation(api_url, "organization_Data", organization_id)
    json = resp.json()

    organization = json['data']['organization']
    if organization is None:
        sys.exit(f"Fatal: no CIViC organization found for organization_id: {organization_id}")

    organization_data = {
        "organization_name": organization['name'],
        "organization_url": organization['url'],
        "organization_member_count": organization['memberCount'],
    }

    return organization_data


def gather_variant_details(variant_id: int) -> dict:
    """execute graphql queries, parse json, return basic variant info"""
    #graphql template name: "variant_VariantDetail"
    resp = run_graphql_operation(api_url, "variant_VariantDetail", variant_id)
    json = resp.json()
    variant_name = json['data']['variant']['name']
    feature_name = json['data']['variant']['feature']['name']
    deprecated = json['data']['variant']['deprecated']
    variant_data = {
        "variant_name": variant_name,
        "feature_name": feature_name,
        "deprecated": deprecated
    }
    #pdb.set_trace()

    return variant_data    


def gather_accepted_variant_data(variant_id: int) -> dict:
    """execute graphql queries, parse json, return detailed info on already accepted variant fields"""
    #graphql template name: "variant_VariantSummary"
    resp = run_graphql_operation(api_url, "variant_VariantSummary", variant_id)
    json = resp.json()

    variant_types = []
    for vts in json['data']['variant']['variantTypes']:
        variant_types.append(vts['name'])

    #pdb.set_trace()
    accepted_variant_data = {
        "allele_registry_id": json['data']['variant']['alleleRegistryId'],
        "name": json['data']['variant']['name'],
        "variant_types": variant_types,
        "variant_aliases": json['data']['variant']['variantAliases'],
        "hgvs_descriptions": json['data']['variant']['hgvsDescriptions'],
        "clinvar_ids": json['data']['variant']['clinvarIds'],
        "reference_build": json['data']['variant']['coordinates']['referenceBuild'],
        "chromosome": json['data']['variant']['coordinates']['chromosome'],
        "start": json['data']['variant']['coordinates']['start'],
        "stop": json['data']['variant']['coordinates']['stop'],
        "reference_bases": json['data']['variant']['coordinates']['referenceBases'],
        "variant_bases": json['data']['variant']['coordinates']['variantBases'],
        "representative_transcript": json['data']['variant']['coordinates']['representativeTranscript'],
        "ensembl_version": json['data']['variant']['coordinates']['ensemblVersion'],
        "coordinate_type": json['data']['variant']['coordinates']['coordinateType']
    }

    return accepted_variant_data


def gather_variant_revisions(variant_id: int, contributor_id: int) -> dict:
    """execute graphql queries, parse json, build a simplied data structure with the variant revision info needed"""

    #Get coordinate ids for variant (takes a variant id)
    #graphql template name: "variant_CoordinateIdsForVariant"
    resp = run_graphql_operation(api_url, "variant_CoordinateIdsForVariant", variant_id)
    json = resp.json()

    open_revision_count_variant = json['data']['variant']['openRevisionCount']
    open_revision_count_coordinates = json['data']['variant']['coordinates']['openRevisionCount']
    variant_coordinates_id = json['data']['variant']['coordinates']['id']

    variant_data = {
        "variant_id": variant_id,
        "open_revision_count_variant": open_revision_count_variant,
        "open_revision_count_coordinates": open_revision_count_coordinates,
        "variant_coordinates_id": variant_coordinates_id,
        "contributor_revisions": 0,
        "variant_revisions": [],
        "coordinate_revisions": []
    }	

    #graphql template name: "variant_VariantDetail"
    resp = run_graphql_operation(api_url, "variant_VariantDetail", variant_id)
    json = resp.json()
    variant_name = json['data']['variant']['name']
    feature_name = json['data']['variant']['feature']['name']
	
    variant_data['variant_name'] = variant_name
    variant_data['feature_name'] = feature_name
    variant_data['name_change'] = False

    #variant_Revisions-Variant (takes a variant id)
    #graphql template name: "variant_Revisions-Variant"
    resp = run_graphql_operation(api_url, 'variant_Revisions-Variant', variant_id)
    json = resp.json()

    #pdb.set_trace()

    revisions = json["data"]["revisions"]["edges"]
    i = 0
    for revision in revisions:
        revision_id = revision['node']['id']
        user_id = revision['node']['creationActivity']['user']['id']
        user_display_name = revision['node']['creationActivity']['user']['displayName']
        if user_id == contributor_id:
            variant_data['contributor_revisions'] += 1
        
        field_name = revision['node']['fieldName']
        revision_values_string = ""
        revision_values_list = []

        #special handling when the revision is the variant "name" itself
        if field_name == 'name':
            current_value = revision['node']['currentValue']
            suggested_value = revision['node']['suggestedValue']
            revision_values_string = f"'{current_value}' -> '{suggested_value}'"
            variant_data['name_change'] = True
        else:
            #diffValue is a GraphQL union: an ObjectFieldDiff (list-type fields, e.g.
            #hgvs_description_ids) has addedObjects; a ScalarFieldDiff (e.g. an unsupported
            #field like chromosome2) has left/right instead. Handle both shapes here so an
            #unexpected field can't crash this parsing step; merge_revision_data() is still
            #responsible for rejecting field names this tool has no comparison logic for.
            diff_value = revision['node']['linkoutData']['diffValue']
            if 'addedObjects' in diff_value:
                #revisions that are lists of things
                for revision_value in diff_value['addedObjects']:
                    revision_display_name = revision_value['displayName']
                    revision_values_list.append(revision_display_name)
                revision_values_string = ",".join(sorted(revision_values_list))
            else:
                #a scalar before/after value pair
                revision_values_string = f"'{diff_value.get('left')}' -> '{diff_value.get('right')}'"

        variant_data["variant_revisions"].append({
            "index": i,
            "revision_id": revision_id,
            "user_id": user_id,
            "user_display_name": user_display_name,
            "field_name": field_name,
            "revision_values_list": revision_values_list,
            "revision_values_string": revision_values_string
        })
        i += 1

    #variant_Revisions-VariantCoordinates (takes a variant _coordinates_ id)
    #graphql template name: "variant_Revisions-VariantCoordinates"
    resp = run_graphql_operation(api_url, "variant_Revisions-VariantCoordinates", variant_coordinates_id)
    json = resp.json()
    coordinate_revisions = json["data"]["revisions"]["edges"]
    i = 0
    for revision in coordinate_revisions:
        revision_id = revision['node']['id']
        user_id = revision['node']['creationActivity']['user']['id']
        if user_id == contributor_id:
            variant_data['contributor_revisions'] += 1
        user_display_name = revision['node']['creationActivity']['user']['displayName']
        field_name = revision['node']['fieldName']
        suggested_value = revision['node']['suggestedValue']

        variant_data["coordinate_revisions"].append({
            "index": i,
            "revision_id": revision_id,
            "user_id": user_id,
            "user_display_name": user_display_name,
            "field_name": field_name,
            "suggested_value": suggested_value
        })
        i += 1
    
    variant_data['open_revisions_non_contributor'] = variant_data['open_revision_count_variant'] - variant_data['contributor_revisions']

    #print(variant_data)

    #To interactively explore json responses that come back from these queryies, place this inline above:
    #pdb.set_trace()
    #json = resp.json()
    #json['data']['revisions']['edges'][0]['node']['creationActivity']['user']['displayName']
    #json['data']['revisions']['edges'][0]['node']['linkoutData']['diffValue']['addedObjects'][0]['displayName']
    #json['data']['revisions']['edges'][0]['node']['fieldName']
    return variant_data


def merge_revision_data(variant_data: dict) -> dict:
    
    #Define the priority order for field names
    FIELD_NAME_PRIORITY = {
        "name": 0,
        "variant_type_ids": 1,
        "variant_alias_ids": 2,
        "hgvs_description_ids": 3,
        "clinvar_entry_ids": 4,
        "reference_build": 5,
        "chromosome": 6,
        "start": 7,
        "stop": 8,
        "reference_bases": 9,
        "variant_bases": 10,
        "representative_transcript": 11,
        "ensembl_version": 12
        # ... add all known field names here
    }
    #First tag each entry with its revision type then combine into "all_revisions"
    all_revisions = [
        {**entry, "revision_type": "variant"}
        for entry in variant_data["variant_revisions"]
    ] + [
        {**entry, "revision_type": "coordinate"}
        for entry in variant_data["coordinate_revisions"]
    ]

    # Validate all field names before sorting
    unknown_fields = {
       entry["field_name"]
       for entry in all_revisions
       if entry["field_name"] not in FIELD_NAME_PRIORITY
    }
    if unknown_fields:
        raise ValueError(f"Fatal: unexpected field_name(s) encountered: {unknown_fields}")

    #Sort revisions by the hard coded priority order of the features
    all_revisions.sort(key=lambda entry: FIELD_NAME_PRIORITY[entry["field_name"]])

    variant_data["all_revisions"] = all_revisions

    return variant_data


#Raw GraphQL enum -> human-readable display form, for values that don't reduce to a plain
#.title() of the enum name (irregular CIViC enum spellings, e.g. no separator between words).
#Sourced from https://docs.civicdb.org/en/latest/model/evidence/{direction,significance}.html
_EVIDENCE_DIRECTION_DISPLAY = {
    "SUPPORTS": "Supports",
    "DOES_NOT_SUPPORT": "Does Not Support",
    "NA": "N/A",
}

_EVIDENCE_SIGNIFICANCE_DISPLAY = {
    "NA": "N/A",
    "SENSITIVITYRESPONSE": "Sensitivity/Response",
    "RESISTANCE": "Resistance",
    "REDUCEDSENSITIVITY": "Reduced Sensitivity",
    "ADVERSERESPONSE": "Adverse Response",
    "POSITIVE": "Positive",
    "NEGATIVE": "Negative",
    "BETTEROUTCOME": "Better Outcome",
    "POOROUTCOME": "Poor Outcome",
    "PREDISPOSITION": "Predisposition",
    "PROTECTIVENESS": "Protectiveness",
    "ONCOGENICITY": "Oncogenicity",
    "GAINOFFUNCTION": "Gain of Function",
    "LOSSOFFUNCTION": "Loss of Function",
    "UNALTEREDFUNCTION": "Unaltered Function",
    "NEOMORPHIC": "Neomorphic",
    "DOMINANTNEGATIVE": "Dominant Negative",
    "UNKNOWN": "Unknown",
}


def humanize_evidence_type(raw_value: str) -> str:
    """CIViC's evidence_type enum values (PREDICTIVE, DIAGNOSTIC, ...) are single words, so a
    plain .title() recovers the display form used in CIViC's own docs (Predictive, Diagnostic)."""
    return raw_value.title() if raw_value else raw_value


def humanize_evidence_direction(raw_value: str) -> str:
    if not raw_value:
        return raw_value
    return _EVIDENCE_DIRECTION_DISPLAY.get(raw_value, raw_value)


def humanize_evidence_significance(raw_value: str) -> str:
    if not raw_value:
        return raw_value
    return _EVIDENCE_SIGNIFICANCE_DISPLAY.get(raw_value, raw_value)


def gather_evidence_summary_data(evidence_id: int) -> dict:
    """execute graphql queries, parse json, return current live details for a CIViC evidence item"""
    #graphql template name: "evidence_EvidenceSummary"
    resp = run_graphql_operation(api_url, "evidence_EvidenceSummary", evidence_id)
    json = resp.json()
    evidence = json['data']['evidenceItem']

    therapy_names = [t['name'] for t in evidence['therapies']]
    phenotype_names = [p['name'] for p in evidence['phenotypes']]

    evidence_data = {
        "evidence_id": evidence['id'],
        "evidence_name": evidence['name'],
        "description": evidence['description'],
        "status": evidence['status'],
        "open_revision_count": evidence['revisions']['totalCount'],
        "open_flag_count": evidence['flags']['totalCount'],
        "evidence_type": evidence['evidenceType'],
        "evidence_direction": evidence['evidenceDirection'],
        "evidence_level": evidence['evidenceLevel'],
        "significance": evidence['significance'],
        "evidence_rating": evidence['evidenceRating'],
        "variant_origin": evidence['variantOrigin'],
        "therapy_interaction_type": evidence['therapyInteractionType'],
        "disease_name": evidence['disease']['name'] if evidence['disease'] else None,
        "molecular_profile_id": evidence['molecularProfile']['id'],
        "molecular_profile_name": evidence['molecularProfile']['name'],
        "therapy_names": therapy_names,
        "phenotype_names": phenotype_names,
        "source_citation": evidence['source']['citation'],
        "source_url": evidence['source']['sourceUrl'],
        "source_type": evidence['source']['sourceType'],
    }

    return evidence_data


#segment __typename values that carry human-readable text, and which key holds it
_COMMENT_TEXT_KEYS = {
    "CommentTextSegment": "text",
    "CommentTagSegment": "displayName",
    "CommentTagSegmentFlagged": "displayName",
    "CommentTagSegmentFlaggedAndWithStatus": "displayName",
    "CommentTagSegmentFlaggedAndDeprecated": "displayName",
    "User": "displayName",
}


def extract_comment_text(parsed_note: list) -> str:
    """Concatenate a revision's creationActivity.parsedNote segments into a plain-text comment.
    In practice almost every real comment is a single CommentTextSegment, but a comment can also
    tag other CIViC entities or mention a user inline; those segments are rendered by their
    displayName so the concatenated text stays readable rather than being silently dropped."""
    parts = []
    for segment in parsed_note:
        key = _COMMENT_TEXT_KEYS.get(segment.get("__typename"))
        if key:
            parts.append(segment.get(key) or "")
    return html.unescape("".join(parts)).strip()


def gather_evidence_revisions(evidence_id: int, contributor_id: int) -> dict:
    """execute graphql queries, parse json, build a simplified data structure with open revision info for an evidence item"""
    #graphql template name: "evidence_Revisions-Evidence"
    resp = run_graphql_operation(api_url, "evidence_Revisions-Evidence", evidence_id)
    json = resp.json()

    revisions = json["data"]["revisions"]["edges"]
    evidence_revisions = []
    contributor_revisions = 0
    for i, revision in enumerate(revisions):
        revision_id = revision['node']['id']
        user_id = revision['node']['creationActivity']['user']['id']
        user_display_name = revision['node']['creationActivity']['user']['displayName']
        if user_id == contributor_id:
            contributor_revisions += 1

        field_name = revision['node']['fieldName']
        revision_values_string = ""
        revision_values_list = []

        #diffValue is a GraphQL union: an ObjectFieldDiff (list/reference-type fields, e.g.
        #molecular_profile_id, therapy_ids) has addedObjects (and currentObjects), whose
        #displayName values are the only human-readable form available (currentValue/
        #suggestedValue are raw id lists for these fields). A ScalarFieldDiff (e.g. description,
        #significance) has left/right, but those are pre-rendered HTML diff markup meant for the
        #CIViC web UI, not plain text - so for scalar fields prefer the node's own
        #currentValue/suggestedValue instead.
        diff_value = revision['node']['linkoutData']['diffValue']
        if 'addedObjects' in diff_value:
            for revision_value in diff_value['addedObjects']:
                revision_values_list.append(revision_value['displayName'])
            revision_values_string = ",".join(sorted(revision_values_list))
            current_object_names = sorted(obj['displayName'] for obj in (diff_value.get('currentObjects') or []))
            current_value = ",".join(current_object_names)
            proposed_value = revision_values_string
        else:
            current_value = revision['node']['currentValue']
            proposed_value = revision['node']['suggestedValue']
            revision_values_string = f"'{current_value}' -> '{proposed_value}'"

        comment = extract_comment_text(revision['node']['creationActivity']['parsedNote'])

        evidence_revisions.append({
            "index": i,
            "revision_id": revision_id,
            "user_id": user_id,
            "user_display_name": user_display_name,
            "field_name": field_name,
            "current_value": current_value,
            "proposed_value": proposed_value,
            "comment": comment,
            "revision_values_list": revision_values_list,
            "revision_values_string": revision_values_string
        })

    evidence_data = {
        "evidence_id": evidence_id,
        "evidence_revisions": evidence_revisions,
        "contributor_revisions": contributor_revisions,
    }

    return evidence_data


def load_blacklisted_variant_ids(filepath: str) -> list:
	"""Load blacklisted variant IDs from a file. One per line. Each line must start with the ID, anything else on the line will be ignored"""
	variant_ids = set()

	with open(filepath, "r") as fh:
		next(fh) #skip header

		for line in fh:
			line = line.strip()
			if not line:
				continue

			# Take the first whitespace-delimited field only
			variant_id_str = line.split()[0]

			try:
				variant_ids.add(int(variant_id_str))
			except ValueError:
				continue

	return variant_ids


def main (variant_id: int, contributor_id: int) -> None:
    """demonstrate functionality of the methods above and variant data retrieved"""
    
    #get user/contributor information from the contributor id
    user_details = gather_user_details(contributor_id)
    print(
        f"\nContributor/user info for contributor id: {contributor_id}\n"
        f"  User name: {user_details['user_name']}\n"
        f"  User display name: {user_details['user_display_name']}\n"
        f"  User role: {user_details['user_role']}"
    )

    #get basic variant ino
    variant_data = gather_variant_details(variant_id)
    print(
        f"\nBasic variant details for variant id: {variant_id}\n"
        f"  Variant name: {variant_data['variant_name']}\n"
        f"  Feature name: {variant_data['feature_name']}\n"
        f"  Deprecated: {variant_data['deprecated']}"
    )

    #get much more detailed info on already accepted variant fields
    accepted_variant_data = gather_accepted_variant_data(variant_id)

    print(
        f"\nAccepted variant details for variant id: {variant_id}\n"
        f"  Allele Registry ID: {accepted_variant_data['allele_registry_id']}\n"
        f"  Variant types: {accepted_variant_data['variant_types']}\n"
        f"  Name: {accepted_variant_data['name']}\n"
        f"  Variant Aliases: {accepted_variant_data['variant_aliases']}\n"
        f"  HGVS Descriptions: {accepted_variant_data['hgvs_descriptions']}\n"
        f"  ClinVar IDs: {accepted_variant_data['clinvar_ids']}\n"
        f"  Reference Build: {accepted_variant_data['reference_build']}\n"
        f"  Chromosome: {accepted_variant_data['chromosome']}\n"
        f"  Start: {accepted_variant_data['start']}\n"
        f"  Stop: {accepted_variant_data['stop']}\n"
        f"  Reference Bases {accepted_variant_data['reference_bases']}\n"
        f"  Variant Bases: {accepted_variant_data['variant_bases']}\n"
        f"  Representative Transcript {accepted_variant_data['representative_transcript']}\n"
        f"  Ensembl Version: {accepted_variant_data['ensembl_version']}\n"
        f"  Coordinate Type: {accepted_variant_data['coordinate_type']}\n"
    )

    #get variant revision summary information
    variant_data = gather_variant_revisions(variant_id, contributor_id)
    print(
        f"\nVariant revision info from gather_variant_revisions()\n"
        f"Variant ID used for graphql query: {variant_data['variant_id']}\n"
        f"  Variant name: {variant_data['variant_name']}\n"
        f"  Feature name: {variant_data['feature_name']}\n"
        f"  Open gene-variant revisions (total): {variant_data['open_revision_count_variant']}\n"
        f"  Open gene-variant revisions from specified contributor: {variant_data['contributor_revisions']}\n"
        f"  Open gene-variant revisions from all others users: {variant_data['open_revisions_non_contributor']}\n"
        f"  Variant coordinates id: {variant_data['variant_coordinates_id']}"
    )
    #iterate through individual variant revisions
    variant_revisions = variant_data['variant_revisions']
    for variant_revision in variant_revisions:
        print(
            f"\nInformation for variant revision: {variant_revision['revision_id']}\n"
            f"  Revision user display name: {variant_revision['user_display_name']} (id: {variant_revision['user_id']})\n"
            f"  Revision field name: {variant_revision['field_name']}\n"
            f"  Revision values(s): {variant_revision['revision_values_list']}\n"
            f"  Revision value(s) string: {variant_revision['revision_values_string']}"
        )

    #iterate through coordinate revisions
    coordinate_revisions = variant_data['coordinate_revisions']
    for coordinate_revision in coordinate_revisions:
        print(
            f"\nInformation for coordinate revision: {coordinate_revision['revision_id']}\n"
            f"  Revision user display name: {coordinate_revision['user_display_name']} (id: {coordinate_revision['user_id']})\n"
            f"  Revision field name: {coordinate_revision['field_name']}\n"
            f"  Revision value(s): {coordinate_revision['suggested_value']}"
        )

    #create a unified revisions object that combines all the revisions together and order them logically
    variant_data = merge_revision_data(variant_data)

    all_revisions = variant_data['all_revisions']

    for revision in all_revisions:
        revision_value = ""
        if revision["revision_type"] == "variant":
            revision_value = revision['revision_values_string']
        elif revision["revision_type"] == "coordinate":
            revision_value = revision['suggested_value']

        print(
            f"\nInformation for combined and ordered revisions: {revision['revision_id']}\n"
            f"  Revision user display name: {revision['user_display_name']} (id: {revision['user_id']})\n"
            f"  Revision field name: {revision['field_name']}\n"
            f"  Revision value(s): {revision_value}"
        )



#only run the main function if this script is being run directly
if __name__ == "__main__":
    #test_variant = 1832 #Example variant POLE S459F (civic.vid: 1832)
    #test_variant = 785 #Example variant giving error
    #test_variant = 4050 #Example with duplicate hgvs expressions
    test_variant = 1686 

    contributor_id = 15 #Example user (Malachi Griffith, user id 15)
    main(test_variant, contributor_id)

