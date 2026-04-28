#!/usr/bin/env python3

from pathlib import Path
import requests
import pdb

api_url = "https://civicdb.org/api/graphql"
base_dir = Path(__file__).resolve().parent

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


def run_graphql_operation(api_url: str, operation_name: str, query_id: int, timeout=(10, 200)) -> str:
    """load graphql query and variable json objects from file, update with a query id, and submit the query to the API"""
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

    #inject the actual query ID into the variables object
    variables_updated = populate_variables_id(variables, query_id)

    resp = requests.post(
        api_url,
        json={
            "query": query,
            "variables": variables_updated
        },
        timeout=timeout
    )

    return resp


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

        #special handling when the revision is the variant "name" itself
        if field_name == 'name':
            current_value = revision['node']['currentValue']
            suggested_value = revision['node']['suggestedValue']
            revision_values_string = f"'{current_value}' -> '{suggested_value}'"
            variant_data['name_change'] = True
        else:
            #revisions that are lists of things
            revision_values = revision['node']['linkoutData']['diffValue']['addedObjects']
            revision_values_list = []
            for revision_value in revision_values:
                revision_display_name = revision_value['displayName']
                revision_values_list.append(revision_display_name)
            revision_values_string = ",".join(sorted(revision_values_list))

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

