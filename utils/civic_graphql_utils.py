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


def run_graphql_operation(api_url, operation_name, query_id, timeout=(10, 200)):
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


def gather_user_details(user_id):
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


def gather_variant_details(variant_id):
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


def gather_accepted_variant_data(variant_id):
    """execute graphql queries, parse json, return detailed info on already accepted variant fields"""
    #graphql template name: "variant_VariantSummary"
    resp = run_graphql_operation(api_url, "variant_VariantSummary", variant_id)
    json = resp.json()

    #json['data']['variant']['variantAliases']

    accepted_variant_data = {
        "variant_aliases": json['data']['variant']['variantAliases']
    }
    pdb.set_trace()

    return accepted_variant_data


def gather_variant_revisions(variant_id, contributor_id):
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


def load_blacklisted_variant_ids(filepath):
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


def main (variant_id, contributor_id):
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
    #iterate through individual revisions
    variant_revisions = variant_data['variant_revisions']
    for variant_revision in variant_revisions:
        print(
            f"\nInformation for variant revision: {variant_revision['revision_id']}\n"
            f"  Revision user display name: {variant_revision['user_display_name']} (id: {variant_revision['user_id']})\n"
            f"  Revision field name: {variant_revision['field_name']}\n"
            f"  Revision value(s): {variant_revision['revision_values_string']}"
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

#only run the main function if this script is being run directly
if __name__ == "__main__":
    #test_variant = 1832 #Example variant POLE S459F (civic.vid: 1832)
    test_variant = 785 #Example variant giving error
    contributor_id = 15 #Example user (Malachi Griffith, user id 15)
    main(test_variant, contributor_id)

