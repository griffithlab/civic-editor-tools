#!/usr/bin/env python3

from pathlib import Path
import requests
import pdb

api_url = "https://civicdb.org/api/graphql"

#update a template graphql object string to inject the query ID to be used
def populate_variables_id(variables_template: str, graphql_id: int) -> str:
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

#load graphql query and variable json objects from file, update with a query id, and submit the query to the API
def run_graphql_operation(api_url, operation_name, query_id, timeout=(10, 200)):
    # Base directory = directory containing this script
    base_dir = Path(__file__).resolve().parent

    query_path = base_dir / f"graphql/{operation_name}_query.json"
    variables_path = base_dir / f"graphql/{operation_name}_variables.json"

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

#execute graphql queries, parse json returns, build a simplied data structure with just the info needed
def gather_variant_revisions(variant_id, contributor_id):
	#Get coordinate ids for variant (takes a variant id)
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

	#variant_VariantDetail
	resp = run_graphql_operation(api_url, "variant_VariantDetail", variant_id)
	json = resp.json()
	variant_name = json['data']['variant']['name']
	feature_name = json['data']['variant']['feature']['name']
	
	variant_data['variant_name'] = variant_name
	variant_data['feature_name'] = feature_name	

	#variant_Revisions-Variant (takes a variant id)
	resp = run_graphql_operation(api_url, 'variant_Revisions-Variant', variant_id)
	json = resp.json()

	revisions = json["data"]["revisions"]["edges"]
	i = 0
	for revision in revisions:
		revision_id = revision['node']['id']
		user_id = revision['node']['creationActivity']['user']['id']
		if user_id == contributor_id:
			ariant_data['contributor_revisions'] += 1
		user_display_name = revision['node']['creationActivity']['user']['displayName']
		revision_values = revision['node']['linkoutData']['diffValue']['addedObjects']
		revision_values_list = []
		for revision_value in revision_values:
			revision_display_name = revision_value['displayName']
			revision_values_list.append(revision_display_name)
		field_name = revision['node']['fieldName']
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

#demonstrate functionality of the methods above and variant data retrieved
def main (variant_id, contributor_id):
	variant_data = gather_variant_revisions(variant_id, contributor_id)
	print(
		f"\n\nVariant revision info from gather_variant_revisions()\n"
		f"Variant ID used for graphql query: {variant_data['variant_id']}\n"
		f"  Variant name: {variant_data['variant_name']}\n"
		f"  Feature name: {variant_data['feature_name']}\n"
		f"  Open gene-variant revisions (total): {variant_data['open_revision_count_variant']}\n"
        f"  Open gene-variant revisions from specified contributor: {variant_data['contributor_revisions']}\n"
		f"  Open gene-variant revisions from all others users: {variant_data['open_revisions_non_contributor']}\n"
		f"  Variant coordinates id: {variant_data['variant_coordinates_id']}"
	)
	variant_revisions = variant_data['variant_revisions']
	for variant_revision in variant_revisions:
		print(
			f"\nInformation for variant revision: {variant_revision['revision_id']}\n"
			f"  Revision user display name: {variant_revision['user_display_name']} (id: {variant_revision['user_id']})\n"
			f"  Revision field name: {variant_revision['field_name']}\n"
			f"  Revision value(s): {variant_revision['revision_values_string']}"
		)

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
	test_variant = 1832 #Example variant POLE S459F (civic.vid: 1832)
	contributor_id = 15 #Example user (Malachi Griffith, user id 15)
	main(test_variant, contributor_id)


