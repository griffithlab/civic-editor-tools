#!/usr/bin/env python3

from pathlib import Path
import requests
import pdb

api_url = "https://civicdb.org/api/graphql"

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

def run_graphql_operation(api_url, operation_name, query_id, timeout=(10, 200)):
    # Base directory = directory containing this script
    base_dir = Path(__file__).resolve().parent

    query_path = base_dir / f"{operation_name}_query.json"
    variables_path = base_dir / f"{operation_name}_variables.json"

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

def main(variant_id):
	#Get coordinate ids for variant (takes a variant id)
	resp = run_graphql_operation(api_url, "variant_CoordinateIdsForVariant", variant_id)
	json = resp.json()

	open_revision_count_variant = json['data']['variant']['openRevisionCount']
	open_revision_count_coordinates = json['data']['variant']['coordinates']['openRevisionCount']
	variant_coordinates_id = json['data']['variant']['coordinates']['id']

	print(
		f"Variant ID used for graphql query: {variant_id}\n"
		f"  Open gene-variant revisions: {open_revision_count_variant}\n"
		f"  Open variant coordinate revisions: {open_revision_count_coordinates}\n"
		f"  Variant coordinates id: {variant_coordinates_id}"
	)

	#variant_Revisions-Variant (takes a variant id)
	resp = run_graphql_operation(api_url, "variant_Revisions-Variant", variant_id)
	json = resp.json()

	#pdb.set_trace()

	revisions = json["data"]["revisions"]["edges"]
	for revision in revisions:
		revision_id = revision['node']['id']
		user_display_name = revision['node']['creationActivity']['user']['displayName']
		revision_values = revision['node']['linkoutData']['diffValue']['addedObjects']
		revision_values_list = []
		for revision_value in revision_values:
			revision_display_name = revision_value['displayName']
			revision_values_list.append(revision_display_name)
		field_name = revision['node']['fieldName']
		revision_values_string = ",".join(sorted(revision_values_list))

		print(
			f"\nInformation for revision: {revision_id}\n"
			f"  Revision user display name: {user_display_name}\n"
			f"  Revision field name: {field_name}\n"
			f"  Revision value(s): {revision_values_string}"
		)

	#variant_Revisions-VariantCoordinates (takes a variant _coordinates_ id)
	resp = run_graphql_operation(api_url, "variant_Revisions-VariantCoordinates", variant_coordinates_id)
	json = resp.json()

	#variant_VariantDetail
	resp = run_graphql_operation(api_url, "variant_VariantDetail", variant_id)
	json = resp.json()

	#To interactively explore json responses that come back from these queryies, place this inline above:
	#pdb.set_trace()
	#json = resp.json()

	#json['data']['revisions']['edges'][0]['node']['creationActivity']['user']['displayName']
	#json['data']['revisions']['edges'][0]['node']['linkoutData']['diffValue']['addedObjects'][0]['displayName']
	#json['data']['revisions']['edges'][0]['node']['fieldName']

if __name__ == "__main__":
	test_variant = 1832 #Example variant POLE S459F (civic.vid: 1832)
	main(test_variant)


