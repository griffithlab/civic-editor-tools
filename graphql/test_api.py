#!/usr/bin/env python3

from pathlib import Path
import requests
import pdb

api_url = "https://civicdb.org/api/graphql"

# Read the JSON file into a string variable
#with open('variant_CoordinateIdsForVariant_query.json', 'r') as file:
#    query = file.read()

#with open('variant_CoordinateIdsForVariant_variables.json', 'r') as file:
#    variables = file.read()

#resp = requests.post(
#   api_url, json={"query": query, "variables": variables}, timeout=(10, 200)
#)

def run_graphql_operation(api_url, operation_name, timeout=(10, 200)):
    # Base directory = directory containing this script
    base_dir = Path(__file__).resolve().parent

    query_path = base_dir / f"{operation_name}_query.json"
    variables_path = base_dir / f"{operation_name}_variables.json"

    with query_path.open("r") as f:
        query = f.read()

    with variables_path.open("r") as f:
        variables = f.read()

    resp = requests.post(
        api_url,
        json={
            "query": query,
            "variables": variables
        },
        timeout=timeout
    )

    return resp

#Get coordinate ids for variant
operation_name = "variant_CoordinateIdsForVariant"
resp = run_graphql_operation(api_url, operation_name)
json = resp.json()
json['data']
open_revision_count_variant = json['data']['variant']['openRevisionCount']
open_revision_count_coordinates = json['data']['variant']['coordinates']['openRevisionCount']

print(open_revision_count_variant)
print(open_revision_count_coordinates)

#pdb.set_trace()
#json = resp.json()



