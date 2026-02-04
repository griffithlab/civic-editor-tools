#!/usr/bin/env python3

import requests
API_URL = "https://civicdb.org/api/graphql"

# Read the JSON file into a string variable
with open('variant_CoordinateIdsForVariant_query.json', 'r') as file:
    payload = file.read()
print(payload)

with open('variant_CoordinateIdsForVariant_variables.json', 'r') as file:
    variables = file.read()
print(variables)

resp = requests.post(
   API_URL, json={"query": payload, "variables": variables}, timeout=(10, 200)
)

print(resp)

