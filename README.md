# CIViC editor tools

A collection of tools that help CIViC editors automate certain complex, time consuming tasks to increase the efficiency of editorial review and reduce errors and inconsistencies.


## Cloning the repository
First get a local copy of this repo

```
cd ~/git
git clone git@github.com:malachig/civic-editor-tools.git

```


## Requirements
This tool has been tested on MacOS only.

Before running for the first time the following python dependencies will need to be installed.

- biopython
- civicpy
- requests
- certifi

Install dependencies with:
```bash
pip3 install -r requirements.txt
```

## Staging of reference data files

Before running for the first time, run these scripts to 

```bash
cd ~/git/civic-editor-tools

./data/ensembl/get_ensembl_fastas.sh
./data/entrez/get_refseq_mappings.sh
./data/refseq/get_mane_summary.sh
./data/refseq/get_refseq_protein.sh

```

## Test each individual module and create pickle caches

Execute the following one at a time and examine the output for any errors
```bash

./civic_graphql_utils.py
./civicpy_utils.py
./clingen_ar_utils.py
./ensembl_utils.py
./entrez_utils.py
./generic_utils.py
./refseq_utils.py

```

## Determine your civic user ID

Go to [CIViC](civicdb.org) and find your user ID from the [Contributors View](https://civicdb.org/users/home). Once you have found your user profile, your user ID will be the number in the URL.  For example: https://civicdb.org/users/15/activity (15 is the ID for user MalachiGriffith)

## Run a test of an editor tool

```bash
# get help docs
./review_variant_coordinates.py --help

# example test command for a single variant
./review_variant_coordinates.py --contributor-id 15 --variant-id 1832

```

If all this has worked smoothly you should be ready to start using the editor tool suite. 

## Troubleshooting

...





