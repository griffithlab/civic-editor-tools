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

Before running for the first time the following python dependencies will need to be installed. If you are unsure it should do no harm to reinstall them.

- biopython
- civicpy
- requests
- certifi

Install dependencies with:
```bash
cd ~/git/civic-editor-tools
pip3 install -r requirements.txt
```

## Staging of reference data files

Before running for the first time, run these scripts to stage needed data files. If you are unsure it should do no harm to run these scripts again.

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

./utils/civic_graphql_utils.py
./utils/civicpy_utils.py
./utils/clingen_ar_utils.py
./utils/ensembl_utils.py
./utils/entrez_utils.py
./utils/generic_utils.py
./utils/refseq_utils.py

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





