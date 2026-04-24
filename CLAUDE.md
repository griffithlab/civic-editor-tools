# CLAUDE.md — civic-editor-tools

## Project overview
A Python CLI toolkit for CIViC database editors to automate variant coordinate review. The main script cross-references ClinGen Allele Registry, NCBI Entrez, Ensembl, and RefSeq to validate genomic coordinates.

- **Main script:** `review_variant_coordinates.py`
- **Utilities:** `utils/` (civic_graphql_utils, civicpy_utils, clingen_ar_utils, ensembl_utils, entrez_utils, refseq_utils, compare_utils, generic_utils)
- **GraphQL queries:** `graphql/`
- **Reference data:** `data/` (staged locally, not committed)

## Python dependencies
Managed in `requirements.txt`:
- `biopython` — Entrez queries, FASTA parsing (`Bio.SeqIO`, `Bio.Entrez`)
- `civicpy` — CIViC API client
- `requests` — HTTP calls (GraphQL, ClinGen AR)
- `certifi` — SSL certificate bundle
- `pymysql` — Ensembl public MySQL (`data/ensembl/get_ensembl_v75_version_numbers.py`)

Install:
```bash
pip3 install -r requirements.txt
```

## Staging reference data
Run once before first use (safe to re-run; skips already-downloaded files):
```bash
./stage_local_data.sh
```

This runs:
- `data/ensembl/get_ensembl_fastas.sh` — Ensembl protein FASTAs
- `data/ensembl/get_ensembl_grch37_info.sh` — Ensembl GRCh37 transcript info
- `data/ensembl/get_ensembl_v75_version_numbers.py` — Ensembl v75 version numbers (via public MySQL)
- `data/entrez/get_refseq_mappings.sh` — Entrez gene-to-RefSeq mappings
- `data/refseq/get_mane_summary.sh` — MANE select transcript info
- `data/refseq/get_refseq_protein.sh` — RefSeq protein FASTA

## Running the main tool
```bash
# Help
./review_variant_coordinates.py --help

# Single variant
./review_variant_coordinates.py --contributor-id 15 --variant-id 1832
```

## Testing utility modules
Each utility module has a `__main__` block and can be run standalone to test and populate pickle caches:
```bash
./utils/civic_graphql_utils.py
./utils/civicpy_utils.py
./utils/clingen_ar_utils.py
./utils/ensembl_utils.py
./utils/entrez_utils.py
./utils/generic_utils.py
./utils/refseq_utils.py
```

## Docker
Build and run a containerized environment (Ubuntu 24.04, Python 3.12):
```bash
docker build -t civic-editor-tools .
docker run -it civic-editor-tools
```

Note: reference data in `data/` is not staged inside the image by default — mount it or run `stage_local_data.sh` inside the container after starting.

## Platform
Tested on macOS only.

## Claude preferences
Do not commit or push to github automatically, wait until told to do so, or let me do it manually after review

