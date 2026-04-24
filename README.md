# CIViC editor tools

A collection of tools that help CIViC editors automate certain complex, time consuming tasks to increase the efficiency of editorial review and reduce errors and inconsistencies.


## What does `review_variant_coordinates.py` do?

CIViC curators frequently submit revisions to variant coordinate fields — genomic position, reference/alternate bases, representative transcript, HGVS expressions, ClinVar IDs, and related annotations. Manually verifying these revisions against authoritative external databases is tedious, error-prone, and time-consuming. `review_variant_coordinates.py` automates that cross-referencing workflow for a CIViC editor, presenting a field-by-field comparison (color-coded green/yellow/red) between what a curator submitted and what authoritative sources say.

### Scope and variant types supported

The tool currently supports **missense single-nucleotide variants (SNVs)** — variants with names following the amino acid substitution convention (e.g., `V600E`, `S459F`). It infers the variant type from the CIViC variant name itself, and skips variants of other types (fusions, indels, copy number events, etc.). Only variants with pending revisions submitted by *other* contributors are processed by default (a user cannot moderate their own submissions).

### Information sources

The tool integrates five external data sources:

| Source | What it provides |
|---|---|
| **CIViC** (GraphQL + CIViCpy) | Accepted variant fields, pending revision history, variant evidence sources |
| **ClinGen Allele Registry (CAR)** | Canonical Allele IDs (CAIDs), GRCh37/GRCh38 genomic coordinates, genomic and transcript HGVS expressions, ClinVar IDs, variant aliases, MANE Select designations |
| **NCBI Entrez / RefSeq** | Gene-to-transcript-to-protein mappings; RefSeq protein FASTA sequences |
| **Ensembl** | Transcript-to-protein mappings, transcript biotypes, protein FASTA sequences, GRCh37 transcript IDs (Ensembl v75 and v87) |
| **MANE** | MANE Select transcript designations for each gene |

### Checks performed

**1. Variant identity and amino acid validation**

The variant name (e.g., `V600E`) is parsed into three components: reference amino acid, codon position, and alternate amino acid. For each candidate transcript supported by the ClinGen Allele Registry, the tool fetches the actual protein sequence (from RefSeq or Ensembl) and confirms that the reference amino acid at the named position is correct. Transcripts where the reference amino acid does not match are excluded. This catches position errors and handles ambiguities such as methionine counting differences (where the initiator methionine is sometimes excluded from numbering).

**2. Transcript compatibility filtering**

ClinGen Allele Registry transcripts for the gene are filtered to retain only protein-coding transcripts (non-coding transcripts with Ensembl IDs are excluded using Ensembl biotype annotations) that have a verifiable protein sequence and a confirmed amino acid match at the variant position.

**3. Ambiguous genomic position detection**

In some cases, the same amino acid substitution name can refer to genuinely distinct genomic events — for example, when two different transcripts both carry the same reference amino acid at the named position but those codons map to different chromosomal locations. The tool checks whether all compatible CAIDs resolve to the same GRCh37 genomic position and warns if they do not.

**4. MANE Select name concordance**

The tool compares the CIViC variant name (converted to three-letter HGVS protein notation, e.g., `p.Val600Glu`) against the protein-level variant name derived from the MANE Select transcript for that gene. A mismatch here may indicate the variant is named according to a non-canonical transcript.

**5. Build37 representative transcript resolution**

CIViC uses GRCh37 coordinates and expects a GRCh37 Ensembl representative transcript (typically from Ensembl v75 or the GRCh37-imported Ensembl v87 annotation). The tool maps each compatible ClinGen transcript (which may be versioned for GRCh38) to its GRCh37 counterpart from both Ensembl v75 and v87, allowing version-number-tolerant matching.

**6. Field-by-field comparison**

For each compatible CAID, the tool compares every relevant CIViC field — both the currently accepted values and each pending revision — against the authoritative ClinGen value. Fields compared include:

- **ClinGen Allele Registry ID** (CAID)
- **Variant type** (inferred from the variant name)
- **Variant aliases** (e.g., HGVS protein short-form names, legacy names)
- **HGVS expressions** (genomic, coding, and protein; compared with transcript version tolerance)
- **ClinVar IDs**
- **Reference genome build** (GRCh37 expected)
- **Chromosome**
- **Genomic start position** (with automatic conversion from ClinGen's 0-based to CIViC's 1-based coordinate system)
- **Genomic stop position**
- **Reference bases**
- **Variant (alternate) bases**
- **Representative transcript** (with both exact and version-tolerant partial matching)
- **Ensembl version** (v75 or v87 expected for GRCh37 annotations)

Results are printed with green (match), yellow (qualified/partial match), and red (mismatch) highlighting so an editor can quickly assess which revision fields are well-supported and which warrant further investigation.

---

## Cloning the repository
First get a local copy of this repo

```
cd ~/git
git clone git@github.com:malachig/civic-editor-tools.git

```

## Requirements
This tool has been tested on MacOS only (but docker images are available, described below).

Before running for the first time the following python dependencies will need to be installed. If you are unsure it should do no harm to reinstall them.

- biopython
- civicpy
- requests
- certifi
- pymysql

Install Python dependencies with:
```bash
cd ~/git/civic-editor-tools
pip3 install -r requirements.txt

```

Install efetch
```bash
#install script
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

#test command
efetch -db nucleotide -id NM_020630.6 -format fasta
```


## Staging of reference data files

Before running for the first time, run these scripts to stage needed data files. If you are unsure it should do no harm to run these scripts again.

```bash
./stage_local_data.sh

```

## Test each individual module (optional)

Execute the following one at a time and examine the output for any errors
```bash
./utils/civic_graphql_utils.py
./utils/civicpy_utils.py
./utils/clingen_ar_utils.py
./utils/compare_utils.py 
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


## Running in Docker

Example command using the latest version of the image
```bash
docker run -it -v /Users/mgriffit/git/civic-editor-tools:/civic-editor-tools griffithlab/civic-editor-tools:latest 
cd /civic-editor-tools
./review_variant_coordinates.py --contributor-id 15 --variant-id 1832

```

Example command using a specific version of the image
```bash
docker run -it -v /Users/mgriffit/git/civic-editor-tools:/civic-editor-tools griffithlab/civic-editor-tools:1.0
cd /civic-editor-tools
./review_variant_coordinates.py --contributor-id 15 --variant-id 1832

```

---

## Developer notes

### Code organisation

```
review_variant_coordinates.py   # main entry point and workflow orchestration
utils/
  generic_utils.py              # amino acid conversions, variant name parsing, connectivity checks
  civic_graphql_utils.py        # CIViC GraphQL queries and response parsing
  civicpy_utils.py              # CIViC Python client (bulk variant retrieval, evidence sources)
  clingen_ar_utils.py           # ClinGen Allele Registry REST API calls and JSON extraction
  entrez_utils.py               # NCBI Entrez/RefSeq transcript-to-protein mappings, MANE data
  ensembl_utils.py              # Ensembl FASTA parsing, transcript maps, build37 transcript IDs
  refseq_utils.py               # RefSeq protein FASTA indexing and sequence retrieval
  compare_utils.py              # ValueComparator class: field-by-field ClinGen vs CIViC comparison
graphql/
  <operation>_query.json        # GraphQL query bodies (one file per operation)
  <operation>_variables.json    # Corresponding variable templates (one file per operation)
data/                           # locally staged reference files (not committed; see stage_local_data.sh)
```

The main script handles workflow logic and user interaction. The `utils/` modules are responsible for all external communication and data wrangling. Neither layer knows about the other's internals — the main script calls utility functions and passes their return values along.

### Utility module conventions

Every utility module has a `__main__` block at the bottom that exercises its key functions against real data. This doubles as a functional test and a way to populate any pickle caches the first time (see *Pickle caching* below). Run any module directly to smoke-test it in isolation:

```bash
./utils/ensembl_utils.py
```

Modules that import each other use a pattern that handles both package and standalone execution:

```python
try:
    from . import generic_utils
except ImportError:
    import generic_utils
```

### GraphQL query pattern

All CIViC GraphQL operations live in `graphql/` as paired JSON files: one for the query body and one for the variables. The naming convention is `<object>_<OperationName>_query.json` / `..._variables.json`. The variable template contains the literal placeholder string `graphql_query_id1`; `civic_graphql_utils.populate_variables_id()` performs a string substitution to inject the actual integer ID before submitting.

To add a new GraphQL operation:
1. Create `graphql/<object>_<OperationName>_query.json` with your query.
2. Create `graphql/<object>_<OperationName>_variables.json` with the variable template, using `graphql_query_id1` as the placeholder wherever the integer ID belongs.
3. Call `civic_graphql_utils.run_graphql_operation(api_url, "<object>_<OperationName>", id)` — it handles file loading, ID injection, and the POST.
4. Add a corresponding gather/parse function in `civic_graphql_utils.py` that extracts the fields you need from `resp.json()`.

### Pickle caching pattern

Several data structures are expensive to build from raw files (Ensembl FASTA parsing, build37 GTF parsing, RefSeq transcript maps). These use a consistent lazy-load-and-cache pattern:

```python
def load_<thing>(...):
    if os.path.exists(pickle_path):
        return load_transcript_map_pickle(pickle_path)
    else:
        data = compile_<thing>(...)           # expensive build step
        save_transcript_map_pickle(data, pickle_path)
        return data
```

The first call builds and saves the pickle; subsequent calls just load it. Pickle files live in `data/` alongside their source files. If you update the underlying source data, delete the corresponding `.pkl` file to force a rebuild on the next run. Follow the same pattern for any new expensive-to-build lookup tables.

### FASTA indexing pattern

Protein sequence lookup uses `Bio.SeqIO.index_db()` (a persistent SQLite-backed index) rather than loading entire FASTA files into memory. Both RefSeq and Ensembl protein FASTAs are pre-indexed into `.idx` files stored next to the FASTA in `data/*/indexed/`. The index is built once (if it doesn't exist) and opened read-only on every subsequent lookup. This allows O(1) retrieval by protein ID from files that are gigabytes in size.

`refseq_utils.get_refseq_protein_indexed()` and `ensembl_utils.get_ensembl_protein_indexed()` follow the same call signature and return a plain string sequence. When adding support for a new sequence type, follow this pattern rather than loading the whole FASTA.

### ClinGen Allele Registry workflow

The ClinGen lookups follow a two-stage pipeline:

**Stage 1 — transcript discovery (per gene, cached across variants)**

1. `get_reference_sequences_by_gene(gene_name)` — fetches all CAR-supported transcripts for the gene.
2. `extract_reference_sequences()` — filters to transcript-type entries, strips predicted transcripts (`XM_`, `XR_`), prefers NCBI IDs over Ensembl when both are present.
3. `keep_latest_transcript_versions()` — deduplicates by base ID, retaining only the highest version number.
4. The gene-level transcript list is cached in the `clingen_transcript_ids` dict in `main()` so the API is only called once per gene, not once per variant.

**Stage 2 — allele resolution (per variant)**

5. Each transcript is filtered to protein-coding only (Ensembl biotype check) and amino acid validated (the reference AA at the named position must match the actual protein sequence).
6. Surviving protein IDs are assembled into protein HGVS expressions (e.g. `NP_006222.2:p.Ser459Phe`) and submitted to `get_allele_by_hgvs()`.
7. The returned protein allele JSON is walked to extract transcript-level CAIDs via `extract_transcript_cas()`.
8. Each CAID is then fetched individually with `get_allele_by_id()` to obtain genomic coordinates, HGVS expressions, aliases, ClinVar IDs, and MANE annotations.

### ValueComparator and adding a new field comparison

`compare_utils.ValueComparator` is the central comparison engine. It uses a dispatch dictionary (`self._dispatch`) that maps CIViC field name strings to handler methods. To add comparison logic for a new field:

1. Add a method `compare_<fieldname>(self, civic_value)` to `ValueComparator`. Use `self.clingen_data[key]` to retrieve the corresponding authoritative value, `self._print_revision_details()` to print the field label and revision ID, and `self._print_match(MatchLevel.MATCH/MISMATCH/QUALIFIED_MATCH, message)` for colour-coded output.
2. Register it in `self._dispatch` inside `__init__`: `"your_field_name": self.compare_<fieldname>`.
3. Add the field name to `FIELD_NAME_PRIORITY` in `civic_graphql_utils.merge_revision_data()` at the position where it should appear in the output. The `merge_revision_data()` function raises a `ValueError` for any field name not listed there, so this step is enforced at runtime.

`MatchLevel.QUALIFIED_MATCH` (yellow) is for partial or version-tolerant matches — use it when the comparison is meaningful but not an exact string match (e.g., same transcript base ID but different version number).

### Adding support for a new variant type

Currently only missense SNVs (`Missense Variant`) are processed end-to-end. The gate is in `main()`:

```python
if variant_type_is_unsupported(..., guessed_gene_variant_type, "Missense Variant"): continue
```

`generic_utils.guess_variant_type()` parses the CIViC variant name with a series of regex patterns (frameshift first, then missense SNV, then splice site). To add a new type:

1. Add a detection branch in `guess_variant_type()` returning a new type string.
2. Add a `parse_<type>_name_components()` function in `generic_utils.py` analogous to `parse_snv_coding_name_components()`.
3. In `main()`, change the `variant_type_is_unsupported()` call (or add a parallel branch) so variants of the new type are not skipped.
4. Implement whatever transcript filtering and HGVS construction logic is needed for the new type before the ClinGen allele lookup.

### Coordinate system convention

ClinGen Allele Registry reports genomic coordinates in **0-based, half-open** intervals (start is 0-based). CIViC stores coordinates in **1-based** positions. The `compare_start()` method in `ValueComparator` adds 1 to the ClinGen start before comparison; `compare_stop()` does not adjust the end coordinate. Keep this in mind whenever you add new genomic coordinate comparisons.

### Debugging tips

Several functions in `civic_graphql_utils.py` contain commented-out `pdb.set_trace()` calls and example JSON path expressions. Uncomment these to drop into an interactive debugger and inspect the raw API responses when developing new GraphQL queries or response parsers. The comments above each `pdb.set_trace()` show example JSON navigation paths for the current query's response structure.

