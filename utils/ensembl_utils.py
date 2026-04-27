#!/usr/bin/env python3

import os
import re
import gzip
import pickle
from Bio import SeqIO
from pathlib import Path

base_dir = Path(__file__).resolve().parent

def build_transcript_to_protein_id_map(fasta_path, transcript_to_protein=None, strip_version=False):
    """
    Parse an Ensembl pep.all FASTA file and return a dictionary mapping
    Ensembl transcript IDs (ENST...) to Ensembl protein IDs (ENSP...).

    Parameters
    ----------
    fasta_path : str
        Path to the .fa.gz protein FASTA file.
    strip_version : bool
        If True, remove version suffix (e.g. .2) from IDs.
    """

    if transcript_to_protein is None:
        transcript_to_protein = {}

    with gzip.open(fasta_path, "rt") as f:
        for line in f:
            if not line.startswith(">"):
                continue

            header = line.strip()

            # Protein ID is first token after ">"
            protein_id = header.split()[0][1:]  # remove ">"
            
            # Find transcript:ENST...
            transcript_id = None
            for field in header.split():
                if field.startswith("transcript:"):
                    transcript_id = field.replace("transcript:", "")
                    break

            if transcript_id is None:
                continue  # skip if no transcript field

            if strip_version:
                protein_id = protein_id.split(".")[0]
                transcript_id = transcript_id.split(".")[0]

            # Only add if transcript not already present
            if transcript_id not in transcript_to_protein:
                transcript_to_protein[transcript_id] = protein_id

    return transcript_to_protein


def build_transcript_biotype_map(fasta_paths, transcript_to_biotype=None):
    """
    Parse one or more gzipped ensembl FASTA files and return a dictionary:
        transcript_id -> transcript_biotype
    """
    if transcript_to_biotype is None:
        transcript_to_biotype = {}

    for fasta_path in fasta_paths:
        with gzip.open(fasta_path, "rt") as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue

                line = line.strip()

                # First token after '>'
                first_token = line[1:].split()[0]
                transcript_id = first_token  # keeps version (e.g. ENST00000389680.2)

                # Find transcript_biotype field
                biotype = None
                for field in line.split():
                    if field.startswith("transcript_biotype:"):
                        biotype = field.split(":", 1)[1]
                        break

                if biotype is None:
                    raise ValueError(
                        f"No transcript_biotype found for {transcript_id} in {fasta_path}"
                    )
                # Only add if transcript not already present
                if transcript_id not in transcript_to_biotype:
                    transcript_to_biotype[transcript_id] = biotype

    return transcript_to_biotype


def load_ensembl_versions(filepath):
    """
    Load Ensembl versions from a file (one version per line).
    Returns a list of strings.
    """
    versions = []

    with open(filepath, "r") as f:
        for line in f:
            v = line.strip()
            if not v:
                continue          # skip blank lines
            if v.startswith("#"):
                continue          # skip comments
            versions.append(v)

    return versions


def compile_transcript_to_protein_map(ensembl_versions_file):
    """
    Use the ensembl peptide parsing method to compile a transcript to protein id map 
    from multiple version of ensembl defined in a file
    """

    print(f"\nLoading ensembl transcript to protein map from multiple versions of ensembl")
    ensembl_transcript_to_protein_map = {}
    versions = load_ensembl_versions(ensembl_versions_file)
    for version in versions:
        pep_fasta_file = base_dir / f"../data/ensembl/version_data/Homo_sapiens.GRCh38.pep.all.v{version}.fa.gz"
        ensembl_transcript_to_protein_map = build_transcript_to_protein_id_map(
            fasta_path=pep_fasta_file,
            transcript_to_protein=ensembl_transcript_to_protein_map,
            strip_version=False
        )
        print(f"  Loaded ensembl version: {version}. Total mappings so far: {len(ensembl_transcript_to_protein_map)}")

    #In addition to the pickle already created, save the ensembl_transcript_to_biotype_map to a TSV to allow inspection
    #note this method only runs if the pickle is being generated, delete it to regenerate both files
    ensembl_transcript_to_protein_map_tsv = base_dir / f"../data/ensembl/ensembl_transcript_to_protein_map.tsv"
    with open(ensembl_transcript_to_protein_map_tsv, "w") as f:
        f.write("transcript_id\tprotein_id\n")
        for transcript_id, protein_id in ensembl_transcript_to_protein_map.items():
            f.write(f"{transcript_id}\t{protein_id}\n")
    print(f"  Saved transcript to protein map to: {ensembl_transcript_to_protein_map_tsv}")

    return ensembl_transcript_to_protein_map


def compile_transcript_to_biotype_map(ensembl_versions_file):
    """
    Use the ensembl biotype parsing method to compile a transcript to biotype map from 
    multiple version of ensembl defined in a file
    """

    print(f"\nLoading ensembl transcript to biotype map from multiple versions of ensembl")
    ensembl_transcript_to_biotype_map = {}
    versions = load_ensembl_versions(ensembl_versions_file)
    for version in versions:
        fasta_files = [base_dir / f"../data/ensembl/version_data/Homo_sapiens.GRCh38.cdna.all.v{version}.fa.gz",
                       base_dir / f"../data/ensembl/version_data/Homo_sapiens.GRCh38.ncrna.v{version}.fa.gz"]
        transcript_to_biotype = build_transcript_biotype_map(
            fasta_files,
            transcript_to_biotype=ensembl_transcript_to_biotype_map
        )
        print(f"  Loaded ensembl version: {version}. Total mappings so far: {len(transcript_to_biotype)}")

    #In addition to the pickle already created, save the ensembl_transcript_to_biotype_map to a TSV to allow inspection
    #note this method only runs if the pickle is being generated, delete it to regenerate both files
    ensembl_transcript_to_biotype_map_tsv = base_dir / f"../data/ensembl/ensembl_transcript_to_biotype_map.tsv"
    with open(ensembl_transcript_to_biotype_map_tsv, "w") as f:
        f.write("transcript_id\tbiotype\n")
        for transcript_id, biotype in ensembl_transcript_to_biotype_map.items():
            f.write(f"{transcript_id}\t{biotype}\n")
    print(f"  Saved transcript to biotype map to: {ensembl_transcript_to_biotype_map_tsv}")

    return ensembl_transcript_to_biotype_map


def build_ensembl_fasta_index(ensembl_versions_file):
    """Creates a persistent index file (.idx) for each ensembl fasta file"""

    versions = load_ensembl_versions(ensembl_versions_file)

    for version in versions:
        ensembl_fasta_path = base_dir / f"../data/ensembl/version_data/indexed/Homo_sapiens.GRCh38.pep.all.v{version}.fa"
        ensembl_fasta_index_path = f"{ensembl_fasta_path}.idx"

        if not os.path.exists(ensembl_fasta_index_path):
            print(f"Building index for version {version}...")
            SeqIO.index_db(str(ensembl_fasta_index_path), str(ensembl_fasta_path), "fasta")


def get_ensembl_protein_indexed(ensembl_protein_id, ensembl_versions_file):
    """
    Given and ensembl protein id and a prebuilt SeqIO index, retrieve the protein sequence
    Search across multiple version of Ensembl but stop searching once the id is found
    """
    protein_seq = None

    versions = load_ensembl_versions(ensembl_versions_file)

    for version in versions:
        index_path = base_dir / f"../data/ensembl/version_data/indexed/Homo_sapiens.GRCh38.pep.all.v{version}.fa.idx"
        index = SeqIO.index_db(str(index_path))

        if ensembl_protein_id not in index:
            continue
        else:
            protein_seq = str(index[ensembl_protein_id].seq)
            return protein_seq

    if protein_seq is None:
        raise ValueError(
            f"No protein sequence found for {ensembl_protein_id} in ensembl versions: {ensembl_versions_file}"
        )

def parse_gtf_attributes(attr_string: str) -> dict:
    """Parse the semicolon-separated key-value attributes in column 9 of a GTF."""
    return {m.group(1): m.group(2) for m in re.finditer(r'(\w+) "([^"]+)"', attr_string)}

def compile_build37_transcripts():
    """
    Parse old ensembl build37 files and build an index of Ensembl transcript IDs (with version numbers) grouped by the Ensembl annotation version
    """

    #First get the v75 transcript ids with versions (parsing a TSV obtained by SQL query of Ensembl database)
    ensembl_v75_input_path = base_dir / f"../data/ensembl/build37/ensembl75_transcripts.tsv"
    annotations_v75 = {}
    with open(ensembl_v75_input_path) as f:
        header = f.readline().strip().split("\t")
        col = {name: i for i, name in enumerate(header)}

        for line in f:
            parts = line.strip().split("\t")
            if not parts or parts == [""]:
                continue

            versioned_id = f"{parts[col['transcript_id']]}.{parts[col['transcript_version']]}"
            annotations_v75[versioned_id] = {
                "gene_id": parts[col["gene_id"]],
                "gene_name": parts[col["gene_name"]],
                "biotype": parts[col["biotype"]],
            }
    
    #Now get the v87 build37 imported transcript ids with verions (parsing a GTF that contains this info)
    ensembl_v87_input_path = base_dir / f"../data/ensembl/build37/Homo_sapiens.GRCh37.87.gtf.gz"
    annotations_v87 = {}
    opener = gzip.open if str(ensembl_v87_input_path).endswith(".gz") else open
    with opener(ensembl_v87_input_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if fields[2] != "transcript":
                continue

            attrs = parse_gtf_attributes(fields[8])

            # Skip if missing required fields
            if not all(k in attrs for k in ("transcript_id", "transcript_version", "gene_id", "gene_name", "transcript_biotype")):
                continue

            versioned_id = f"{attrs['transcript_id']}.{attrs['transcript_version']}"
            annotations_v87[versioned_id] = {
                "gene_id": attrs["gene_id"],
                "gene_name": attrs["gene_name"],
                "biotype": attrs["transcript_biotype"],
            }

    build37_ensembl_transcripts = {
        "v75": annotations_v75,
        "v87": annotations_v87,
    }

    return build37_ensembl_transcripts


def save_transcript_map_pickle(transcript_map, output_path):
    """Given a dictionary and an output path, save a pickle file"""
    with open(output_path, "wb") as f:
        pickle.dump(transcript_map, f, protocol=pickle.HIGHEST_PROTOCOL)


def load_transcript_map_pickle(input_path):
    """Given a path to a pickle file, load it"""
    with open(input_path, "rb") as f:
        return pickle.load(f)


def load_ensembl_transcript_to_protein_map(ensembl_versions_file):
    """Load and return an ensembl transcript to protein map from a pickle, create it if needed"""
    transcript_to_protein_map_path = base_dir / f"../data/ensembl/ensembl_transcript_to_protein.pkl"
    ensembl_transcript_to_protein_map = {}
    if os.path.exists(transcript_to_protein_map_path):
        print(f"Transcript to protein map pickle exists, loading directly from: {transcript_to_protein_map_path}")
        ensembl_transcript_to_protein_map = load_transcript_map_pickle(transcript_to_protein_map_path)
    else:
        print(f"Transcript to protein map pickle does NOT exist, creating and saving to: {transcript_to_protein_map_path}")
        ensembl_transcript_to_protein_map = compile_transcript_to_protein_map(ensembl_versions_file)
        save_transcript_map_pickle(ensembl_transcript_to_protein_map, transcript_to_protein_map_path)

    return ensembl_transcript_to_protein_map

def load_ensembl_transcript_to_biotype_map(ensembl_versions_file):
    """Load and return an ensembl transcript to biotype map from a pickle, create it if needed"""
    transcript_to_biotype_map_path = base_dir / f"../data/ensembl/ensembl_transcript_to_biotype.pkl"
    ensembl_transcript_to_biotype_map = {}
    if os.path.exists(transcript_to_biotype_map_path):
        print(f"Transcript to biotype map pickle exists, loading directly from: {transcript_to_biotype_map_path}")
        ensembl_transcript_to_biotype_map = load_transcript_map_pickle(transcript_to_biotype_map_path)
    else:
        print(f"Transcript to biotype map pickle does NOT exist, creating and saving to: {transcript_to_biotype_map_path}")
        ensembl_transcript_to_biotype_map = compile_transcript_to_biotype_map(ensembl_versions_file)
        save_transcript_map_pickle(ensembl_transcript_to_biotype_map, transcript_to_biotype_map_path)

    return ensembl_transcript_to_biotype_map

def load_build37_ensembl_transcripts():
    """Load and return a build37 ensembl transcripts with version numbers man from a pickle, create it if needed"""
    build37_ensembl_transcripts = {}
    build37_ensembl_transcripts_path = base_dir / f"../data/ensembl/build37/build37_ensembl_transcripts.pkl"
    
    if os.path.exists(build37_ensembl_transcripts_path):
        print(f"Build 37 ensembl transcripts pickle exists, loading directly from: {build37_ensembl_transcripts_path}")
        build37_ensembl_transcripts = load_transcript_map_pickle(build37_ensembl_transcripts_path)
    else:
        print(f"Build37 ensembl transcripts pickle does NOT exist, creating and saving to: {build37_ensembl_transcripts_path}")
        build37_ensembl_transcripts = compile_build37_transcripts()
        save_transcript_map_pickle(build37_ensembl_transcripts, build37_ensembl_transcripts_path)

    return build37_ensembl_transcripts

def main():

    #load ensembl transcript to protein map from a single file (keep version numbers)
    fasta_path = base_dir / f"../data/ensembl/version_data/Homo_sapiens.GRCh38.pep.all.v115.fa.gz"
    ensembl_transcript_to_protein_map = build_transcript_to_protein_id_map(
         fasta_path,
         strip_version=False
    )
    print(f"\nWith versions")
    print(f"Imported {len(ensembl_transcript_to_protein_map)} Ensembl transcript to protein mappings")
    enst1 = "ENST00000641515.2"
    ensp1 = ensembl_transcript_to_protein_map.get(enst1)
    print(f"Test ensembl transcript to protein map for: {enst1}: {ensp1}\n")

    #load ensembl transcript to protein map from a single file (remove version numbers)
    ensembl_transcript_to_protein_map = build_transcript_to_protein_id_map(
         fasta_path,
         strip_version=True
    )
    print(f"\nWithout versions")
    print(f"Imported {len(ensembl_transcript_to_protein_map)} Ensembl transcript to protein mappings")
    enst2 = "ENST00000641515"
    ensp2 = ensembl_transcript_to_protein_map.get(enst2)
    print(f"Test ensembl transcript to protein map for: {enst2}: {ensp2}\n")

    ensembl_versions_file = base_dir / f"../data/ensembl/ensembl_versions.txt"

    #Now create a more comprehensive transcript to protein map that incorporates multiple Ensembl versions
    ensembl_transcript_to_protein_map = {}
    transcript_to_protein_map_path = base_dir / f"../data/ensembl/ensembl_transcript_to_protein.pkl"
    if os.path.exists(transcript_to_protein_map_path):
        print(f"Transcript to protein map pickle exists, loading directly from: {transcript_to_protein_map_path}")
        ensembl_transcript_to_protein_map = load_transcript_map_pickle(transcript_to_protein_map_path)
    else:
        print(f"Transcript to protein map pickle does NOT exist, creating and saving to: {transcript_to_protein_map_path}")
        ensembl_transcript_to_protein_map = compile_transcript_to_protein_map(ensembl_versions_file)
        save_transcript_map_pickle(ensembl_transcript_to_protein_map, transcript_to_protein_map_path)

    ensp1 = ensembl_transcript_to_protein_map.get(enst1)
    print(f"Test ensembl transcript to protein map using multi-version object for: {enst1}: {ensp1}\n")

    #load ensembl transcript to biotype map from a single file
    fasta_files = [base_dir / f"../data/ensembl/version_data/Homo_sapiens.GRCh38.cdna.all.v115.fa.gz", 
                   base_dir / f"../data/ensembl/version_data/Homo_sapiens.GRCh38.ncrna.v115.fa.gz"]
    ensembl_transcript_to_biotype_map = build_transcript_biotype_map(fasta_files)

    print(f"\nImported {len(ensembl_transcript_to_biotype_map)} Ensembl transcript to biotype mappings")
    print(f"{enst1} -> {ensembl_transcript_to_biotype_map[enst1]}\n")

    #Now create a more comprehensive transcript to biotype map that incorporated multiple Ensembl versions
    ensembl_transcript_to_biotype_map = {}
    transcript_to_biotype_map_path = base_dir / f"../data/ensembl/ensembl_transcript_to_biotype.pkl"
    if os.path.exists(transcript_to_biotype_map_path):
        print(f"Transcript to biotype map pickle exists, loading directly from: {transcript_to_biotype_map_path}")
        ensembl_transcript_to_biotype_map = load_transcript_map_pickle(transcript_to_biotype_map_path)
    else:
        print(f"Transcript to biotype map pickle does NOT exist, creating and saving to: {transcript_to_biotype_map_path}")
        ensembl_transcript_to_biotype_map = compile_transcript_to_biotype_map(ensembl_versions_file)
        save_transcript_map_pickle(ensembl_transcript_to_biotype_map, transcript_to_biotype_map_path)

    print(f"\nImported {len(ensembl_transcript_to_biotype_map)} Ensembl transcript to biotype mappings")
    print(f"{enst1} -> {ensembl_transcript_to_biotype_map[enst1]}")

    #Build fasta indexes for each ensembl peptide fasta
    build_ensembl_fasta_index(ensembl_versions_file)

    #Extract a protein sequence from an ensembl protein id
    ensp_seq = get_ensembl_protein_indexed(ensp1, ensembl_versions_file)

    print(f"\nProtein sequence for ensembl protein: {ensp1}")
    print(f"{ensp_seq}")


    #Extract build37 transcript IDs with versions from ensembl v75 and imported ensembl v87
    print(f"\nGet build 37 Ensembl transcripts with versions:")
    build37_ensembl_transcripts = {}
    build37_ensembl_transcripts_path = base_dir / f"../data/ensembl/build37/build37_ensembl_transcripts.pkl"

    if os.path.exists(build37_ensembl_transcripts_path):
        print(f"Build 37 ensembl transcripts pickle exists, loading directly from: {build37_ensembl_transcripts_path}")
        build37_ensembl_transcripts = load_transcript_map_pickle(build37_ensembl_transcripts_path)
    else:
        print(f"Build37 ensembl transcripts pickle does NOT exist, creating and saving to: {build37_ensembl_transcripts_path}")
        build37_ensembl_transcripts = compile_build37_transcripts()
        save_transcript_map_pickle(build37_ensembl_transcripts, build37_ensembl_transcripts_path)

    build37_ensembl_transcripts_v75 = build37_ensembl_transcripts["v75"]
    build37_ensembl_transcripts_v87 = build37_ensembl_transcripts["v87"]

    print(f"\nEnsembl v75:")
    for versioned_id, entry in list(build37_ensembl_transcripts_v75.items())[:5]:
        print(versioned_id, entry)
    print(f"\nEnsembl v87:")
    for versioned_id, entry in list(build37_ensembl_transcripts_v87.items())[:5]:
        print(versioned_id, entry)



if __name__ == "__main__":
    main()




