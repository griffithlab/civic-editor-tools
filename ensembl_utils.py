#!/usr/bin/env python3

import gzip

def parse_ensembl_pep_fasta(fasta_path, strip_version=False):
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

            transcript_to_protein[transcript_id] = protein_id

    return transcript_to_protein


def build_transcript_biotype_map(fasta_paths):
    """
    Parse one or more gzipped ensembl FASTA files and return a dictionary:
        transcript_id -> transcript_biotype
    """
    transcript_map = {}

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

                transcript_map[transcript_id] = biotype

    return transcript_map


def main():

    mapping = parse_ensembl_pep_fasta(
        "data/ensembl/Homo_sapiens.GRCh38.pep.all.v115.fa.gz",
         strip_version=False
    )
    print(f"\nWith versions")
    print(f"Imported {len(mapping)} Ensembl transcript to protein mappings")
    enst1 = "ENST00000641515.2"
    ensp = mapping.get(enst1)
    print(f"Test mapping for: {enst1}: {ensp}\n")

    mapping = parse_ensembl_pep_fasta(
        "data/ensembl/Homo_sapiens.GRCh38.pep.all.v115.fa.gz",
         strip_version=True
    )
    print(f"\nWithout versions")
    print(f"Imported {len(mapping)} Ensembl transcript to protein mappings")
    enst2 = "ENST00000641515"
    ensp = mapping.get(enst2)
    print(f"Test mapping for: {enst2}: {ensp}\n")

    fasta_files = ["data/ensembl/Homo_sapiens.GRCh38.cdna.all.v115.fa.gz", 
                   "data/ensembl/Homo_sapiens.GRCh38.ncrna.v115.fa.gz"]
    transcript_to_biotype = build_transcript_biotype_map(fasta_files)

    print(f"Imported {len(transcript_to_biotype)} Ensembl transcript to biotype mappings")
    print(f"{enst1} -> {transcript_to_biotype[enst1]}")


if __name__ == "__main__":
    main()




