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

def main():

    mapping = parse_ensembl_pep_fasta(
        "data/ensembl/Homo_sapiens.GRCh38.pep.all.Ensembl.v115.fa.gz",
         strip_version=False
    )
    print(f"\nWith versions")
    print(f"Imported {len(mapping)} Ensembl transcript to protein mappings")
    enst = "ENST00000641515.2"
    ensp = mapping.get(enst)
    print(f"Test mapping for: {enst}: {ensp}\n")

    mapping = parse_ensembl_pep_fasta(
        "data/ensembl/Homo_sapiens.GRCh38.pep.all.Ensembl.v115.fa.gz",
         strip_version=True
    )
    print(f"\nWithout versions")
    print(f"Imported {len(mapping)} Ensembl transcript to protein mappings")
    enst = "ENST00000641515"
    ensp = mapping.get(enst)
    print(f"Test mapping for: {enst}: {ensp}\n")


if __name__ == "__main__":
    main()




