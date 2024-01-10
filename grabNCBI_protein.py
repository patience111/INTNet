from Bio import Entrez, SeqIO

def fetch_protein_info(protein_name):
    Entrez.email = "xiangyu270208@gmail.com"  # Set your email for NCBI API usage

    # Search for the protein of interest in NCBI
    handle = Entrez.esearch(db="protein", term=protein_name)
    record = Entrez.read(handle)
    handle.close()

    if len(record["IdList"]) == 0:
        print(f"No records found for protein: {protein_name}")
        return

    # Retrieve the first protein record
    protein_id = record["IdList"][0]
    handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text")
    protein_record = SeqIO.read(handle, "genbank")
    handle.close()

    # Print the information
    print(dir(protein_record))
    print("LOCUS       {}".format(protein_record.id))
    print("DEFINITION  {}".format(protein_record.description))
    print("ACCESSION   {}".format(protein_record.id))
    print("DBLINK      {}".format(protein_record.annotations.get("dblink", "")))
    print("SOURCE      {}".format(protein_record.annotations.get("source", "")))

    # Additional information can be printed based on your requirements
    # For example, you can print features, references, etc.

    # Print the sequence
    print("ORIGIN")
    sequence = str(protein_record.seq)
    for i in range(0, len(sequence), 60):
        print(" " + sequence[i:i + 60])

# Example usage with the provided protein name
protein_name = "integron [Burkholderia pseudomallei]"
fetch_protein_info(protein_name)
