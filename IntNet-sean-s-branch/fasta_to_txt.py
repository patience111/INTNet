from Bio import SeqIO

def read_fasta_file(file_path):
    """
    Reads a FASTA file and returns a list of sequences.
    """
    sequences = []

    with open(file_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):

            sequences.append(record)

    return sequences

# Specify the path to your FASTA file
fasta_file_path = r"C:\Users\seanji\Desktop\IntNet-sean-s-branch\TFILE1.fasta"

# Read the FASTA file
fasta_sequences = read_fasta_file(fasta_file_path)

# Print the sequences
data = []
for sequence in fasta_sequences:
    data.append(sequence.id)
data_clean = list(range(len(data)))

print(data[5000][:3]=='NC_' or 'NZ_')
data_clean[5000]=data[5000].split("_")[0]+ data[5000].split("_")[1]
print(data_clean[5000])
print(data[20330].split("|"))
print(data[20330][:10]=='Smyshlyaev' or data[20330][:5]=='Zhang' or len(data[20330].split("|"))==2 or len(data[20330].split("|"))==3 )


# Data Clean with selection rule specified
for i in range(len(data)):
    
    if data[i][:3]=='NC_' or  data[i][:3]=='NZ_':
        data_clean[i]=data[i].split("_")[0]+ "_"+ data[i].split("_")[1] # Starts with NZ/NC: require format "NC_015945.1"
    if data[i][:3]=='GCF':
        data_clean[i]=data[i].split("_")[2]+ "_"+ data[i].split("_")[3]  # Starts with GCF: require convert "GCF_020309905.1_NZ_CP084323.1_4621" to "NZ_CP084323.1"
    if data[i][:3]=='gnl':
        data_clean[i]=data[i].split("|")[1] # Starts with gnl: require convert "gnl|AB257723|idbi84_uncultured_bacterium|intI|p100|q100|s100" to "AB257723"
    if data[i][:10]=='Smyshlyaev' or data[i][:5]=='Zhang' or len(data[i].split("|"))==2 or len(data[i].split("|"))==3:
        data_clean[i]='0'   # The rest : to be solved

for i in range(len(data)):
    if data_clean[i]==i:   # For the rest, all of them have a clear accession number, simply fetch them
        data_clean[i]=data[i].split("_")[:1][0]

print(data_clean)

# Open the file in write mode
output_file_path = r"C:\Users\seanji\Desktop\IntNet-sean-s-branch\output3.txt"
with open(output_file_path, 'w') as output_file:
    # Write each line of data to the file
    for line in data_clean:
        output_file.write(str(line) + '\n')