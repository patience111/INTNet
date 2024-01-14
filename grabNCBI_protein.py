from Bio import Entrez, SeqIO
import pandas as pd
import json 


import fasta_to_txt
gene_id_list=fasta_to_txt.data_clean

def fetch_protein_info(gene_id):
    Entrez.email = "xiangyu270208@gmail.com"  # Set your email for NCBI API usage

    # Retrieve the gene record
    handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")
    gene_record = SeqIO.read(handle, "genbank")
    handle.close()
    

    list_of_titles=[]
    for i in range(len(gene_record.annotations["references"])):
        list_of_titles.append(gene_record.annotations["references"][i].title)
    

    if 'organism' not in list(gene_record.features[0].qualifiers.keys()):
        gene_record.features[0].qualifiers['organism'] = '0'
    if 'isolation_source' not in list(gene_record.features[0].qualifiers.keys()):
        gene_record.features[0].qualifiers['isolation_source'] = '0'
    if 'collection_date' not in list(gene_record.features[0].qualifiers.keys()):
        gene_record.features[0].qualifiers['collection_date'] = '0'
    if 'country' not in list(gene_record.features[0].qualifiers.keys()):
        gene_record.features[0].qualifiers['country'] = '0'
    if 'db_xref' not in list(gene_record.features[0].qualifiers.keys()):
        gene_record.features[0].qualifiers['db_xref'] = '0'


    data=[gene_record.id,gene_record.description,gene_record.annotations["source"],
        gene_record.features[0].qualifiers["organism"][0],gene_record.features[0].qualifiers["isolation_source"][0],
        gene_record.features[0].qualifiers["collection_date"][0],gene_record.features[0].qualifiers["country"][0],
        gene_record.features[0].qualifiers["db_xref"][0]]
    

    for i in range(len(gene_record.dbxrefs)):
        data.append(gene_record.dbxrefs[i])
    for j in range(len(list_of_titles)):
        data.append(list_of_titles[j])


    
    ## Print the information
    #print("NCBI Accession Number  {}".format(gene_record.id))
    #print("Definition             {}".format(gene_record.description))
    #print("DBLink                 {}".format(gene_record.dbxrefs))
    #print("Source                 {}".format(gene_record.annotations["source"]))
    #print("References             {}".format(list_of_titles))
    #print("Organism               {}".format(gene_record.features[0].qualifiers["organism"][0]))
    #print("Isolation Source       {}".format(gene_record.features[0].qualifiers["isolation_source"][0]))
    #print("Collection Date        {}".format(gene_record.features[0].qualifiers["collection_date"][0]))
    #print("Country                {}".format(gene_record.features[0].qualifiers["country"][0]))
    #print("Taxon ID               {}".format(gene_record.features[0].qualifiers["db_xref"][0]))
    
 
    
    return data
    
    # Additional information can be printed based on your requirements
    # For example, you can print features, references, etc.
         

# Example usage with the provided gene id
data_file="stored_data.json"
for i in range(3372,len(gene_id_list)):
  print(i)
  gene_id = gene_id_list[i]
  print(gene_id)
  data = fetch_protein_info(gene_id)

  try:
        with open(data_file, 'r') as file:
            existing_data = json.load(file)
  except FileNotFoundError:
        existing_data = []
  
  existing_data.append(data)

  with open(data_file, 'w') as file:
        json.dump(existing_data, file)

with open(data_file, 'r') as file:
    ds = json.load(file) 

df = pd.DataFrame(ds,columns=['NCBI Accession Number','Definition','Source','Organism', 'Isolation Source', 'Coolection Date', 'Country', 'Taxon ID', 'BioProject', 'Biosample', 'Ref'  ])


#note: in some truncations, might be some repeated items, need to clean 
#note: i=1814 failed


#def unique(lists):
#    unique_list = pd.Series(lists).drop_duplicates().tolist()
#    for x in unique_list:
#        print(x)
#unique(lists)