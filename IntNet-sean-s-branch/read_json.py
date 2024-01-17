import json
import pandas as pd


data_file=r"C:\Users\seanji\Desktop\IntNet-sean-s-branch\stored_data.json"
with open(data_file, 'r') as file:
    ds = json.load(file) 

df = pd.DataFrame(ds)
display(df)
#df = pd.DataFrame(ds,columns=['NCBI Accession Number','Definition','Source','Organism', 'Isolation Source', 'Coolection Date', 'Country', 'Taxon ID', 'BioProject', 'Biosample', 'Ref'  ])


#note: in some truncations, might be some repeated items, need to clean 
#note: i=1814 failed


#def unique(lists):
#    unique_list = pd.Series(lists).drop_duplicates().tolist()
#    for x in unique_list:
#        print(x)
#unique(lists)