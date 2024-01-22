import pandas as pd
import numpy as np
import json 

from Bio import Entrez 

import xml.etree.ElementTree as ET 

  

# Always provide your email when using NCBI Entrez 

Entrez.email = "XIANGYU270208@GMAIL.COM"  # Replace with your email 

  

def fetch_biosample_details(biosample_id): 

    # Fetch the biosample record using efetch 

    handle = Entrez.efetch(db="biosample", id=biosample_id, retmode="xml") 

    record = handle.read() 

  

    # Parse the XML record 

    root = ET.fromstring(record) 

  

    # Extracting information 

    biosample_data = {} 

    for biosample in root.findall('BioSample'): 

        # Get various attributes 

        biosample_data['BioSampleID'] = biosample.get('id') 

        biosample_data['Name'] = biosample.find('Description/Title').text 

        # Add other attributes as needed 

  

        # Extract attributes from the Attributes section 

        for attribute in biosample.find('Attributes').findall('Attribute'): 

            attr_name = attribute.get('attribute_name') 

            biosample_data[attr_name] = attribute.text 

  

    return biosample_data 


file_path =  r"C:\Users\seanji\Desktop\IntNet-sean-s-branch\Biosamples.csv"
df = pd.read_csv(file_path, dtype='object')



data_file=r"C:\Users\seanji\Desktop\IntNet-sean-s-branch\biosample_last.json"
data=[]
with open(data_file, 'w') as file:
        json.dump(data, file)


for i in range(6600,len(df)):
    print(i)
    biosample_id = df.iloc[i]['BioSample']
    if pd.isnull(biosample_id):
        details = {}
    else:
      details =  fetch_biosample_details(biosample_id)
    
    with open(data_file, 'r') as file:
     existing_data = json.load(file)
    

    if len(existing_data)==i-6600:
     existing_data.append(details)

    with open(data_file, 'w') as file:
        json.dump(existing_data, file)


