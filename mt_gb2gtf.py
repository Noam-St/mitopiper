#%%
import pandas as pd
import numpy as np
import os, sys
PATH = os.getcwd()
from Bio import SeqIO
from Bio import Entrez

if sys.platform != 'win32':
    SLASH = '/'
else:
    SLASH = '\\'
try:
  org_df = pd.read_csv(f'{PATH}{SLASH}src{SLASH}final.csv')
  org_df.organism = org_df.organism.str.replace(' ','_')
  org_df = org_df.set_index('organism')
except FileNotFoundError: raise FileNotFoundError('final.csv is required in the src folder for the program to operate.')




#%%

def main(ID, org):
  """
  """
  gb_path = f'{PATH}{SLASH}src{SLASH}genbank{SLASH}{ID}.gbk'
  record = SeqIO.read(gb_path, 'genbank')
  features = record.features
  gff_path = f'{PATH}{SLASH}src{SLASH}mtDNA_ref{SLASH}{org}{SLASH}{org}.gtf'
  wanted_types = ['gene', 'tRNA', 'rRNA', 'CDS']
  double = ['transcript', 'exon']
  with open(gff_path, 'w') as gff:
    gff.write(f'#!RefSeq-ID {ID}\n#!Organism {org}\n#!chrM-only\n')
    for feat in features:
      cur_line = 'chrM\tfrom_gb\t'
      if feat.type in wanted_types:
        gene_biotype = 'Mt_' + feat.type
        feat_type = feat.type.replace('tRNA','transcript').replace('rRNA','transcript').replace('CDS','transcript')
        temp_frame = feature.qualifiers.get('codon_start')
        frame = temp_frame if temp_frame != '' else '.'      
        gene_id = feature.qualifiers.get('db_xref')[0].replace('GeneID:')
        gene_name = feature.qualifiers.get('gene')
        transcript_id = feature.qualifiers.get('db_xref')[0].replace('GeneID:')
        if feat_type == 'transcript':
          for i in double:
            if len(feat.location.parts) > 1:
              c=0
              for loc in feat.location.parts:
                start = loc.start
                end = loc.end  
                strand = '+' if loc.strand == 1 else '-'
                gff_line = f'chrM\tmt_gb2gtf\t{feat_type}\t{start}\t{end}\t.\t{strand}\t{frame}\tgene_id {gene_id}; gene_version \"1\"; gene_name {gene_name}; transcript_id {transcript_id + str(c)}; 


#%%
net_handle = Entrez.efetch(db ='nucleotide', id = 'NC_024511', rettype = 'gff3', retmode = 'text') 