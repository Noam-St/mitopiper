import pandas as pd
import numpy as np
import os, sys
PATH = os.getcwd()
from Bio import SeqIO


if sys.platform != 'win32':
    SLASH = '/'
else:
    SLASH = '\\'
try:
  org_df = pd.read_csv(f'{PATH}{SLASH}src{SLASH}final.csv')
  org_df.organism = org_df.organism.str.replace(' ','_')
  org_df = org_df.set_index('organism')
except FileNotFoundError: raise FileNotFoundError('final.csv is required in the src folder for the program to operate.')






def main(ID, org):
  """
  """
  gb_path = f'{PATH}{SLASH}src{SLASH}genbank{SLASH}{ID}.gbk'
  record = SeqIO.read(gb_path, 'genbank')
  features = record.features
  gff_path = f'{PATH}{SLASH}src{SLASH}mtDNA_ref{SLASH}{org}{SLASH}{org}.gtf'
  wanted_types = ['gene', 'tRNA', 'rRNA', 'CDS']
  with open(gff_path, 'w') as gff:
    gff.write(f'#!RefSeq-ID {ID}\n#!Organism {org}\n#!chrM-only\n')
    for feat in features:
      cur_line = 'chrM\tfrom_gb\t'
      if feat.type in wanted_types:
        for location in feat.location.parts:
          start = location.start
          end = location.end
          strand = '+' if location.strand == 1 else '-'
        cur_line += feat.type + '\t' #TODO finish this quickly

