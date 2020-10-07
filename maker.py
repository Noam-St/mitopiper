#%%
import pandas as pd
import sys, os
path = os.getcwd()
from Bio import SeqIO
from Bio import Entrez



if sys.platform != 'win32':
    slash = '/'
else:
    slash = '\\'
try: org_df = pd.read_csv(f'{path}{slash}final.csv', index_col = 'organism')
except FileNotFoundError: raise FileNotFoundError('final.csv is required in the same folder for the program to operate.')

def file_reader(filename):
  """
  Read files in the format: X=Y to a dictionary

  Parameters
  ----------
  filepath : str
    path to txt file in the expected format
  
  Returns
  -------
  format_dict : dict
    Dictionary with the left side of the equal sign as key and right side as value
  """
  format_dict = {}
  with open(f'{path}{slash}{filename}.txt', 'r') as file:
    for line in file:
      line = line.rstrip('\n')
      line = line.split('=')
      format_dict[line[0]] = line[1]
    return format_dict
def file_creator(filename, format_dict):
  """
  Create text files in the format: X=Y, place file in cwd

  Parameters
  ----------
  filename : str
    Name of file to be created
  format_dict : dict
    Dictionary by-which to write the file in this format: key=value
  """
  with open(f'{path}{slash}{filename}.txt','w') as file:
    for k,v in format_dict.items():
      file.write(f'{k}={v}\n')
  return

def settings():
  """
  """
  finished = False
  while not finished:
    print("Welcome to the settings menu\nTo change ")
    break

def pipeline():
  """
  """
  pass
def dir_checker(dirname):
  """
  """
  if not os.path.isdir(f'{path}{slash}{dirname}'):
    os.mkdir(f'{path}{slash}{dirname}')


if not os.path.exists(f'{path}{slash}credentials.txt'):
  EMAIL = input('Please enter your NCBI account email:')
  API = input('Please enter your NCBI API (can be found in NCBI account page)')
  file_creator('credentials', {'Email':EMAIL, 'API':API})

creds = file_reader('credentials')

Entrez.email = creds['Email']
Entrez.api_key = creds['API']


default_pipe = {'wget' : 1, 'fastqc1' : 1}
default_settings = {}
dir_checker('fastq')
dir_checker('mtDNA_ref')
dir_checker('error_output')
dir_checker('std_output')
 
if not os.path.exists(f'{path}{slash}pipeline.txt'):
  file_creator('pipeline', default_pipe)
if not os.path.exists(f'{path}{slash}settings.txt'):
  file_creator('settings', default_settings)



def ref_maker(org_name):
  """
  """
  ID = org_df.loc[org_name, 'RefSeq']
  filename = f'{path}{slash}{ID}.gbk'
  if not os.path.isfile(filename):
    net_handle = Entrez.efetch(db = 'nucleotide', id = ID, rettype = 'gb', retmode = 'text')
    out_handle = open(filename, 'w')
    out_handle.write(net_handle.read())
    net_handle.close()
    out_handle.close()
  record = SeqIO.read(filename, 'genbank')
  seq = str(record.seq)
  with open(f'{path}{slash}mtDNA_ref{slash}{org_name}.fa', 'w') as ref:
    ref.write('>MT\n')
    ref.write(seq)

#%%

tsv_name = input(""" Please enter TSV file name (without extension).
IMPORTANT: The tsv must be in the same folder as this code, and in the following format:\n
study_accession\tsample_accession\trun_accession\ttax_id\tscientific_name\tfastq_ftp\tsubmitted_ftp\tsra_ftp\n
The format is available for download from the EBI website\n""")
tsv_path = path + slash + tsv_name + '.tsv'
tsv_df = pd.read_csv(tsv_path, delimiter = '\t')
org_name = tsv_df['scientific_name'][0]

change = int(input('Enter 1 to change pipeline\n Enter 2 to change settings\n Enter 0 to continue'))
work_on = input('Choose which lines of the tsv to process with the current pipeline.\nenter a single row number (i.e 2), a range (i.e 2-3) or 0 for all rows')

if not os.path.exists(f'{path}{slash}mtDNA_ref{slash}{org_name}'):
  ref_maker(org_name)

#%%
if change == 2:
  settings()
elif change == 1:
  pipeline()

cur_pipe = file_reader('pipeline')
cur_settings = file_reader('settings')

def csh_maker(cur_pipe, cur_settings, work_on, tsv_df):
  """
  """
  with open('run.csh', 'w') as file:
    file.write('#!/bin/csh\n')
    file.write(f'#$ -o {path}{slash}std_output\n')
    file.write(f'#$ -e {path}{slash}error_output\n')
    file.write(f'#$ -cwd\n')
    file.write(f'#$ -N {tsv_df.study_accession[0]}\n')
    file.write(f'#$ -q bioinfo.q\n')
    if '-' in work_on:
      work_on = [int(i) for i in work_on.split('-')]
      start = work_on[0]
      end = work_on[1]
    else:
      work_on = int(work_on)
      start = work_on
      if work_on == 0: end = len(tsv_df)
      else: end = start+1
    for i,row in tsv_df.iloc[start:end,:].iterrows():
      if len(row.fastq_ftp.split(';')) > 1: PAIRED = True
      else : PAIRED = False
      if PAIRED:
        ftps = row.fastq_ftp.split(';')
        for ind, ftp in enumerate(ftps):
          if cur_pipe['wget']:
            file.write(f'echo \"Running wget on {row.sample_accession}\"\n')
            file.write(f'wget -nH -nd -P {path}{slash}fastq{slash} -O {row.sample_accession}_{ind + 1}.fastq.gz {ftp}\n')
          if cur_pipe['fastqc1']:
            file.write(f'echo \"Running fastqc on {row.sample_accession}\"\n')
            file.write(f'fastqc {path}{slash}fastq{slash}{row.sample_accession}_{ind + 1}.fastq.gz\n')
      else:
        ftp = row.fastq_ftp
        if cur_pipe['wget']:
          file.write(f'echo \"Running wget on {row.sample_accession}\"\n')
          file.write(f'wget -nH -nd -P {path}{slash}fastq{slash} -O {row.sample_accession}.fastq.gz {ftp}\n')
        if cur_pipe['fastqc1']:
          file.write(f'echo \"Running fastqc on {row.sample_accession}\"\n')
          file.write(f'fastqc {path}{slash}fastq{slash}{row.sample_accession}.fastq.gz\n')
#%%
csh_maker(cur_pipe, cur_settings, work_on, tsv_df)