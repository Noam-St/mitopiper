__author__ = 'Noam Shtolz'
"""
#########################################################################

This program creates a custom csh file designed to be executed with qsub.
The program is designed to be configurable, the user can customize the processing pipeline,
the settings for each of the pipeline programs used, and the paths of the programs utallized.

#########################################################################
"""
#%%
import pandas as pd
import sys, os
PATH = os.getcwd()
from Bio import SeqIO
from Bio import Entrez
from textwrap import wrap


if sys.platform != 'win32':
    slash = '/'
else:
    slash = '\\'
try:
  org_df = pd.read_csv(f'{PATH}{slash}final.csv')
  org_df.organism = org_df.organism.str.replace(' ','_')
  org_df = org_df.set_index('organism')
except FileNotFoundError: raise FileNotFoundError('final.csv is required in the same folder for the program to operate.')
def file_reader(filename):
  """
  Read files in the format: X=Y to a dictionary

  Parameters
  ----------
  filePATH : str
    PATH to txt file in the expected format
  
  Returns
  -------
  format_dict : dict
    Dictionary with the left side of the equal sign as key and right side as value
  """
  format_dict = {}
  with open(f'{PATH}{slash}{filename}.txt', 'r') as file:
    for line in file:
      line = line.rstrip('\n')
      line = line.split('=')
      try: line[1] = int(line[1])
      except ValueError: pass
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
  with open(f'{PATH}{slash}{filename}.txt','w') as file:
    for k,v in format_dict.items():
      file.write(f'{k}={v}\n')
  return

def dir_checker(dirname):
  """
  check if a directory exists, if it does not, create it.

  Parameters
  ----------
  dirname : str
    Name of the directory to be checked
  """
  if not os.path.isdir(f'{PATH}{slash}{dirname}'):
    os.mkdir(f'{PATH}{slash}{dirname}')


if not os.path.exists(f'{PATH}{slash}credentials.txt'):
  EMAIL = input('Please enter your NCBI account email:')
  API = input('Please enter your NCBI API (can be found in NCBI account page)')
  file_creator('credentials', {'Email':EMAIL, 'API':API})

creds = file_reader('credentials')




default_pipe =\
{
'wget' : 1,
'fastqc_pre' : 1,
'trim' : 1,
'align+index' : 1,
'mpileup' : 1,
'depth' : 0
}
default_settings =\
{
  'trim':{
    'quality':20,
    'adapter':''
  }
}
default_paths =\
{
  'fastqc' : r'/gpfs0/biores/apps/FastQC_v0.11.7/fastqc',
  'trim' : r'/gpfs0/biores/apps/TrimGalore/TrimGalore-0.4.5/trim_galore',
  'bwa' : r'/gpfs0/biores/apps/bwa/bwa-v0.7.16',
  'samtools': r'/gpfs0/biores/apps/samtools/samtools-1.3/bin/samtools'
}
defaults = (default_pipe, default_settings, default_paths)
dir_checker('fastq')
dir_checker('mtDNA_ref')
dir_checker('error_output')
dir_checker('std_output')
dir_checker('Genbank_DB')
dir_checker('fastqc_pre')
dir_checker('fastqc_post')
dir_checker('post_trim')
dir_checker('alignment')

configs = ['pipeline', 'settings', 'paths']
for i,config in enumerate(configs):
  if not os.path.exists(f'{PATH}{slash}{config}.txt'):
    file_creator(config, defaults[i])

def fasta_wrap(seq):
  """
  """
  return '\n'.join(wrap(seq, 60))
def ref_maker(org_name):
  """
  Create a fasta file for indexing purposes.

  Parameters
  ----------
  org_name : str
    Name of the organism to name the file.
  
  Returns
  -------
    Absolute path to the fasta file created.
  """
  ID = org_df.loc[org_name, 'RefSeq']
  filename = f'{PATH}{slash}genbank_DB{slash}{ID}.gbk'
  if not os.path.isfile(filename):
    Entrez.email = creds['Email']
    Entrez.api_key = creds['API']
    net_handle = Entrez.efetch(db = 'nucleotide', id = ID, rettype = 'gb', retmode = 'text')
    out_handle = open(filename, 'w')
    out_handle.write(net_handle.read())
    net_handle.close()
    out_handle.close()
  record = SeqIO.read(filename, 'genbank')
  seq = str(record.seq)
  fasta = f'{PATH}{slash}mtDNA_ref{slash}{org_name}{slash}{org_name}.fasta'
  with open(fasta, 'w') as ref:
    ref.write('>MT\n')
    ref.write(fasta_wrap(seq))
  return fasta

def ref_main(org_name):
  """
  """
  choice = True
  while choice:
    print('It seems an reference fasta does not exist in the default location\nTo create one automatically, press 1\nTo enter an index path, press 2\n')
    choice = input()
    try: choice = int(choice)
    except ValueError: print('Please enter 1 or 2 only')
    if choice == 1:
      dir_checker(f'mtDNA_ref{slash}{org_name}')
      fasta = ref_maker(org_name)
      add_to_csh = f'samtools faidx {fasta}\nbwa index {fasta}\n'
      return fasta, add_to_csh
    elif choice == 2:
      print('Please enter the full path to the reference fasta\n')
      fasta = input()
      return fasta, ''
    else:
      print('Please enter 1 or 2 only')




def settings():
  """
  """
  finished = False
  settings = file_reader('settings')
  while not finished:
    print(f"Welcome to the settings menu\nThe current settings are:\n")
    break

def pipeline():
  """
  """
  running = True
  pipeline = file_reader('pipeline')
  while running:
    print("Welcome to the pipeline customizer\nThe current pipeline is:\n")
    [print(f'{k}\n') for i,k,v in zip(range(len(pipeline)),pipeline.keys(), pipeline.values()) if v != 0]
    print('To change the pipeline enter 1\nTo exit pipeline enter 0\n')
    running = int(input())
    if running:
      pipeline_changes = {i:k for i,k in enumerate(pipeline.keys())}
      print('To remove or add a feature, simply enter its number\n')
      [print(f'{i} {k} - current setting: {v}\n') for i,k,v in zip(pipeline_changes.keys(), pipeline_changes.values(), pipeline.values())]
      to_change = pipeline_changes[int(input())]
      if pipeline[to_change]:
        pipeline[to_change] = 0
      else:
        pipeline[to_change] = 1
  print(pipeline)
  file_creator('pipeline', pipeline)
  return

def paths():
  """
  """
  pass

def menu():
  """
  """
  choice = True
  while choice:
    print('#####################\n\nWelcome to the reconfiguration menu\nTo change the pipeline enter 1\nTo change the settings enter 2\nTo exit enter 0\n\n#####################\n')
    try:
      choice = int(input())
    except ValueError:
      print('Please enter either 1, 2, 3 or 0!')
    if choice == 1:
      pipeline()
    elif choice == 2:
      settings()
    elif choice == 3:
      paths()
    elif not choice:
      pass
    else: print("Please enter either 1, 2, 3 or 0!")
  return




tsv_name = input('Please enter TSV file name (without extension).\nIMPORTANT: The tsv must be in the same folder as this code, and in the following format:\nstudy_accession\tsample_accession\trun_accession\ttax_id\tscientific_name\tfastq_ftp\tsubmitted_ftp\tsra_ftp\nThe format is available for download from the EBI website\n')
tsv_path = PATH + slash + tsv_name + '.tsv'
tsv_df = pd.read_csv(tsv_path, delimiter = '\t')
tsv_df.scientific_name = tsv_df.scientific_name.str.replace(' ','_')
org_name = tsv_df['scientific_name'][0]

change = int(input('Please enter 1 to change the settings or pipeline.\nEnter 0 to continue\n'))
if change:
  menu()
print(f'Choose which lines of the tsv to process with the current pipeline.\nEnter a single row number (i.e 2) or a range (i.e 2-3)\nEnter all to work on all rows\n')
[print(f'{ind}\t{samp}') for ind,samp in zip(tsv_df.index,tsv_df.sample_accession)]
work_on = input()




cur_pipe = file_reader('pipeline')
cur_settings = file_reader('settings')
cur_path = file_reader('paths')

def csh_maker(cur_pipe, cur_settings, work_on, tsv_df):
  """
  """
  with open('run.csh', 'w') as file:
    file.write('#!/bin/csh\n')
    file.write(f'#$ -o {PATH}{slash}std_output\n')
    file.write(f'#$ -e {PATH}{slash}error_output\n')
    file.write(f'#$ -cwd\n')
    file.write(f'#$ -N {tsv_df.study_accession[0]}\n')
    file.write(f'#$ -q bioinfo.q\n')
  
    if '-' in work_on:
      work_on = [int(i) for i in work_on.split('-')]
      start = work_on[0]
      end = work_on[1]
    elif work_on == 'all':
      end = len(tsv_df)
      start = 0
    else:
      work_on = int(work_on)
      start = work_on
      end = start+1
    for _,row in tsv_df.iloc[start:end,:].iterrows():
      if len(row.fastq_ftp.split(';')) > 1: PAIRED = True
      else : PAIRED = False
      ftps = row.fastq_ftp.split(';')
      sample = row.sample_accession
      file.write(f'echo \"Current organism: {row.scientific_name}\"\n')
      for ind, ftp in enumerate(ftps):
        if PAIRED:
          ind+=1
          ind = '_' + str(ind)
        else: ind = ''
        if cur_pipe['wget']:
          file.write(f'echo \"Running wget on {sample}{ind}\"\n')
          file.write(f'wget -nH -nd -O {PATH}{slash}fastq{slash}{sample}{ind}.fastq.gz {ftp}\n')
        if cur_pipe['fastqc_pre']:
          file.write(f'echo \"Running fastqc on {sample}{ind}\"\n')
          file.write(f'{cur_path["fastqc"]} --outdir {PATH}{slash}fastqc_pre {PATH}{slash}fastq{slash}{sample}{ind}.fastq.gz\n')
    for _,row in tsv_df.iloc[start:end,:].iterrows():
      org_name = row.scientific_name
      sample = row.sample_accession
      fasta = f'{PATH}{slash}mtDNA_ref{slash}{org_name}{slash}{org_name}.fasta'
      file.write(f'echo \"Running TrimGalore on {sample}\"\n')
      if cur_pipe['trim']:
        if PAIRED:
          file.write(f'trim_galore -o {PATH}{slash}post_trim --paired --fastqc_args \"--outdir {PATH}{slash}fastqc_post\" {PATH}{slash}fastq{slash}{sample}_1.fastq.gz {PATH}{slash}fastq{slash}{sample}_2.fastq.gz\n')
        else:
          file.write(f'trim_galore -o {PATH}{slash}post_trim --fastqc_args \"--outdir {PATH}{slash}fastqc_post\" {PATH}{slash}fastq{slash}{sample}.fastq.gz\n')
      if not os.path.exists(fasta):
        fasta, add_to_csh = ref_main(org_name)
        file.write(add_to_csh)
      if cur_pipe['align+index']:
        if PAIRED:
          file.write(f'bwa mem {fasta} {PATH}{slash}post_trim{slash}{sample}_1_val_1.fq.gz {PATH}{slash}post_trim{slash}{sample}_2_val_2.fq.gz > {PATH}{slash}alignment{slash}{sample}.sam\n')
        else:
          file.write(f'bwa mem {fasta} {PATH}{slash}post_trim{slash}{sample}_trimmed.fq.gz < {PATH}{slash}alignment{slash}{sample}.sam\n')
        file.write(f'samtools view -h -b {PATH}{slash}alignment{slash}{sample}.sam > {PATH}{slash}alignment{slash}{sample}.bam\n')
        file.write(f'samtools index {PATH}{slash}alignment{slash}{sample}.bam\n')
      if cur_pipe['mpileup']:
        file.write(f'samtools mpileup -f {fasta} -Q 30 -q 0 -a > {PATH}{slash}pileup{slash}{sample}.pileup')
      #CONTINUE FROM HERE add some settings.



    

csh_maker(cur_pipe, cur_settings, work_on, tsv_df)

# %%

# %%
