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
from ast import literal_eval
import configparser as cp
import json

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
  with open(f'{PATH}{slash}{filename}.json', 'r') as file:
    format_dict = json.load(file)
  return format_dict

def file_creator(filename, format_dict):
  """
  Create json files of the format_dict, place file in cwd

  Parameters
  ----------
  filename : str
    Name of file to be created
  format_dict : dict
    Dictionary by-which to write the file in this format: key=value
  """
  filepath = f'{PATH}{slash}{filename}.json'
  with open(filepath, 'w') as file:
    json.dump(format_dict, sort_keys = True, fp = file, indent = 1)

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


if not os.path.exists(f'{PATH}{slash}credentials.json'):
  EMAIL = input('Please enter your NCBI account email:')
  API = input('Please enter your NCBI API (can be found in NCBI account page)')
  file_creator('credentials', {'Email':EMAIL, 'API':API})

creds = file_reader('credentials')

default_pipe =\
{
'wget' : {
  'key':'wget',
  'description':'Download command',
  'switch':True,
  'path':''
},
'trim' : {
  'key':'trim_galore',
  'description':'Trims adapters',
  'switch':True,
  'path':r'/gpfs0/biores/apps/TrimGalore/TrimGalore-0.4.5/trim_galore'
},
'alignment' : {
  'key':'bwa index',
  'description':'Align reads to reference',
  'switch':True,
  'path':r'/gpfs0/biores/apps/bwa/bwa-v0.7.16'
},
'view':{
  'key':'samtools view',
  'description':'Convert sam to bam',
  'switch':True,
  'path':r'/gpfs0/biores/apps/samtools/samtools-1.3/bin/samtools view'
},
'index':{
  'key':'samtools index',
  'description':'Index bam for fast access via samtools',
  'switch':True,
  'path':r'/gpfs0/biores/apps/samtools/samtools-1.3/bin/samtools index'
},
'mpileup':{
  'key':'samtools mpileup',
  'description':'Create a readable mpileup format',
  'switch':True,
  'path':r'/gpfs0/biores/apps/samtools/samtools-1.3/bin/samtools mpileup'
},
'fastqc':{
  'key':'fastqc',
  'description':'Create a readable quality report in html format',
  'switch':True,
  'path':r'/gpfs0/biores/apps/FastQC_v0.11.7/fastqc'
}
}

default_settings =\
{
  'trim':{
    'quality':{'description':'Trim low-quality ends from reads', 'key':'-q', 'value':20, 'switch':True}
  },
  'alignment':{
    'seed':{'description':'Minimum seed length', 'key':'-k', 'value':19, 'switch':True}
  },
  'view':{
    'header':{'description':'Include the header in the output', 'key':'-h', 'value':'', 'switch':True},
    'bam':{'description':'Output in BAM format', 'key':'-b', 'value':'', 'switch':1}
  },
  'index':{},
  'mpileup':{
    'quality':{'description':'Minimum base quality for a base to be considered', 'key':'-Q', 'value':30, 'switch':True}
  },
  'fastqc':{},

}

defaults = (default_pipe, default_settings)
dir_checker('fastq')
dir_checker('mtDNA_ref')
dir_checker('error_output')
dir_checker('std_output')
dir_checker('Genbank_DB')
dir_checker('fastqc_pre')
dir_checker('fastqc_post')
dir_checker('post_trim')
dir_checker('alignment')
dir_checker('samples')
dir_checker('pileup')

configs = ['pipeline', 'settings', 'paths']
for i,config in enumerate(configs):
  if not os.path.exists(f'{PATH}{slash}{config}.json'):
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




def change_settings():
  """
  """
  settings = file_reader('settings')
  pipeline = file_reader('pipeline')
  options = list(settings.keys())
  print(f"Welcome to the settings menu\n")
  while True:
    [print(f'{i}\t{v}\n') for i,v in enumerate(options) if pipeline[v]['switch']]
    print('To change settings for a specific step, enter the number next to it\n')
    print('Type exit to exit\n')
    choice = input()
    if choice == 'exit': break
    choice = int(choice)
    chosen_stage = options[choice]
    configs = list(settings[chosen_stage].keys())
    print(f'Index\tName\tCurrent_setting\tDescription\n')
    [
      print(f'{i}\t{v}\t{s}\t{d}\n') for i,v,d,s in zip(range(len(configs)),configs,
      [i['description'] for i in settings[chosen_stage].values()],[i['value'] for i in settings[chosen_stage].values()])
    ]
    if settings[chosen_stage] == {}:
      print(f'It seems no settings were added to {chosen_stage} yet\nto add some return first to the main menu by typing exit\n')
    else:
      print(f'What would you like to change in {chosen_stage}?\nEnter the number next to the setting\n')
      change = int(input())
      chosen_setting = configs[change]
      setting_value = settings[chosen_stage][chosen_setting]['value']
      setting_switch = settings[chosen_stage][chosen_setting]['switch']
      if settings[chosen_stage][chosen_setting]['value'] == '':
        settings[chosen_stage][chosen_setting]['switch'] = not setting_switch
        print(f'\nSetting switched to {not bool(setting_switch)}\n')

      else:
        print(f'You chose: {chosen_setting}\nPlease choose a new value, current value {setting_value}:')
        new_val = input()
        settings[chosen_stage][chosen_setting]['value'] = new_val
        print('\nSetting changed successfuly!\n')
    
  file_creator('settings', settings)

def change_pipeline():
  """
  """
  running = True
  pipeline = file_reader('pipeline')
  while running:
    print("Welcome to the pipeline customizer\nThe current pipeline is:\n")
    [print(f'{k}\n') for i,k in zip(range(len(pipeline)),pipeline.keys()) if pipeline[k]['switch']]
    print('To change the pipeline enter 1\nTo exit pipeline enter 0\n')
    running = int(input())
    if running:
      pipeline_changes = {i:k for i,k in enumerate(pipeline.keys())}
      print('To remove or add a feature, simply enter its number\n')
      [print(f'{i} {k} - current setting: {pipeline[k]["switch"]}\n') for i,k in zip(pipeline_changes.keys(), pipeline_changes.values())]
      to_change = pipeline_changes[int(input())]
      pipeline[to_change]['switch'] = not pipeline[to_change]['switch']
  file_creator('pipeline', pipeline)

def change_paths():
  """
  """
  pass

def add_new_setting():
  """
  """
  settings = file_reader('settings')
  setting_keys = list(settings.keys())
  [print(f'{i}\t{k}\n') for i,k in enumerate(setting_keys)]
  print('\nPlease type the pipeline step to add a setting to\n')
  step = int(input())
  chosen_step = setting_keys[step]
  print(f'You chose {chosen_step}\nPlease enter the setting name\n')
  new_name = input()
  settings[chosen_step][new_name] = {}
  print(f'Please enter a description for {new_name}\n')
  settings[chosen_step][new_name]['description'] = input()
  print(f'Please enter a key for {new_name} (a key is the name that needs to be added to the command, must be identical to the one defined in {chosen_step} docs\n')
  key = input()
  try: int(key)
  except ValueError: pass
  settings[chosen_step][new_name]['key'] = key
  print(f'Please enter a value for {new_name} (if the setting is an on/off switch, enter nothing!)\n')
  value = input()
  settings[chosen_step][new_name]['value']  = value
  settings[chosen_step][new_name]['switch'] = 1
  file_creator('settings', settings) 
  print(settings)
  print('Setting added successfuly!\n')

def add_new_pipeline():
  """
  """
  pipeline = file_reader('pipeline')
  print('To begin defining a new pipeline step, first enter its step\nIMPORTANT: the name must be different from the current pipeline steps\n')
  name = input()
  pipeline[name] = {k:'' for k,v in pipeline['trim'].items() if k != 'switch'}
  for k in pipeline[name].keys():
    print(f'Please enter a {k}\n')
    v = input()
    pipeline[name][k] = v
  print('\nPipeline step added successfuly!\n')

def change_menu():
  """
  """
  choice = True
  print('####################\n\nWelcome to the configuration menu\n')
  menu = [lambda :0,change_settings,add_new_setting, change_pipeline, add_new_pipeline, change_paths]
  while choice:
    [print(f"Enter {i+1}\tto {x.__name__}\n") for i,x in enumerate(menu[1:])]
    print("Enter 0\t to exit\n")
    try:
      choice = int(input())
    except ValueError:
      print('Please only choose one of the listed options!\n')
    try: menu[choice]()
    except IndexError:
      print("Please only choose one of the listed options!\n")




tsv_name = input('Please enter TSV file name (without extension).\nIMPORTANT: The tsv must be in a folder named samples (the folder was created if it doesnt exist)\nIMPORTANT: the tsv msut be in the following format:\nstudy_accession\tsample_accession\trun_accession\ttax_id\tscientific_name\tfastq_ftp\tsubmitted_ftp\tsra_ftp\nThe format is available for download from the EBI website\n')
tsv_path = f'{PATH}{slash}samples{slash}{tsv_name}.tsv'
tsv_df = pd.read_csv(tsv_path, delimiter = '\t')
tsv_df.scientific_name = tsv_df.scientific_name.str.replace(' ','_')
org_name = tsv_df['scientific_name'][0]

change = int(input('Please enter 1 to enter the configuration menu.\nEnter 0 to continue\n'))
if change:
  change_menu()
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
      sample = row.run_accession
      file.write(f'echo \"Current organism: {row.scientific_name}\"\n\n')
      for ind, ftp in enumerate(ftps):
        if PAIRED:
          ind+=1
          ind = '_' + str(ind)
        else: ind = ''
        if cur_pipe['wget']:
          file.write(f'echo \"Downloading sample: {sample}{ind}\"\n\n')
          file.write(f'wget -nH -nd -O {PATH}{slash}fastq{slash}{sample}{ind}.fastq.gz {ftp}\n')
        if cur_pipe['fastqc_pre']:
          file.write(f'echo \"Checking quality of {sample}{ind}\"\n\n')
          file.write(f'{cur_path["fastqc"]} --outdir {PATH}{slash}fastqc_pre {PATH}{slash}fastq{slash}{sample}{ind}.fastq.gz\n')
    for _,row in tsv_df.iloc[start:end,:].iterrows():
      org_name = row.scientific_name
      sample = row.run_accession
      fasta = f'{PATH}{slash}mtDNA_ref{slash}{org_name}{slash}{org_name}.fasta'
      file.write(f'echo \"Current organism: {row.scientific_name}\"\n\n')
      if cur_pipe['trim']:
        file.write(f'echo \"Trimming {sample}\"\n\n')
        if PAIRED:
          file.write(f'trim_galore -o {PATH}{slash}post_trim --paired {PATH}{slash}fastq{slash}{sample}_1.fastq.gz {PATH}{slash}fastq{slash}{sample}_2.fastq.gz\n')
          file.write(f'{cur_path["fastqc"]} --outdir {PATH}{slash}fastqc_post {PATH}{slash}post_trim{slash}{sample}_1_val_1.fq.gz\n')
          file.write(f'{cur_path["fastqc"]} --outdir {PATH}{slash}fastqc_post {PATH}{slash}post_trim{slash}{sample}_2_val_2.fq.gz\n')

        else:
          file.write(f'trim_galore -o {PATH}{slash}post_trim {PATH}{slash}fastq{slash}{sample}.fastq.gz\n')
          file.write(f'{cur_path["fastqc"]} --outdir {PATH}{slash}fastqc_post {PATH}{slash}post_trim{slash}{sample}_trimmed.fq.gz\n')
      if not os.path.exists(fasta):
        fasta, add_to_csh = ref_main(org_name)
        file.write(add_to_csh)
      if cur_pipe['align+index']:
        file.write(f'echo \"Alignment stage start for {sample}\"\n\n')
        if PAIRED:
          file.write(f'bwa mem {fasta} {PATH}{slash}post_trim{slash}{sample}_1_val_1.fq.gz {PATH}{slash}post_trim{slash}{sample}_2_val_2.fq.gz > {PATH}{slash}alignment{slash}{sample}.sam\n')
        else:
          file.write(f'bwa mem {fasta} {PATH}{slash}post_trim{slash}{sample}_trimmed.fq.gz > {PATH}{slash}alignment{slash}{sample}.sam\n')
        file.write(f'echo \"Mpileup stage start for {sample}\"\n\n')
        file.write(f'samtools view -h -b {PATH}{slash}alignment{slash}{sample}.sam > {PATH}{slash}alignment{slash}{sample}.bam\n')
        file.write(f'samtools sort {PATH}{slash}alignment{slash}{sample}.bam -o {PATH}{slash}alignment{slash}{sample}.bam\n')
        file.write(f'samtools index {PATH}{slash}alignment{slash}{sample}.bam\n')
      if cur_pipe['mpileup']:
        file.write(f'samtools mpileup -f {fasta} -Q 30 -q 0 {PATH}{slash}alignment{slash}{sample}.bam > {PATH}{slash}pileup{slash}{sample}.pileup\n')
      #CONTINUE FROM HERE add some settings.



    

csh_maker(cur_pipe, cur_settings, work_on, tsv_df)

# %%

# %%
