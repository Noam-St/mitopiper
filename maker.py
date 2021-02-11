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
import sys, os, json
PATH = os.getcwd()
from Bio import SeqIO
from Bio import Entrez
from textwrap import wrap
from ast import literal_eval

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
  with open(f'{PATH}{slash}configs{slash}{filename}.json', 'r') as file:
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
  filepath = f'{PATH}{slash}configs{slash}{filename}.json'
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

def fasta_wrap(string, per = 60):
  """
  Wraps a string to a default of 60 characters per line

  Parameters
  ----------
  string : str
    The string to be wrapped
  
  int : per
    Amount of characters per line

  Returns
  -------
  str
    The wrapped string  
  """
  return '\n'.join(wrap(string, per))

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
    if creds['Email'] != '':
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
  Create a reference sequence .fasta file and return the path to it

    Parameters
    ----------
    org_name : str
      Name of the organism to make a reference file for

    Returns
    -------
    fasta : str
      Path to the created fasta file
    
    add_to_csh : str
      string to add to the csh (incase the user chooses to automatically create a fasta)
  """
  refs = file_reader('refpaths')
  try:
    refs[org_name]['refpath']
    refs[org_name]['gtfpath']
    return ''
  except KeyError:
    pass
  choice = True
  while choice:
    print('It seems a reference fasta does not exist in the default location\nTo create one automatically, press 1\nTo enter an index path, press 2\n')
    choice = input()
    try: choice = int(choice)
    except ValueError: print('Please enter 1 or 2 only')
    if choice == 1:
      dir_checker(f'mtDNA_ref{slash}{org_name}')
      fasta = ref_maker(org_name)
      add_to_csh = f'samtools faidx {fasta}\nbwa index {fasta}\n'
    elif choice == 2:
      print('Please enter the full path to the reference fasta\n')
      fasta = input()
      add_to_csh = ''
    else:
      print('Please enter 1 or 2 only')
    print('Please enter path to gtf file for this organism\nPress enter if none needed\n')
    gtf = input()
    refs[org_name] = {}
    refs[org_name]['refpath'] = fasta
    refs[org_name]['gtfpath'] = gtf
    file_creator('refpaths', refs)
    return add_to_csh

def change_general_settings():
  """
  Change general settings defined on a per-run basis.
  """
  general_settings = file_reader('general_settings')
  options = list(general_settings.keys())
  while True:
    [print(f'{i}\t{v}\n') for i,v in enumerate(options)]
    print('Which setting would you like to change?\n')
    print('Type exit to exit\n')
    choice = input()
    if choice == 'exit': break
    try:
      choice = int(choice)
      chosen = general_settings[options[choice]]
    except (ValueError, IndexError):
      print('Please enter either a number or exit!\n')
      continue
    print(f"Please enter a new value for {options[choice]}\nCurrent value: {chosen[0]}\nDescription: {chosen[1]}")
    new_value = input()
    if options[choice] == 'free_form':
      if new_value == '':
        chosen[0] = ''
      else:
        chosen[0] += new_value + '\n'
    else:
      chosen[0] = new_value
    general_settings[options[choice]] = chosen

  file_creator(filename = 'general_settings', format_dict = general_settings)

def change_settings():
  """
  Settings change menu, in this menu you can reconfigure any of the setting options for any
  of the pipeline steps
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
    try: choice = int(choice)
    except ValueError:
      print('Please enter either a number or exit!\n')
      continue
    try: chosen_stage = options[choice]
    except IndexError:
      print('Please enter a number within range!')
      continue
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
      try: chosen_setting = configs[change]
      except IndexError:
        print('Please enter a number within range!')
        continue
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
  Pipeline menu, in this menu you can remove/add commands to the pipeline or reorder it
  """
  running = True
  pipeline = file_reader('pipeline')
  order = file_reader('run_order')
  while running:
    print("Welcome to the pipeline customizer\nThe current pipeline is:\n")
    [print(f'{i}\n') for i in order]
    print('To change the pipeline enter 1\nTo change the pipeline order enter 2\nTo exit pipeline enter 0\n')
    try:
      running = int(input())
      if running == 1:
        pipeline_changes = {i:k for i,k in enumerate(pipeline.keys())}
        print('To remove or add a feature, simply enter its number\nIMPORTANT: if you remove a command and add it later, it will be in the last position!\n')
        [print(f'{i} {k} - current setting: {pipeline[k]["switch"]}\n') for i,k in zip(pipeline_changes.keys(), pipeline_changes.values())]
        to_change = pipeline_changes[int(input())]
        pipeline[to_change]['switch'] = not pipeline[to_change]['switch']
        if to_change not in order:
          order.append(to_change)
      elif running == 2:
        print('The current run order is:\n')
        [print(f'{i}\t{v}') for i,v in enumerate(order)]
        print('To move a command enter the number next to it\nNOTE: you can also reorder by entering a sequence of numbers seperated by space (i.e: \"0 3 5 ...\")\n')
        to_move = input()
        if not len(to_move) > 1:
          to_move = int(to_move)
          print(f'Who would you like to replace {order[to_move]} with?\n')
          where = int(input())
          temp = order[where]
          order[where] = order[to_move]
          order[to_move] = temp
        else:
          order = [order[int(i)] for i in to_move.split(' ')]
          for i in pipeline.keys():
            if i not in order:
              pipeline[i]['switch'] = False
    except ValueError:
      print('Please enter a number!\n')
    order = [i for i in order if pipeline[i]['switch']]
    file_creator('pipeline', pipeline)
    file_creator('run_order', order)

def add_new_setting():
  """
  add a new setting menu, currently the user has to reconfigure every menu option.
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
  print(f'Please enter a value for {new_name} (if the setting is an on/off switch, enter nothing!)\nKeywords gtf, fasta fastq and fastq_pre will be replaced with appropriate paths\n')
  value = input()
  settings[chosen_step][new_name]['value']  = value
  settings[chosen_step][new_name]['switch'] = 1
  file_creator('settings', settings) 
  print('Setting added successfuly!\n')

def add_new_pipeline():
  """
  Add new pipeline step menu, currently the user has to configure every option
  """
  pipeline = file_reader('pipeline')
  order = file_reader('run_order')
  settings = file_reader('settings')
  print('To begin defining a new pipeline step, first enter its name\nIMPORTANT: the name must be different from the current pipeline steps\n')
  name = input()
  pipeline[name] = {k:'' for k,v in pipeline['trim'].items() if k != 'switch'}
  for k in pipeline[name].keys():
    print(f'Please enter a {k}\n')
    try: print(explanation_pipe[k] + '\n')
    except KeyError:pass
    v = input()
    pipeline[name][k] = v
  pipeline[name]['switch'] = True
  order.append(name)
  settings[name] = {}
  file_creator('run_order', order)
  file_creator('pipeline', pipeline)
  file_creator('settings', settings)
  print('\nPipeline step added successfuly!\n')

def change_menu():
  """
  Main reconfiguration menu, sends to (currently 4) sub-menus, after exiting sub-menus, returns
  the user here.
  """
  choice = True
  print('####################\n\nWelcome to the configuration menu\n')
  menu = [lambda :0, change_general_settings, change_settings, change_pipeline, add_new_setting, add_new_pipeline]
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

def csh_commands(i, row, settings, pipeline, order, fasta, gtf, PAIRED, org):
  """
  Create a csh fomratted string designed to be run by qsub, works differnetly for PAIRED and unPAIRED samples.

  Paremeters
  ----------
  row : pd.Series
    A pandas series which corresponds to a single row in the sample tsv
  settings : dict
    A nested dictionary which contains the defined settings (added to the command) for each pipeline step
  pipeline : dict
    A dictionary which contains the necessary parameters for each pipeline step
  order : list
    The run-order of all commands that are turned on
  fasta : str
    The path to the ref sequence fasta file
  PAIRED : bool
    True if the sample is paired, false otherwise
  
  Returns
  -------
  csh_str : str
    The string that contains all the commands - to be added to the csh file in main.
  """
  csh_str = ''
  org = org.split(' ')
  try:
    org = org[0][0] + org[1][0]
  except IndexError:
    org = org[0].split('_')
    org = org[0][0] + org[1][0]
  org = org.upper()
  sample = f'{row.run_accession}_{i}_{org}'
  for command in order:
    to_f_string = ['fasta','gtf', 'command_input', 'fastq', 'fastq_pre', 'sample']
    command_key = pipeline[command]['key'] if pipeline[command]['path'] == '' else pipeline[command]['path']
    command_outdir = pipeline[command]['outdir']
    command_inputdir = pipeline[command]['inputdir']
    input_ext = pipeline[command]['input_ext']
    output_ext = pipeline[command]['output_ext']
    if PAIRED:
      fastq = f'{PATH}{slash}{pipeline["trim"]["outdir"]}{slash}{sample}_1_val_1.fq.gz {PATH}{slash}{pipeline["trim"]["outdir"]}{slash}{sample}_2_val_2.fq.gz'
      fastq_pre = f'{PATH}{slash}{pipeline["wget"]["outdir"]}{slash}{sample}_1.fastq.gz {PATH}{slash}{pipeline["wget"]["outdir"]}{slash}{sample}_2.fastq.gz'
    else:
      fastq = f'{PATH}{slash}{pipeline["trim"]["outdir"]}{slash}{sample}_trimmed.fq.gz'
      fastq_pre = f'{PATH}{slash}{pipeline["wget"]["outdir"]}{slash}{sample}.fastq.gz'
    for i,v in settings[command].items():
      if v['value'] == 'fasta': v['value'] = fasta
      if 'output_file' in i:
        v['value'] += slash + sample + output_ext
    if settings[command] != {}:
      for k in settings[command].keys():
        if settings[command][k]['value'] in to_f_string:
          settings[command][k]['value'] = eval(settings[command][k]["value"])
      command_settings = ' '.join([\
        f'{v["key"]} {v["value"]} ' for v in settings[command].values() if v['switch'] == True])
    else:
      command_settings=''
    command_input = f'{PATH}{slash}{command_inputdir}{slash}{sample}{input_ext}'
    command_output = f'{PATH}{slash}{command_outdir}{slash}{sample}{output_ext}'
    structure = pipeline[command]['structure']
    default_structure = f'{command_key} {command_settings} {command_input}'
    if structure == '':
      structure = default_structure
    else:
      struc_list = structure.split(' ')
      struc_list = ['command_' + i if i in ['input','key','settings'] else i for i in struc_list]
      for k,v in enumerate(struc_list):
        struc_list[k] = eval(v)
      try: pass
      except NameError:
        print(f'You entered an invalid command structure! reset to default, to change please re-create the pipeline step: {struc_list}') 
        structure = default_structure
      structure = ' '.join(struc_list)
    csh_str+= f'echo \"Running {command}: {pipeline[command]["description"]} on {sample}\"\n'
    if PAIRED:
      ftps = row.fastq_ftp.split(';')
      if command == 'wget':
        for i in range(2):
          csh_str+=f'{command_key} {command_settings} -O {PATH}{slash}{command_outdir}{slash}{sample}_{i+1}.fastq.gz ftp://{ftps[i]}\n'
          csh_str+='sleep 30s\n'
      elif command == 'trim':
        csh_str+=f'{command_key} {command_settings} --paired {fastq_pre}\n'
      elif command == 'fastqc':
        for i in range(2):
          csh_str+=f'{command_key} {command_settings} {PATH}{slash}{pipeline["trim"]["outdir"]}{slash}{sample}_{i+1}_val_{i+1}.fq.gz\n'
      elif command == 'fastqc_pre':
        for i in range(2):
          csh_str+=f'{command_key} {command_settings} {PATH}{slash}{pipeline["wget"]["outdir"]}{slash}{sample}_{i+1}.fastq.gz\n'
      elif command == 'alignment':
        csh_str+=f'{command_key} {command_settings} {fastq} > {PATH}{slash}{command_outdir}{slash}{sample}.sam\n'
      else:
        csh_str+=structure
        csh_str+=f' > {command_output}\n' if command_outdir != '' else '\n'
    else:
      if command == 'wget':
        csh_str+=f'{command_key} {command_settings} -O {PATH}{slash}{command_outdir}{slash}{sample}.fastq.gz ftp://{row.fastq_ftp}\n'
      elif command == 'trim':
        csh_str+=f'{command_key} {command_settings} {fastq_pre}\n'
      elif command == 'fastqc':
        csh_str+=f'{command_key} {command_settings} {PATH}{slash}{pipeline["trim"]["outdir"]}{slash}{sample}_trimmed.fq.gz\n'
      elif command == 'fastqc_pre':
        csh_str+=f'{command_key} {command_settings} {PATH}{slash}{pipeline["wget"]["outdir"]}{slash}{sample}.fastq.gz\n'
      elif command == 'alignment':
        csh_str+=f'{command_key} {command_settings} {fastq} > {PATH}{slash}{command_outdir}{slash}{sample}.sam\n'
      else:
        csh_str+=structure
        csh_str+=f' > {command_output}\n' if command_outdir != '' else '\n'
  return csh_str

def csh_maker(work_on, tsv_df):
  """
  Makes a single .csh file for each line of .tsv, and a master.csh to be run by qsub

  Parameters
  ----------
  work_on : int/str
    The lines of .tsv that the user wishes to process
  tsv_df
    The df version of samples.tsv
  """
  master = open(f'{PATH}{slash}master.csh', 'w')
  master.write('#!/bin/csh\n')
  master.write('#$ -cwd\n')
  master.write('#$ -q bioinfo.q\n')
  master.write('#$ -N master.csh\n')
  pipeline = file_reader('pipeline')
  general_settings = file_reader('general_settings')
  [dir_checker(i['outdir']) for i in pipeline.values() if i['outdir'] != '']
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
  for i,row in tsv_df.iloc[start:end,:].iterrows():
    settings = file_reader('settings')
    pipeline = file_reader('pipeline')
    order = file_reader('run_order')
    org_name = row.scientific_name

    if len(row.fastq_ftp.split(';')) > 1: PAIRED = True
    else : PAIRED = False
    with open(f'run{i}.csh', 'w') as file:
      file.write(f'#!/bin/{general_settings["environment"][0]}\n')
      file.write(f'#$ -o {PATH}{slash}std_output\n')
      file.write(f'#$ -e {PATH}{slash}error_output\n')
      file.write(f'#$ -cwd\n')
      file.write(f'#$ -N {row.run_accession}\n')
      file.write(f'#$ -q {general_settings["qsub_queue"][0]}\n')
      file.write(general_settings["free_form"][0])
      file.write(f'echo \"Current organism: {row.scientific_name}\"\n\n')
      csh_add_ref = ref_main(org_name)
      if csh_add_ref != '': file.write(csh_add_ref)
      refs = file_reader('refpaths')
      fasta = refs[org_name]['refpath']
      gtf = refs[org_name]['gtfpath']

      to_csh = csh_commands(i, row, settings, pipeline, order, fasta, gtf, PAIRED, org_name)
      file.write(to_csh)
      master.write(f'sleep 1s\n')
      master.write(f'qsub run{i}.csh\n')
  master.close()
def rows_of_links(tsv_path):
  """
  """
  tsv_df = pd.DataFrame(columns = ['study_accession', 'sample_accession', 'experiment_accession', 'run_accession', 'scientific_name', 'fastq_ftp'])
  with open(tsv_path, 'r') as file:
    lines = file.readlines()
    if 'fastq.gz' not in lines[0]:
      lines = lines[1:]
    print("What is the organism's name?")
    tsv_df.fastq_ftp = lines
    tsv_df.scientific_name = input()
    tsv_df.run_accession = 0 #PLACEHOLDER TODO define run accession (FINSIH THIS FUNCTION)

  
def main():
  """
  """
  print('#########\nWelcome to the semi-automatic .csh maker\n#########')
  while True:
    print('Please enter TSV file name (without extension).\nIMPORTANT: The tsv must be in a folder named samples (the folder was created if it doesnt exist)\nIMPORTANT: the tsv must be in the following format:\nstudy_accession\tsample_accession\trun_accession\ttax_id\tscientific_name\tfastq_ftp\tsubmitted_ftp\tsra_ftp\nThe format is available for download from the EBI website\n')
    tsv_name = input()
    tsv_path = f'{PATH}{slash}samples{slash}{tsv_name}.tsv'
    try: 
      tsv_df = pd.read_csv(tsv_path, delimiter = '\t')
      break
    except FileNotFoundError:
      try:
        tsv_path = tsv_path.replace('.tsv','.txt')
        tsv_df = pd.read_csv(tsv_path, delimiter = '\t')
      except FileNotFoundError:
        print('\nthe tsv file does not exist or is not inside the samples folder!\n')
  try: tsv_df.scientific_name = tsv_df.scientific_name.str.replace(' ','_')
  except AttributeError:
    formats = [rows_of_links]
    print("It appears your file is not in EBI format.\n")
    [print(f'Enter {i}\tto {v.__name__}' for i,v in enumerate(formats))]
    choice = input()
    tsv_df = formats[choice](tsv_path)
    
  print('Please enter 1 to enter the configuration menu.\nEnter 0 to continue\n')
  change = int(input())
  if change:
    change_menu()
  print(f'Choose which lines of the tsv to process with the current pipeline.\nEnter a single row number (i.e 2) or a range (i.e 2-3)\nEnter all to work on all rows\n')
  [print(f'{ind}\t{samp}') for ind,samp in zip(tsv_df.index,tsv_df.run_accession)]
  work_on = input()
  csh_maker(work_on, tsv_df)

  
dir_checker('mtDNA_ref')
dir_checker('error_output')
dir_checker('std_output')
dir_checker('Genbank_DB')
dir_checker('samples')
dir_checker('configs')
explanation_pipe = {
  'key':'The command needed to run the pipeline step (i.e samtools view)',
  'description':'A description of what the command does',
  'path':'Path to the command in the cluster (or wherever the code is being run)',
  'outdir':'Output directory to pipe the resulting file into\nNOTICE: Leave empty if you wish to define output as a setting.',
  'output_ext': 'What extension is the command supposed to output? (i.e .bam)',
  'inputdir': 'If the command requires a file to run, what is the name of the directory it is in? (i.e alignment)',
  'input_ext': 'What is the extension of the file the command runs on? (i.e .sam)',
  'structure':'Command structure, enter nothing for default format: key settings input\IMPORTANT: A non-default structure can only include the following keywords: sample, key, input, settings, gtf, fasta, fastq and fastq_pre seperated by a SINGLE space'
}
default_pipe =\
{
'wget' : {
  'key':'wget',
  'description':'Download command',
  'switch':True,
  'path':r'',
  'outdir':'fastq',
  'output_ext':'',
  'inputdir':'',
  'input_ext':'',
  'structure':''
},
'trim' : {
  'key':'trim_galore',
  'description':'Trims adapters',
  'switch':True,
  'path':r'',
  'outdir':'post_trim',
  'output_ext':'.fq.gz',
  'inputdir':'fastq',
  'input_ext':'.fastq',
  'structure':'key settings fastq_pre'
},
'fastqc':{
  'key':'fastqc',
  'description':'Create a readable quality report in html format',
  'switch':True,
  'path':r'fastqc',
  'outdir':'post_trim',
  'output_ext':'',
  'inputdir':'',
  'input_ext':'',
  'structure':''
},
'fastqc_pre':{
  'key':'fastqc',
  'switch':True,
  'description':'Create a readable quality report in html format',
  'path':r'fastqc',
  'outdir':'fastq',
  'output_ext':'',
  'inputdir':'',
  'input_ext':'',
  'structure':''
},
'bwa-align' : {
  'key':'bwa mem',
  'description':'Align reads to reference',
  'switch':True,
  'path':r'',
  'outdir':'alignment',
  'output_ext':'.sam',
  'inputdir':'post_trim',
  'input_ext':'.fq.gz',
  'structure':'key settings fasta fastq'
},
'view':{
  'key':'samtools view',
  'description':'Convert sam to bam',
  'switch':True,
  'path':r'',
  'outdir':'alignment',
  'output_ext':'.bam',
  'inputdir':'alignment',
  'input_ext':'.sam',
  'structure':''
},
'index':{
  'key':'samtools index',
  'description':'Index bam for fast access via samtools',
  'switch':True,
  'path':r'',
  'outdir':'alignment',
  'output_ext':'.bam.bai',
  'inputdir':'alignment',
  'input_ext':'.bam',
  'structure':''
},
'mpileup':{
  'key':'samtools mpileup',
  'description':'Create a readable mpileup format',
  'switch':True,
  'path':r'',
  'outdir':'pileup',
  'output_ext':'.pileup',
  'inputdir':'alignment',
  'input_ext':'.bam',
  'structure':''
},
'sort':{
  'key': 'samtools sort',
  'description':'Sort reads',
  'switch':True,
  'path':r'',
  'outdir':'',
  'output_ext':'.bam',
  'inputdir':'alignment',
  'input_ext':'.bam',
  'structure':''
},
'htseq-count':{
  'key': 'htseq-count',
  'description': 'Measure expression',
  'switch': False,
  'path': r'',
  'outdir':'count',
  'output_ext':'.txt',
  'inputdir':'alignment',
  'input_ext':'.bam',
  'structure': 'key settings input gtf'
},
'STAR-alignment':{
  'key':'STAR',
  'description':'Align reads to ref',
  'switch': False,
  'path': 'STAR',
  'outdir':'alignment',
  'output_ext':'.sam',
  'inputdir':'post_trim',
  'input_ext':'.fq.gz',
  'structure':'key settings'
}
}

default_settings =\
{
  'trim':{
    'quality':{'description':'Trim low-quality ends from reads', 'key':'-q', 'value':20, 'switch':True},
    'outdir':{'description':'Choose output directory', 'key':'-o', 'value':f'{PATH}{slash}post_trim', 'switch':True}
  },
  'bwa-align':{
    'seed':{'description':'Minimum seed length', 'key':'-k', 'value':19, 'switch':True},
  },
  'view':{
    'header':{'description':'Include the header in the output', 'key':'-h', 'value':'', 'switch':True},
    'bam':{'description':'Output in BAM format', 'key':'-b', 'value':'', 'switch':True}
  },
  'index':{},
  'mpileup':{
    'quality':{'description':'Minimum base quality for a base to be considered', 'key':'-Q', 'value':30, 'switch':True},
    'reference':{'description':'Fasta type ref file','key':'-f', 'value':'fasta', 'switch':True}
  },
  'fastqc':{
    'outdir':{'description':'Output directory', 'key' : '--outdir', 'value':f'{PATH}{slash}post_trim', 'switch':True}
  },
   'fastqc_pre':{
    'outdir':{'description':'Output directory', 'key' : '--outdir', 'value':f'{PATH}{slash}fastq', 'switch':True}
  },
  'wget':{},
  'sort':{
    'output_file':{'description':'Output file', 'key' : '-o', 'value':f'{PATH}{slash}alignment','switch':True}
  },
  'htseq-count':{},
  'STAR-alignment':{
    'genomeDir':{'description':'Input ref', 'key' : '--genomeDir', 'value':'fasta', 'switch':True},
    'readFilesIn':{'description':'Input fastq', 'key': '--readFilesIn', 'value':'fastq', 'switch':True}
  }
}
default_order = ['wget','fastqc_pre','trim','fastqc','bwa-align','view','sort','index','mpileup']
default_main = {
  'qsub_queue' : ['bioinfo.q', 'choose which cluster queue to send the jobs into'],
  'environment':['bash', 'choose the linux environment to use'],
  'free_form':['', 'any command you type here will be added to every single .csh file created, enter nothing to wipe setting']
}

defaults = (default_pipe, default_settings, default_order, default_main)
if not os.path.exists(f'{PATH}{slash}configs{slash}refpaths.json'):
  file_creator('refpaths', {})
if not os.path.exists(f'{PATH}{slash}configs{slash}credentials.json'):
  EMAIL = input('Please enter your NCBI account email:')
  API = input('Please enter your NCBI API (can be found in NCBI account page)')
  file_creator('credentials', {'Email':EMAIL, 'API':API})

creds = file_reader('credentials')
configs = ['pipeline', 'settings', 'run_order', 'general_settings']
for i,config in enumerate(configs):
  if not os.path.exists(f'{PATH}{slash}configs{slash}{config}.json'):
    file_creator(config, defaults[i])
main()
print('Enter \"qsub master.csh\" into the console to run the commands\nIMPORTANT: Each sample will be executed simultaneously in qsub you run this command, make sure you dont run too many at once!\n')
