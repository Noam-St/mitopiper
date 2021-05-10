"""
#########################################################################

This program creates a custom csh file designed to be executed with qsub.
The program is designed to be configurable, the user can customize the processing pipeline,
the settings for each of the pipeline programs used, and the paths of the programs utallized.

#########################################################################
"""
__author__ = 'Noam Shtolz'


import sys, shutil, os, math, json
from textwrap import wrap

#%%
import pandas as pd
from Bio import Entrez, SeqIO
from tabulate import tabulate

DASH = '--------------------------------------------'
PATH = os.getcwd()
if sys.platform != 'win32':
    SLASH = '/'
else:
    SLASH = '\\'
try:
  org_df = pd.read_csv(f'{PATH}{SLASH}src{SLASH}final.csv')
  org_df.organism = org_df.organism.str.replace(' ','_')
  org_df = org_df.set_index('organism')
except FileNotFoundError: raise FileNotFoundError('final.csv is required in the src folder for the program to operate.')


def file_reader(filename : str, directory : str = 'current') -> dict: 
  """
  Read .json files in the format key:value to a dictionary

  Parameters
  ----------
  filename : str
    PATH to txt file in the expected format
  Returns
  -------
  format_dict : dict
    Dictionary with the left side of the equal sign as key and right side as value
  """
  if directory == '':
    with open(f'{PATH}{SLASH}config{SLASH}{filename}.json', 'r') as file:
      format_dict = json.load(file)
  else:
    with open(f'{PATH}{SLASH}config{SLASH}{directory}{SLASH}{filename}.json', 'r') as file:
      format_dict = json.load(file)
  return format_dict


def file_creator(filename : str, format_dict : dict, directory : str = 'current') -> None:
  """
  Create json files of the format_dict, place file in cwd

  Parameters
  ----------
  filename : str
    Name of file to be created
  format_dict : dict
    Dictionary by-which to write the file in this format: key=value
  """
  if directory == '':
    filepath = f'{PATH}{SLASH}config{SLASH}{filename}.json'
  else:
    filepath = f'{PATH}{SLASH}config{SLASH}{directory}{SLASH}{filename}.json'
  with open(filepath, 'w') as file:
    json.dump(format_dict, sort_keys = True, fp = file, indent = 1)


def dir_checker(dirname : str) -> None:
  """
  check if a directory exists, if it does not, create it.

  Parameters
  ----------
  dirname : str
    Name of the directory to be checked
  """
  dirname = dirname.replace('/', SLASH)
  if not os.path.isdir(f'{PATH}{SLASH}{dirname}'):
    os.mkdir(f'{PATH}{SLASH}{dirname}')


def fasta_wrap(string : str, per : int = 60) -> str:
  """
  Wraps a string to a default of 60 characters per line

  Parameters
  ----------
  string : str
    The string to be wrapped
  per : int
    Amount of characters per line
  Returns
  -------
  str
    The wrapped string  
  """
  return '\n'.join(wrap(string, per))


def ref_maker(org_name : str) -> str:
  """
  CURRENTLY UNUSED

  Create a fasta file for indexing purposes.

  Parameters
  ----------
  org_name : str
    Name of the organism to name the file.
  Returns
  -------
  fasta : str
    Absolute path to the fasta file created.
  """
  ID = org_df.loc[org_name, 'RefSeq']
  filename = f'{PATH}{SLASH}src{SLASH}genbank{SLASH}{ID}.gbk'
  creds = file_reader('credentials')
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
  fasta = f'{PATH}{SLASH}src{SLASH}mtDNA_ref{SLASH}{org_name}{SLASH}{org_name}.fasta'
  with open(fasta, 'w') as ref:
    ref.write('>MT\n')
    ref.write(fasta_wrap(seq))
  return fasta, ID


def tabulate_dict(dicty : dict) -> tabulate:
  """
  """
  temp_df = pd.DataFrame(dicty)
  return tabulate(temp_df, headers = 'keys', tablefmt = 'psql')


def gff3_maker(org_name : str) -> str:
  """
  CURRENTLY UNUSED

  Create a gtf(gff3) file if needed, return the path to the gtf file
  """
  ID = org_df.loc[org_name, 'RefSeq']
  gtf = f'{PATH}{SLASH}src{SLASH}mtDNA_ref{SLASH}{org_name}{SLASH}{org_name}.gff3'
  creds = file_reader('credentials')
  if not os.path.isfile(gtf):
    if creds['Email'] != '':
      Entrez.email = creds['Email']
      Entrez.api_key = creds['API']
    net_handle = Entrez.efetch(db = 'nucleotide', id = ID, rettype = 'gtf', retmode = 'text')
    out_handle = open(gtf, 'w')
    out_handle.write(net_handle.read())
    net_handle.close()
    out_handle.close()
  return gtf


def save_settings(source : str, output : str) -> None:
  """
  Transfer settings from source to output
  """
  dir_checker('config/' + output)
  for jsonfile in os.listdir(f'{PATH}{SLASH}config{SLASH}{source}'):
    if '.json' in jsonfile:
      jsonfile = jsonfile.replace('.json','')
      file_creator(filename = jsonfile, format_dict = file_reader(filename = jsonfile, directory = source), directory = output)


def sa_index_calc(glen : int) -> float:
  """
  Calculate the SA index according to STAR manual
  """
  return round(min(14, math.log2(glen)/2 - 1),0)


def save_current_settings() -> None:
  """
  Save the current settings preset
  """
  print("Please name your preset\n")
  name = input()
  save_settings(source = 'current', output = name)
  print('\nSettings saved successfuly!\n')
  return


def ref_main(org_name : str) -> str:
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
  refs = file_reader('refpaths', 'global')
  pipeline = file_reader('pipeline')
  try:
    refs[org_name]['refpath']
    refs[org_name]['gtfpath']
    return ''
  except KeyError:
    pass
  choice = True
  while choice:
    print('It seems a reference fasta+gtf do not exist in the default location\nTo create them automatically, press 1\nNOTICE: The reference fasta+gtf will be for the MT only and in NCBI format!\nTo enter pre-indexed fasta and gff pathes, press 2.\n')
    choice = input()
    try: choice = int(choice)
    except ValueError: print('Please enter 1 or 2 only')
    if choice == 1:
      dir_checker(f'src{SLASH}mtDNA_ref{SLASH}{org_name}')
      try: ID = org_df.loc[org_name, 'RefSeq']
      except KeyError:
        print(f"Warning:\tIndex can not be created for {org_name}!\n")
        return False
      fasta = f'{PATH}{SLASH}src{SLASH}mtDNA_ref{SLASH}{org_name}{SLASH}{org_name}.fasta'
      os.system(f'wget -O {fasta} \"https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id={ID}\"\n')
      temp_fasta = open(fasta, 'r')
      glen = len(''.join([i.rstrip('\n') for i in temp_fasta.readlines()[1:]]))
      sa_index = sa_index_calc(glen)
      temp_fasta.close()
      fasta_dir = os.path.dirname(fasta)
      try: 
        STAR_command = pipeline['STAR-alignment']['path'] if pipeline['STAR-alignment']['path'] != '' else  pipeline['STAR-alignment']['key']
        STAR_command += f' --runMode genomeGenerate --genomeDir {fasta_dir} --runThreadN 4 --genomeFastaFiles {fasta} --genomeSAindexNbases {sa_index}\n'
      except KeyError: STAR_command = ''
      try:
        BWA_command = f'samtools faidx {fasta}\n'
        BWA_command += pipeline['bwa-align']['path'] if pipeline['bwa-align']['path'] != '' else  pipeline['bwa-align']['key'].split(' ')[0]
        BWA_command += f' index {fasta}\n'
      except KeyError: BWA_command = ''
      gtf = f'{PATH}{SLASH}src{SLASH}mtDNA_ref{SLASH}{org_name}{SLASH}{org_name}.gff3'
      GTF_command = f'wget -O {gtf} \"https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id={ID}\"\n'
      add_to_csh = GTF_command + STAR_command + BWA_command 
    elif choice == 2:
      print('Please enter the full path to the reference fasta\n')
      fasta = input()
      add_to_csh = ''
      print('Please enter path to gtf file for this organism\nEnter nothing if annotations are not needed.\n')
      gtf = input()
    else:
      print('Please enter 1 or 2 only')
 
    refs[org_name] = {}
    refs[org_name]['refpath'] = fasta
    refs[org_name]['gtfpath'] = gtf
    file_creator('refpaths', refs, 'global') #NOTE(Noam) Currently set to global, think about changing to local pre-configurations to enable multiple ref files for the same organism depending on preset.
    return add_to_csh


def change_general_settings() -> None:
  """
  Change general settings defined on a per-run basis.
  """
  general_settings = file_reader('general_settings')
  options = list(general_settings.keys())
  while True:
    temp_setting = {'Setting name' : options, 'Current setting' : [general_settings[v][0] for v in options]}
    print(tabulate_dict(temp_setting))
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
    print(f"Please enter a new value for {options[choice]}\nCurrent value: {chosen[0]}\nDescription: {chosen[1]}\n")
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
  print("Setting changed successfuly!\n")


def change_pipeline_settings() -> None:
  """
  Settings change menu, in this menu you can reconfigure any of the setting options for any of the pipeline steps
  """
  settings = file_reader('settings')
  pipeline = file_reader('pipeline')
  options = list(settings.keys())
  order = file_reader('run_order')
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
    temp_config_df = {
      'Setting name':configs,
      'Description':[i['description'] for i in settings[chosen_stage].values()],
      'Current setting':[i['value'] for i in settings[chosen_stage].values()]
      }
    print(tabulate_dict(temp_config_df))
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
        print(f'You chose: {chosen_setting}\nPlease choose a new value, current value: {setting_value}\n\nThe value may include several Keywords: _gtf, _fasta, _fastq, _fastq_pre, _fastq_dir, _fastq_dir_pre, _bam and _sam\nThese keywords will be replaced with appropriate paths\n')
        print(f'Additionally, you may add the combination of output_command_name or dir_command_name to include the output of another command in the setting value\nNOTE: you are adding a setting to {chosen_stage}. You can grab the output of the following commands with the current run order\n')
        try: print('\t'.join(order[:order.index(chosen_stage)]))
        except ValueError: print('Step isnt active!')
        new_val = input()
        settings[chosen_stage][chosen_setting]['value'] = new_val
        print('\nSetting changed successfuly!\n')
    
  file_creator('settings', settings)


def change_pipeline() -> None:
  """
  Pipeline menu, in this menu you can remove/add commands to the pipeline or reorder it
  """
  running = True
  pipeline = file_reader('pipeline')
  order = file_reader('run_order')
  while running:
    print("Welcome to the pipeline customizer\nThe current pipeline is:\n")
    [print(f'{i}\n') for i in order]
    print('To change the pipeline enter 1\nTo change the pipeline order enter 2\nTo exit this menu enter 0\n')
    try:
      running = int(input())
      if running == 1:
        pipeline_changes = {i:k for i,k in enumerate(pipeline.keys())}
        print('To remove or add a feature, simply enter its number\nIMPORTANT: if you remove a command and add it later, it will appended in the last position!\n')
        [print(f'{i} {k} - current setting: {pipeline[k]["switch"]}\n') for i,k in zip(pipeline_changes.keys(), pipeline_changes.values())]
        to_change = pipeline_changes[int(input())]
        pipeline[to_change]['switch'] = not pipeline[to_change]['switch']
        if to_change not in order:
          order.append(to_change)
      elif running == 2:
        print('The current run order is:\n')
        [print(f'{i}\t{v}') for i,v in enumerate(order)]
        print('To move a command enter the number next to it\nNOTE: you can also reorder quickly by entering a sequence of numbers seperated by space (i.e: \"0 3 5 ...\")\n')
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
    print('Pipeline changed successfuly!\n')


def key_fixer(key : str) -> str:
  """
  Fix a key that does not contain dashes or double dashes according to unix rules

  Parameters
  ----------
  key : str
    The key to fix
  Returns
  -------
  key : str
    The fixed key
  """
  if len(key) > 1:
    key = '--' + key.replace(' ','-')
  else:
    key = '-' + key
  return key


def add_new_pipeline_setting() -> None:
  """
  add a new setting menu, currently the user has to reconfigure every menu option.
  """
  settings = file_reader('settings')
  order = file_reader('run_order')
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
  if '-' not in key:
    key = key_fixer(key)
  try: int(key)
  except ValueError: pass
  settings[chosen_step][new_name]['key'] = key
  print(f'Please enter a value for {new_name} (if the setting is an on/off switch, enter nothing!)\nYou can also enter a keyword which will be replaced with the appropriate file path, Keywords: _gtf, _fasta, _fastq, _fastq_pre, _fastq_dir, _fastq_dir_pre, _bam and _sam, _command_output and _command_input will be replaced with appropriate paths\n')
  print(f'Additionaly, you may add the combination of output_command_name or dir_command_name to include the output or directory of another command in the setting value\nNOTE: you are adding a setting to {chosen_step}. You can grab the output of the following commands with the current run order\n')
  try: print('\t'.join(order[:order.index(chosen_step)]))
  except ValueError: print('Current command isnt active!')
  value = input()
  settings[chosen_step][new_name]['value']  = value
  settings[chosen_step][new_name]['switch'] = 1
  file_creator('settings', settings) 
  print(f'Setting added successfuly!\nName\tKey\tValue\n{new_name}\t{key}\t{value}\n')


def add_new_pipeline() -> None:
  """
  Add new pipeline step menu, currently the user has to configure every option
  """
  explanation_pipe = {
  'key':'The command needed to run the pipeline step (i.e samtools view)',
  'description':'A description of what the command does',
  'path':'Path to the command in the cluster (or wherever the code is being run)',
  'outdir':'Output directory to pipe the resulting file into\nNOTICE: Leave empty if you wish to define output as a setting.',
  'output_ext': 'What extension is the command supposed to output? (i.e .bam)',
  'inputdir': 'If the command requires a file to run, what is the name of the directory it is in? (i.e alignment)',
  'input_ext': 'What is the extension of the file the command runs on? (i.e .sam)',
  'structure':'Command structure, enter nothing for default format: key settings input\IMPORTANT: A non-default structure can only include the following keywords: sample, key, settings, gtf, fasta, fastq, fastq_pre fastq_dir, fastq_pre_dir, sam, bam, command_output, command_input and output_pipe seperated by a SINGLE space'
}

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


def manually_change_ref_or_gtf() -> None:
  """
  Manually change fasta or gtf path after assigning them through the regular process.
  """
  ref = file_reader('refpaths', 'global')
  if ref == {}:
    print('Please run the program once before attempting to change ref or gtf manually!')
    return
  keys = list(ref.keys())
  while True:
    print("Please choose an organism to add reference fasta or gtf path to\nEnter exit to return to main menu")
    [print(f'{i}\t{v}\n') for i,v in enumerate(keys)]
    org = input()
    if org == 'exit': break
    try: org = keys[int(org)]
    except IndexError:
      print("Please enter an organism from the list!\n")
      continue
    print('Would you like to change gtf, fasta or both?\n')
    work_on = ['ref', 'gtf', 'both']
    [print(f'{i}\t{v}\n') for i,v in enumerate(work_on)]
    choice = input()
    try: work_on = work_on[int(choice)]
    except IndexError:
      print('Please enter one of the possible options!\n')
      continue
    if work_on == 'both':
      work_on = ['ref', 'gtf']
    else:
      work_on = [work_on]
    for i in work_on:
      print(f"Please enter the {i} absolute path!\n")
      path = input()
      i += 'path'
      ref[org][i] = path
    break
  file_creator(filename = 'refpaths', format_dict = ref, directory = 'global') 


def change_menu():
  """
  Main reconfiguration menu, sends to (currently 4) sub-menus, after exiting sub-menus, returns
  the user here.
  """
  choice = True
  print('####################\n\nWelcome to the configuration menu\n')
  menu = [lambda :0, change_general_settings, change_pipeline_settings, change_pipeline, add_new_pipeline_setting, add_new_pipeline,manually_change_ref_or_gtf, save_current_settings]
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


def csh_commands(i, row, settings, pipeline, order, fasta, gtf, PAIRED, org, tsv_dir):
  """
  Create a csh fomratted string designed to be run by qsub, works differnetly for PAIRED and unPAIRED samples.

  Paremeters
  ----------
  i : int
    The index of the current row
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
  gtf : str
    The path to the GTF file (usually GFF3)
  PAIRED : bool
    True if the sample is paired, false otherwise
  org : str
    The name of the current organism
  tsv_dif : str
    The path the sample tsv is located in
  Returns
  -------
  csh_str : str
    The string that contains all the commands - to be added to the csh file in main.
  """
  csh_str = ''
  org_name = org
  org = org.split(' ')
  try:
    org = org[0][0] + org[1][0]
  except IndexError:
    org = org[0].split('_')
    org = org[0][0] + org[1][0]
  org = org.upper()
  sample = f'{row.run_accession}_{i}_{org}'
  command_outputs = {f'output_{command}':f'{tsv_dir}{SLASH}{pipeline[command]["outdir"]}{SLASH}{sample}{pipeline[command]["output_ext"]}' for command in pipeline.keys()}
  command_dirs = {f'dir_{command}':f'{tsv_dir}{SLASH}{pipeline[command]["outdir"]}' for command in pipeline.keys()}
  to_f_string = ['_fasta', '_fasta_dir', '_gtf', '_command_input', '_fastq', '_fastq_pre', '_sample', '_fastq_dir', '_fastq_pre_dir', '_bam', '_sam', '_STAR_out', '_fastq1', '_fastq2', '_org_name', '_command_output']
  for command in order:
    command_key = pipeline[command]['key'] if pipeline[command]['path'] == '' else (pipeline[command]['path'] +' ' + ' '.join(pipeline[command]['key'].split(' ')[1:]))
    command_outdir = pipeline[command]['outdir']
    command_inputdir = pipeline[command]['inputdir']
    input_ext = pipeline[command]['input_ext']
    output_ext = pipeline[command]['output_ext']
    if PAIRED:
      fastq = f'{tsv_dir}{SLASH}{pipeline["trim"]["outdir"]}{SLASH}{sample}_1_val_1.fq.gz {tsv_dir}{SLASH}{pipeline["trim"]["outdir"]}{SLASH}{sample}_2_val_2.fq.gz'
      fastq1, fastq2 = fastq.split(' ')
      fastq_dir =  f'{tsv_dir}{SLASH}{pipeline["trim"]["outdir"]}'
      fastq_pre = f'{tsv_dir}{SLASH}{pipeline["wget"]["outdir"]}{SLASH}{sample}_1.fastq.gz {tsv_dir}{SLASH}{pipeline["wget"]["outdir"]}{SLASH}{sample}_2.fastq.gz'
      fastq_pre_dir = f'{tsv_dir}{SLASH}{pipeline["wget"]["outdir"]}'
    else:
      fastq = f'{tsv_dir}{SLASH}{pipeline["trim"]["outdir"]}{SLASH}{sample}_trimmed.fq.gz'
      fastq_pre = f'{tsv_dir}{SLASH}{pipeline["wget"]["outdir"]}{SLASH}{sample}.fastq.gz'
      fastq_dir =  f'{tsv_dir}{SLASH}{pipeline["trim"]["outdir"]}'
      fastq_pre_dir = f'{tsv_dir}{SLASH}{pipeline["wget"]["outdir"]}'
      fastq1 = fastq
      fastq2 = fastq

    bam =  f'{tsv_dir}{SLASH}{pipeline["bwa-align"]["outdir"]}{SLASH}{sample}*.bam'
    sam =  f'{tsv_dir}{SLASH}{pipeline["bwa-align"]["outdir"]}{SLASH}{sample}*.sam'
    STAR_out = f'{tsv_dir}{SLASH}{pipeline["bwa-align"]["outdir"]}{SLASH}{sample}'
    fasta_dir = os.path.dirname(fasta)
    command_input = f'{tsv_dir}{SLASH}{command_inputdir}{SLASH}{sample}*{input_ext}'
    command_output = f'{tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}{output_ext}'
    output_pipe = f' > {command_output}\n'
    for i,v in settings[command].items():
      if v['value'] == 'fasta': v['value'] = fasta
      if 'output_file' in i:
        v['value'] += SLASH + sample + output_ext
    if settings[command] != {}:
      for k in settings[command].keys():
        cur_val = settings[command][k]['value']
        if cur_val in to_f_string:
          cur_val = eval(cur_val.replace('_','',1))
        elif cur_val in command_outputs.keys():
          cur_val = command_outputs[cur_val]
        elif cur_val in command_dirs.keys():
          cur_val = command_dirs[cur_val]

        settings[command][k]['value'] = cur_val
      command_settings = ' '.join([\
        f'{v["key"]} {v["value"]} ' for v in settings[command].values() if v['switch'] == True])
    else:
      command_settings=''

    structure = pipeline[command]['structure']
    default_structure = f'{command_key} {command_settings} {command_input}'
    if structure == '':
      structure = default_structure
    else:
      struc_list = structure.split(' ')
      struc_list = ['command_' + i if i in ['input','key','settings', 'outdir', 'output'] else i for i in struc_list]
      for k,v in enumerate(struc_list):
        if v in command_outputs.keys():
          struc_list[k] = command_outputs[v]
        elif v in command_dirs.keys():
          struc_list[k] = command_dirs[v]
        elif v in to_f_string:
          struc_list[k] = eval(v.replace('_','',1))
        try: struc_list[k] = eval(v)
        except NameError: pass
      structure = ' '.join(struc_list)
    csh_str+= f'echo \"Running {command}: {pipeline[command]["description"]} on {sample}\"\n'
    if PAIRED:
      ftps = row.fastq_ftp.split(';')
      if command == 'wget':
        for i in range(2):
          csh_str+=f'{command_key} {command_settings} -O {tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}_{i+1}.fastq.gz ftp://{ftps[i]}\n'
          csh_str+='sleep 30s\n'
      elif command == 'trim':
        csh_str+=f'{command_key} {command_settings} --paired {fastq_pre}\n'
      elif command == 'fastqc':
        for i in range(2):
          csh_str+=f'{command_key} {command_settings} {tsv_dir}{SLASH}{pipeline["trim"]["outdir"]}{SLASH}{sample}_{i+1}_val_{i+1}.fq.gz\n'
      elif command == 'fastqc_pre':
        for i in range(2):
          csh_str+=f'{command_key} {command_settings} {tsv_dir}{SLASH}{pipeline["wget"]["outdir"]}{SLASH}{sample}_{i+1}.fastq.gz\n'
      elif command == 'alignment':
        csh_str+=f'{command_key} {command_settings} {fastq} > {tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}.sam\n'
      else:
        csh_str+=structure
        csh_str+='\n'
    else:
      if command == 'wget':
        csh_str+=f'{command_key} {command_settings} -O {tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}.fastq.gz ftp://{row.fastq_ftp}\n'
      elif command == 'trim':
        csh_str+=f'{command_key} {command_settings} {fastq_pre}\n'
      elif command == 'fastqc':
        csh_str+=f'{command_key} {command_settings} {tsv_dir}{SLASH}{pipeline["trim"]["outdir"]}{SLASH}{sample}_trimmed.fq.gz\n'
      elif command == 'fastqc_pre':
        csh_str+=f'{command_key} {command_settings} {tsv_dir}{SLASH}{pipeline["wget"]["outdir"]}{SLASH}{sample}.fastq.gz\n'
      elif command == 'alignment':
        csh_str+=f'{command_key} {command_settings} {fastq} > {tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}.sam\n'
      else:
        csh_str+=structure
        csh_str += '\n'
  return csh_str


def csh_maker(work_on, tsv_df, tsv_name):
  """
  Makes a single .csh file for each line of .tsv, and a master.csh to be run by qsub

  Parameters
  ----------
  work_on : int/str/list
    The lines of .tsv that the user wishes to process
  tsv_df
    The df version of samples.tsv
  tsv_name : str
    The name of the tsv, given as input by the user
  """
  tsv_dir = f"{PATH}{SLASH}output{SLASH}{tsv_name}"
  general_settings = file_reader('general_settings')
  master = open(f'{PATH}{SLASH}master.csh', 'w')
  master.write(f'#!/bin/bash\n')
  pipeline = file_reader('pipeline')
  [dir_checker('output' + '/' + tsv_name + '/' + i['outdir']) for i in pipeline.values() if i['outdir'] != '']
  tsvs = []
  if '-' in work_on:
    work_on = [int(i) for i in work_on.split('-')]
    start = work_on[0]
    end = work_on[1]
    work_on = range(start, end)
  elif work_on == 'all':
    end = len(tsv_df)
    start = 0
    work_on = range(start, end)
  elif type(work_on) == list:
    pass
  else:
    work_on = int(work_on)
    if work_on > len(tsv_df):
      print("Error:\tRow entered is out of bounds! Running on last row\n\n")
      work_on = len(tsv_df) - 1
    work_on = range(work_on, work_on + 1)
  bad_orgs = []
  for i,row in tsv_df.iloc[work_on, :].iterrows():
    settings = file_reader('settings')
    pipeline = file_reader('pipeline')
    order = file_reader('run_order')
    org_name = row.scientific_name
    if org_name in bad_orgs:
      print(f'Warning:\tSkipped organism {org_name} because it does not exist in final.csv!\n\n')
      continue
    csh_add_ref = ref_main(org_name)
    if csh_add_ref == False:
      bad_orgs.append(org_name)
      continue
    if len(row.fastq_ftp.split(';')) > 1: PAIRED = True
    else : PAIRED = False
    with open(f'{tsv_dir}{SLASH}run{i}.csh', 'w') as file:
      file.write(f'#!/bin/{general_settings["environment"][0]}\n')
      file.write(f'#$ -o {tsv_dir}{SLASH}stdout\n')
      file.write(f'#$ -e {tsv_dir}{SLASH}stderr\n')
      file.write(f'#$ -cwd\n')
      file.write(f'#$ -N {row.run_accession}_{i}\n')
      file.write(f'#$ -q {general_settings["qsub_queue"][0]}\n')
      file.write(general_settings["free_form"][0])
      file.write(f'echo \"Current organism: {row.scientific_name}\"\n\n')
      if csh_add_ref != '': file.write(csh_add_ref)
      refs = file_reader('refpaths', 'global')
      fasta = refs[org_name]['refpath']
      gtf = refs[org_name]['gtfpath']

      to_csh = csh_commands(i, row, settings, pipeline, order, fasta, gtf, PAIRED, org_name, tsv_dir)
      file.write(to_csh)
      tsvs.append(f'{tsv_dir}{SLASH}run{i}.csh')

  master.write(f'for file in {" ".join(tsvs)}\ndo\nqsub $file\ndone')
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


def defaults():
  """
  Set the defaults - hard coded.
  """

  dir_checker('config')
  dir_checker('config/default')
  dir_checker('config/global')

  #if not os.path.exists(f"{PATH}{SLASH}config{SLASH}default{SLASH}creds.json"):
  #  print("Please enter your NCBI account email:\n")
  #  EMAIL = input()
  #  print("Please enter your NCBI API key: (can be found in the NCBI account page)\n")
  #  API = input()
  #  file_creator('credentials', {'Email':EMAIL, 'API':API}, directory = 'default')
  
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
  'structure':'key settings command_input output_pipe'
},
'sort':{
  'key': 'samtools sort',
  'description':'Sort reads',
  'switch':True,
  'path':r'',
  'outdir':'',
  'output_ext':'',
  'inputdir':'alignment',
  'input_ext':'.bam',
  'structure':'key settings command_input'
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
  'structure': 'key settings input gtf output_pipe'
},
'STAR-alignment':{
  'key':'STAR',
  'description':'Align reads to ref',
  'switch': False,
  'path': 'STAR',
  'outdir':'',
  'output_ext':'',
  'inputdir':'post_trim',
  'input_ext':'.fq.gz',
  'structure':'key settings'
}
}

  default_settings =\
{
  'trim':{
    'quality':{'description':'Trim low-quality ends from reads', 'key':'-q', 'value':20, 'switch':True},
    'outdir':{'description':'Choose output directory', 'key':'-o', 'value':'_fastq_dir', 'switch':True}
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
    'reference':{'description':'Fasta type ref file','key':'-f', 'value':'_fasta', 'switch':True}
  },
  'fastqc':{
    'outdir':{'description':'Output directory', 'key' : '--outdir', 'value':'_fastq_dir', 'switch':True}
  },
   'fastqc_pre':{
    'outdir':{'description':'Output directory', 'key' : '--outdir', 'value':f'_fastq_pre_dir', 'switch':True}
  },
  'wget':{},
  'sort':{
    'output':{'description':'Output directory', 'key' : '-o', 'value': '_bam', 'switch':True}
  },
  'htseq-count':{
    'ID':{'description': 'ID attribute to use', 'key' : '-i', 'value' : 'ID', 'switch':True},
    'Additional attribute':{'description': 'Additional attribute to use', 'key':'--additional-attr', 'value':'gene','switch':True},
    'Feature column name':{'description': 'Feature column name to use','key':'-t', 'value':'gene', 'switch':True},
    'stranded':{'description': 'Whether the reads are stranded', 'key':'-s','value':'no','switch':True},
    'format':{'description':'Format to read input in','key':'-f','value':'bam','switch':True},
    'ambiguous':{'description':'How to handle ambiguous reads', 'key':'--nonunique','value':'all','switch':True}
  },
  'STAR-alignment':{
    'genomeDir':
    {
    'description':'Input ref',
    'key' : '--genomeDir',
    'value':'_fasta_dir',
    'switch':True
    },
    'readFilesIn':
    {
    'description':'Input fastq',
    'key': '--readFilesIn',
    'value':'_fastq',
    'switch':True
    },
    "outSAMtype":
    {
    "description": "output format",
    "key": "--outSAMtype",
    "switch": True,
    "value": "BAM SortedByCoordinate"
    },
    "readFilesCommand":
    {
    "description": "Input format",
    "key": "--readFilesCommand",
    "switch": True,
    "value": "gunzip -c"
    },
    "runThreadN":
    {
    "description": "The amount of threads to use",
    "key": "--runThreadN",
    "switch": True,
    "value": 4
    },
    "outFileNamePrefix":
    {
    "description": "Output prefix",
    "key": '--outFileNamePrefix',
    'switch': True,
    'value': '_STAR_out'
    },
    'limit sort RAM':
    {
    'description':'whatever',
    'key':'--limitBAMsortRAM',
    'value':1205130643,
    'switch':True
    }}
}
  default_order = ['wget','fastqc_pre','trim','fastqc','bwa-align','view','sort','index','mpileup']
  default_main = {
  'qsub_queue' : ['bioinfo.q', 'choose which cluster queue to send the jobs into'],
  'environment':['bash', 'choose the linux environment to use'],
  'free_form':['', 'any command you type here will be added to every single .csh file created, enter nothing to wipe setting']
}
  defaults = (default_pipe, default_settings, default_order, default_main)
  
  file_creator('refpaths', {}, directory = 'global')
  configs = ['pipeline', 'settings', 'run_order', 'general_settings']
  for i,config in enumerate(configs):
    file_creator(config, defaults[i], directory = 'default')


def tsv_handler():
  """
  TODO HERE, This function will set up a sample directory tree for each sample, based on the tsv name.
  """

  while True:
    print('Please enter TSV file name (without extension) or the number next to a pre-existing TSV file.\nIMPORTANT: The tsv must be in a folder named samples (the folder was created if it doesnt exist)\nIMPORTANT: the tsv must be in the following format:\nstudy_accession\tsample_accession\trun_accession\ttax_id\tscientific_name\tfastq_ftp\tsubmitted_ftp\tsra_ftp\nThe format is available for download from the EBI website\n')
    tsvs = os.listdir(f'{PATH}{SLASH}output')
    [print(f'{i}\t{v}') for i,v in enumerate(tsvs) if os.path.isdir(f'{PATH}{SLASH}output{SLASH}{v}')]
    tsv = input()
    try:
      tsv = tsvs[int(tsv)]
      break
    except (IndexError, AttributeError, TypeError, ValueError):
      tsv_path = f"{PATH}{SLASH}{tsv}.tsv"
      try:
        assert os.path.exists(tsv_path)
        break
      except AssertionError:
        try:
          assert os.path.exists(tsv_path.replace('.tsv','.txt'))
          break
        except AssertionError: print("Please enter a valid TSV file name!\n")
      
  dir_checker('output/' + tsv)
  dir_checker('output/' + tsv + '/' + 'stderr')
  dir_checker('output/' + tsv + '/' + 'stdout')
  new_tsv_path = f'{PATH}{SLASH}output{SLASH}{tsv}{SLASH}{tsv}.tsv'
  try:
    shutil.move(src = tsv_path, dst = new_tsv_path)
  except (FileNotFoundError, NameError):
    pass
  assert os.path.exists(new_tsv_path), "TSV file does not exist!"
  return new_tsv_path, tsv


def main():
  """
  Manage the program
  
  * Asks the user to choose workflow.

  * Asks for the tsv name or a choice of previously used tsvs.

  * Users can either access the config menu, handled by change_menu function or continue.

  * If they continue, the users have to choose which TSV lines they wish to process.

  * main runs csh_maker on each line and saves the settings at the end of runtime.

  """
  dir_checker('output')
  dir_checker('src')
  dir_checker('src/mtDNA_ref')
  dir_checker('src/genbank')
  dir_checker('src/samples')

  if not os.path.exists(f'{PATH}{SLASH}config{SLASH}default'):
    defaults()
  print(f'{DASH}\nWelcome to the semi-automatic .csh maker\n{DASH}')
  workflows = os.listdir(f'{PATH}{SLASH}config')
  while True:
    print('Please choose a pre-configured workflow\n')
    [print(f'{i}\t{v}') for i,v in enumerate(workflows) if (v not in  ('global','current')) and (os.path.isdir(f'{PATH}{SLASH}config{SLASH}{v}'))]
    workflow = input()
    try:
      workflow = workflows[int(workflow)] 
      save_settings(source = workflow, output = 'current')
      break
    except (IndexError, ValueError, TypeError):
      print('Please enter the number next to the configuration name\n\n')

  tsv_path, tsv_name = tsv_handler()
  tsv_df = pd.read_csv(tsv_path, delimiter = '\t')
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
  print(f'Choose which lines of the tsv to process with the current pipeline.\nEnter a single row number (i.e 2) or a range (i.e 2-3)\nEnter all to work on all rows\nTo work on specific rows, enter a list of rows seperated by a SINGLE space (1 2 3 4)\n')
  wanted_cols = ['index', 'scientific_name', 'run_accession', 'sample_title', 'library_strategy', 'library_selection', 'read_count']
  colnames = ['Index', 'ID', 'Sample accession', 'Sample', 'Library Strategy', 'Library Selection', 'Read Count']
  actual_cols = [i for i in wanted_cols if i in tsv_df.columns]
  sliced_df = tsv_df.loc[:, actual_cols]
  print(tabulate(sliced_df, headers = 'keys', tablefmt = 'psql'))
  while True:
    work_on = input()
    if ' ' in work_on:
      try:
        work_on = [int(i) for i in work_on.split(' ')]
      except ValueError:
        print('The list of rows to work on must only consist of integers seperated by a SINGLE SPACE!\nTry again\n')
        continue
    break
  csh_maker(work_on, tsv_df, tsv_name)
  save_settings(source = 'current', output = 'last_used')
  print('Enter \"bash master.csh\" into the console to run the commands\nIMPORTANT: Each sample will be executed simultaneously in qsub you run this command, make sure you dont run too many at once!\n')


if __name__ == '__main__':
  main()
