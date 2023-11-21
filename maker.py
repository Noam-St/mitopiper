"""
#########################################################################

This program creates a custom csh file designed to be executed with qsub.
The program is designed to be configurable, the user can customize the processing pipeline,
the settings for each of the pipeline programs used, and the paths of the programs utallized.

#########################################################################
"""
__author__ = 'Noam Shtolz'


import sys, shutil, os, math, json, re
from textwrap import wrap
import consts
import pandas as pd
from Bio import Entrez, SeqIO
from tabulate import tabulate
from importlib import reload
import random
from ast import literal_eval

reload(consts)

DASH = '##############################################'
PATH = os.getcwd()
if sys.platform != 'win32':
    SLASH = '/'
else:
    SLASH = '\\'
try:
  # Read in the final.csv file, which contains the organism names and their corresponding mt reference
  org_df = pd.read_csv(os.path.join(PATH, 'src', 'final.csv'))
  org_df.organism = org_df.organism.str.replace(' ','_')
  org_df = org_df.set_index('organism')
  org_df_exists = True
except FileNotFoundError:
  print('final.csv is required in the src folder for the script to automatically download mt references')
  org_df_exists = False

def check_input(n, options, special = None, allow = None):
  """
  Recieve a user input that is supposed to be a number within a certain range.
  Make sure the input is convertable to int and within the range.

  Parameters
  ----------
  n : str
    The user input
  options : list
    The list of options to choose from
  special : list
    A list of special options that are not in the options list
  allow : list
    A list of allowed inputs that are not in the options list
  
  Returns
  -------
  bool
    True if the input is valid, False otherwise
  """
  # If the user enters an allowed input, the program will return True
  if allow != None:
    if n in allow:
      return True
  try:
    n = int(n)
  # If the user enters a non-integer, the program will ask them to enter a valid option
  except ValueError:
    print(f'Error: {n} is not a number!\nPlease only choose one of the listed options!\n')
    return False
  try: options[n]
  except IndexError:
    print(f"Error: the number {n} is out of range!\nPlease only choose one of the listed options!\n")
    return False
  if special != None:
    if options[n] in special:
      print(f'Error: {n} is not a valid option!\nPlease only choose one of the listed options!\n')
      return False
  return True

def input_to_list(user_in : str, list_size : int) -> list:
  """
  Converts user input to a list of numbers
  """
  try:
    user_in = int(user_in)
    if user_in > list_size:
      print("Error:\tRow entered is out of bounds! Running on last row\n\n")
      user_in = list_size - 1
    user_in = list(range(user_in, user_in + 1))
  except ValueError:
    if '-' in user_in: #Input is a range of numbers (x-y)
      user_in = [int(i) for i in user_in.split('-')]
      start = user_in[0]
      end = user_in[1]
      user_in = list(range(start, end))
    elif user_in == 'all': #Input is all of the list
      start = 0
      end = list_size
      user_in = list(range(start, end))
    elif ' ' in user_in: #Input is a list of numbers (indicating specific indices)
      try:
        user_in = [int(i) for i in user_in.split(' ')]
      except ValueError:
        print('The list of rows to work on must only consist of integers seperated by a SINGLE SPACE!\nTry again\n')
        return []
    elif ',' in user_in:
      try:
        user_in = [int(i) for i in user_in.split(',')]
      except ValueError:
        print('The list of rows to work on must only consist of integers seperated by a SINGLE SPACE!\nTry again\n')
        return []
    else:
      print('The list of rows to work on must only consist of integers seperated by a SINGLE SPACE!\nTry again\n')
      return []
  # If the user enters a list longer than the total length of the list, it will be truncated
  if len(user_in) > list_size:
    user_in = user_in[:list_size]
    print('Warning: list length is longer than the total number of rows in the tsv. Truncating list\n\n')
  return user_in

def pprint(string):
  """
  Pretty print a string by adding newlines and sorrounding it with dashes
  """
  # If the string contains newlines, find the longest line and use that as the length of the string
  if '\n' in string:
    lens = [len(i) for i in string.split('\n')]
    string_len = max(lens)
  else:
    string_len = len(string)
  # Add dashes to the edges of the string
  dashes = ''.join(['#' for i in range(string_len)])
  print(f'{dashes}\n{string}\n{dashes}\n')

def file_reader(filename : str, directory : str = 'current') -> dict: 
  """
  Read .json files in the format key:value from a dictionary

  Parameters
  ----------
  filename : str
    PATH to txt file in the expected format
  Returns
  -------
  format_dict : dict
    Dictionary with the left side of the equal sign as key and right side as value
  """
  # If the directory argument is empty, read the file directly from the config folder TODO(Noam): Check if this is necessary...
  if directory == '':
    with open(os.path.join(PATH, 'config', filename + '.json'), 'r') as file:
      format_dict = json.load(file)
  # If the directory argument is not empty, read the file from the config folder, then the directory folder
  else:
    with open(os.path.join(PATH, 'config', directory, filename + '.json'), 'r') as file:
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
  # If the directory argument is empty, create the file directly in the config folder #TODO(Noam): Again, check if this is necessary...
  if directory == '':
    filepath = os.path.join(PATH, 'config', filename + '.json')
  # If the directory argument is not empty, create the file in the config folder, then the directory folder
  else:
    filepath = os.path.join(PATH, 'config', directory, filename + '.json')
  with open(filepath, 'w') as file:
    json.dump(format_dict, sort_keys = True, fp = file, indent = 4)

def validate_configs():
  """
  Iterate over all config folders, make sure that all the keys in each .json file match the keys within the default config, if it doesn't add the new key and re-save the config
  """
  resave = False
  # List of config files to check
  CONFIG_NAMES = ['settings', 'pipeline', 'general_settings']
  # Iterate over all config folders
  for dir in os.listdir(os.path.join(PATH, 'config')):
    # Ignore the default and global folders
    if dir not in ['default', 'global']:
      # Iterate over all files in the current config folder
      for file in os.listdir(os.path.join(PATH, 'config', dir)):
        # Remove the .json extension from the file name
        file = file.replace('.json','')
        # If the file is one of the config files to check, read it and the default config, then compare the keys with the default config
        if file in CONFIG_NAMES:
          file_dict = file_reader(file, dir)
          def_dict = file_reader(file, 'default')
          for key in def_dict:
            # If the key is not in the file, add it and re-save the file
            if key not in file_dict:
              file_dict[key] = def_dict[key]
              resave = True
          # TODO(Noam): Add a check for the "order" parameter in each config file. 
          if file == 'settings':
            rnd_stage = random.choice([i for i in file_dict.keys() if file_dict[i] != {}])
            rnd_stage_param = random.choice(list(file_dict[rnd_stage].keys()))
            if 'order' not in file_dict[rnd_stage][rnd_stage_param].keys():
              resave = True
              for stage in file_dict.keys():
                file_dict[stage] = assign_order(file_dict[stage])
        if resave:
          print('Resaved successfully')
          file_creator(file, file_dict, dir)
          resave = False

def validate_workflow(tsv_dir):
  """
  Validate the workflow by checking that all the stages in the workflow are valid and make sense
  This will be invoked before the user is asked whether to change the workflow  
  """
  settings = file_reader('settings')
  pipeline = file_reader('pipeline')
  run_order = file_reader('run_order')
  # Check whether there are files in the tsv_dir/fastq directory
  if (not os.path.exists(os.path.join(tsv_dir, 'fastq'))) or (len(os.listdir(os.path.join(tsv_dir, 'fastq'))) == 0):
    fastq_empty = True
  else:
    fastq_empty = False
  # Check whether the run_order contains wget
  if run_order[0] != 'wget' and fastq_empty:
    print('Warning: No fastq files found in the fastq directory and wget is not in the run order. This will result in the workflow not working properly. Please add wget to the run order or add fastq files to the fastq directory')
    print('Wget must be the first stage in the run order!')

def guess_threads(settings, pipeline):
  """
  Attempt to guess the number of threads to use for the current run based on the maximum number of threads specified for any of the active steps in the settings file, default to 4
  """
  threads = 4
  thread_terms = ['threads', 'cores', 'cpus', 'cpu', 'thread', 'core']
  thread_values = []
  active_stages = [i for i in pipeline.keys() if pipeline[i]['switch']]

  for stage in active_stages:
    for param in settings[stage].keys():
      if settings[stage][param]['switch'] and any(term in param.lower() for term in thread_terms):
        value = settings[stage][param]['value']
        if isinstance(value, int):
          thread_values.append(value)
        elif isinstance(value, str):
          if value.isdigit():
            thread_values.append(int(value))
          elif value.lower() == 'max':
            thread_values.append(8)
  return threads

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

def tabulate_dict(dicty: dict, settings_mode = False) -> tabulate:
  """
  Recieve a dictionary and return a pretty-printed table where each column is a key and the values are spread into rows.
  """
  if not settings_mode: temp_df = pd.DataFrame.from_dict(dicty, orient = 'index', columns = ['Value'])
  else: temp_df = pd.DataFrame.from_dict(dicty)
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

def save_settings(source: str, output: str) -> None:
  """
  Transfer settings from source to output
  """
  dir_checker('config/' + output)
  for jsonfile in os.listdir(os.path.join(PATH, 'config', source)):
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
  dir_checker('config/global')
  if not os.path.exists(os.path.join(PATH, 'config', 'global', 'refpaths.json')):
    file_creator('refpaths', {}, directory = 'global')
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
      fasta = os.path.join(PATH, 'src', 'mtDNA_ref', org_name, f'{org_name}.fasta')
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
      chrominfo = ''
      nuc_fasta = ''
    elif choice == 2:
      print('Please enter the full path to the main reference fasta\n')
      fasta = input()
      print('Please enter the full path to the nuclear reference fasta\nNOTE: This is only required if you wish to map to both the nuclear and mitochondrial genomes separately\n')
      nuc_fasta = input()
      add_to_csh = ''
      print('Please enter path to gtf file for this organism\nEnter nothing if annotations are not required.\n')
      gtf = input()
      print('Please enter path to ChromInfo file for this organism\nEnter nothing if chromosome info is not required.\n')
      chrominfo = input()
    else:
      print('Please enter 1 or 2 only')
 
    refs[org_name] = {}
    refs[org_name]['refpath'] = fasta.replace(r'"', '')
    refs[org_name]['nuc_refpath'] = nuc_fasta.replace(r'"', '')
    refs[org_name]['gtfpath'] = gtf.replace(r'"', '')
    refs[org_name]['chrominfopath'] = chrominfo.replace(r'"', '')
    file_creator('refpaths', refs, 'global') #NOTE(Noam) Currently set to global, think about changing to local pre-configurations to enable multiple ref files for the same organism depending on preset.
    return add_to_csh


def change_general_settings() -> None:
  """
  Change general settings defined on a per-run basis.
  """
  # Read the general settings file
  general_settings = file_reader('general_settings')
  # Create a list of the options
  options = list(general_settings.keys())
  while True:
    temp_setting = {'Name' : options, 'Value' : [general_settings[v][0] for v in options]}
    print(tabulate_dict(temp_setting, settings_mode = True))
    print('Which setting would you like to change?\n')
    print('Type exit to exit\n')
    choice = input()
    if choice == 'exit': break
    valid = check_input(choice, options)
    if not valid: continue
    choice = int(choice)
    chosen = general_settings[options[choice]]
    # If the setting is a boolean, change it to the opposite
    if options[choice] == 'ignore_paired_end':
      chosen[0] = not chosen[0]
    else:
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

def change_order(stage):
  """
  Change the order of the parameters for a given stage.
  Save the new order to the settings file.

  Parameters
  ----------
  stage : str
    The stage to change the order of.
  settings : dict
    The settings dictionary.
  """
  # Read the settings file
  settings = file_reader('settings')
  stage_params = settings[stage]
  print(f'To change the order in which the parameters appear for {stage}, enter the number next to it.\nNOTE: you can also reorder quickly by entering a sequence of numbers seperated by space (i.e: \"0 3 5 ...\")\n')
  stage_df = pd.DataFrame.from_dict(stage_params, orient = 'index').sort_values('order')[['order', 'description']]
  # Create a sorted list of parameters according to their order
  order = list(stage_df.index)
  print(tabulate(stage_df, headers = 'keys', tablefmt = 'psql'))
  while True:
    choice = input()
    if choice == 'exit': break
    choice = input_to_list(choice, len(order))
    if choice == []: continue
    else: break
  if not len(choice) > 1:
    choice = int(choice[0])
    print(f'Who would you like to replace {order[choice]} with?\n')
    while True:
      where = input()
      valid = check_input(where, order)
      if not valid: continue
      else: break
    where = int(where)
    temp = order[where]
    order[where] = order[choice]
    order[choice] = temp
  else:
    order = [order[i] for i in choice] + [v for v in order if v not in [order[i] for i in choice]]
  # Apply the order to the settings
  for i,v in enumerate(order):
    stage_params[v]['order'] = i
  for param in stage_params:
    if param not in order:
      stage_params[param]['switch'] = False
  settings[stage] = stage_params
  file_creator('settings', settings)

def assign_order(stage_dict):
  """
  Assign a random order to the parameters in a given stage.
  """
  for i,v in enumerate(stage_dict.keys()):
    stage_dict[v]['order'] = i
  return stage_dict

def change_pipeline_settings() -> None:
  """
  Settings change menu, in this menu you can reconfigure any of the setting options for any of the pipeline steps
  """
  # Grab the current settings
  settings = file_reader('settings')
  pipeline = file_reader('pipeline')
  options = list(settings.keys())
  order = file_reader('run_order')
  print(f"Welcome to the settings menu\n")
  while True:
    # Print the current settings that are turned on
    [print(f'{i}\t{v}') for i,v in enumerate(options) if pipeline[v]['switch']]
    print('To change settings for a specific step, or their order, enter the number next to it\n')
    print('Type exit to exit\n')
    choice = input()
    valid = check_input(choice, options, special = [v for v in options if pipeline[v]['switch'] == False], allow = ['exit'])
    if not valid: continue
    if choice == 'exit': break
    choice = int(choice)
    chosen_stage = options[choice]
    configs = list(settings[chosen_stage].keys())
    # Create a temporary dataframe to print the settings
    temp_config_df = {
      'Name':configs,
      'Description':[i['description'] for i in settings[chosen_stage].values()],
      'Value':[i['value'] if i['value'] != '' else i['switch'] for i in settings[chosen_stage].values()]
      }
    print(tabulate_dict(temp_config_df, settings_mode=True))
    # If there are no settings for this stage, print a message and return to the main menu
    if settings[chosen_stage] == {}:
      print(f'It seems no settings were added to {chosen_stage} yet\nto add some return first to the main menu by typing exit\n')
    # If there are settings, print them and ask the user which one they would like to change
    else:
      print(f'Which setting would you like to modify in {chosen_stage}?\nEnter the number next to the setting\nNOTE: To reorder the setting\'s parameters, enter r\n')
      while True:
        change = input()
        valid = check_input(change, configs, allow = ['r'])
        if not valid: continue
        # If the user chose to reorder the settings, call the change_order function and return to the main menu
        if change == 'r': change_order(chosen_stage); return 
        break
      change = int(change)
      # Grab the chosen setting, its current value and its switch (on/off)
      chosen_setting = configs[change]
      setting_value = settings[chosen_stage][chosen_setting]['value']
      setting_switch = settings[chosen_stage][chosen_setting]['switch']
      # If the setting is a boolean, change it to the opposite
      if setting_value == '':
        settings[chosen_stage][chosen_setting]['switch'] = not setting_switch
        print(f'\nBoolean setting switched to {not bool(setting_switch)}\n')
      else:
        print(f'You chose: {chosen_setting}')
        print(f'Please choose a new value, current value: {setting_value}\n')
        print('The value may include several Keywords:\n_gtf, _chrominfo, _fasta, _fastq, _fastq_pre, _fastq_sample_path_no_ext, _fastq_dir, _fastq_pre_dir, _bam, _sam, _org_name\nThese keywords will be replaced with appropriate paths\n')
        print(f'Additionally, you may add the combination of output_command_name or dir_command_name to include the output of another command in the setting value\nNOTE: you are adding a setting to {chosen_stage}. You can grab the output of the following commands with the current run order')
        try: print('\t'.join(order[:order.index(chosen_stage)]))
        except ValueError:
          print('This step is not active!')
          continue
        print(f'EXAMPLE: dir_{order[order.index(chosen_stage)-1]}\n')
        new_val = input()
        settings[chosen_stage][chosen_setting]['value'] = new_val
        print('\nSetting changed successfuly!\n')
    
  file_creator('settings', settings)


def modify_pipeline() -> None:
  """
  Pipeline menu, in this menu you can remove/add commands to the pipeline or reorder it
  """
  running = True
  # Grab the current pipeline and run order
  pipeline = file_reader('pipeline')
  order = file_reader('run_order')
  while running:
    print("Welcome to the pipeline customizer\nThe current pipeline order is:\n")
    # Print the current pipeline order
    [print(f'\t{v}') for v in order]
    print('To change the pipeline enter 1\nTo change the pipeline order enter 2\nTo exit this menu enter 0\n')
    running = input()
    # Check if the input is valid
    valid = check_input(running, [0,1,2])
    if not valid: continue
    running = int(running)
    # IF the user chose to activate/inactivate specific pipeline steps
    if running == 1:
      pipeline_changes = {i:k for i,k in enumerate(pipeline.keys())}
      print('To remove or add a feature, simply enter its number\nIMPORTANT: if you remove a command and add it later, it will appended in the last position!\n')
      print('ACTIVE COMMANDS:\n')
      [print(f'{i} {k}') for i,k in zip(pipeline_changes.keys(), pipeline_changes.values()) if pipeline[k]['switch']]
      print('INACTIVE COMMANDS:\n')
      [print(f'{i} {k}') for i,k in zip(pipeline_changes.keys(), pipeline_changes.values()) if not pipeline[k]['switch']]
      while True:
        to_change = input()
        to_change = input_to_list(to_change, len(pipeline_changes))
        if to_change == []: continue
        else: break
      for i in to_change:
        step = pipeline_changes[i]
        pipeline[step]['switch'] = not pipeline[step]['switch']
        if i not in order:
          order.append(step)
    # If the user chose to reorder the pipeline
    elif running == 2:
      print('The current run order is:\n')
      [print(f'{i}\t{v}') for i,v in enumerate(order)]
      print('To move a command enter the number next to it\nNOTE: you can also reorder quickly by entering a sequence of numbers seperated by space (i.e: \"0 3 5 ...\")\n')
      while True:
        to_move = input()
        to_move = input_to_list(to_move, len(order))
        if to_move == []: continue
        else: break
      if not len(to_move) > 1:
        to_move = int(to_move[0])
        print(f'Who would you like to replace {order[to_move]} with?\n')
        while True:
          where = input()
          valid = check_input(where, order)
          if not valid: continue
          else: break
        where = int(where)
        temp = order[where]
        order[where] = order[to_move]
        order[to_move] = temp
      else:
        order = [order[i] for i in to_move]
        for i in pipeline.keys():
          if i not in order:
            pipeline[i]['switch'] = False
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
  [print(f'{i}\t{k}') for i,k in enumerate(setting_keys)]
  print(f'{DASH}\nPlease type the pipeline step to add a setting to\n{DASH}')
  step = input()
  valid = check_input(step, setting_keys)
  if not valid: return
  step = int(step)
  chosen_step = setting_keys[step]
  print(f"You chose {chosen_step}\nPlease enter the setting name\n")
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
  print(f'Please enter a value for {new_name} (if the setting is an on/off switch, enter nothing!)\nYou can also enter a keyword which will be replaced with the appropriate file path, Keywords: _gtf, _fasta, _fastq, _fastq_pre, _fastq_dir, _fastq_dir_pre, _bam and _sam, _chrominfo, _command_output and _command_input will be replaced with appropriate paths\nYou can use _if_SE or _if_PE if the setting is an on/off flag that only appears when SE or PE respectively\n')
  print(f'Additionaly, you may add the combination of output_command_name or dir_command_name to include the output or directory of another command in the setting value\nNOTE: you are adding a setting to {chosen_step}. You can grab the output of the following commands with the current run order\n')
  try: print('\t'.join(order[:order.index(chosen_step)]))
  except ValueError: print('Current command isnt active!')
  value = input()
  settings[chosen_step][new_name]['value']  = value
  settings[chosen_step][new_name]['switch'] = 1
  # Define the order of the setting to be last in the list
  settings[chosen_step][new_name]['order'] = len(settings[chosen_step].keys())
  file_creator('settings', settings) 
  print(f'Setting added successfuly!\nName\tKey\tValue\n{new_name}\t{key}\t{value}\n')

def add_new_pipeline_step() -> None:
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
  'structure':'Command structure, enter nothing for default format: key settings input\nIMPORTANT: A non-default structure can only include the following keywords: sample, key, settings, gtf, fasta, fastq, fastq_pre fastq_dir, fastq_pre_dir, sam, bam, command_output, command_input and output_pipe seperated by a SINGLE space'
}

  pipeline = file_reader('pipeline')
  order = file_reader('run_order')
  settings = file_reader('settings')
  print('To begin defining a new pipeline step, first enter its name\nIMPORTANT: the name must be different from the current pipeline steps\n')
  name = input()
  # Remove spaces from input
  name = name.replace(' ','_')
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

def modify_chosen_step(step_name):
  """
  Change the chosen parameter in the chosen step to the chosen value

  Parameters
  ----------
  step_name : str
    The name of the step to modify
  """
  pipeline = file_reader('pipeline')
  step = pipeline[step_name]
  step_key_list = [i for i in list(step.keys()) if i not in ('switch', 'path')]
  while True:
    print(f'Which parameter would you like to change in {step_name}?\nType \"exit\" to return to main menu\n')
    temp_df = pd.DataFrame([step_key_list, [step[k] for k in step_key_list]]).T
    print(tabulate(temp_df, tablefmt= 'psql'))
    chosen_param = input()
    if check_input(chosen_param, options = list(step.keys()), allow='exit'):
      if chosen_param == 'exit': break
      chosen_param_name = step_key_list[int(chosen_param)]
      print(f'Please enter the new value for {chosen_param_name}')
      new_value = input()
      pipeline[step_name][chosen_param_name] = new_value
      file_creator('pipeline', pipeline)
      print(f'{chosen_param_name} in {step_name} changed to {new_value} successfuly!\n')
      break
    else:
      print('Please enter a valid input!\n')
      continue

def modify_pipeline_step():
  """
  Modify pipeline step menu
  The user chooses a pipeline step to modify, the chosen pipeline step is passed to modify_chosen_step function.
  """
  pipeline = file_reader('pipeline')
  pipe_keys = list(pipeline.keys())
  while True:
    print("Pleast choose a pipeline step to modify")
    [print(f'{i}\t{k}') for i, k in enumerate(pipe_keys)]
    print(f'Type \"exit\" to return to main menu')
    step_ind = input()
    if check_input(step_ind, options = pipe_keys, allow='exit'):
      if step_ind == 'exit': break
      modify_chosen_step(pipe_keys[int(step_ind)])
      break
    else:
      print('Please enter a valid input!\n')
      continue

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
    pprint("Please choose an organism to add reference fasta or gtf path to\nEnter exit to return to main menu")
    [print(f'{i}\t{v}') for i,v in enumerate(keys)]
    org = input()
    if org == 'exit': break
    try: org = keys[int(org)]
    except IndexError:
      print("Please enter an organism from the list!\n")
      continue
    print(f'Current reference pathes for {org}:\n')
    [print(f'{k.replace("path"," path")}\n{v}\n') for k,v in ref[org].items()]
    print('Would you like to change gtf, fasta, chromInfo all?\n')
    work_on = ['ref', 'nuc_ref', 'gtf', 'chrominfo', 'all']
    [print(f'{i}\t{v}') for i,v in enumerate(work_on)]
    choice = int(input())
    
    try: work_on = work_on[int(choice)]
    except IndexError:
      print('Please enter one of the possible options!\n')
      continue
    if work_on == 'all':
      work_on = ['ref', 'nuc_ref', 'gtf', 'chrominfo']
    else:
      work_on = [work_on]
    for i in work_on:
      print(f"Please enter the {i} absolute path!\n")
      path = input()
      i += 'path'
      ref[org][i] = path
    break
  file_creator(filename = 'refpaths', format_dict = ref, directory = 'global') 

def copy_settings_and_pipeline_from_one_config_to_another() -> None:
  """
  Copy all existing pipeline steps and settings from an existing pipeline to another existing pipeline without running over them.
  """
  all_configs = [i for i in os.listdir(os.path.join(PATH, 'config')) if i not in ['default', 'current', 'global', 'last_used']]

  print('You chose to copy settings and pipeline steps from one config to another.\nNOTE: this will not overwrite existing settings and steps!\n\nEnter "exit" to return to main menu')
  while True:
    print('Please enter the name of the pipeline you want to copy the steps and settings from!\n')
    [print(f'{i}\t{v}') for i,v in enumerate(all_configs)]
    try:
      copy_from = input()
      if copy_from == 'exit': return
      print('Great! Now enter the name of the config you want to copy the steps and settings to!\n')
      copy_to = input()
      copy_from = int(copy_from)
      copy_to = int(copy_to)
      copy_from = all_configs[copy_from]
      copy_to = all_configs[copy_to]
      break
    except (IndexError, ValueError):
      print('Please enter valid configs from the list!\n')
      continue

  from_pipeline = file_reader('pipeline', copy_from)
  from_settings = file_reader('settings', copy_from)
  to_pipeline = file_reader('pipeline', copy_to)
  to_settings = file_reader('settings', copy_to)

  # Only copy pipeline steps that are not already in the pipeline and settings that are not already in the settings
  for k,v in from_pipeline.items():
    if k not in to_pipeline.keys():
      to_pipeline[k] = v
  for k,v in from_settings.items():
    if k not in to_settings.keys():
      to_settings[k] = v
  file_creator('pipeline', to_pipeline, copy_to)
  file_creator('settings', to_settings, copy_to)
  print(f'\nSettings and pipeline steps copied to {copy_to} from {copy_from} successfuly!\n')


def change_menu():
  """
  Main reconfiguration menu, sends to sub-menus, after exiting sub-menus, returns
  the user here.
  """
  choice = True
  pprint('Welcome to the configuration menu!')
  # Menu options as a list of functions, the user chooses a number to choose a function
  menu = [lambda : 0, change_general_settings, change_pipeline_settings, modify_pipeline, add_new_pipeline_setting, add_new_pipeline_step, modify_pipeline_step, manually_change_ref_or_gtf, save_current_settings, copy_settings_and_pipeline_from_one_config_to_another]
  while choice:
    [print(f"Enter {i+1}\tto {x.__name__}\n") for i,x in enumerate(menu[1:])]
    print("Enter 0\tto exit\n")
    choice = input()
    valid = check_input(choice, menu)
    if not valid: continue
    else: choice = int(choice)
    menu[choice]()

def csh_commands(i, row, s_id, settings, pipeline, order, fasta, gtf, chrominfo, fasta_nuc, PAIRED, org, tsv_dir, SAMPLE_LEVEL):
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
    try:
      org = org[0][0] + org[1][0]
    except IndexError:
      org = org[0][0]
  org = org.upper()
  #if SAMPLE_LEVEL: s_id = row.sample_accession
  #else: s_id = row.run_accession

  sample = f'{s_id}_{i}_{org}'
  command_outputs = {f'output_{command}':f'{tsv_dir}{SLASH}{pipeline[command]["outdir"]}{SLASH}{sample}{pipeline[command]["output_ext"]}' for command in pipeline.keys()}
  command_dirs = {f'dir_{command}':f'{tsv_dir}{SLASH}{pipeline[command]["outdir"]}' for command in pipeline.keys()}
  to_f_string = [
  '_fasta', '_fasta_dir', '_fasta_nuc', '_fasta_nuc_dir', '_gtf', '_command_input', '_command_output', '_command_output_dir_path', '_command_input_dir_path','_fastq', '_fastq_pre', '_sample', '_fastq_sample_path_no_ext','_fastq_dir', '_fastq_pre_dir', '_bam', '_sam', '_STAR_out', '_fastq1', '_fastq2', '_org_name','_STAR_final_log', '_chrominfo'
  ]
  if not PAIRED and SAMPLE_LEVEL:
    row.fastq_ftp = [i.split(';')[0] for i in row.fastq_ftp]
  elif not PAIRED and not SAMPLE_LEVEL:
    row.fastq_ftp = row.fastq_ftp.split(';')[0]
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
    fastq_sample_path_no_ext = f'{fastq_pre_dir}{SLASH}{sample}'
    if_SE = False if PAIRED else True
    if_PE = True if PAIRED else False

    bam =  f'{tsv_dir}{SLASH}{pipeline["bwa-align"]["outdir"]}{SLASH}{sample}*.bam'
    sam =  f'{tsv_dir}{SLASH}{pipeline["bwa-align"]["outdir"]}{SLASH}{sample}*.sam'
    STAR_out = f'{tsv_dir}{SLASH}{pipeline["bwa-align"]["outdir"]}{SLASH}{sample}'
    STAR_final_log = f'{tsv_dir}{SLASH}{pipeline["bwa-align"]["outdir"]}{SLASH}{sample}Log.final.out'
    fasta_dir = os.path.dirname(fasta).replace(r'"', '')
    fasta_nuc_dir = os.path.dirname(fasta_nuc).replace(r'"', '')
    command_input = f'{tsv_dir}{SLASH}{command_inputdir}{SLASH}{sample}*{input_ext}'
    command_output = f'{tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}{output_ext}'
    command_output_dir_path = f'{tsv_dir}{SLASH}{command_outdir}'
    command_input_dir_path = f'{tsv_dir}{SLASH}{command_outdir}'

    output_pipe = f' > {command_output}\n'
    for i,v in settings[command].items():
      if v['value'] == 'fasta': v['value'] = fasta
      if v['value'] == 'fasta_nuc': v['value'] = fasta_nuc
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
        elif cur_val == '_if_SE':
          settings[command][k]['switch'] = False if PAIRED else True
          cur_val = ''
        elif cur_val == '_if_PE':
          settings[command][k]['switch'] = True if PAIRED else False
          cur_val = ''
        settings[command][k]['value'] = cur_val
      cur_keys = sorted(list(settings[command].keys()), key = lambda i: settings[command][i]['order'],reverse=False)
      command_settings = ' '.join([\
        f'{settings[command][k]["key"]} {settings[command][k]["value"]} ' for k in cur_keys if settings[command][k]['switch'] == True])
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
      if command == 'wget':
        if not SAMPLE_LEVEL:
          ftps = row.fastq_ftp.split(';')
          for i in range(2):
            csh_str+=f'{command_key} {command_settings} -O {tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}_{i+1}.fastq.gz {ftps[i]}\n'
            # Wait 30s between each wget command to avoid overloading the server
            csh_str+=f'sleep 30\n'
        else: # If SAMPLE_LEVEL
          for i in range(2):
            for run_ind, ftp in enumerate(row.fastq_ftp):
              ftps = ftp.split(';')
              csh_str+=f'{command_key} {command_settings} -O {tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}_{i+1}_{run_ind}.fastq.gz {ftps[i]}\n'
              # Wait 30s between each wget command to avoid overloading the server
              csh_str+=f'sleep 30\n'
            # Concat all fastq files for each sample for forward and reverse reads and remove the individual files afterwards
            csh_str+=f'cat {tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}_{i+1}_*.fastq.gz > {tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}_{i+1}.fastq.gz\n'
            csh_str+=f'rm {tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}_{i+1}_*.fastq.gz\n'
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
        if not SAMPLE_LEVEL:
          csh_str+=f'{command_key} {command_settings} -O {tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}.fastq.gz {row.fastq_ftp}\n'
          # Wait 30s between each wget command to avoid overloading the server
          csh_str+=f'sleep 30\n'
        else: # If SAMPLE_LEVEL
          for run_ind, ftp in enumerate(row.fastq_ftp):
            csh_str+=f'{command_key} {command_settings} -O {tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}_{run_ind}.fastq.gz {ftp}\n'
            # Wait 30s between each wget command to avoid overloading the server
            csh_str+=f'sleep 30\n'
          # Concat all fastq files for each sample and remove the individual files afterwards
          csh_str+=f'cat {tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}_*.fastq.gz > {tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}.fastq.gz\n'
          csh_str+=f'rm {tsv_dir}{SLASH}{command_outdir}{SLASH}{sample}_*.fastq.gz\n'
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

def sample_lvl_csh_maker(work_on, tsv_df, tsv_name, queue_mode):
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
  settings = file_reader('settings')
  pipeline = file_reader('pipeline')
  order = file_reader('run_order')
  [dir_checker('output' + '/' + tsv_name + '/' + i['outdir']) for i in pipeline.values() if i['outdir'] != '']
  tsvs = []
  bad_orgs = []
  samples_df = tsv_df.iloc[work_on, :]
  threads = guess_threads(settings, pipeline)
  for i,row in samples_df.iterrows():
    settings = file_reader('settings')
    pipeline = file_reader('pipeline')
    # Check if sample level metadata file
    if type(row.run_accession) == list: SAMPLE_LEVEL = True
    else: SAMPLE_LEVEL = False
    org_name = row.scientific_name
    if org_name in bad_orgs:
      print(f'Warning:\tSkipped organism {org_name} because it does not exist in final.csv!\n\n')
      continue
    csh_add_ref = ref_main(org_name)
    if csh_add_ref == False:
      bad_orgs.append(org_name)
      continue
    
    if SAMPLE_LEVEL:
      fastq_ftp = row.fastq_ftp[0]
      accs = row.sample_accession
    else:
      fastq_ftp = row.fastq_ftp
      accs = row.run_accession
    if len(fastq_ftp.split(';')) > 1:
      PAIRED = True
    else:
      PAIRED = False
    try:
      if general_settings['ignore_paired_end'][0]:
        # If ignore_paired_end is True, then ignore the paired end status of the sample
        PAIRED = False
    except KeyError:
      pass
    with open(f'{tsv_dir}{SLASH}run{i}.csh', 'w') as file:
      if queue_mode == 'sge':
        file.write(f'#!/bin/{general_settings["environment"][0]}\n')
        file.write(f'#$ -o {tsv_dir}{SLASH}stdout\n')
        file.write(f'#$ -e {tsv_dir}{SLASH}stderr\n')
        file.write(f'#$ -cwd\n')
        file.write(f'#$ -N s_{accs}_{i}\n')
        file.write(f'#$ -q {general_settings["qsub_queue"][0]}\n')
        file.write(f'#$ -S /bin/bash\n')
        file.write(f'#$ -V\n')
        file.write(f'#$ -pe shared {general_settings["qsub_threads"][0]}\n')
        file.write(f'#$ -l mem_free={general_settings["qsub_memory"][0]}\n')
      elif queue_mode == 'sbatch':
        file.write(f'#!/bin/{general_settings["environment"][0]}\n')
        file.write(f'#SBATCH --nodes=1\n')
        file.write(f'#SBATCH --ntasks=1\n')
        file.write(f'#SBATCH --job-name=s_{accs}_{i}\n')
        file.write(f'#SBATCH --output={tsv_dir}{SLASH}stdout\n')
        file.write(f'#SBATCH --error={tsv_dir}{SLASH}stderr\n')
        file.write(f'#SBATCH --partition={general_settings["qsub_queue"][0]}\n')
        file.write(f'#SBATCH --cpus-per-task={general_settings["qsub_threads"][0]}\n')
        file.write(f'#SBATCH --mem={general_settings["qsub_memory"][0]}\n')
        file.write(f'#SBATCH --time=1-20:00:00\n')
        file.write(f'#SBATCH --account=shtolz\n')
      file.write(general_settings["free_form"][0])
      file.write(f'echo \"Current organism: {row.scientific_name}\"\n\n')
      if csh_add_ref != '': file.write(csh_add_ref)
      refs = file_reader('refpaths', 'global')
      fasta = refs[org_name]['refpath']
      gtf = refs[org_name]['gtfpath']
      try:
        fasta_nuc = refs[org_name]['nuc_refpath']
      except KeyError:
        fasta_nuc = ''
      try:
        chrominfo = refs[org_name]['chrominfopath']
      except KeyError:
        chrominfo = ''
      to_csh = csh_commands(i, row, accs, settings, pipeline, order, fasta, gtf, chrominfo, fasta_nuc, PAIRED, org_name, tsv_dir, SAMPLE_LEVEL)
      file.write(to_csh)
      tsvs.append(f'{tsv_dir}{SLASH}run{i}.csh')
  return tsvs

def make_master_file(tsvs, queue_mode, qsub_queue):
  """
  """
  tsv_dirs = set([i.split('run')[0] for i in tsvs])
  with open(os.path.join(PATH, 'master.csh'), 'w') as master:
    master.write(f'#!/bin/bash\n')
    if queue_mode == 'sge':
      master.write(f'for file in {" ".join(tsvs)}\ndo\nqsub $file\nsleep 0.5\ndone\n')
    elif queue_mode == 'sbatch':
      master.write(f'for file in {" ".join(tsvs)}\ndo\nsbatch $file\nsleep 0.5\ndone\n')
    #for tsv_dir in tsv_dirs:
    #  master.write(f'mkdir {os.path.join(PATH, tsv_dir, "multiqc")}\n')
    #  master.write(f'multiqc {os.path.join(PATH, tsv_dir)} -o {os.path.join(PATH, tsv_dir, "multiqc")}\n')


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
  dir_checker(os.path.join('config', 'default'))
  dir_checker(os.path.join('config','global'))

  DEFAULT_PIPE = consts.DEFAULT_PIPE
  DEFAULT_SETTINGS = consts.DEFAULT_SETTINGS
  DEFAULT_ORDER = consts.DEFAULT_ORDER
  DEFAULT_MAIN = consts.DEFAULT_MAIN
  DEFAULTS = (DEFAULT_PIPE, DEFAULT_SETTINGS, DEFAULT_ORDER, DEFAULT_MAIN)
  file_creator('refpaths', {}, directory = 'global')
  configs = ['pipeline', 'settings', 'run_order', 'general_settings']
  for i,config in enumerate(configs):
    file_creator(config, DEFAULTS[i], directory = 'default')

def metadata_input() -> list:
  """
  Recieve the requested tsvs from the user and return them as a list of tsv names
  """
  while True:
    print('\nPlease enter the ENA metadata file name (without extension) or the number next to a pre-existing metadata file.\nIf you want to work on several metadata files at once you can do this by entering either a range of numbers separated by a dash (x-y),\nseveral numbers separated by space (x y z), or the word \"all\".\n')
    print(f'IMPORTANT: The metadata file must be in a the same folder as this script!\nIMPORTANT: all the metadata files must be in the following format:\n\nstudy_accession\tsample_accession\trun_accession\ttax_id\tscientific_name\tfastq_ftp\tsubmitted_ftp\tsra_ftp\n\nThe format is available for download from the ENA website\n')
    # List all previously used experiments
    tsvs = sorted(os.listdir(f'{PATH}{SLASH}output'))
    [print(f'{i}\t{v}') for i,v in enumerate(tsvs) if os.path.isdir(f'{PATH}{SLASH}output{SLASH}{v}')]
    tsv = input()
    #If the input is the name of a new TSV file read it
    if re.match(r'\D+', tsv) and tsv != 'all':
      # If file does not include extension, add it
      tsv = tsv.split('.')[0]
      tsv_path = f"{PATH}{SLASH}{tsv}.tsv"
      try:
        assert os.path.exists(tsv_path)
      except AssertionError:
        try: assert os.path.exists(tsv_path.replace('.tsv','.txt'))
        except AssertionError:
          print("Error: File name does not exist!\nPlease enter a valid TSV file name!\n")
          continue
      tsvs.append(tsv)
      tsv = tsvs.index(tsv)
      break
    else:
      break
  tsv_inds = input_to_list(tsv, len(tsvs))
  tsvs = [v for i, v in enumerate(tsvs) if i in tsv_inds]
  return tsvs

def metadata_to_sample_level(tsv_path):
  """
  Receive path to ENA metadata tsv file and convert it from run level to sample level
  Save the new tsv file in the same folder as the original tsv file
  """
  tsv_df = pd.read_csv(tsv_path, sep = '\t', header = 0)
  tsv_df_sample_level = []
  sample_to_runs = {}
  sample_to_links = {}
  for sample in tsv_df.sample_accession.unique():
    tsv_df_sample = {}
    cur_runs = tsv_df[tsv_df.sample_accession == sample].run_accession.values.tolist()
    cur_links = tsv_df[tsv_df.sample_accession == sample].fastq_ftp.values.tolist()
    sample_to_runs[sample] = cur_runs
    sample_to_links[sample] = cur_links
    tsv_df_sample['sample_accession'] = sample
    tsv_df_sample['run_accession'] = cur_runs
    tsv_df_sample['fastq_ftp'] = cur_links
    for col in tsv_df.columns:
      if col not in ['sample_accession', 'run_accession', 'fastq_ftp']:
        tsv_df_sample[col] = tsv_df[tsv_df.sample_accession == sample][col].values[0]
    tsv_df_sample_level.append(tsv_df_sample)
  tsv_df_sample_level = pd.DataFrame.from_dict(tsv_df_sample_level)
  # Replace the tsv file with the sample level tsv file
  tsv_df_sample_level.to_csv(tsv_path, sep = '\t', index = False)
  return tsv_df_sample_level


def metadata_handler(tsv):
  """
  Invoked in cases where the user wants to work on a new metadata file.
  Move the tsv file to the appropriate input folder and rename it to the name of the tsv file.

  """
  new_tsv = tsv.replace(' ', '_')
  tsv_path = os.path.join(PATH, tsv + '.tsv')
  if not os.path.exists(tsv_path):
    tsv_path = tsv_path.replace('.tsv', '.txt')
  dir_checker('output/' + new_tsv)
  dir_checker('output/' + new_tsv + '/' + 'stderr')
  dir_checker('output/' + new_tsv + '/' + 'stdout')
  new_tsv_path = f'{PATH}{SLASH}output{SLASH}{new_tsv}{SLASH}{new_tsv}.tsv'
  try:
    shutil.move(src = tsv_path, dst = new_tsv_path)
    metadata_to_sample_level(new_tsv_path)
  except (FileNotFoundError, NameError):
    # This is a workaround for the case where the user enters a tsv name that does not exist or a tsv that is already in the output folder
    pass
  assert os.path.exists(new_tsv_path), "TSV file does not exist!"
  return new_tsv_path, tsv

def main():
  """
  Run the program
  
  * Asks the user to choose workflow.

  * Asks for the tsv name or a choice of previously used tsvs.

  * Users can either access the config menu, handled by change_menu function or continue.

  * If they continue, the users have to choose which TSV lines they wish to process.

  * main runs sample_lvl_csh_maker on each line and saves the settings at the end of runtime.

  """
  # Check if the required directories exist, if not create them
  dir_checker('output')
  dir_checker('src')
  dir_checker(os.path.join('src', 'mtDNA_ref'))
  dir_checker(os.path.join('src', 'genbank'))

  # Check if the default config exists, if not create it
  if not os.path.exists(f'{PATH}{SLASH}config{SLASH}default'):
    defaults()
  validate_configs()
  pprint(f'Welcome to the semi-automatic .csh maker for pipelines!')
  # Choose a queue mode (sge or sbatch)
  while True:
    print('Please choose a queue mode:\n')
    print('Enter 0 for sge\nEnter 1 for sbatch\n')
    queue_mode = input()
    valid = check_input(queue_mode, [0,1])
    if not valid:
      continue
    queue_mode = int(queue_mode)
    if queue_mode == 0:
      queue_mode = 'sge'
    elif queue_mode == 1:
      queue_mode = 'sbatch'
    break
  tsvs = metadata_input()
  # Check if the user has a pre-configured workflow and if so, ask if they want to use it
  workflows = sorted(os.listdir(f'{PATH}{SLASH}config'))
  while True:
    print('Please choose a pre-configured workflow from the list below:\n')
    # Print the list of all workflows a part for the global and current directories which are used internally by the program
    [print(f'{i}\t{v}') for i,v in enumerate(workflows) if (v not in  ('global','current')) and (os.path.isdir(f'{PATH}{SLASH}config{SLASH}{v}'))]
    workflow = input()
    # Test the input
    valid = check_input(workflow, workflows, special = ('global','current'))
    if not valid:
      continue
    workflow = workflows[int(workflow)] 
    # Save the chosen workflow in the "current" directory
    save_settings(source = workflow, output = 'current')
    break
  run_list = []
  while True:
    print('\nMain menu\n')    
    pprint('Enter 1 to change the configuration\nEnter 0 to continue to sample selection\n')
    change = input()
    valid = check_input(change, [0,1])
    if not valid:
      continue
    change = int(change)
    break
  # If the user wants to change the configuration, run the change_menu function
  if change:
    change_menu()
  for tsv in tsvs:
    tsv_path, tsv_name = metadata_handler(tsv)
    tsv_df = pd.read_csv(tsv_path, delimiter = '\t')

    try: tsv_df.scientific_name = tsv_df.scientific_name.str.replace(' ','_')
    except AttributeError:
      formats = [rows_of_links]
      print("It appears your file is not in EBI format.\n")
      [print(f'Enter {i}\tto {v.__name__}' for i,v in enumerate(formats))]
      choice = input()
      tsv_df = formats[choice](tsv_path)
    # Make sure the fastq_ftp col exists
    if 'fastq_ftp' not in tsv_df.columns:
      print('The fastq_ftp column does not exist in the tsv file. Please add it and try again.')
      print('Skipping to next tsv file...')
      continue
    for i in ['run_accession', 'fastq_ftp']:
      if i in tsv_df.columns:
        try:
          tsv_df[i] = tsv_df[i].apply(literal_eval)
        except (ValueError, SyntaxError):
          pass
    tsv_dir = os.path.dirname(tsv_path)
    pprint(f'Choose which lines of the tsv to process with the current pipeline.\nEnter a single row number (i.e 2) or a range (i.e 2-3)\nEnter all to work on all rows\nTo work on specific rows, enter a list of rows seperated by a SINGLE space (1 2 3 4)')
    wanted_cols = ['index', 'scientific_name', 'run_accession', 'sample_accession', 'sample_title', 'library_strategy', 'library_selection', 'read_count']
    colnames = ['Index', 'ID', 'Run accession', 'Sample accession', 'Sample', 'Library Strategy', 'Library Selection', 'Read Count']
    actual_cols = [i for i in wanted_cols if i in tsv_df.columns]
    sliced_df = tsv_df.loc[:, actual_cols]
    # if the type of run_accession col values is a list, convert to string for printing
    if isinstance(sliced_df.run_accession[0], list):
      sliced_df.run_accession = sliced_df.run_accession.apply(lambda x: ' '.join(x))
    print(tabulate(sliced_df, headers = 'keys', tablefmt = 'psql'))
    while True:
      work_on = input()
      work_on = input_to_list(work_on, len(tsv_df))
      if work_on == []:
        continue
      break
    run_list += sample_lvl_csh_maker(work_on, tsv_df, tsv_name, queue_mode)
  make_master_file(run_list, qsub_queue = queue_mode)
  save_settings(source = 'current', output = 'last_used')
  print('Enter \"bash master.csh\" into the console to run the commands\nIMPORTANT: Each sample will be executed simultaneously in qsub you run this command, make sure you dont run too many at once!\n')


if __name__ == '__main__':
  main()