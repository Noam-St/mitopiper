"""
This program cleans up temporary files left behind by maker.py\n
usage: python cleanup.py [--forced] folder1 folder2 folder3...\n
The program will clean up all folders listed and automatically remove any run*.csh or master*.csh files.
"""

import os, sys, argparse

PATH = os.getcwd()

if sys.platform != 'win32':
    SLASH = '/'
else:
    SLASH = '\\'

def yn(c, forced):
  """
  """
  if forced:
    return True
  print(f'Clear {c} in ALL sample folders? (y/n)\n')
  choice = input()
  if type(choice) != str:
    return False
  elif 'y' in choice or 'Y' in choice:
    return True

def remove_files(clears, forced):
  for file in os.listdir(PATH):
    if 'run' in file or 'master' in file:
      os.remove(f'{PATH}{SLASH}{file}')
  for c in clears:
    yn_result = yn(c, forced)
    for f in os.listdir(os.path.join(PATH, 'output')):
      f = os.path.join(PATH, 'output', f)
      if os.path.isdir(f):
        direc = os.path.join(PATH, 'output', f)
        try:
          [os.remove(os.path.join(direc, c, file)) for file in os.listdir(os.path.join(direc, c)) if yn_result]
        except FileNotFoundError:
          print(f'Folder {c} does not exist in {os.path.basename(f)}!\n')
          continue
        [os.remove(os.path.join(direc, temp)) for temp in os.listdir(direc) if 'run' in temp]

def main(args):
  remove_files(args.clears, args.forced) 
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description = __doc__)
  parser.add_argument('clears', nargs = '+', default = ['stdout', 'stderr'], help = 'Choose which folders inside the sample directories to remove')
  parser.add_argument('--forced', action = 'store_true', default = False, help = 'Add this flag if you wish the program to run with no confirmations')

  args = parser.parse_args()
  main(args)


