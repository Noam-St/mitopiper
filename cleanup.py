import os, sys

PATH = os.getcwd()

if sys.platform != 'win32':
    SLASH = '/'
else:
    SLASH = '\\'

def yn(choice):
  """
  """
  if type(choice) != str:
    return False
  if 'y' in choice or 'Y' in choice:
    return True

def remove_files():
  for file in os.listdir(PATH):
    if 'run' in file or 'master' in file:
      os.remove(f'{PATH}{SLASH}{file}')
  clears = ['std_output', 'error_output']
  for c in clears:
    print(f'Clear {c}? (y/n)')
    choice = input()
    for f in os.listdir(f'{PATH}{SLASH}output'):
      if os.path.isdir(f):
        direc = f'{PATH}{SLASH}output{SLASH}{f}'
        [os.remove(f'{direc}{SLASH}{c}{SLASH}{file}') for file in os.listdir(f'{direc}{SLASH}{c}') if yn(choice)]
        [os.remove(f'{direc}{SLASH}{temp}') for temp in os.listdir(direc) if 'run' in temp]
    
if __name__ == '__main__':
  remove_files()

