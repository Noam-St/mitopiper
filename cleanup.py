import os, sys

PATH = os.getcwd()

if sys.platform != 'win32':
    slash = '/'
else:
    slash = '\\'

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
      os.remove(f'{PATH}{slash}{file}')
  clears = ['std_output', 'error_output']
  for c in clears:
    print(f'Clear {c}? (y/n)')
    choice = input()
    [os.remove(f'{PATH}{slash}{c}{slash}{file}') for file in os.listdir(f'{PATH}{slash}{c}') if yn(choice)]
    
if __name__ == '__main__':
  remove_files()

