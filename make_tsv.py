import os, sys, argparse, re
import pandas as pd

def make_tsv(path, org, libstrat):
    """
    """
    COLS = ['scientific_name', 'fastq_ftp', 'run_accession', 'sample_title']
    if libstrat != '':
        COLS.append('library_strategy')
    tsv_df = pd.DataFrame(columns = COLS)
    if libstrat != '':
        tsv_df['library_strategy'] = libstrat
    tsv_df['scientific_name'] = org
    file_list = []
    for file in os.listdir(path):
        file_list.append(file.replace('fastq.gz',''))
    files = ''.join(file_list)
    if 'R1' in files and 'R2' in files:
        PAIRED = True
    else: 
        PAIRED = False
    for file in file_list:
        to_append = []
        file = file.replace('.fastq.gz', '')
        to_append.append(org)
        file = re.sub(pattern = r'_R\d', repl = '')







def path_check(path : str) -> str:
    """
    Check whether the path is valid
    """
    if os.path.isfile(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f'{path} is not a valid path!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = __doc__, prog = 'Maker.py TSV creator')
    parser.add_argument('--org', help = 'Full taxonomic species name', metavar = 'org_name', type = str)
    parser.add_argument('-s','--sample-folder', help = 'full path to the sample folder where all the fastq files are located', metavar = 'sample_folder', type = path_check)
    parser.add_argument('-l, --libstrat', help = 'Library protocol', metavar = 'library_strategy', type = str, default = '')
    parser.add_argument('-o', '--output', help = 'Output folder', metavar = 'path', type = path_check, default = os.getcwd())
    parser.add_argument('-n', '--name', help = 'Choose how to name the TSV file', metavar = 'tsv_name', type = str, default = 'test')