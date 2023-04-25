"""
Automatically create a TSV file to be used by maker.py script
NOTICE : This script requires a sample folder to run on, the samples inside the sample folder will be copied to the output location inside maker.py's folder structure and renamed to fit maker.py (_i_org_initials_1/2 if paired) will be added to the end of the file.
After running this script, the samples should be ready for initial steps (fastqc, trimming, alignment). 
IMPORTANT : Wget will not work because an ftp link cannot be generated!
"""

import os, sys, argparse, re, random, shutil
import pandas as pd
from maker import dir_checker
import logging

PATH = os.getcwd()

def gen_id():
    """
    """
    MAX = 10000000
    MIN = 1000000
    num = random.randint(MIN, MAX)
    return 'SRR' + str(num)

def make_tsv(sample_files, org, libstrat, ignore_paired, force_paired):
    """
    """
    if sample_files.endswith('/') or sample_files.endswith('\\'):
        sample_files = os.path.dirname(sample_files)
    NAME = os.path.split(sample_files)[-1]    
    dir_checker('output')
    dir_checker(os.path.join('output', NAME))
    dir_checker(os.path.join('output', NAME, 'fastq'))
    COLS = ['scientific_name', 'fastq_ftp', 'run_accession', 'sample_title']
    if libstrat != '':
        COLS.append('library_strategy')
        lib = True
    tsv_dict = {k:[] for k in COLS}
    file_list = []
    for file in os.listdir(sample_files):
        # Grab all the file format suffixes
        if (file.endswith('.fastq.gz') or file.endswith('.fq.gz') or file.endswith('.fastq') or file.endswith('.fq')):
            file_name = os.path.split(file)[-1]
            print(f'Copying {file_name} to {os.path.join(PATH, "output", NAME, "fastq")}')
            shutil.copy(src = os.path.join(sample_files, file), dst = os.path.join(PATH, 'output', NAME, 'fastq'))            
            print(f'Copied {file} to {os.path.join(PATH, "output", NAME, "fastq")} successfully!')
            for suffix in ['.fastq.gz', '.fq.gz', '.fastq', '.fq', '.gz']:
                file = file.replace(suffix, '')
            file_list.append(file)
    files = ''.join(file_list)
    if ('R1' in files and 'R2' in files and not ignore_paired) or force_paired:
        # Report that paired library detected
        print('Paired library detected!')
        #PAIRED
        file_list.sort()
        for i in range(0, len(file_list), 2):
            #file = re.sub(pattern = r'_R\d', repl = '')
            file1 = file_list[i]
            file2 = file_list[i + 1]
            print(f'Paired files detected: {file1} and {file2}')
            tsv_dict['scientific_name'].append(org)
            tsv_dict['fastq_ftp'].append(f'{file1};{file2}')
            file1 = re.sub(pattern = r'(_R\d|_\d{1}_)', repl = '', string = file1)
            tsv_dict['run_accession'].append(file1)
            tsv_dict['sample_title'].append('_'.join(file1.split('_')[0:2]))
            if lib:
                tsv_dict['library_strategy'].append(libstrat)
    else: 
        for file in file_list:
            tsv_dict['scientific_name'].append(org)
            tsv_dict['fastq_ftp'].append(file)
            tsv_dict['run_accession'].append(re.sub(pattern = r'_R\d', repl = '', string = file))
            tsv_dict['sample_title'].append('_'.join(file.split('_')[0:2]))
            if lib:
                tsv_dict['library_strategy'].append(libstrat)
    tsv_df = pd.DataFrame(tsv_dict)
    tsv_df.to_csv(os.path.join(PATH, 'output', NAME, f'{NAME}.tsv'), sep = '\t', header = True, index = False)
    print(f'Program complete. TSV file created at {os.path.join(PATH, "output", NAME, f"{NAME}.tsv")}.')
    print('To analyze the samples, please run maker.py and find the samples in the experiment selection menu')
    return NAME

def fix_names(org, NAME):
    """
    """
    initials = org.split(' ')
    try:
        initials = initials[0][0] + initials[1][0]
    except IndexError:
        initials = initials[0].split('_')
        try: initials = initials[0][0] + initials[1][0]
        except IndexError: initials = initials[0][0]
    initials = initials.upper()

    tsv_df = pd.read_csv(os.path.join(PATH, 'output', NAME, f'{NAME}.tsv'), sep = '\t')
    suffixes = ['.fastq.gz', '.fq.gz', '.fastq', '.fq', '.gz']
    c = 0
    for i, row in tsv_df.iterrows():
        if len(row.fastq_ftp.split(';')) > 1: PAIRED = True
        else: PAIRED = False
        if PAIRED:
            i = 1
            for f in row.fastq_ftp.split(';'):
                # Try the old file with several suffixes
                for suffix in suffixes:
                    old = os.path.join(PATH, 'output', NAME, 'fastq', f'{f}{suffix}')
                    if os.path.exists(old): break
                f = re.sub(pattern = r'(_R\d|_\d{1}_)', repl = '', string = f)
                if '.gz' in suffix:
                    new = os.path.join(PATH, 'output', NAME, 'fastq', f'{f}_{c}_{initials}_{i}.fastq.gz')
                else:
                    new = os.path.join(PATH, 'output', NAME, 'fastq', f'{f}_{c}_{initials}_{i}.fastq')
                os.rename(src = old, dst = new)
                i+=1
        else:
            f = row.fastq_ftp
            for suffix in suffixes:
                old = os.path.join(PATH, 'output', NAME, 'fastq', f'{f}{suffix}')
                if os.path.exists(old): break
            f = re.sub(pattern = r'_R\d', repl = '', string = f)
            if '.gz' in suffix:
                new = os.path.join(PATH, 'output', NAME, 'fastq', f'{f}_{c}_{initials}.fastq.gz')
            else:
                new = os.path.join(PATH, 'output', NAME, 'fastq', f'{f}_{c}_{initials}_{i}.fastq')
            os.rename(old, new)
        c+=1

def main(args):
    """
    """
    NAME = make_tsv(args.sample_folder, args.org, args.lib_strat, args.ignore_paired, args.force_paired)
    fix_names(args.org, NAME)
    print(f'TSV file created successfuly!\n')
    print(f'File name\t{NAME}\n')
    print('To run pipeline on the samples, launch maker.py and select the appropriate name.\nThe created TSV file is already located within the sample folder, and the sample files were transferred into the sample file location.')

def path_check(path : str) -> str:
    """
    Check whether the path is valid
    """
    if os.path.exists(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f'{path} is not a valid path!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = __doc__, prog = 'Maker.py TSV creator')
    parser.add_argument('--org', help = 'Full taxonomic species name', metavar = 'org_name', type = str)
    parser.add_argument('-s','--sample-folder', help = 'full path to the sample folder where all the fastq files are located', metavar = 'sample_folder', type = path_check)
    parser.add_argument('-l','--lib-strat', help = 'Library protocol', metavar = 'library_strategy', type = str, default = '')
    # Ignore paired arguement
    parser.add_argument('--ignore-paired', help = 'Ignore paired library detection', action = 'store_true', default = False)
    # Force paired arguement
    parser.add_argument('--force-paired', help = 'Force paired library detection', action = 'store_true', default = False)


    args = parser.parse_args()
    main(args)