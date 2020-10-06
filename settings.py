"""
Created in 10/06/2020 by Noam Shtolz
This file defines the settings to run RNA-seq_pipeline.py by.
Must be included for RNA-seq_pipeline.py to run properly.
"""




import os, sys
import subprocess
path = os.getcwd()

if sys.platform != 'win32':
    slash = '/'
else:
    slash = '\\'
# %%


tsv_name = input(""" Please enter TSV file name.
                  IMPORTANT: The tsv must be in the same folder as this code, and in the following format:\n
                  study_accession\tsample_accession\trun_accession\ttax_id\tscientific_name\tfastq_ftp\tsubmitted_ftp\tsra_ftp\n
                  The format is available for download from the EBI website""")
settings = bool(int(input('If you would like to change the run-settings, please enter 1, enter 0 otherwise.')))
tsv_path = path + slash + tsv_name
subprocess.run(f'qsub -q bioinfo.q -cwd -N {tsv_name} -b y python3.6 RNA-seq_pipeline.py')
