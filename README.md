# Semi automatic csh file creator for genomics pipeline using qsub

This program's main function is creating csh files which can be run using qsub and execute a pipeline sequentially on every sample chosen.
After cloning the repository, to run the script type `python maker.py` into the console.
The script will then request the name of the input tsv file, the tsv file must be in the same folder as the code and in this format:

study_accession   sample_accession    experiment_accession    run_accession   tax_id    scientific_name   fastq_ftp   submitted_ftp   sra_ftp

This inout file in this format can be downloaded from ENA.

Usage:
`python maker.py`

Immediately after execution, you will be prompted to pick a pre-configured workflow.
Workflows can be saved in the configuration menu, saved workflows will be shown here as well as a *last_used* workflow.
After choosing your preset, you will need to enter a formatted TSV file. previously used tsv files are saved as seperate instances and can will show up on the option list, if you want to enter a new sample file, you will have to enter its full name. Afterwards it will be saved permenantly.

Having chosen both the workflow preset and the sample TSV file, you can either enter the configuration menu or continue to choose which samples to run the selected pipeline on.

The *default* pipeline is:
- download with wget
- trimming with trim_galore
- quality check with fastq
- alignment with bwa mem
- sam to bam conversion with samtools view
- read sorting with samtools sort
- indexing with samtools index
- pileup file creation with samtools pileup

**Customization menus**
- change_general_settings:
  Change settings that are not specific to a single pipeline step
    - environment: the linux environment to run the pipeline algorithms on.
    - free_form: an option to enter any free-from BASH code to each run.csh file created by this program (USE CAREFULLY).
    - qsub_queue: select which qsub queue to send the pipeline algorithm steps to.
- change_pipeline_settings:
  Change the currently defined settings for a chosen pipeline step.
- change_pipeline:
  Change the pipeline structure (run order and steps).
- add_new_pipeline_setting
  Permenantly add a new setting to one of the defined pipeline steps, this setting can later be changed using *change_pipeline_setting*.
- add_new_pipeline
  Add a new pipeline step, this pipeline step can later be change using *change_pipeline*.
- manually_change_ref_to_gtf
  Freely changed the paths to reference fastas or gtf files for a certain organism.
- save_current_settings
  Save the current setup as a pre-defined workflow for later use.
- exit
  Exit customization and continue to sample-selection.

Following customization, you enter *sample-selection*.
All the samples available in the TSV file will show up as a table in the command window, to choose which samples to run the program on, you must enter either a single `index` from the table (i.e 0), a range (i.e 1-4) or `all` in order to run the pipeline on every single sample

If a reference fasta is not found automatically (first run on this organism), you will have 2 options:
- Create a reference automatically
  The program will attempt to contact the NCBI servers in order to download the appropriate mtDNA reference for this organism
- Enter a reference fasta manually **RECOMMENDED**
  Enter the *full* path to the reference fasta.

If a gtf file is not found automatically (first run on this organism), you will be asked to proved a gtf path, this is not necessary if your pipeline does not require gtf to run.

After the script run is completed, a `master.csh` file will be created in the main program directory, type `qsub master.csh` into the console to run the file and execute the commands.

***Adding new pipeline steps***
If you choose to enter a new custom pipeline step, you must do it carefully.
You will be prompted to set every single setting sequentially:
- name
  Enter a name for this new pipeline step which will be shown in the run-order.
- description
  Enter an informative description for fututre reference.
- input extension
  Enter the file format the pipeline step runs on, for example, `bwa mem` runs on `.fastq.gz` files, file suffixes **MUST** start with a `.`.
- input directory
  Enter the directory the input file is expected to be in, for example, `bwa mem` picks up inputs from post_trim folder usually, that is the folder `trim_galore` outputs to by default.
- command
  This is the command which calls that pipeline step algorithm, for example, the command for the default alignment step is `bwa mem`.
- output directory
  This is the directory where the output files should be created, for example, `bwa mem` outputs files into the alignment folder by default.
- output extension
  The format of the files this pipeline step should output natively, this will determine how the stdout will be piped `>` into a file.
- cluster location
  If this pipeline step is not runnable via the environment, and requires a specific path, enter it here. Otherwise, leave empty.
- **command structure**
  This step is very important and must be used carefully. most commands run in the following format *key settings inputfile*, however, if your new defined step requires something different, this is the place to change it accordingly, the structure can be defined using several keywords. Such as key (command name), settings (all the defined settings) and more. The full list is available in the program.


