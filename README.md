# Semi Automatic csh file creator for RNA-seq pipeline using qsub

This project script creates csh file that can be run easily using qsub.

After cloning the repository, to run the script type `python3.6 maker.py` into the console.
The script will then request the name of the input tsv file, the tsv file must be in this format:

study_accession   sample_accession    experiment_accession    run_accession   tax_id    scientific_name   fastq_ftp   submitted_ftp   sra_ftp

and inside **samples** folder.
The default pipeline is:
- download with wget
- trimming with trim_galore
- quality check with fastq
- alignment with bwa mem
- sam to bam conversion with samtools view
- read sorting with samtools sort
- indexing with samtools index
- pileup file creation with samtools pileup

The pipeline is run with the default parameters for each program, it is possible to customize the parameters or add new parameters in the configuration menu.
Also in the configuration menu, it is possible to add new pipeline steps, reorder the pipeline and remove unwanted steps.

After the script run is completed, a `master.csh` file will be created, type `qsub master.csh` into the console to run the file and execute the commands
