"""pipemaker: generate HPC job scripts for sequencing-data pipelines.

pipemaker turns a sample metadata table (ENA or SRA) and a configurable
workflow into one job script per sample for SGE, SLURM or LSF, plus a master
script that submits them with the right dependencies.
"""

__version__ = "2.0.0.dev0"
__author__ = "Noam Shtolz"
