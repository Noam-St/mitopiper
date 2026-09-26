"""Backwards-compatible wrapper for ``pipemaker make-tsv``.

usage: python make_tsv.py --org "Homo sapiens" -s /path/to/fastq_folder [-l RNA-Seq]
"""

import sys

from pipemaker.cli import main

if __name__ == "__main__":
    argv = sys.argv[1:]
    folder = None
    for flag in ("-s", "--sample-folder"):
        if flag in argv:
            i = argv.index(flag)
            folder = argv[i + 1]
            del argv[i:i + 2]
    if folder is None:
        sys.exit("usage: python make_tsv.py --org ORG -s SAMPLE_FOLDER")
    sys.exit(main(["make-tsv", folder] + argv))
