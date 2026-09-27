"""Backwards-compatible wrapper for ``mitopiper cleanup``.

usage: python cleanup.py [--forced] [--keep-csh] folder1 folder2 ...
"""

import sys

from mitopiper.cli import main

if __name__ == "__main__":
    argv = ["--yes" if a == "--forced" else "--keep-scripts" if a == "--keep-csh" else a for a in sys.argv[1:]]
    sys.exit(main(["cleanup"] + argv))
