"""Follow the jobs launched by master.sh and track their errors.

usage: python monitor.py [--watch [SECONDS]] [--errors] [--job NAME] [experiment ...]

Equivalent to ``mitopiper monitor``; see ``python monitor.py --help``.
"""

import sys

from mitopiper.cli import main

if __name__ == "__main__":
    sys.exit(main(["monitor"] + sys.argv[1:]))
