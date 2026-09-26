"""Backwards-compatible entry point.

The original single-file maker.py has been refactored into the ``pipemaker``
package.  ``python maker.py`` still starts the interactive menus; see
``pipemaker --help`` for the scriptable sub-commands.
"""

import sys

from pipemaker.cli import main

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["interactive"]))
