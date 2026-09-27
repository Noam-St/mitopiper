"""Backwards-compatible entry point.

The original single-file maker.py has been refactored into the ``mitopiper``
package.  ``python maker.py`` still starts the interactive menus; see
``mitopiper --help`` for the scriptable sub-commands.
"""

import sys

from mitopiper.cli import main

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["interactive"]))
