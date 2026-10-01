"""`kira-spliceqc-py ...` forwards the command line to the binary."""

import subprocess
import sys

from . import binary


def main() -> int:
    return subprocess.call([binary(), *sys.argv[1:]])


if __name__ == "__main__":
    sys.exit(main())
