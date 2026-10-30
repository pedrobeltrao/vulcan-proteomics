#!/usr/bin/env python3
"""Write the node fingerprint (vulcan.hw.fingerprint) as JSON: `python bench/hwinfo.py [out.json]`."""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from vulcan.hw import fingerprint  # noqa: E402


def main(argv):
    s = json.dumps(fingerprint(), indent=2, sort_keys=True)
    if len(argv) > 1 and argv[1] != "-":
        with open(argv[1], "w") as fh:
            fh.write(s + "\n")
    else:
        print(s)


if __name__ == "__main__":
    main(sys.argv)
