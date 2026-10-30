#!/usr/bin/env python3
"""3g step 3: run the E9 stream driver (scripts/e9/stream.py, unchanged) on D8 runs against the 3g index, with the
download sha256 log redirected to the 3g work dir (the tracked E9 manifest file is not appended to).
    python3 scripts/decoy_realism/stream_3g.py --start 0 --end 100 --work /tmp/vulcan_3g --no-c \
        --index-p /tmp/vulcan_3g/index/human_contam_3g
Then: conda_envs/pred/bin/python scripts/e9/score_worker.py --work /tmp/vulcan_3g --index-p <3g index> \
        --predcache /tmp/vulcan_3g/predcache --exit-when-idle 900
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "e9"))
import stream  # noqa: E402

stream.SHA = os.environ.get("VULCAN_3G_SHA", "/tmp/vulcan_3g/d8_sha256.txt")

if __name__ == "__main__":
    stream.main()
