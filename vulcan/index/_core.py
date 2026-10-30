"""Import the vulcan_core PyO3 extension (sage-core 99407db6; built by scripts/m2/build_vulcan_core.sh)."""
from __future__ import annotations

import os
import sys

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_DIR = os.environ.get("VULCAN_CORE_DIR", os.path.join(_ROOT, "tools", "vulcan_core"))
if _DIR not in sys.path:
    sys.path.insert(0, _DIR)
try:
    import vulcan_core as core  # noqa: E402
except ImportError as e:  # pragma: no cover
    raise ImportError(f"vulcan_core not found in {_DIR}; run scripts/m2/build_vulcan_core.sh") from e

__all__ = ["core"]
