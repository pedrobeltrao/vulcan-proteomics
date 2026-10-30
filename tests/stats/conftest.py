import os

import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from vulcan import hw  # noqa: E402

# tests must be quick and not oversubscribe: cap at 4 unless VULCAN_THREADS says otherwise
os.environ.setdefault("VULCAN_THREADS", str(min(4, hw.n_threads())))
hw.pin_thread_pools()

import pytest  # noqa: E402


def pytest_addoption(parser):
    try:
        parser.addoption("--runslow", action="store_true", default=False,
                         help="run slow (>=1e4-replicate) Monte-Carlo tests")
    except ValueError:  # already registered by another conftest
        pass


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: >=1e4-replicate Monte-Carlo tests (run with --runslow)")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        return
    skip = pytest.mark.skip(reason="slow Monte-Carlo test; use --runslow")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip)
