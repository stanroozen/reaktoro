"""
conftest.py for parity tests.
Ensures the parity script directory is on sys.path so wrapper imports work.
"""

import sys
import os

_here = os.path.dirname(os.path.abspath(__file__))
if _here not in sys.path:
    sys.path.insert(0, _here)
