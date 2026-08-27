"""Runner: Quartz solubility in H2O-CO2, PerplexDEW backend, ExtendedDH model."""

import sys
import os

sys.argv = [sys.argv[0], "--backend", "PerplexDEW", "--dh-model", "ExtendedDH"]

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import quartz_H2OCO2_solubility as _q

_q.main()
