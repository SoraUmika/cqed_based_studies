from __future__ import annotations

import sys
from pathlib import Path


WORKSPACE_ROOT = Path(__file__).resolve().parents[3]
SIM_ROOT_CANDIDATES = (
    Path("C:/Users/jl82323/Box/Shyam Shankar Quantum Circuits Group/Users/Users_JianJun/cQED_simulation"),
    Path("C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group/Users/Users_JianJun/cQED_simulation"),
)
SIM_ROOT = next((path for path in SIM_ROOT_CANDIDATES if (path / "cqed_sim").exists()), SIM_ROOT_CANDIDATES[0])

for _path in (WORKSPACE_ROOT, SIM_ROOT):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))
