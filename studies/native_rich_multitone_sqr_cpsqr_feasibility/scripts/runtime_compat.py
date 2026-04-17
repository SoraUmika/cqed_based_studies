"""Runtime compatibility helpers for the Windows study environment."""

from __future__ import annotations

import os
import platform
import sys
from pathlib import Path


def patch_windows_qutip_import() -> None:
    """Avoid a blocking WMI call during qutip/scipy import on Windows."""
    if os.name != "nt":
        return
    if getattr(platform, "_cqed_qutip_wmi_patched", False):
        return

    original_query = getattr(platform, "_wmi_query", None)
    if callable(original_query):

        def _raise_oserror(*args, **kwargs):
            raise OSError("WMI disabled for qutip import compatibility")

        platform._wmi_query = _raise_oserror

    platform._cqed_qutip_wmi_patched = True


def ensure_cqed_sim_on_path() -> Path:
    """Add the patched local ``cQED_simulation`` checkout to ``sys.path``."""
    script_dir = Path(__file__).resolve().parent
    candidates: list[Path] = []

    env_value = os.environ.get("CQED_SIM_ROOT", "").strip()
    if env_value:
        candidates.append(Path(env_value).expanduser())

    for parent in script_dir.parents:
        candidates.append(parent / "cQED_simulation")

    for candidate in candidates:
        if candidate.exists():
            resolved = candidate.resolve()
            candidate_str = str(resolved)
            if candidate_str not in sys.path:
                sys.path.insert(0, candidate_str)
            return resolved

    raise FileNotFoundError(
        "Could not locate the patched cQED_simulation checkout. "
        "Set CQED_SIM_ROOT or place the checkout beside the workspace root."
    )


patch_windows_qutip_import()
ensure_cqed_sim_on_path()
