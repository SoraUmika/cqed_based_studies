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


def ensure_repo_import_path() -> None:
    study_scripts = Path(__file__).resolve().parent
    workspace_root = study_scripts.parents[2]
    cqed_repo = workspace_root.parent / "cQED_simulation"
    if cqed_repo.exists():
        repo_str = str(cqed_repo)
        if repo_str not in sys.path:
            sys.path.insert(0, repo_str)


patch_windows_qutip_import()
ensure_repo_import_path()
