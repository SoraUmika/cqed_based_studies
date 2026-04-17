"""Runtime compatibility helpers for the Windows study environment."""

from __future__ import annotations

import os
import platform


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


patch_windows_qutip_import()