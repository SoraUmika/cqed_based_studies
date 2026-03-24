"""Runtime compatibility helpers for the Windows study environment.

QuTiP transitively imports SciPy, which on this machine triggers a blocking WMI
query through ``platform.machine()`` during module import. The simulation stack
does not require that WMI result, so we force the standard library to take its
non-WMI fallback path before any ``qutip`` or ``cqed_sim`` import occurs.
"""

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
