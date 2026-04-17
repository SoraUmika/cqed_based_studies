from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common as c


STYLE_PATH = c.STUDY_ROOT.parents[1] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
if STYLE_PATH.exists():
    plt.style.use(str(STYLE_PATH))


def load_frontier() -> dict[str, Any]:
    path = c.DATA_DIR / "grape_frontier_extension.json"
    return json.loads(path.read_text(encoding="utf-8"))


def complete_rows(frontier: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for payload in frontier.values():
        if str(payload.get("status", "")) != "complete":
            continue
        rows.append(payload)
    rows.sort(key=lambda row: float(row["duration_ns"]))
    return rows


def save_figure(fig: plt.Figure, stem: str) -> None:
    png_path = c.FIG_DIR / f"{stem}.png"
    pdf_path = c.FIG_DIR / f"{stem}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    rows = complete_rows(load_frontier())
    if not rows:
        raise RuntimeError("No complete frontier rows available.")

    durations = [float(row["duration_ns"]) for row in rows]
    replay = [float(row["best_replay_fidelity"]) for row in rows]
    leakage = [float(row["best_replay_leakage_worst"]) for row in rows]
    open_process = [float(row["open_process"]["process_fidelity"]) for row in rows]
    basis_fid = [float(row["open_process"]["mean_basis_state_fidelity"]) for row in rows]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

    ax = axes[0]
    ax.plot(durations, replay, "o-", color="#4477AA", label="Best replay fidelity")
    ax.plot(durations, open_process, "s--", color="#AA3377", label="Open-process fidelity")
    ax.axhline(0.99, color="#444444", linestyle=":", linewidth=1.0, label="0.99 target")
    ax.set_xlabel("Duration (ns)")
    ax.set_ylabel("Fidelity")
    ax.set_ylim(0.85, 1.005)
    ax.set_title("Validated GRAPE frontier")
    ax.legend(loc="lower right", fontsize=8)

    ax = axes[1]
    ax.plot(durations, leakage, "o-", color="#CC6677", label="Worst replay leakage")
    ax.plot(durations, basis_fid, "s--", color="#228833", label="Mean open basis fidelity")
    ax.set_xlabel("Duration (ns)")
    ax.set_ylabel("Leakage / fidelity")
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Replay leakage and open-basis fidelity")
    ax.legend(loc="lower right", fontsize=8)

    fig.suptitle("N_cav=12 GRAPE extension", y=1.02)
    save_figure(fig, "grape_frontier_extension")
    print("[grape-frontier-plot] wrote grape_frontier_extension.{png,pdf}", flush=True)


if __name__ == "__main__":
    main()