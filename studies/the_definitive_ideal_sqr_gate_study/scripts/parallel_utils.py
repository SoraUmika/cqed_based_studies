from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import os
import traceback
from typing import Any, Callable, Sequence

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - optional dependency
    def tqdm(iterable, **_: Any):  # type: ignore[misc]
        return iterable


@dataclass
class TaskResult:
    item: Any
    ok: bool
    value: Any | None = None
    error: str | None = None


def default_worker_count() -> int:
    cpu_count = os.cpu_count() or 2
    return max(1, int(cpu_count) - 1)


def add_parallel_cli(parser: Any) -> None:
    parser.add_argument(
        "--n-workers",
        type=int,
        default=default_worker_count(),
        help="Parallel worker count. Defaults to os.cpu_count() - 1.",
    )
    parser.add_argument(
        "--sequential",
        action="store_true",
        help="Disable multiprocessing and run tasks serially for debugging.",
    )


def _run_wrapped(fn: Callable[[Any], Any], item: Any) -> TaskResult:
    try:
        return TaskResult(item=item, ok=True, value=fn(item))
    except Exception:
        return TaskResult(item=item, ok=False, error=traceback.format_exc())


def run_tasks(
    items: Sequence[Any],
    fn: Callable[[Any], Any],
    *,
    n_workers: int,
    sequential: bool,
    desc: str,
) -> list[TaskResult]:
    results: list[TaskResult] = []
    if sequential or n_workers <= 1:
        for item in tqdm(items, desc=desc):
            results.append(_run_wrapped(fn, item))
        return results

    with ProcessPoolExecutor(max_workers=int(n_workers)) as executor:
        future_map = {executor.submit(_run_wrapped, fn, item): item for item in items}
        for future in tqdm(as_completed(future_map), total=len(future_map), desc=desc):
            results.append(future.result())
    return results
