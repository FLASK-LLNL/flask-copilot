#!/usr/bin/env python3
"""
Summarize checkpoint/session data stored in the database.
"""

import argparse
import json
from datetime import datetime
from typing import Any

from sqlalchemy import select
from sqlalchemy.orm import sessionmaker

from db_backend.database.engine import sync_engine
from db_backend.database import models


def mask_url(url: str) -> str:
    if "@" not in url or "://" not in url:
        return url
    scheme, rest = url.split("://", 1)
    if "@" not in rest or ":" not in rest.split("@", 1)[0]:
        return url
    creds, host = rest.split("@", 1)
    user = creds.split(":", 1)[0]
    return f"{scheme}://{user}:***@{host}"


def safe_len(value: Any) -> int:
    if value is None:
        return 0
    if isinstance(value, (list, tuple, dict)):
        return len(value)
    if isinstance(value, str):
        try:
            parsed = json.loads(value)
            if isinstance(parsed, (list, tuple, dict)):
                return len(parsed)
        except Exception:
            return len(value)
    return 1


def short(value: Any, max_len: int = 28) -> str:
    if value is None:
        return ""
    text = str(value)
    if len(text) <= max_len:
        return text
    return text[: max_len - 3] + "..."


def fmt_dt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, datetime):
        return value.strftime("%Y-%m-%d %H:%M:%S")
    return str(value)


def print_table(headers: list[str], rows: list[list[Any]]) -> None:
    if not rows:
        print("(no rows)")
        return
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(str(cell)))
    fmt = " | ".join([f"{{:{w}}}" for w in widths])
    sep = "-+-".join(["-" * w for w in widths])
    print(fmt.format(*headers))
    print(sep)
    for row in rows:
        print(fmt.format(*row))


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize checkpoint data.")
    parser.add_argument("--limit", type=int, default=0, help="Limit number of experiments shown (0 = all)")
    args = parser.parse_args()

    if sync_engine is None:
        print("Database engine not initialized (MariaDB is not reachable).")
        print("\nOptions:")
        print("  1. Start MariaDB:  brew services start mariadb")
        print("  2. Open SSH tunnel: ssh -L 32636:cz-marathe1-mymariadb1.apps.czapps.llnl.gov:32636 oslic.llnl.gov")
        print("  3. Use SQLite:     USE_SQLITE_FALLBACK=1 python checkpoint_report.py --limit 20")
        return 1

    SessionLocal = sessionmaker(bind=sync_engine)

    try:
        with SessionLocal() as session:
            projects = session.execute(select(models.Project)).scalars().all()
            experiments = session.execute(select(models.Experiment)).scalars().all()
    except Exception as exc:
        print(f"ERROR: Could not query database: {exc}")
        print("\nPossible fixes:")
        print("  1. Start MariaDB:  brew services start mariadb")
        print("  2. Open SSH tunnel: ssh -L 32636:cz-marathe1-mymariadb1.apps.czapps.llnl.gov:32636 oslic.llnl.gov")
        print("  3. Use SQLite:     USE_SQLITE_FALLBACK=1 python checkpoint_report.py --limit 20")
        return 1

    # Summary
    try:
        db_url = mask_url(str(sync_engine.url))
    except Exception:
        db_url = "(unknown)"

    running_count = sum(1 for exp in experiments if exp.is_running)

    print("=== Database Summary ===")
    print(f"Database URL: {db_url}")
    print(f"Projects: {len(projects)}")
    print(f"Experiments: {len(experiments)}")
    print(f"Running experiments: {running_count}")
    print()

    # Project summary table
    proj_rows = []
    exp_by_project = {}
    for exp in experiments:
        exp_by_project.setdefault(exp.project_id, []).append(exp)

    for proj in projects:
        exps = exp_by_project.get(proj.id, [])
        last_modified = max([e.last_modified for e in exps], default=proj.last_modified)
        proj_rows.append([
            proj.id,
            short(proj.name, 32),
            proj.user,
            len(exps),
            fmt_dt(last_modified),
        ])

    print("=== Projects ===")
    print_table(["project_id", "name", "user", "experiments", "last_modified"], proj_rows)
    print()

    # Experiment details table
    experiments_sorted = sorted(experiments, key=lambda e: e.last_modified or e.created_at or datetime.min, reverse=True)
    if args.limit and args.limit > 0:
        experiments_sorted = experiments_sorted[: args.limit]

    exp_rows = []
    for exp in experiments_sorted:
        exp_rows.append([
            exp.id,
            short(exp.name, 24),
            exp.user,
            short(exp.smiles, 24),
            short(exp.problem_type, 16),
            "yes" if exp.is_running else "no",
            fmt_dt(exp.last_modified),
            safe_len(exp.tree_nodes),
            safe_len(exp.edges),
            safe_len(exp.metrics_history),
            safe_len(exp.visible_metrics),
            "yes" if exp.graph_state else "no",
            "yes" if exp.experiment_context else "no",
        ])

    print("=== Experiments (Checkpoint Data) ===")
    print_table(
        [
            "experiment_id",
            "name",
            "user",
            "smiles",
            "problem",
            "running",
            "last_modified",
            "nodes",
            "edges",
            "metrics",
            "visible",
            "graph",
            "context",
        ],
        exp_rows,
    )

    if not experiments:
        print("No experiments found. Try running a checkpoint or saving a session in the UI.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
