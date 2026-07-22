"""Input/output helpers for the PFAS technical note workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import pandas as pd

from .schema import REQUIRED_COLUMNS, REQUIRED_FILES


def read_archive(input_dir: str | Path) -> Dict[str, pd.DataFrame]:
    """Read and validate the synthetic case CSV archive.

    Parameters
    ----------
    input_dir:
        Directory containing the technical note CSV files.

    Returns
    -------
    dict
        Mapping from file stem to pandas DataFrame.

    Raises
    ------
    FileNotFoundError
        If a required CSV file is missing.
    ValueError
        If required columns are absent.
    """
    input_path = Path(input_dir)
    if not input_path.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_path}")

    missing_files = [name for name in REQUIRED_FILES if not (input_path / name).exists()]
    if missing_files:
        raise FileNotFoundError(
            "Missing required technical note archive files: " + ", ".join(missing_files)
        )

    tables: Dict[str, pd.DataFrame] = {}
    for file_name in REQUIRED_FILES:
        key = Path(file_name).stem
        tables[key] = pd.read_csv(input_path / file_name)

    validate_archive(tables)
    return tables


def validate_archive(tables: Dict[str, pd.DataFrame]) -> None:
    """Validate required columns in the technical note archive."""
    for file_name, required_cols in REQUIRED_COLUMNS.items():
        key = Path(file_name).stem
        if key not in tables:
            raise ValueError(f"Table {key} was not loaded")
        missing_cols = [c for c in required_cols if c not in tables[key].columns]
        if missing_cols:
            raise ValueError(
                f"Table {file_name} is missing required columns: {', '.join(missing_cols)}"
            )


def ensure_output_dirs(output_dir: str | Path) -> Dict[str, Path]:
    """Create output directories and return their paths."""
    out = Path(output_dir)
    tables = out / "tables"
    figures = out / "figures"
    logs = out / "logs"
    for p in (out, tables, figures, logs):
        p.mkdir(parents=True, exist_ok=True)
    return {"root": out, "tables": tables, "figures": figures, "logs": logs}


def write_table(df: pd.DataFrame, output_dir: str | Path, file_name: str) -> Path:
    """Write a DataFrame to CSV and return the file path."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / file_name
    df.to_csv(path, index=False)
    return path
