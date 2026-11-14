"""Reproduce every V3 dataset, figure and output table:  python run_all.py"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src" / "python"))

from stcool.pipeline import main  # noqa: E402

if __name__ == "__main__":
    main(ROOT)
