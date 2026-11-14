"""Shared figure style. Colour follows the entity in every figure."""
from __future__ import annotations

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = dict(blue="#2a78d6", orange="#eb6834", aqua="#1baf7a", yellow="#eda100", magenta="#e87ba4",
         green="#008300", violet="#6250d6", red="#e34948", ink="#0b0b0b", ink2="#52514e",
         muted="#8a8985", grid="#e6e5e1", surface="#ffffff")

ENTITY = {
    "original": C["ink2"],
    "baseline": C["blue"],
    "nominal": C["orange"],
    "recommended": C["aqua"],
    "v1 robust": C["violet"],
    "v2 recommended": C["magenta"],
    "adaptive": C["green"],
    "target": C["red"],
}


def apply_style():
    plt.rcParams.update({
        "figure.facecolor": C["surface"], "axes.facecolor": C["surface"], "savefig.facecolor": C["surface"],
        "axes.edgecolor": C["muted"], "axes.labelcolor": C["ink"], "xtick.color": C["ink2"], "ytick.color": C["ink2"],
        "axes.grid": True, "grid.color": C["grid"], "grid.linewidth": 0.8, "axes.spines.top": False,
        "axes.spines.right": False, "font.size": 9.5, "axes.titlesize": 10.5, "axes.titleweight": "bold",
        "axes.titlelocation": "left", "legend.frameon": False, "lines.linewidth": 1.8,
        "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
        # Pinned so that a local matplotlibrc cannot change the rendering (values used for the published figures).
        "figure.dpi": 200, "text.hinting": "no_hinting",
    })


def save(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)
