"""Minimalist, print-safe figure style (matplotlib only).

Palette (validated for colour-vision deficiency, OKLab dE >= 12 between neighbours):
  BLUE   #3556a8  - primary series      (dark in greyscale)
  TEAL   #16a085  - secondary series    (mid grey)
  OCHRE  #e6a519  - tertiary series     (light grey; always paired with a dark marker edge)
  Neutrals: ink #1f2328, mid grey #6b7280, light grey #c9ced6 (null bands, reference lines)
Every series is also distinguished by marker shape and line style or hatch,
so that the figures remain readable in black-and-white print.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, TEAL, OCHRE = "#3556a8", "#16a085", "#e6a519"
INK, GREY, LIGHT = "#1f2328", "#6b7280", "#c9ced6"
SERIES = [
    dict(color=BLUE, marker="o", ls="-"),
    dict(color=OCHRE, marker="s", ls="--"),
    dict(color=TEAL, marker="^", ls="-."),
]
MM = 1 / 25.4
W1, W2 = 85 * MM, 180 * MM   # single- and double-column widths


def apply():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "Carlito", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
        "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK, "ytick.color": INK,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "axes.grid.axis": "y", "grid.color": "#e5e7eb", "grid.linewidth": 0.5,
        "axes.axisbelow": True,
        "lines.linewidth": 1.2, "lines.markersize": 4.5,
        "legend.frameon": False,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.dpi": 300, "figure.dpi": 120,
        "hatch.linewidth": 0.5,
    })


def panel_label(ax, letter):
    ax.text(-0.02, 1.08, letter, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="bottom", ha="right")


def save(fig, name, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"{name}.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
