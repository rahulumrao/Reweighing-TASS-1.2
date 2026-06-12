#!/usr/bin/env python3
"""
Plot TASS example output using styles from:
https://github.com/rahulumrao/Matplotlib_pyplot_script

Adapted from 1_plot.py and 3_error_plot.py; 2D map follows example/splot.gpi.
Probability plot uses PROB.dat_* from REWEIGHTING TOOL = prob runs.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.ticker import FixedLocator, FormatStrFormatter
from scipy.interpolate import RectBivariateSpline

ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT / "example"
OUT_DIR = ROOT / "docs" / "_static" / "images" / "example"

DPI = 200
X2D_LO, X2D_HI = 2.0, 4.0
Y2D_LO, Y2D_HI = 1.5, 4.5


def _style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "cm",
            "figure.dpi": DPI,
            "savefig.dpi": DPI,
            "savefig.bbox": "tight",
        }
    )


def _load_delta_g() -> tuple[np.ndarray, np.ndarray]:
    umb, err = [], []
    with (EXAMPLE / "delta_G.dat").open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            umb.append(float(parts[1]))
            err.append(float(parts[3]))
    return np.asarray(umb), np.asarray(err)


def _umbrella_centers() -> list[float]:
    centers: list[float] = []
    with (EXAMPLE / "input.inp").open() as fh:
        for line in fh:
            parts = line.split()
            if len(parts) == 2:
                try:
                    centers.append(float(parts[0]))
                except ValueError:
                    continue
    return centers


def plot_1d_free_energy(out: Path) -> None:
    """Adapted from Matplotlib_pyplot_script/1_plot.py."""
    data = np.loadtxt(EXAMPLE / "free_energy.dat")
    data1 = np.loadtxt(EXAMPLE / "interp_free_energy.dat")

    x, fes = data[:, 0], data[:, 1]
    x1, fes1 = data1[:, 0], data1[:, 1]
    z, z1 = fes.min(), fes1.min()

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.scatter(x, fes - z, color="tab:blue", alpha=0.6, s=45, label="Mean-force PMF")
    ax.plot(x1, fes1 - z1, color="tab:red", lw=2.0, label="B-spline interpolation")

    ax.set_xlim(x.min() - 0.05, x.max() + 0.05)
    ax.set_ylim(min(0.0, (fes - z).min()) - 0.5, (fes - z).max() + 1.0)
    ax.set_xlabel(r"CV$_1$ (umbrella coordinate)", fontsize=14, fontweight="bold")
    ax.set_ylabel(r"$\Delta F$ (kcal mol$^{-1}$)", fontsize=14, fontweight="bold")
    ax.grid(color="grey", linestyle="-.", linewidth=0.8)
    ax.legend(frameon=True)
    ax.tick_params(labelsize=12)
    fig.savefig(out / "free_energy_1d.png")
    plt.close(fig)


def plot_1d_with_error(out: Path) -> None:
    """Adapted from Matplotlib_pyplot_script/3_error_plot.py."""
    if not (EXAMPLE / "delta_G.dat").exists():
        return

    data1 = np.loadtxt(EXAMPLE / "interp_free_energy.dat")
    x1, fes1 = data1[:, 0], data1[:, 1]
    z1 = fes1.min()

    umb, sqrt_dg = _load_delta_g()
    yerr = np.interp(x1, umb, sqrt_dg)

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(x1, fes1 - z1, color="tab:orange", lw=2.0, label="Interpolated PMF")
    ax.fill_between(
        x1,
        (fes1 - z1) - yerr,
        (fes1 - z1) + yerr,
        color="tab:orange",
        alpha=0.35,
        label=r"$\pm\sqrt{\Delta G}$",
    )

    ax.set_xlim(x1.min() - 0.05, x1.max() + 0.05)
    ax.set_ylabel(r"$\Delta F$ (kcal mol$^{-1}$)", fontsize=14, fontweight="bold")
    ax.set_xlabel(r"CV$_1$ (umbrella coordinate)", fontsize=14, fontweight="bold")
    ax.grid(color="grey", linestyle="-.", linewidth=0.8)
    ax.legend(frameon=True)
    ax.tick_params(labelsize=12)
    fig.savefig(out / "free_energy_1d_error.png")
    plt.close(fig)


def plot_1d_probability(out: Path) -> None:
    """1D unbiased probability along umbrella CV (PROB.dat_*)."""
    prob_files = sorted(EXAMPLE.glob("PROB.dat_*"), key=lambda p: int(p.name.split("_")[-1]))
    if not prob_files:
        raise FileNotFoundError("No PROB.dat_* files — run with input_prob.tass first")

    centers = _umbrella_centers()
    nr = len(prob_files)
    cmap = plt.colormaps["coolwarm"]

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    for idx, path in enumerate(prob_files):
        ir = int(path.name.split("_")[-1])
        data = np.loadtxt(path)
        x, p = data[:, 0], data[:, 1]
        t = idx / max(nr - 1, 1)
        color = cmap(t)
        s0 = centers[ir - 1] if ir - 1 < len(centers) else ir
        ax.fill_between(x, 0.0, p, color=color, alpha=0.4, linewidth=0)
        ax.plot(x, p, color=color, lw=1.3, alpha=1.0, label=f"$s_0={s0:.1f}$")

    ax.set_xlim(1.8, 4.2)
    ax.set_ylim(bottom=0.0)
    ax.xaxis.set_major_locator(FixedLocator(np.arange(2.0, 4.25, 0.2)))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.set_xlabel(r"CV$_1$ (umbrella coordinate)", fontsize=14, fontweight="bold")
    ax.set_ylabel(r"$P(s)$ (arb. units)", fontsize=14, fontweight="bold")
    ax.set_title("1D probability along umbrella CV per replica window")
    ax.grid(color="grey", linestyle="-.", linewidth=0.8, alpha=0.4)
    ax.tick_params(labelsize=12)
    ax.legend(fontsize=8, ncol=2, loc="upper right", frameon=True)

    fig.savefig(out / "prob_1d.png")
    plt.close(fig)



def _smooth_grid(z: np.ndarray, sigma: float = 1.0) -> np.ndarray:
    from scipy.ndimage import gaussian_filter

    if np.ma.isMaskedArray(z):
        filled = z.filled(np.nan)
        mask = z.mask
    else:
        filled = np.asarray(z, dtype=float)
        mask = ~np.isfinite(filled)

    smoothed = gaussian_filter(np.nan_to_num(filled, nan=0.0), sigma=sigma)
    if mask.any():
        return np.ma.array(smoothed, mask=mask)
    return smoothed


def _load_interp_2d_grid(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load interp_free_energy_2D.dat on its regular grid within splot.gpi axis limits."""
    data = np.loadtxt(path)
    mask = (
        (data[:, 0] >= X2D_LO)
        & (data[:, 0] <= X2D_HI)
        & (data[:, 1] >= Y2D_LO)
        & (data[:, 1] <= Y2D_HI)
    )
    data = data[mask]
    x_vals = np.unique(data[:, 0])
    y_vals = np.unique(data[:, 1])
    z_grid = data[:, 2].reshape(len(x_vals), len(y_vals))
    return x_vals, y_vals, z_grid


def _fes_blue_white_cmap():
    """Deep blue at low ΔF (minima), fading to white at high ΔF."""
    from matplotlib.colors import LinearSegmentedColormap

    colors = [
        (0.00, "#0a0040"),
        (0.18, "#1a237e"),
        (0.38, "#1565c0"),
        (0.58, "#42a5f5"),
        (0.72, "#90caf9"),
        (0.84, "#cce5ff"),
        (0.94, "#eef6ff"),
        (1.00, "#ffffff"),
    ]
    return LinearSegmentedColormap.from_list("fes_blue_white", colors)


def plot_2d_pm3d(out: Path) -> None:
    """2D FES from interp_free_energy_2D.dat — same source and axes as example/splot.gpi."""
    path = EXAMPLE / "interp_free_energy_2D.dat"
    x_vals, y_vals, z_grid = _load_interp_2d_grid(path)
    delta_f = z_grid - z_grid.min()

    # Upsample onto a uniform grid so contour lines are smooth (raw grid is 28×1201).
    n_fine = 200
    xi = np.linspace(x_vals.min(), x_vals.max(), n_fine)
    yi = np.linspace(y_vals.min(), y_vals.max(), n_fine)
    spline = RectBivariateSpline(x_vals, y_vals, delta_f, kx=3, ky=3)
    delta_fine = _smooth_grid(spline(xi, yi, grid=True), sigma=1.0)
    x_plot, y_plot = np.meshgrid(xi, yi, indexing="ij")

    df_vmax = 25.0
    white_cutoff = 26.0
    z_plot = np.ma.masked_where(delta_fine > white_cutoff, delta_fine)

    cmap = _fes_blue_white_cmap()
    cmap.set_bad(color="white")

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.set_facecolor("white")
    surf = ax.pcolormesh(
        x_plot,
        y_plot,
        z_plot,
        cmap=cmap,
        shading="gouraud",
        vmin=0.0,
        vmax=df_vmax,
    )
    fig.colorbar(surf, ax=ax, label=r"$\Delta F$ (kcal mol$^{-1}$)")

    contour_z = np.where(z_plot.mask, np.nan, z_plot)
    ax.contour(
        x_plot,
        y_plot,
        contour_z,
        levels=np.arange(2.0, df_vmax, 2.0),
        colors="black",
        linewidths=0.55,
        alpha=0.65,
        corner_mask=False,
    )

    ax.set_xlim(X2D_LO, X2D_HI)
    ax.set_ylim(Y2D_LO, Y2D_HI)
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.set_xlabel(r"CV$_1$ (umbrella coordinate)", fontsize=14, fontweight="bold")
    ax.set_ylabel(r"CV$_3$ (TASS coordinate)", fontsize=14, fontweight="bold")
    ax.tick_params(labelsize=12)
    fig.savefig(out / "free_energy_2d.png", facecolor="white")
    plt.close(fig)


def main() -> None:
    _style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if (EXAMPLE / "free_energy.dat").exists():
        plot_1d_free_energy(OUT_DIR)
        plot_1d_with_error(OUT_DIR)
        plot_2d_pm3d(OUT_DIR)

    plot_1d_probability(OUT_DIR)

    print("Wrote plots to", OUT_DIR)
    for name in sorted(p.name for p in OUT_DIR.glob("*.png")):
        print(" ", name)


if __name__ == "__main__":
    main()
