"""Debug figures for the golden cases.

One PNG per case. One row per engine. Columns:
  1. Section (XZ at y = 0.5): solved field, truth horizons dashed black, solved horizons in colour,
     data near the section as dots (value) and arrows (normal).
  2. Section: gradient direction error in degrees (0 good, 180 reversed). Same horizons.
  3. Section: unit gradient arrows, solved (blue) over truth (grey).
  4. Plan (XY at z = 0.5): signed value error. Blue = solved too low, red = too high.
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from . import geomodels  # noqa: E402
from .harness import OUT_DIR, angle_deg  # noqa: E402

N = 61  # section resolution
SECTION_Y = 0.5
PLAN_Z = 0.5
SLAB = 0.08  # data within this distance of a section are drawn


def _plane_grid(kind: str):
    u = np.linspace(0, 1, N)
    a, b = np.meshgrid(u, u)
    if kind == "xz":
        pts = np.c_[a.ravel(), np.full(a.size, SECTION_Y), b.ravel()]
    else:
        pts = np.c_[a.ravel(), b.ravel(), np.full(a.size, PLAN_Z)]
    return a, b, pts


def _horizons(ax, a, b, field, colors, styles, widths=1.2):
    if np.isfinite(field).sum() < 4:
        return
    ax.contour(a, b, field.reshape(a.shape), levels=geomodels.HORIZONS, colors=colors,
               linestyles=styles, linewidths=widths)


def _data(ax, case, kind):
    v, n = case.values, case.normals
    if kind == "xz":
        mv, mn = np.abs(v[:, 1] - SECTION_Y) < SLAB, np.abs(n[:, 1] - SECTION_Y) < SLAB
        ax.plot(v[mv, 0], v[mv, 2], "o", ms=3.5, mfc="white", mec="k", mew=0.8, zorder=5)
        ax.plot(v[~mv, 0], v[~mv, 2], "o", ms=2.5, mfc="none", mec="0.5", mew=0.6, zorder=5)
        if mn.any():
            nn = n[mn]
            u = nn[:, 3:6] / np.linalg.norm(nn[:, 3:6], axis=1, keepdims=True)
            ax.quiver(nn[:, 0], nn[:, 2], u[:, 0], u[:, 2], color="k", scale=22, width=0.006, zorder=6)
        ax.plot(n[~mn, 0], n[~mn, 2], "^", ms=3, mfc="none", mec="0.5", mew=0.6, zorder=5)
    else:
        mv, mn = np.abs(v[:, 2] - PLAN_Z) < SLAB, np.abs(n[:, 2] - PLAN_Z) < SLAB
        ax.plot(v[mv, 0], v[mv, 1], "o", ms=3.5, mfc="white", mec="k", mew=0.8, zorder=5)
        ax.plot(n[mn, 0], n[mn, 1], "^", ms=4, mfc="white", mec="k", mew=0.8, zorder=5)


def render_case(case_name: str, results: dict[str, dict], out_dir: Path | None = None,
                suffix: str = "") -> Path:
    """`results` maps engine name -> harness.compute() result."""
    case = geomodels.build(case_name)
    out_dir = Path(out_dir or OUT_DIR / "figures")
    out_dir.mkdir(parents=True, exist_ok=True)
    engines = list(results)
    fig, axes = plt.subplots(len(engines), 4, figsize=(15, 3.55 * len(engines)), squeeze=False,
                             constrained_layout=True)
    xa, za, xz = _plane_grid("xz")
    xb, yb, xy = _plane_grid("xy")
    truth_xz, tgrad_xz = case.truth(xz), case.truth_grad(xz)
    truth_xy = case.truth(xy)
    fmin, fmax = float(truth_xz.min()), float(truth_xz.max())
    arrow = slice(None, None, 6)

    for row, eng in enumerate(engines):
        interp, m = results[eng]["interp"], results[eng]["metrics"]
        v_xz = np.asarray(interp.evaluate_value(xz), float)
        g_xz = np.asarray(interp.evaluate_gradient(xz), float)
        v_xy = np.asarray(interp.evaluate_value(xy), float)
        ax1, ax2, ax3, ax4 = axes[row]

        im = ax1.contourf(xa, za, v_xz.reshape(xa.shape), levels=np.linspace(fmin, fmax, 21),
                          cmap="cividis", extend="both")
        _horizons(ax1, xa, za, truth_xz, "k", "--", 1.4)
        _horizons(ax1, xa, za, v_xz, "#ff7f0e", "-", 1.6)
        _data(ax1, case, "xz")
        ax1.set_title(f"{eng}: field (y={SECTION_Y})   value rms {m['value_rms']:.3g}", fontsize=9)

        ang = angle_deg(g_xz, tgrad_xz).reshape(xa.shape)
        im2 = ax2.pcolormesh(xa, za, ang, cmap="magma_r", vmin=0, vmax=45, shading="auto")
        _horizons(ax2, xa, za, truth_xz, "k", "--", 1.0)
        ax2.set_title(f"gradient error  median {m['angle_median']:.1f}  p90 {m['angle_p90']:.1f} deg", fontsize=9)

        gu = g_xz / np.linalg.norm(g_xz, axis=1, keepdims=True)
        tu = tgrad_xz / np.linalg.norm(tgrad_xz, axis=1, keepdims=True)
        _horizons(ax3, xa, za, truth_xz, "0.6", "-", 0.8)
        sx, sz = xa[arrow, arrow], za[arrow, arrow]
        shape = xa.shape
        ax3.quiver(sx, sz, tu[:, 0].reshape(shape)[arrow, arrow], tu[:, 2].reshape(shape)[arrow, arrow],
                   color="0.55", scale=18, width=0.007)
        ax3.quiver(sx, sz, gu[:, 0].reshape(shape)[arrow, arrow], gu[:, 2].reshape(shape)[arrow, arrow],
                   color="#1f77b4", scale=18, width=0.004)
        ax3.set_title("unit gradient: solved (blue) over truth (grey)", fontsize=9)

        e = (v_xy - truth_xy).reshape(xb.shape)
        lim = max(0.02, float(np.nanpercentile(np.abs(e), 98)))
        im4 = ax4.pcolormesh(xb, yb, e, cmap="RdBu_r", vmin=-lim, vmax=lim, shading="auto")
        _data(ax4, case, "xy")
        ax4.set_title(f"plan z={PLAN_Z}: value error (solved - truth)", fontsize=9)

        for ax in (ax1, ax2, ax3):
            ax.set_xlabel("x")
            ax.set_ylabel("z")
        ax4.set_xlabel("x")
        ax4.set_ylabel("y")
        for ax in axes[row]:
            ax.set_aspect("equal")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
        fig.colorbar(im, ax=ax1, shrink=0.8, label="s")
        fig.colorbar(im2, ax=ax2, shrink=0.8, label="deg")
        fig.colorbar(im4, ax=ax4, shrink=0.8, label="error")

    fig.suptitle(f"{case.title}\n{case.why}", fontsize=10, ha="left", x=0.01)
    path = out_dir / f"{case_name}{suffix}.png"
    fig.savefig(path, dpi=80)
    plt.close(fig)
    return path


def render_diff(case_name: str, engine: str, points: np.ndarray, golden: dict, current: dict,
                out_dir: Path | None = None) -> Path:
    """Where did the field move? Maps the |current - golden| value and angle change on two sections."""
    out_dir = Path(out_dir or OUT_DIR / "failures")
    out_dir.mkdir(parents=True, exist_ok=True)
    n = round(len(points) ** (1 / 3))
    shape = (n, n, n)  # [z, y, x]
    dv = (current["value"] - golden["value"]).reshape(shape)
    da = angle_deg(current["gradient"], golden["gradient"]).reshape(shape)
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    i = n // 2
    views = {"XZ (y mid)": (dv[:, i, :], da[:, i, :]), "XY (z mid)": (dv[i, :, :], da[i, :, :]),
             "YZ (x mid)": (dv[:, :, i], da[:, :, i])}
    for col, (name, (v2, a2)) in enumerate(views.items()):
        lim = max(1e-12, float(np.nanmax(np.abs(v2))))
        p = axes[0, col].imshow(v2, origin="lower", extent=(0, 1, 0, 1), cmap="RdBu_r", vmin=-lim, vmax=lim)
        axes[0, col].set_title(f"{name}: value change (current - golden)")
        fig.colorbar(p, ax=axes[0, col], shrink=0.8)
        q = axes[1, col].imshow(a2, origin="lower", extent=(0, 1, 0, 1), cmap="magma_r", vmin=0,
                                vmax=max(1e-6, float(np.nanmax(a2))))
        axes[1, col].set_title(f"{name}: gradient direction change (deg)")
        fig.colorbar(q, ax=axes[1, col], shrink=0.8)
    fig.suptitle(f"{case_name} / {engine}: current vs golden", fontsize=11)
    path = out_dir / f"{case_name}__{engine}__diff.png"
    fig.savefig(path, dpi=80)
    plt.close(fig)
    return path


def render_surfaces(case_name: str, engine: str, interp, n: int = 41) -> bytes:
    """3D view of the solved horizons (colour per horizon) over the true horizons (black mesh).

    White spheres are value data, cones are normals. Returns PNG bytes.
    """
    import pyvista as pv

    case = geomodels.build(case_name)
    eps = 1e-3  # P1/DFI return NaN exactly on the lower faces (finding F2)
    u = np.linspace(eps, 1 - eps, n)
    zz, yy, xx = np.meshgrid(u, u, u, indexing="ij")
    pts = np.c_[xx.ravel(), yy.ravel(), zz.ravel()]
    grid = pv.ImageData(dimensions=(n, n, n), spacing=((1 - 2 * eps) / (n - 1),) * 3, origin=(eps,) * 3)
    grid["solved"] = np.asarray(interp.evaluate_value(pts), float)
    grid["truth"] = case.truth(pts)

    pl = pv.Plotter(off_screen=True, window_size=(760, 560))
    pl.set_background("white")
    pl.add_mesh(pv.Box(bounds=(0, 1, 0, 1, 0, 1)), style="wireframe", color="#999999", line_width=1)
    colors = ("#1b9e77", "#d95f02", "#7570b3")
    for level, color in zip(geomodels.HORIZONS, colors):
        solved = grid.contour([level], scalars="solved")
        if solved.n_points:
            pl.add_mesh(solved, color=color, opacity=0.85, smooth_shading=True)
        truth = grid.contour([level], scalars="truth")
        if truth.n_points:
            pl.add_mesh(truth, color="black", style="wireframe", opacity=0.3, line_width=1)
    pl.add_mesh(pv.PolyData(case.values[:, :3]), color="white", point_size=9,
                render_points_as_spheres=True)
    pl.add_mesh(pv.PolyData(case.values[:, :3]), color="black", point_size=11,
                render_points_as_spheres=True, opacity=0.35)
    pl.add_arrows(case.normals[:, :3], case.normals[:, 3:6] / np.linalg.norm(case.normals[:, 3:6], axis=1, keepdims=True),
                  mag=0.12, color="red")
    pl.camera_position = [(2.1, -1.7, 1.5), (0.5, 0.5, 0.5), (0, 0, 1)]
    img = pl.screenshot(return_img=True)
    pl.close()
    import io

    from matplotlib import image as mpimg

    buf = io.BytesIO()
    mpimg.imsave(buf, img, format="png")
    return buf.getvalue()
