"""The 2-D sheets of a phase-atlas run as PDF files, in the style of the paper's Fig. 5.

    python3 sa_plot_sheets.py <sheets_dir>

seeded_atlas_plot.m writes the inputs into <sheets_dir> and runs this script; it can also be run by hand on a
sheets folder (for example after installing the Python packages). Inputs:
  sheets_2d.json       axis labels and lattice, run label, colours, and one example profile per sheet (inset)
  sheets_2d_roots.csv  one row per quasi-steady state: sheet, grid indices, coordinates, dn2 and the flux
                       (vx, vy: the rate of change of the two plane coordinates, from the state's sources)
Outputs, in <sheets_dir>:
  <sheet>.pdf          one per sheet (file name: '+' -> p, '-' -> m, '/' -> _, so 1+/0 is 1p_0.pdf)
  all_families.pdf     every sheet on one plane
Each sheet is drawn as in Fig. 5: one translucent layer per state at its grid point (darker where several
states of the sheet share a point), unit flux arrows on every few grid points (at a point with several states,
the state with the highest dn2), the nullclines (solid: vy = 0, dashed: vx = 0; zero contours over the
sheet's grid points) and one example profile of the main variable as an inset, its grid point circled.
Needs numpy, matplotlib and contourpy (requirements_sheets.txt).
"""
import csv
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
import contourpy  # noqa: E402

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 6, "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                     "xtick.major.size": 2.2, "ytick.major.size": 2.2, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})
INK, GRID, ARROW = "#1f2430", "#9aa0aa", "#8c8c8c"
FS_TICK, FS_LABEL, FS_SUB, FS_NOTE = 6.0, 8.3, 8.3, 4.5
NULL_DASH = (4, 2.5)
N_ARROWS = 34                      # about this many arrows across each axis (every 3rd grid point of 101)


def tex(s):
    """a MATLAB TeX label such as '\\Delta n_R' as matplotlib mathtext"""
    return "$" + s.replace(" ", "\\ ") + "$" if "\\" in s or "_" in s else s


def load(sd):
    meta = json.load(open(os.path.join(sd, "sheets_2d.json")))
    for k in ("insets", "colours", "files"):          # MATLAB's jsonencode writes an empty map as []
        if not isinstance(meta.get(k), dict):
            meta[k] = {}
    if isinstance(meta.get("order"), str):
        meta["order"] = [meta["order"]]
    R = []
    with open(os.path.join(sd, "sheets_2d_roots.csv")) as fh:
        for q in csv.DictReader(fh):
            R.append(dict(root_id=q["root_id"], fam=q["family"], ix=int(q["ix"]), iy=int(q["iy"]),
                          x=float(q["x"]), y=float(q["y"]), dn2=float(q["dn2"]),
                          vx=float(q["vx"]), vy=float(q["vy"])))
    return meta, R


class Plane:
    """the lattice of the plane: grid point (ix, iy) is at (ix * hx, iy * hy) around the anchor"""

    def __init__(self, meta):
        self.hx, self.hy = meta["step"]
        self.x0, self.y0 = meta["anchor"]
        (self.ix0, self.iy0), (self.ix1, self.iy1) = meta["index_min"], meta["index_max"]
        self.AX = self.x0 + np.arange(self.ix0, self.ix1 + 1) * self.hx
        self.AY = self.y0 + np.arange(self.iy0, self.iy1 + 1) * self.hy
        self.xlim = (self.AX[0] - self.hx / 2, self.AX[-1] + self.hx / 2)
        self.ylim = (self.AY[0] - self.hy / 2, self.AY[-1] + self.hy / 2)
        n = max(len(self.AX), len(self.AY))
        self.every = max(1, int(round(n / N_ARROWS)))

    def grid(self, roots):
        """per grid point, the state with the highest dn2 -> arrays [iy, ix] of vx, vy; states per point"""
        shape = (len(self.AY), len(self.AX))
        best, cnt = {}, np.zeros(shape, int)
        for r in roots:
            k = (r["iy"] - self.iy0, r["ix"] - self.ix0)
            cnt[k] += 1
            if k not in best or r["dn2"] > best[k]["dn2"]:
                best[k] = r
        VX, VY = np.full(shape, np.nan), np.full(shape, np.nan)
        for k, r in best.items():
            VX[k], VY[k] = r["vx"], r["vy"]
        return dict(vx=VX, vy=VY, n=cnt, null_x=self.zero_lines(VX), null_y=self.zero_lines(VY))

    def zero_lines(self, F):
        if np.sum(np.isfinite(F)) < 4:
            return []
        m = np.ma.masked_invalid(F)
        g = contourpy.contour_generator(self.AX, self.AY, m, line_type="Separate")
        return [s for s in g.lines(0.0) if len(s) > 1]

    def exist_mesh(self, ax, n, col, alpha):
        """one translucent layer per state: a grid point with k states is covered k times"""
        ex, ey = np.r_[self.AX - self.hx / 2, self.AX[-1] + self.hx / 2], np.r_[self.AY - self.hy / 2,
                                                                                 self.AY[-1] + self.hy / 2]
        for k in range(1, int(n.max()) + 1 if n.size else 1):
            m = np.ma.masked_where(n < k, np.ones_like(n, float))
            ax.pcolormesh(ex, ey, m, cmap=ListedColormap([col]), alpha=alpha, shading="flat", zorder=1,
                          linewidth=0, rasterized=True)

    def arrows(self, ax, g, w_in, h_in, frac=0.45, shaft_pt=0.4):
        """unit arrows, direction true on the page, length frac x the drawn spacing of the arrows"""
        sx = (np.arange(self.ix0, self.ix1 + 1) % self.every) == 0
        sy = (np.arange(self.iy0, self.iy1 + 1) % self.every) == 0
        VX, VY = g["vx"][np.ix_(sy, sx)], g["vy"][np.ix_(sy, sx)]
        XX, YY = np.meshgrid(self.AX[sx], self.AY[sy])
        ok = np.isfinite(VX) & np.isfinite(VY)
        if not ok.any():
            return
        kx = w_in / (self.xlim[1] - self.xlim[0]); ky = h_in / (self.ylim[1] - self.ylim[0])   # inch per unit
        ux, uy = VX[ok] * kx, VY[ok] * ky
        m = np.hypot(ux, uy); m[m == 0] = 1
        L = frac * self.every * min(self.hx * kx, self.hy * ky)                                 # arrow length (in)
        ax.quiver(XX[ok], YY[ok], L * ux / m / kx, L * uy / m / ky, angles="xy", scale_units="xy", scale=1.0,
                  pivot="middle", color=ARROW, width=shaft_pt / 72 / w_in, headwidth=3.5, headlength=4,
                  headaxislength=3.5, zorder=2)

    def axes(self, ax, meta):
        ax.set_xlim(*self.xlim); ax.set_ylim(*self.ylim)
        ax.axhline(0, color=GRID, lw=0.36, zorder=0); ax.axvline(0, color=GRID, lw=0.36, zorder=0)
        ax.tick_params(labelsize=FS_TICK, pad=1.5)
        ax.set_xlabel(tex(meta["labels"][0]), fontsize=FS_LABEL, labelpad=1)
        ax.set_ylabel(tex(meta["labels"][1]), fontsize=FS_LABEL, labelpad=1)


def nulls(ax, g, col, lw):
    for s in g["null_y"]:
        ax.plot(s[:, 0], s[:, 1], color=col, lw=lw, ls="-", zorder=4, solid_capstyle="round")
    for s in g["null_x"]:
        ax.plot(s[:, 0], s[:, 1], color=col, lw=lw, ls=(0, NULL_DASH), zorder=4)


def inset(ax, P, meta, fam, col, rs):
    """the example profile of the sheet (chosen by seeded_atlas_plot), its grid point circled; placed in the
    corner of the plot farthest from the sheet's states"""
    e = meta["insets"].get(fam)
    if not e:
        return
    xs, ms = np.asarray(e["x"], float), np.asarray(e["m"], float)
    ax.plot(e["at"][0], e["at"][1], marker="o", ms=3.4, mfc="none", mec=INK, mew=0.7, zorder=9)
    sx, sy = np.array([r["x"] for r in rs]), np.array([r["y"] for r in rs])
    u = (sx - P.xlim[0]) / (P.xlim[1] - P.xlim[0]); v = (sy - P.ylim[0]) / (P.ylim[1] - P.ylim[0])
    corners = {(0.08, 0.62): (0.24, 0.77), (0.64, 0.62): (0.80, 0.77), (0.08, 0.12): (0.24, 0.27),
               (0.64, 0.12): (0.80, 0.27)}
    pos = max(corners, key=lambda c: np.min(np.hypot(u - corners[c][0], v - corners[c][1])) if len(u) else 0)
    ia = ax.inset_axes([pos[0], pos[1], 0.32, 0.30])
    Lc = meta["L"]
    ia.plot(xs, ms, color=col, lw=0.8)
    ia.axvline(Lc / 2, color=GRID, lw=0.4, ls=(0, (2, 2)))
    ia.set_xlim(0, Lc); ia.set_ylim(0, max(1e-12, 1.15 * ms.max()))
    ia.set_xticks([0, Lc / 2, Lc]); ia.set_xticklabels(["0", "", f"{Lc:g}"])
    ia.tick_params(labelsize=FS_NOTE, length=1.5, pad=1, width=0.4)
    ia.set_xlabel("$x$", fontsize=FS_NOTE, labelpad=0)
    ia.set_ylabel(tex(meta["main_variable"]), fontsize=FS_NOTE, labelpad=4, rotation=0, va="center")
    for sp in ("top", "right"):
        ia.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ia.spines[sp].set_linewidth(0.4)


def note(fig, meta, what):
    fig.text(0.005, 0.005, f"{meta['title']} | {what}", fontsize=FS_NOTE - 1, color=GRID, ha="left",
             va="bottom", wrap=True)


def sheet_pdf(sd, meta, P, R, fam):
    rs = [r for r in R if r["fam"] == fam]
    g = P.grid(rs)
    col = meta["colours"][fam]
    aw, ah, LM, BOT, TOP, RM = 3.0, 2.2, 0.5, 0.42, 0.28, 0.1
    W, H = LM + aw + RM, BOT + ah + TOP
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([LM / W, BOT / H, aw / W, ah / H])
    P.axes(ax, meta)
    ax.set_title(f"{fam}  ({len(rs)} states)", fontsize=FS_SUB, color=INK, pad=3)
    P.exist_mesh(ax, g["n"], col, 0.30)
    P.arrows(ax, g, aw, ah)
    nulls(ax, g, col, 0.9)
    inset(ax, P, meta, fam, col, rs)
    note(fig, meta, f"sheet {fam}")
    fig.savefig(os.path.join(sd, meta["files"][fam]))
    plt.close(fig)


def all_pdf(sd, meta, P, R):
    order = meta["order"]
    aw, ah, LM, BOT, TOP, RM = 3.6, 2.6, 0.5, 0.42, 0.28, 2.0
    W, H = LM + aw + RM, BOT + ah + TOP
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([LM / W, BOT / H, aw / W, ah / H])
    P.axes(ax, meta)
    G = {f: P.grid([r for r in R if r["fam"] == f]) for f in order}
    for f in order:
        P.exist_mesh(ax, G[f]["n"], meta["colours"][f], 0.07)
    P.arrows(ax, P.grid(R), aw, ah)
    for f in order:
        nulls(ax, G[f], meta["colours"][f], 1.17)
    h = [plt.Rectangle((0, 0), 1, 1, fc=meta["colours"][f], alpha=0.35, ec="none") for f in order]
    lab = [f"{f} ({int(G[f]['n'].sum())})" for f in order]
    lx, ly = tex(meta["labels"][0]), tex(meta["labels"][1])
    h += [plt.Line2D([], [], color=INK, lw=1.17), plt.Line2D([], [], color=INK, lw=1.17, ls=(0, NULL_DASH)),
          plt.Line2D([], [], color=ARROW, marker=r"$\rightarrow$", ls="", ms=7)]
    lab += [f"d{ly}/dt = 0", f"d{lx}/dt = 0", "flux"]
    ax.legend(h, lab, fontsize=FS_TICK, frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              handlelength=1.8, handletextpad=0.5, labelspacing=0.4, ncol=1 if len(order) <= 20 else 2)
    note(fig, meta, "all sheets")
    fig.savefig(os.path.join(sd, meta["files"]["all_families"]))
    plt.close(fig)


def main(sd):
    meta, R = load(sd)
    P = Plane(meta)
    for f in meta["order"]:
        sheet_pdf(sd, meta, P, R, f)
    all_pdf(sd, meta, P, R)
    print(f"sheets {sd}: {len(meta['order']) + 1} PDF files")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python3 sa_plot_sheets.py <sheets_dir>")
    main(sys.argv[1])
