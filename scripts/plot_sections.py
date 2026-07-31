#!/usr/bin/env python3
"""Transverse cross-sections (cuts perpendicular to the axis) of the
axisymmetric body-of-revolution mesh. Revolving the meridional mesh, a cut at
station x is a set of concentric shells (fine BL rings near the body -> poly
core -> farfield), divided azimuthally into ~square cells."""
import sys, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


def load(path):
    L = open(path).read().split("\n"); i = 0
    while not L[i].startswith("POINTS"): i += 1
    npts = int(L[i].split()[1]); i += 1
    nums = []
    while len(nums) < 3 * npts:
        nums += [float(x) for x in L[i].split()]; i += 1
    pts = [(nums[3*k], nums[3*k+1]) for k in range(npts)]
    while not L[i].startswith("CELLS"): i += 1
    nc = int(L[i].split()[1]); i += 1
    cells = []
    for _ in range(nc):
        v = [int(x) for x in L[i].split()]; i += 1
        cells.append(v[1:])
    return pts, cells


def shells_at(pts, cells, x0):
    """Radii where vertical line x=x0 crosses meridional cell edges."""
    seen = set(); rs = []
    for c in cells:
        m = len(c)
        for k in range(m):
            a, b = pts[c[k]], pts[c[(k+1) % m]]
            if (a[0] < x0) != (b[0] < x0):
                t = (x0 - a[0]) / (b[0] - a[0])
                r = a[1] + t * (b[1] - a[1])
                key = round(r, 5)
                if r >= -1e-9 and key not in seen:
                    seen.add(key); rs.append(max(r, 0.0))
    rs.sort()
    return rs


def section(ax, pts, cells, x0, band=0.015, title="", rmax=None):
    rs = shells_at(pts, cells, x0)
    if len(rs) < 2:
        ax.set_title(title + " (mimo těleso)"); ax.set_axis_off(); return
    rbody = rs[0]
    view = rmax if rmax else rs[-1]
    patches, colors = [], []
    for i in range(len(rs) - 1):
        r0, r1 = rs[i], rs[i + 1]
        if r0 > view:
            break
        rm, dr = 0.5 * (r0 + r1), (r1 - r0)
        if dr <= 1e-6:
            continue
        n = int(min(220, max(6, round(2 * math.pi * rm / dr))))
        bl = (rm - rbody) < band
        for j in range(n):
            t0 = 2 * math.pi * j / n
            t1 = 2 * math.pi * (j + 1) / n
            poly = [(r0*math.cos(t0), r0*math.sin(t0)),
                    (r1*math.cos(t0), r1*math.sin(t0)),
                    (r1*math.cos(t1), r1*math.sin(t1)),
                    (r0*math.cos(t1), r0*math.sin(t1))]
            patches.append(Polygon(poly, closed=True))
            colors.append("#ffb066" if bl else "#8fbfe8")
    pc = PatchCollection(patches, facecolors=colors, edgecolors="#333",
                         linewidths=0.25)
    ax.add_collection(pc)
    # body cross-section (solid disc)
    th = [2*math.pi*k/200 for k in range(201)]
    ax.fill([rbody*math.cos(t) for t in th], [rbody*math.sin(t) for t in th],
            color="white", zorder=3)
    ax.plot([rbody*math.cos(t) for t in th], [rbody*math.sin(t) for t in th],
            color="#333", lw=0.6, zorder=4)
    ax.set_xlim(-view, view); ax.set_ylim(-view, view); ax.set_aspect("equal")
    ax.set_title(title); ax.set_axis_off()


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "out/sphere/mesh.vtk"
    out = sys.argv[2] if len(sys.argv) > 2 else "/tmp/sph_sections.png"
    pts, cells = load(src)
    stations = [(-0.35, "x = -0.35 (u nosu)"), (-0.1, "x = -0.1 (max. tloušťka)"),
                (0.2, "x = 0.2 (u zádi)")]
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    for ax, (x0, t) in zip(axes, stations):
        section(ax, pts, cells, x0, title=t, rmax=0.16)
    fig.suptitle("Příčné řezy osově symetrickou sítí (kolmo na osu, zoom u tělesa)",
                 fontsize=14)
    plt.tight_layout()
    plt.savefig(out, dpi=140); print("wrote", out)
