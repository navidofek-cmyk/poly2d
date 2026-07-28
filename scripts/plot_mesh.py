#!/usr/bin/env python3
"""Render a poly2d VTK polygon mesh to PNG (prism = orange, core = blue)."""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection


def load(path):
    L = open(path).read().split("\n")
    i = 0
    while not L[i].startswith("POINTS"):
        i += 1
    npts = int(L[i].split()[1]); i += 1
    nums = []
    while len(nums) < 3 * npts:
        nums += [float(x) for x in L[i].split()]; i += 1
    pts = [(nums[3 * k], nums[3 * k + 1]) for k in range(npts)]
    while not L[i].startswith("CELLS"):
        i += 1
    nc = int(L[i].split()[1]); i += 1
    cells = []
    for _ in range(nc):
        v = [int(x) for x in L[i].split()]; i += 1
        cells.append(v[1:])
    region = []
    while i < len(L) and not L[i].startswith("SCALARS region"):
        i += 1
    if i < len(L):
        i += 2  # skip SCALARS + LOOKUP_TABLE
        while len(region) < nc:
            region += [int(x) for x in L[i].split()]; i += 1
    return pts, cells, region


def plot(path, out, xlim=None, ylim=None, lw=0.3):
    pts, cells, region = load(path)
    polys = [[pts[j] for j in c] for c in cells]
    colors = ["#ffb066" if (region and region[k] == 1) else "#8fbfe8"
              for k in range(len(cells))]
    fig, ax = plt.subplots(figsize=(14, 6))
    pc = PolyCollection(polys, facecolors=colors, edgecolors="#222", linewidths=lw)
    ax.add_collection(pc)
    if xlim:
        ax.set_xlim(*xlim)
    else:
        ax.set_xlim(min(p[0] for p in pts), max(p[0] for p in pts))
    if ylim:
        ax.set_ylim(*ylim)
    else:
        ax.set_ylim(min(p[1] for p in pts), max(p[1] for p in pts))
    ax.set_aspect("equal")
    ax.set_title(f"{path}   cells={len(cells)}")
    plt.tight_layout()
    plt.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    path = sys.argv[1]
    out = sys.argv[2]
    xlim = ylim = None
    if len(sys.argv) > 6:
        xlim = (float(sys.argv[3]), float(sys.argv[4]))
        ylim = (float(sys.argv[5]), float(sys.argv[6]))
    lw = float(sys.argv[7]) if len(sys.argv) > 7 else 0.3
    plot(path, out, xlim, ylim, lw)
