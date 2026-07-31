#!/usr/bin/env python3
"""3D view of the extruded duct: the 2D cross-section mesh extruded along z into
a block pierced by a square through-hole."""
import sys, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from collections import defaultdict


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


def main(src, out, Lz=3.0):
    pts, cells = load(src)
    ec = defaultdict(int)
    for c in cells:
        m = len(c)
        for k in range(m):
            u, v = c[k], c[(k+1) % m]
            ec[(min(u, v), max(u, v))] += 1
    boundary = [e for e, n in ec.items() if n == 1]

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection="3d")

    # front cross-section (z=0): filled cells
    front = [[(pts[j][0], pts[j][1], 0.0) for j in c] for c in cells]
    ax.add_collection3d(Poly3DCollection(front, facecolors="#cfe3f7",
                                         edgecolors="#37506b", linewidths=0.2))
    # back cross-section (z=Lz): outline only
    back = [[(pts[j][0], pts[j][1], Lz) for j in c] for c in cells]
    ax.add_collection3d(Poly3DCollection(back, facecolors="#eef4fb",
                                         edgecolors="#9bb3cc", linewidths=0.15))
    # walls (outer box + hole tunnel): extruded boundary edges
    walls = []
    for u, v in boundary:
        a, b = pts[u], pts[v]
        walls.append([(a[0], a[1], 0), (b[0], b[1], 0), (b[0], b[1], Lz), (a[0], a[1], Lz)])
    ax.add_collection3d(Poly3DCollection(walls, facecolors="#f0c48a",
                                         edgecolors="#7a5a2a", linewidths=0.12, alpha=0.35))

    ax.set_xlim(0, 2); ax.set_ylim(0, 2); ax.set_zlim(0, Lz)
    try: ax.set_box_aspect((2, 2, Lz))
    except Exception: pass
    ax.view_init(elev=32, azim=-52)
    ax.set_title("Kvádr se čtvercovou dírou skrz — 3D extruze průřezu")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z (proud)")
    plt.tight_layout(); plt.savefig(out, dpi=140); print("wrote", out)


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "out/duct/mesh.vtk"
    out = sys.argv[2] if len(sys.argv) > 2 else "/tmp/duct3d.png"
    main(src, out)
