#!/usr/bin/env python3
"""3D views of the axisymmetric case: revolve the NACA 0012 meridional mesh
around the x-axis to show the body of revolution inside the sphere."""
import sys, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np


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


def naca_yt(x):
    t = 0.12
    return 5*t*(0.2969*math.sqrt(x) - 0.126*x - 0.3516*x*x + 0.2843*x**3 - 0.1036*x**4)


def body_curve(n=160, chord=1.0, x0=-0.5):
    xs, rs = [], []
    for i in range(n + 1):
        xl = 0.5 * (1 - math.cos(math.pi * i / n))
        xs.append(x0 + chord * xl); rs.append(chord * naca_yt(xl))
    return np.array(xs), np.array(rs)


def revolve(x, r, nth, th0=0.0, th1=2*math.pi):
    th = np.linspace(th0, th1, nth)
    X = np.repeat(x[:, None], nth, axis=1)
    Y = r[:, None] * np.cos(th)[None, :]
    Z = r[:, None] * np.sin(th)[None, :]
    return X, Y, Z


def mesh_slice(pts, cells, sign=1.0):
    segs = []
    for c in cells:
        m = len(c)
        for k in range(m):
            a, b = pts[c[k]], pts[c[(k+1) % m]]
            segs.append([(a[0], sign*a[1], 0.0), (b[0], sign*b[1], 0.0)])
    return segs


def set_axes_equal(ax, R):
    ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_zlim(-R, R)
    try: ax.set_box_aspect((1, 1, 1))
    except Exception: pass


def overview(path, out, R=5.0):
    pts, cells = load(path)
    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection="3d")

    # translucent sphere (farfield), cut open 90 deg so the interior is visible
    u = np.linspace(0, 1.5*math.pi, 60); v = np.linspace(0, math.pi, 30)
    sx = R*np.outer(np.cos(u), np.sin(v))
    sy = R*np.outer(np.sin(u), np.sin(v))
    sz = R*np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(sz, sx, sy, color="#4a90d9", alpha=0.10, linewidth=0)

    # body of revolution
    x, r = body_curve()
    X, Y, Z = revolve(x, r, 80)
    ax.plot_surface(X, Y, Z, color="#d98a3a", alpha=1.0, linewidth=0, antialiased=True)

    # meridional mesh on the cut plane (z=0, y=r>=0)
    lc = Line3DCollection(mesh_slice(pts, cells), colors="#333", linewidths=0.25)
    ax.add_collection3d(lc)

    set_axes_equal(ax, R)
    ax.set_title("NACA 0012 body of revolution in a sphere (axisymmetric)")
    ax.view_init(elev=22, azim=-60)
    ax.set_axis_off()
    plt.tight_layout(); plt.savefig(out, dpi=140); print("wrote", out)


def bodyview(path, out):
    pts, cells = load(path)
    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_subplot(111, projection="3d")

    x, r = body_curve()
    # body surface (270 deg so the BL/mesh cut is visible)
    X, Y, Z = revolve(x, r, 90, 0, 1.5*math.pi)
    ax.plot_surface(X, Y, Z, color="#d98a3a", alpha=1.0, linewidth=0)

    # boundary-layer shell (body offset by the prism-band thickness), translucent
    band = 0.015
    Xb, Yb, Zb = revolve(x, np.where(r > 1e-6, r + band, r), 90, 0, 1.5*math.pi)
    ax.plot_surface(Xb, Yb, Zb, color="#8fbfe8", alpha=0.18, linewidth=0)

    # meridional mesh slice near the body, on the cut plane
    lc = Line3DCollection(mesh_slice(pts, cells), colors="#333", linewidths=0.35)
    ax.add_collection3d(lc)

    B = 0.9
    ax.set_xlim(-B, B); ax.set_ylim(-B, B); ax.set_zlim(-B, B)
    try: ax.set_box_aspect((1, 1, 1))
    except Exception: pass
    ax.set_title("Body of revolution + boundary-layer shell (cut view)")
    ax.view_init(elev=18, azim=-70)
    ax.set_axis_off()
    plt.tight_layout(); plt.savefig(out, dpi=140); print("wrote", out)


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "out/sphere/mesh.vtk"
    overview(src, sys.argv[2] if len(sys.argv) > 2 else "/tmp/sph3d_over.png")
    bodyview(src, sys.argv[3] if len(sys.argv) > 3 else "/tmp/sph3d_body.png")
