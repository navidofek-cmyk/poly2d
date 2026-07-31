#!/usr/bin/env python3
"""Dependency-free sanity check of an OpenFOAM polyMesh.

Verifies the invariants OpenFOAM's checkMesh cares about most:
  * label / list sizes are consistent,
  * owner < neighbour on every internal face,
  * boundary patches tile [nInternal, nFaces) contiguously,
  * every cell is geometrically closed  ->  sum of outward face-area
    vectors per cell is ~0  (this is what catches owner/neighbour
    orientation bugs).
"""
import re
import sys


def read_foam_list(path):
    txt = open(path).read()
    txt = re.sub(r'/\*.*?\*/', '', txt, flags=re.S)
    txt = re.sub(r'//.*', '', txt)
    # drop FoamFile header block
    b = txt.find('}')
    txt = txt[b + 1:]
    return txt


def read_points(path):
    txt = read_foam_list(path)
    pts = re.findall(r'\(([^()]*)\)', txt)
    out = []
    for p in pts:
        c = p.split()
        if len(c) == 3:
            out.append((float(c[0]), float(c[1]), float(c[2])))
    return out


def read_faces(path):
    txt = read_foam_list(path)
    faces = []
    for m in re.finditer(r'(\d+)\s*\(([\d\s]+)\)', txt):
        faces.append([int(x) for x in m.group(2).split()])
    return faces


def read_labels(path):
    txt = read_foam_list(path)
    m = re.search(r'\((.*)\)', txt, flags=re.S)
    return [int(x) for x in m.group(1).split()]


def read_boundary(path):
    txt = read_foam_list(path)
    patches = []
    for m in re.finditer(r'(\w+)\s*\{([^}]*)\}', txt, flags=re.S):
        body = m.group(2)
        t = re.search(r'type\s+(\w+)', body)
        nf = re.search(r'nFaces\s+(\d+)', body)
        sf = re.search(r'startFace\s+(\d+)', body)
        if nf and sf:
            patches.append((m.group(1), t.group(1) if t else '?',
                            int(nf.group(1)), int(sf.group(1))))
    return patches


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def main(d):
    pts = read_points(d + '/points')
    faces = read_faces(d + '/faces')
    owner = read_labels(d + '/owner')
    neigh = read_labels(d + '/neighbour')
    bnd = read_boundary(d + '/boundary')
    nF = len(faces)
    nInt = len(neigh)
    nCells = max(max(owner), max(neigh) if neigh else 0) + 1
    ok = True

    def chk(cond, msg):
        nonlocal ok
        print(('  OK  ' if cond else ' FAIL ') + msg)
        ok = ok and cond

    print(f"[{d}]")
    print(f"  points={len(pts)} faces={nF} internal={nInt} cells={nCells} patches={len(bnd)}")
    chk(len(owner) == nF, "len(owner) == nFaces")
    chk(all(0 <= min(f) and max(f) < len(pts) for f in faces), "face point ids in range")
    chk(all(len(f) >= 3 for f in faces), "faces have >= 3 points")
    chk(all(owner[i] < neigh[i] for i in range(nInt)), "owner < neighbour (internal)")

    # boundary patches tile [nInt, nF)
    exp = nInt
    tile = True
    for name, t, nf, sf in bnd:
        tile = tile and (sf == exp)
        exp += nf
    chk(tile and exp == nF, "boundary patches tile [nInternal, nFaces) contiguously")

    # face use count: internal 2, boundary 1
    use = [0] * nCells
    for i in range(nF):
        use[owner[i]] += 1
        if i < nInt:
            use[neigh[i]] += 1
    # closedness: sum of signed area vectors per cell ~ 0
    cellVec = [[0.0, 0.0, 0.0] for _ in range(nCells)]
    maxbad = 0.0
    for i, f in enumerate(faces):
        p0 = pts[f[0]]
        n = [0.0, 0.0, 0.0]
        for k in range(1, len(f) - 1):
            a = tuple(pts[f[k]][j] - p0[j] for j in range(3))
            b = tuple(pts[f[k + 1]][j] - p0[j] for j in range(3))
            c = cross(a, b)
            n = [n[j] + 0.5 * c[j] for j in range(3)]
        for j in range(3):
            cellVec[owner[i]][j] += n[j]      # owner: outward
        if i < nInt:
            for j in range(3):
                cellVec[neigh[i]][j] -= n[j]   # neighbour: inward
    for c in range(nCells):
        m = max(abs(v) for v in cellVec[c])
        maxbad = max(maxbad, m)
    chk(maxbad < 1e-9, f"cells geometrically closed (max |sum area| = {maxbad:.2e})")

    for name, t, nf, sf in bnd:
        print(f"    patch {name:<14} type={t:<8} nFaces={nf}")
    print("  RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else 'out/rect/constant/polyMesh'))
