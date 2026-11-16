"""Checks (and optionally regenerates) the (n, T) triangulation used for the spectral-function interpolation.

TriangFiles/nTs.dat      180 points, columns (T [MeV], n); 15 temperatures (40, 50, ..., 180) with 12 densities each, row after row
TriangFiles/triangles.dat 308 triangles (three point indices each)

The shipped triangulation has a strip structure: between two neighbouring temperature rows the 11 quadrilaterals
(a_j, a_j+1, b_j+1, b_j) are each cut into two triangles along one of the two diagonals (14 strips x 11 quads x 2 = 308).
The original notebook that produced `triangles.dat` is not part of this repository and the diagonal choice cannot be reproduced by a
simple rule (shorter diagonal, Delaunay), so the shipped file should be kept. This script
  python3 tools/triangulation.py            verifies the shipped file (every quad covered exactly twice, no overlaps, total area)
  python3 tools/triangulation.py --write f  writes a new strip triangulation to f, cutting each quad along its shorter diagonal
                                            (in coordinates T/10 MeV and n/n_max-of-row-pair; change `scale` below if needed)
"""
import sys
import numpy as np

NT, NN = 15, 12


def load(path="TriangFiles/"):
    return np.loadtxt(path + "nTs.dat"), np.loadtxt(path + "triangles.dat").astype(int)


def area(P, t):
    a, b, c = P[t[0]], P[t[1]], P[t[2]]
    return 0.5 * abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def quads():
    for s in range(NT - 1):
        for j in range(NN - 1):
            yield s, j, 12 * s + j, 12 * s + j + 1, 12 * (s + 1) + j, 12 * (s + 1) + j + 1


def check(P, tri):
    ref = {tuple(sorted(t)) for t in tri}
    bad = 0
    for s, j, a0, a1, b0, b1 in quads():
        d0 = {tuple(sorted((a0, a1, b1))), tuple(sorted((a0, b0, b1)))}
        d1 = {tuple(sorted((a0, a1, b0))), tuple(sorted((a1, b0, b1)))}
        if not (d0 <= ref or d1 <= ref):
            bad += 1
            print("quad not covered by two triangles of one diagonal: strip %d, column %d" % (s, j))
    tot = sum(area(P, t) for t in tri)
    quad_area = sum(area(P, (a0, a1, b1)) + area(P, (a0, b0, b1)) for _, _, a0, a1, b0, b1 in quads())
    print("triangles: %d (expected %d), duplicates: %d, bad quads: %d" % (len(tri), (NT - 1) * (NN - 1) * 2, len(tri) - len(ref), bad))
    print("sum of triangle areas %.6g, sum of quad areas %.6g" % (tot, quad_area))
    return bad == 0 and len(ref) == len(tri) == (NT - 1) * (NN - 1) * 2 and abs(tot - quad_area) < 1e-6 * tot


def build(P, scale=(0.1, None)):
    out = []
    for s, j, a0, a1, b0, b1 in quads():
        nmax = max(P[[a0, a1, b0, b1], 1])
        X = lambda i: np.array([P[i, 0] * scale[0], P[i, 1] / nmax])
        d0, d1 = np.linalg.norm(X(a0) - X(b1)), np.linalg.norm(X(a1) - X(b0))
        out += [(a0, a1, b1), (a0, b0, b1)] if d0 <= d1 else [(a0, a1, b0), (a1, b0, b1)]
    return np.array(out)


if __name__ == "__main__":
    P, tri = load()
    ok = check(P, tri)
    print("shipped triangulation consistent:", ok)
    if "--write" in sys.argv:
        f = sys.argv[sys.argv.index("--write") + 1]
        new = build(P)
        np.savetxt(f, new, fmt="%d")
        print("wrote", f, "consistent:", check(P, new))
