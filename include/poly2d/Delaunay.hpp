#pragma once
// -----------------------------------------------------------------------------
// Delaunay -- Bowyer-Watson incremental triangulation.
//
// Simple O(n^2) implementation: robust enough for a few thousand points, which
// is plenty for meshing this kind of domain. Each triangle stores its
// circumcenter and squared circumradius; the circumcenters are the Voronoi
// vertices used by the mesher.
// -----------------------------------------------------------------------------
#include "Vec2.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <vector>

namespace poly2d {

struct Triangle {
    int a, b, c;      // indices into the point list
    Vec2 cc;          // circumcenter (a Voronoi vertex)
    double cr2;       // squared circumradius
};

struct Triangulation {
    std::vector<Vec2> pts;
    std::vector<Triangle> tris;
};

inline bool circumcircle(const Vec2& A, const Vec2& B, const Vec2& C,
                         Vec2& cc, double& cr2) {
    const double d = 2.0 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
    if (std::abs(d) < 1e-20) return false; // collinear
    const double A2 = A.x * A.x + A.y * A.y;
    const double B2 = B.x * B.x + B.y * B.y;
    const double C2 = C.x * C.x + C.y * C.y;
    cc.x = (A2 * (B.y - C.y) + B2 * (C.y - A.y) + C2 * (A.y - B.y)) / d;
    cc.y = (A2 * (C.x - B.x) + B2 * (A.x - C.x) + C2 * (B.x - A.x)) / d;
    cr2 = norm2(A - cc);
    return true;
}

inline Triangulation triangulate(std::vector<Vec2> points) {
    Triangulation T;
    T.pts = points;
    const int n = (int)points.size();
    if (n < 3) return T;

    // Bounding box -> super triangle.
    Vec2 lo{std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    Vec2 hi{-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};
    for (const auto& p : points) {
        lo.x = std::min(lo.x, p.x); lo.y = std::min(lo.y, p.y);
        hi.x = std::max(hi.x, p.x); hi.y = std::max(hi.y, p.y);
    }
    const double dx = hi.x - lo.x, dy = hi.y - lo.y;
    const double dmax = std::max(dx, dy);
    const Vec2 mid{(lo.x + hi.x) * 0.5, (lo.y + hi.y) * 0.5};
    const double R = 20.0 * (dmax + 1.0);

    std::vector<Vec2>& pts = T.pts;
    const int s0 = (int)pts.size();
    pts.push_back({mid.x - R, mid.y - R});
    pts.push_back({mid.x + R, mid.y - R});
    pts.push_back({mid.x,     mid.y + R});

    auto makeTri = [&](int a, int b, int c) -> Triangle {
        Triangle t{a, b, c, {}, 0.0};
        circumcircle(pts[a], pts[b], pts[c], t.cc, t.cr2);
        return t;
    };

    std::vector<Triangle> tris;
    tris.push_back(makeTri(s0, s0 + 1, s0 + 2));

    std::vector<char> bad;
    for (int ip = 0; ip < n; ++ip) {
        const Vec2& p = pts[ip];
        bad.assign(tris.size(), 0);

        // Edge -> count for the cavity boundary.
        std::map<std::pair<int, int>, int> edgeCount;
        auto addEdge = [&](int u, int v) {
            if (u > v) std::swap(u, v);
            edgeCount[{u, v}]++;
        };

        for (size_t it = 0; it < tris.size(); ++it) {
            const Triangle& t = tris[it];
            if (norm2(p - t.cc) <= t.cr2 * (1.0 + 1e-12)) {
                bad[it] = 1;
                addEdge(t.a, t.b);
                addEdge(t.b, t.c);
                addEdge(t.c, t.a);
            }
        }

        // Remove bad triangles (swap-erase, order irrelevant).
        std::vector<Triangle> keep;
        keep.reserve(tris.size());
        for (size_t it = 0; it < tris.size(); ++it)
            if (!bad[it]) keep.push_back(tris[it]);
        tris.swap(keep);

        // Re-triangulate the cavity: boundary edges appear exactly once.
        for (const auto& [e, cnt] : edgeCount) {
            if (cnt == 1) tris.push_back(makeTri(e.first, e.second, ip));
        }
    }

    // Drop triangles that touch the super-triangle vertices.
    T.tris.clear();
    T.tris.reserve(tris.size());
    for (const auto& t : tris) {
        if (t.a >= s0 || t.b >= s0 || t.c >= s0) continue;
        T.tris.push_back(t);
    }
    T.pts.resize(s0); // drop super vertices from the stored point list
    return T;
}

} // namespace poly2d
