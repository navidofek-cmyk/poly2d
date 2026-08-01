#pragma once
// -----------------------------------------------------------------------------
// Airfoils -- contour generators for the mesher.
//   * naca4()      : cambered NACA 4-digit (e.g. 0012, 4412), sharp or blunt TE
//   * loadSeligDat(): read a Selig .dat coordinate file (TE->top->LE->bottom->TE)
//   * placeAirfoil(): scale by chord, rotate by angle of attack, translate
//   * truncateTE() : clip at x=teCut -> small blunt trailing edge
//
// All generators return an AirfoilShape whose points form a closed loop; blunt
// shapes also carry the two trailing-edge corner nodes so the mesher can skip
// prism on the (short) TE base.
// -----------------------------------------------------------------------------
#include "Vec2.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

namespace poly2d {

struct AirfoilShape {
    std::vector<Vec2> pts;     // closed contour
    bool bluntTE = false;
    Vec2 teU{}, teL{};         // trailing-edge corners (blunt TE only)
};

inline Vec2 lerp(const Vec2& a, const Vec2& b, double t) {
    return {a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t};
}

// Interior points of a rounded trailing edge: a downstream semicircle from teU
// to teL (n-1 points, ordered teU->teL). A rounded (corner-free) TE lets the
// boundary layer wrap it smoothly -- like the leading edge -- instead of two
// sharp base corners where the layers collide and distort.
inline std::vector<Vec2> baseArc(const Vec2& teU, const Vec2& teL, int n) {
    std::vector<Vec2> out;
    const Vec2 C{(teU.x + teL.x) * 0.5, (teU.y + teL.y) * 0.5};
    const Vec2 a{teU.x - C.x, teU.y - C.y};        // |a| = radius
    Vec2 perp{-a.y, a.x};                          // a rotated +90, same length
    if (perp.x < 0.0) { perp.x = -perp.x; perp.y = -perp.y; }  // bulge downstream (+x)
    for (int k = 1; k < n; ++k) {
        const double th = std::numbers::pi * k / n;
        out.push_back({C.x + std::cos(th) * a.x + std::sin(th) * perp.x,
                       C.y + std::cos(th) * a.y + std::sin(th) * perp.y});
    }
    return out;
}

// Cambered NACA 4-digit. camber m [%], position p [tenths], thickness t [%].
// baseCells: number of segments across a blunt TE base (>=1).
inline AirfoilShape naca4(int mI, int pI, int tI, int nPerSide, double teCut = 1.0,
                          int baseCells = 1) {
    const double m = mI / 100.0, p = pI / 10.0, t = tI / 100.0;
    auto yt = [&](double x) {
        return 5.0 * t * (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x +
                          0.2843 * x * x * x - 0.1036 * x * x * x * x);
    };
    auto camber = [&](double x, double& yc, double& dy) {
        if (m <= 0.0 || p <= 0.0) { yc = 0.0; dy = 0.0; return; }
        if (x < p) { yc = m / (p * p) * (2 * p * x - x * x); dy = 2 * m / (p * p) * (p - x); }
        else { yc = m / ((1 - p) * (1 - p)) * ((1 - 2 * p) + 2 * p * x - x * x);
               dy = 2 * m / ((1 - p) * (1 - p)) * (p - x); }
    };
    std::vector<Vec2> up, lo;
    for (int i = 0; i <= nPerSide; ++i) {
        const double xl = teCut * 0.5 * (1.0 - std::cos(std::numbers::pi * i / nPerSide));
        double yc, dy; camber(xl, yc, dy);
        const double th = std::atan(dy), yy = yt(xl);
        up.push_back({xl - yy * std::sin(th), yc + yy * std::cos(th)});
        lo.push_back({xl + yy * std::sin(th), yc - yy * std::cos(th)});
    }
    const bool blunt = teCut < 1.0 - 1e-9;
    AirfoilShape s;
    s.bluntTE = blunt;
    for (int i = 0; i <= nPerSide; ++i) s.pts.push_back(up[i]);            // LE -> TE (upper)
    if (blunt)                                                            // rounded TE cap
        for (const auto& q : baseArc(up[nPerSide], lo[nPerSide], baseCells))
            s.pts.push_back(q);
    for (int i = (blunt ? nPerSide : nPerSide - 1); i >= 1; --i) s.pts.push_back(lo[i]);
    if (blunt) { s.teU = up[nPerSide]; s.teL = lo[nPerSide]; }
    return s;
}

// Read a Selig .dat file: a closed loop, TE -> upper -> LE -> lower -> TE.
inline AirfoilShape loadSeligDat(const std::string& path) {
    AirfoilShape s;
    std::ifstream f(path);
    std::string line;
    std::vector<Vec2> pts;
    while (std::getline(f, line)) {
        std::istringstream is(line);
        double x, y;
        if (is >> x >> y) {
            if (x >= -0.5 && x <= 2.0 && y >= -1.0 && y <= 1.0) pts.push_back({x, y});
        }
    }
    // drop a duplicated closing TE point
    if (pts.size() > 2 && dist(pts.front(), pts.back()) < 1e-9) pts.pop_back();
    s.pts = std::move(pts);
    s.bluntTE = false;
    return s;
}

// Clip a local-coordinate contour at x = teCut -> small blunt trailing edge.
// baseCells: number of segments across the base (>=1).
inline AirfoilShape truncateTE(const AirfoilShape& in, double teCut, int baseCells = 1) {
    AirfoilShape s;
    for (const auto& p : in.pts) if (p.x <= teCut) s.pts.push_back(p);
    if (s.pts.size() < 3) return in;
    s.bluntTE = true;
    s.teU = s.pts.front();   // first kept (upper, x~teCut)
    s.teL = s.pts.back();    // last kept (lower, x~teCut)
    const auto arc = baseArc(s.teU, s.teL, baseCells);   // teU->teL order
    for (int k = (int)arc.size() - 1; k >= 0; --k) s.pts.push_back(arc[k]);  // append teL->teU
    return s;
}

// Scale by chord, rotate by angle of attack (deg), translate leading edge to le.
inline AirfoilShape placeAirfoil(AirfoilShape s, double chord, Vec2 le, double aoaDeg) {
    const double c = std::cos(-aoaDeg * std::numbers::pi / 180.0);
    const double sn = std::sin(-aoaDeg * std::numbers::pi / 180.0);
    auto tf = [&](Vec2 q) {
        const double x = q.x * chord, y = q.y * chord;
        return Vec2{le.x + x * c - y * sn, le.y + x * sn + y * c};
    };
    for (auto& p : s.pts) p = tf(p);
    if (s.bluntTE) { s.teU = tf(s.teU); s.teL = tf(s.teL); }
    return s;
}

} // namespace poly2d
