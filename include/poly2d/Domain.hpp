#pragma once
// -----------------------------------------------------------------------------
// Domain -- a general 2D region described by oriented boundary loops.
//
//   * The outer loop must be oriented counter-clockwise (CCW).
//   * Every hole loop must be oriented clockwise (CW).
//
// With that convention the domain interior always lies to the LEFT of every
// directed boundary segment, so the inward normal is a uniform leftNormal(dir).
// Orientation is enforced automatically in build().
//
// Each loop carries:
//   * a patch name per segment  -> drives OpenFOAM boundary patches,
//   * a PrismSpec               -> number / height / growth of prism layers.
// This keeps the mesher itself completely geometry-agnostic.
// -----------------------------------------------------------------------------
#include "Vec2.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

namespace poly2d {

// Boundary-layer (prism) specification for one loop.
struct PrismSpec {
    int nLayers = 0;
    double firstHeight = 0.0;   // thickness of the first (wall) layer
    double growth = 1.2;        // geometric growth ratio between layers

    // Distance from the wall to the outer edge of layer i (i in [0..nLayers]).
    double edge(int i) const {
        double b = 0.0, h = firstHeight;
        for (int k = 0; k < i; ++k) { b += h; h *= growth; }
        return b;
    }
    // Center offset of layer i (where the prism seed row is placed).
    double center(int i) const { return 0.5 * (edge(i) + edge(i + 1)); }
    double totalThickness() const { return edge(nLayers); }
};

struct Loop {
    std::vector<Vec2> nodes;          // polygon vertices (last connects to first)
    std::vector<std::string> patch;   // patch[i] labels segment nodes[i]->nodes[i+1]
    bool hole = false;                // true for interior holes
    double hBnd = 1.0;                // target boundary sampling spacing
    PrismSpec prism;                  // boundary-layer settings for this loop
};

// A single oriented boundary segment with a precomputed inward normal.
struct Segment {
    Vec2 a, b;
    Vec2 inward;      // unit inward normal (into the domain)
    std::string patch;
    int loop = -1;    // owning loop index
};

class Domain {
public:
    std::vector<Loop> loops;

    // ---- Convenience builders ----------------------------------------------

    void addRectangle(double x0, double y0, double x1, double y1,
                      const std::string& left, const std::string& right,
                      const std::string& bottom, const std::string& top,
                      double hBnd, PrismSpec prism = {}) {
        Loop L;
        L.hole = false;
        L.hBnd = hBnd;
        L.prism = prism;
        L.nodes = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
        L.patch = {bottom, right, top, left};
        loops.push_back(std::move(L));
    }

    void addCircle(Vec2 center, double radius, const std::string& patch,
                   double hBnd, bool hole, PrismSpec prism = {}) {
        Loop L;
        L.hole = hole;
        L.hBnd = hBnd;
        L.prism = prism;
        const int n = std::max(16, (int)std::lround(2.0 * std::numbers::pi * radius / hBnd));
        L.nodes.reserve(n);
        L.patch.assign(n, patch);
        for (int i = 0; i < n; ++i) {
            const double t = 2.0 * std::numbers::pi * i / n;
            L.nodes.push_back({center.x + radius * std::cos(t),
                               center.y + radius * std::sin(t)});
        }
        loops.push_back(std::move(L));
    }

    // Add an arbitrary closed polyline as a loop (patch name shared by all segs).
    void addPolyLoop(std::vector<Vec2> pts, const std::string& patch, bool hole,
                     double hBnd, PrismSpec prism = {}) {
        Loop L;
        L.hole = hole;
        L.hBnd = hBnd;
        L.prism = prism;
        L.nodes = std::move(pts);
        L.patch.assign(L.nodes.size(), patch);
        loops.push_back(std::move(L));
    }

    // ---- Finalise ----------------------------------------------------------
    void build() {
        segments_.clear();
        loopRange_.assign(loops.size(), {0, 0});
        for (int li = 0; li < (int)loops.size(); ++li) {
            Loop& L = loops[li];
            const double area = signedArea(L.nodes);
            const bool wantCCW = !L.hole;          // outer CCW, holes CW
            if ((area > 0.0) != wantCCW) {
                std::reverse(L.nodes.begin(), L.nodes.end());
                std::reverse(L.patch.begin(), L.patch.end());
                std::rotate(L.patch.begin(), L.patch.begin() + (L.patch.size() - 1),
                            L.patch.end());
            }
            const int start = (int)segments_.size();
            const int n = (int)L.nodes.size();
            for (int i = 0; i < n; ++i) {
                Segment s;
                s.a = L.nodes[i];
                s.b = L.nodes[(i + 1) % n];
                s.inward = normalized(leftNormal(s.b - s.a));
                s.patch = L.patch[i];
                s.loop = li;
                segments_.push_back(s);
            }
            loopRange_[li] = {start, (int)segments_.size()};
        }
    }

    const std::vector<Segment>& segments() const { return segments_; }

    // ---- Queries -----------------------------------------------------------

    int nearestSegment(const Vec2& p, double* outDist = nullptr) const {
        int best = -1;
        double bestD = 1e300;
        for (int i = 0; i < (int)segments_.size(); ++i) {
            const double d = segDistance(p, segments_[i]);
            if (d < bestD) { bestD = d; best = i; }
        }
        if (outDist) *outDist = bestD;
        return best;
    }

    double distanceToLoop(const Vec2& p, int li) const {
        double best = 1e300;
        for (int i = loopRange_[li].first; i < loopRange_[li].second; ++i)
            best = std::min(best, segDistance(p, segments_[i]));
        return best;
    }

    std::string patchAt(const Vec2& p) const {
        const int s = nearestSegment(p);
        return (s >= 0) ? segments_[s].patch : std::string("boundary");
    }

    Vec2 snapToBoundary(const Vec2& p) const {
        const int s = nearestSegment(p);
        return (s >= 0) ? projectOnSeg(p, segments_[s]) : p;
    }

    static Vec2 reflectLine(const Segment& s, const Vec2& p) {
        const Vec2 n = s.inward; // unit
        const double d = dot(p - s.a, n);
        return p - 2.0 * d * n;
    }

    // Point-inside test (even-odd rule) over all loops.
    bool inside(const Vec2& p) const {
        bool in = false;
        for (const auto& L : loops) {
            const int n = (int)L.nodes.size();
            for (int i = 0, j = n - 1; i < n; j = i++) {
                const Vec2& a = L.nodes[i];
                const Vec2& b = L.nodes[j];
                if (((a.y > p.y) != (b.y > p.y)) &&
                    (p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x))
                    in = !in;
            }
        }
        return in;
    }

    void bbox(Vec2& lo, Vec2& hi) const {
        lo = {1e300, 1e300};
        hi = {-1e300, -1e300};
        for (const auto& L : loops)
            for (const auto& p : L.nodes) {
                lo.x = std::min(lo.x, p.x); lo.y = std::min(lo.y, p.y);
                hi.x = std::max(hi.x, p.x); hi.y = std::max(hi.y, p.y);
            }
    }

    // ---- Static geometry helpers -------------------------------------------

    static double segDistance(const Vec2& p, const Segment& s) {
        return dist(p, projectOnSeg(p, s));
    }

    static Vec2 projectOnSeg(const Vec2& p, const Segment& s) {
        const Vec2 ab = s.b - s.a;
        const double L2 = norm2(ab);
        double t = (L2 > 0.0) ? dot(p - s.a, ab) / L2 : 0.0;
        t = std::clamp(t, 0.0, 1.0);
        return s.a + t * ab;
    }

    static double signedArea(const std::vector<Vec2>& poly) {
        double a = 0.0;
        const int n = (int)poly.size();
        for (int i = 0; i < n; ++i)
            a += cross(poly[i], poly[(i + 1) % n]);
        return 0.5 * a;
    }

private:
    std::vector<Segment> segments_;
    std::vector<std::pair<int, int>> loopRange_; // [begin,end) into segments_
};

} // namespace poly2d
