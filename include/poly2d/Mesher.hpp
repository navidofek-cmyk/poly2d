#pragma once
// -----------------------------------------------------------------------------
// Mesher -- builds a 2D polygonal (polyhedral) mesh with prism boundary layers.
//
// Pipeline:
//   1. Core seeds        : graded Poisson-disk sampling driven by a size field,
//                          restricted to the domain minus the prism bands.
//   2. Prism seeds       : structured rings offset inward from every loop that
//                          requests boundary layers (geometric growth).
//   3. Mirror seeds      : near-boundary seeds are reflected across the wall so
//                          Voronoi edges land exactly on the boundary.
//   4. Delaunay + Voronoi: the polygonal cells are the Voronoi duals of the real
//                          seeds; Voronoi vertices are triangle circumcenters.
//   5. Snap              : boundary Voronoi nodes are projected onto the wall.
//
// The result is a clean, boundary-conforming polygonal mesh whose near-wall
// cells form a structured prism layer.
// -----------------------------------------------------------------------------
#include "Delaunay.hpp"
#include "Domain.hpp"
#include "Vec2.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <numbers>
#include <random>
#include <unordered_map>
#include <vector>

namespace poly2d {

// Final polygonal mesh (2D). Cells reference welded node indices, CCW.
struct PolyMesh2D {
    std::vector<Vec2> nodes;
    std::vector<std::vector<int>> cells;
    std::vector<int> cellType;   // 0 = polyhedral core, 1 = prism layer
};

class Mesher {
public:
    struct Options {
        std::function<double(Vec2)> sizeField;  // target edge length h(x); required
        unsigned seed = 1234567u;                // RNG seed (deterministic)
        int poissonTries = 30;                   // Bridson candidate attempts
        int lloydIters = 0;                      // Lloyd relaxation passes on core
        bool transitionRing = false;             // add a 1:1 prism->poly ring
    };

    Mesher(const Domain& dom, Options opt) : dom_(dom), opt_(std::move(opt)) {}

    PolyMesh2D generate() {
        buildSeeds();
        return buildVoronoi();
    }

private:
    // ------------------------------------------------------------------ seeds
    struct Seed { Vec2 p; int type; };

    const Domain& dom_;
    Options opt_;
    std::vector<Seed> real_;     // core + prism seeds (output cells)
    std::vector<Vec2> ghost_;    // mirror seeds (not output)
    std::mt19937 rng_{};
    Vec2 lo_, hi_;
    double diag_ = 1.0;

    // Segments that actually grow prism, with their band thickness. Core seeds
    // are excluded near these (but not near prism-free walls such as a blunt TE).
    std::vector<std::pair<const Segment*, double>> prismSegs_;

    double sizeAt(const Vec2& p) const { return opt_.sizeField(p); }

    double prismBand(int li) const {
        const auto& pr = dom_.loops[li].prism;
        return (pr.nLayers > 0) ? pr.totalThickness() : 0.0;
    }

    void collectPrismSegments() {
        prismSegs_.clear();
        for (const auto& s : dom_.segments()) {
            const Loop& L = dom_.loops[s.loop];
            if (L.prism.nLayers <= 0) continue;
            if (L.prismSkip && L.prismSkip(s.a, s.b)) continue;
            prismSegs_.emplace_back(&s, prismBand(s.loop));
        }
    }

    bool inPrismBand(const Vec2& p) const {
        // Keep the free (relaxed) core clear of the prism band (and of the
        // aligned transition ring, if enabled, which sits ~0.5 h beyond it).
        const double margin = (opt_.transitionRing ? 1.5 : 0.6) * sizeAt(p);
        for (const auto& [seg, band] : prismSegs_)
            if (Domain::segDistance(p, *seg) < band + margin) return true;
        return false;
    }

    bool validCore(const Vec2& p) const {
        return dom_.inside(p) && !inPrismBand(p);
    }

    void buildSeeds() {
        rng_.seed(opt_.seed);
        dom_.bbox(lo_, hi_);
        diag_ = dist(lo_, hi_);

        collectPrismSegments();
        poissonCore();
        prismRings();
        dedupReal();
        lloydRelax(opt_.lloydIters);
        ghost_ = computeMirrors(real_);
    }

    // Lloyd relaxation: move each *core* seed to its Voronoi-cell centroid a few
    // times. Prism seeds stay fixed (structured layers), so only the polyhedral
    // core is smoothed into regular, rounded cells (ANSYS-like honeycomb).
    void lloydRelax(int iters) {
        for (int it = 0; it < iters; ++it) {
            const int nReal = (int)real_.size();
            std::vector<Vec2> ghosts = computeMirrors(real_);
            std::vector<Vec2> all;
            all.reserve(nReal + ghosts.size());
            for (const auto& s : real_) all.push_back(s.p);
            for (const auto& g : ghosts) all.push_back(g);

            Triangulation T = triangulate(all);
            std::vector<std::vector<int>> inc(all.size());
            for (int ti = 0; ti < (int)T.tris.size(); ++ti) {
                inc[T.tris[ti].a].push_back(ti);
                inc[T.tris[ti].b].push_back(ti);
                inc[T.tris[ti].c].push_back(ti);
            }
            for (int v = 0; v < nReal; ++v) {
                if (real_[v].type != 0) continue;         // move core only
                auto& tris = inc[v];
                if ((int)tris.size() < 3) continue;
                const Vec2 c = real_[v].p;
                std::sort(tris.begin(), tris.end(), [&](int i, int k) {
                    const Vec2 pi = T.tris[i].cc - c, pk = T.tris[k].cc - c;
                    return std::atan2(pi.y, pi.x) < std::atan2(pk.y, pk.x);
                });
                std::vector<Vec2> poly;
                poly.reserve(tris.size());
                for (int ti : tris) poly.push_back(T.tris[ti].cc);
                const Vec2 g = polygonCentroid(poly);
                if (validCore(g)) real_[v].p = g;          // full Lloyd step
            }
        }
    }

    static Vec2 polygonCentroid(const std::vector<Vec2>& p) {
        double A = 0.0;
        Vec2 c{0.0, 0.0};
        const int n = (int)p.size();
        for (int i = 0; i < n; ++i) {
            const Vec2& a = p[i];
            const Vec2& b = p[(i + 1) % n];
            const double cr = cross(a, b);
            A += cr;
            c.x += (a.x + b.x) * cr;
            c.y += (a.y + b.y) * cr;
        }
        if (std::abs(A) < 1e-30) {                        // degenerate -> average
            Vec2 m{0, 0};
            for (const auto& q : p) { m.x += q.x; m.y += q.y; }
            return {m.x / n, m.y / n};
        }
        return {c.x / (3.0 * A), c.y / (3.0 * A)};
    }

    // --- graded Poisson-disk (Bridson) over the core region -----------------
    void poissonCore() {
        // Spatial hash sized to the smallest expected spacing.
        double hmin = 1e300;
        for (double fx = 0.05; fx <= 0.95; fx += 0.1)
            for (double fy = 0.05; fy <= 0.95; fy += 0.1)
                hmin = std::min(hmin, sizeAt({lo_.x + fx * (hi_.x - lo_.x),
                                             lo_.y + fy * (hi_.y - lo_.y)}));
        hmin = std::max(hmin, diag_ * 1e-4);
        const double cell = hmin / std::numbers::sqrt2;
        const int gx = std::max(1, (int)std::ceil((hi_.x - lo_.x) / cell) + 1);
        const int gy = std::max(1, (int)std::ceil((hi_.y - lo_.y) / cell) + 1);

        std::vector<std::vector<int>> grid(gx * gy);
        std::vector<Vec2> pts;
        auto gidx = [&](const Vec2& p) {
            int ix = std::clamp((int)((p.x - lo_.x) / cell), 0, gx - 1);
            int iy = std::clamp((int)((p.y - lo_.y) / cell), 0, gy - 1);
            return std::pair<int, int>{ix, iy};
        };
        auto fits = [&](const Vec2& p, double r) {
            auto [ix, iy] = gidx(p);
            const int rad = (int)std::ceil(r / cell) + 1;
            for (int j = std::max(0, iy - rad); j <= std::min(gy - 1, iy + rad); ++j)
                for (int i = std::max(0, ix - rad); i <= std::min(gx - 1, ix + rad); ++i)
                    for (int id : grid[j * gx + i]) {
                        const double need = std::max(r, sizeAt(pts[id]));
                        if (dist(p, pts[id]) < need) return false;
                    }
            return true;
        };
        auto push = [&](const Vec2& p) {
            auto [ix, iy] = gidx(p);
            grid[iy * gx + ix].push_back((int)pts.size());
            pts.push_back(p);
        };

        // Find a starting point somewhere in the core region.
        std::uniform_real_distribution<double> ux(lo_.x, hi_.x), uy(lo_.y, hi_.y);
        Vec2 start{};
        bool haveStart = false;
        for (int k = 0; k < 20000 && !haveStart; ++k) {
            Vec2 c{ux(rng_), uy(rng_)};
            if (validCore(c)) { start = c; haveStart = true; }
        }
        if (!haveStart) return;
        push(start);
        std::vector<int> active{0};

        std::uniform_real_distribution<double> u01(0.0, 1.0);
        while (!active.empty()) {
            std::uniform_int_distribution<int> pick(0, (int)active.size() - 1);
            const int ai = pick(rng_);
            const Vec2 p = pts[active[ai]];
            const double r = sizeAt(p);
            bool found = false;
            for (int t = 0; t < opt_.poissonTries; ++t) {
                const double ang = 2.0 * std::numbers::pi * u01(rng_);
                const double rr = r * (1.0 + u01(rng_)); // annulus [r, 2r]
                const Vec2 c{p.x + rr * std::cos(ang), p.y + rr * std::sin(ang)};
                if (!validCore(c)) continue;
                if (!fits(c, sizeAt(c))) continue;
                push(c);
                active.push_back((int)pts.size() - 1);
                found = true;
                break;
            }
            if (!found) {
                active[ai] = active.back();
                active.pop_back();
            }
        }
        for (const auto& p : pts) real_.push_back({p, 0});
    }

    // --- structured prism rings around every loop that requests layers ------
    void prismRings() {
        for (int li = 0; li < (int)dom_.loops.size(); ++li) {
            const Loop& L = dom_.loops[li];
            if (L.prism.nLayers <= 0) continue;
            const int n = (int)L.nodes.size();
            for (int i = 0; i < n; ++i) {
                const Vec2 a = L.nodes[i];
                const Vec2 b = L.nodes[(i + 1) % n];
                if (L.prismSkip && L.prismSkip(a, b)) continue;  // e.g. blunt TE base
                const Vec2 dir = b - a;
                const double len = norm(dir);
                if (len <= 0.0) continue;
                const Vec2 inward = normalized(leftNormal(dir));
                const int sub = std::max(1, (int)std::lround(len / L.hBnd));
                // One transition polyhedron per prism column, aligned 1:1 with
                // the layers, sitting just outside the last prism layer. It is
                // fixed (type 2, not relaxed) so the prism->poly interface stays
                // one-prism-to-one-polyhedron and continuous.
                const double toff = L.prism.totalThickness() + 0.5 * L.hBnd;
                for (int k = 0; k < sub; ++k) {
                    const double t = (k + 0.5) / sub;         // cell-centered
                    const Vec2 base = a + t * dir;
                    for (int layer = 0; layer < L.prism.nLayers; ++layer)
                        real_.push_back({base + L.prism.center(layer) * inward, 1});
                    if (opt_.transitionRing)
                        real_.push_back({base + toff * inward, 2}); // transition ring
                }
            }
        }
    }

    // --- remove coincident real seeds (corner overlaps etc.) ----------------
    void dedupReal() {
        double tol = 1e300;
        for (const auto& L : dom_.loops)
            if (L.prism.nLayers > 0) tol = std::min(tol, L.prism.firstHeight);
        if (tol > 1e299) tol = diag_ * 1e-3;
        tol *= 0.3;
        real_ = weld(real_, tol);
    }

    static std::vector<Seed> weld(const std::vector<Seed>& in, double tol) {
        std::vector<Seed> out;
        std::unordered_map<std::int64_t, std::vector<int>> grid;
        const double inv = 1.0 / std::max(tol, 1e-300);
        auto key = [&](const Vec2& p) {
            auto qx = (std::int64_t)std::floor(p.x * inv);
            auto qy = (std::int64_t)std::floor(p.y * inv);
            return (qx << 32) ^ (qy & 0xffffffff);
        };
        for (const auto& s : in) {
            bool dup = false;
            for (int dx = -1; dx <= 1 && !dup; ++dx)
                for (int dy = -1; dy <= 1 && !dup; ++dy) {
                    auto it = grid.find(key({s.p.x + dx * tol, s.p.y + dy * tol}));
                    if (it == grid.end()) continue;
                    for (int id : it->second)
                        if (dist(out[id].p, s.p) < tol) { dup = true; break; }
                }
            if (!dup) {
                grid[key(s.p)].push_back((int)out.size());
                out.push_back(s);
            }
        }
        return out;
    }

    // --- reflect near-boundary seeds so cells conform to the walls ----------
    std::vector<Vec2> computeMirrors(const std::vector<Seed>& seeds) const {
        const auto& segs = dom_.segments();
        std::vector<Vec2> ghosts;
        for (const auto& s : seeds) {
            const double band = 0.9 * sizeAt(s.p);
            std::vector<std::pair<double, int>> near;
            for (int i = 0; i < (int)segs.size(); ++i) {
                const double d = Domain::segDistance(s.p, segs[i]);
                if (d < band) near.emplace_back(d, i);
            }
            if (near.empty()) continue;
            std::sort(near.begin(), near.end());
            for (auto& [d, i] : near)
                ghosts.push_back(Domain::reflectLine(segs[i], s.p));
            // corner: double reflection across the two closest segments
            if (near.size() >= 2) {
                const Segment& s0 = segs[near[0].second];
                const Segment& s1 = segs[near[1].second];
                ghosts.push_back(Domain::reflectLine(s1, Domain::reflectLine(s0, s.p)));
                ghosts.push_back(Domain::reflectLine(s0, Domain::reflectLine(s1, s.p)));
            }
        }
        // weld ghosts against each other, keep those actually outside the domain
        std::vector<Seed> tmp;
        tmp.reserve(ghosts.size());
        for (const auto& g : ghosts)
            if (!dom_.inside(g)) tmp.push_back({g, -1});
        tmp = weld(tmp, diag_ * 1e-4);
        std::vector<Vec2> out;
        out.reserve(tmp.size());
        for (const auto& t : tmp) out.push_back(t.p);
        return out;
    }

    // ------------------------------------------------------------- Voronoi
    PolyMesh2D buildVoronoi() {
        const int nReal = (int)real_.size();
        std::vector<Vec2> all;
        all.reserve(nReal + ghost_.size());
        // tiny deterministic jitter breaks exact cocircularity of the rings
        std::uniform_real_distribution<double> j(-1.0, 1.0);
        const double jitter = diag_ * 1e-6;
        for (const auto& s : real_)
            all.push_back({s.p.x + jitter * j(rng_), s.p.y + jitter * j(rng_)});
        for (const auto& g : ghost_)
            all.push_back({g.x + jitter * j(rng_), g.y + jitter * j(rng_)});

        Triangulation T = triangulate(all);

        // vertex -> incident triangles
        std::vector<std::vector<int>> inc(all.size());
        for (int ti = 0; ti < (int)T.tris.size(); ++ti) {
            const Triangle& t = T.tris[ti];
            inc[t.a].push_back(ti);
            inc[t.b].push_back(ti);
            inc[t.c].push_back(ti);
        }

        PolyMesh2D mesh;
        // Voronoi vertex id == triangle id; remap used ones to compact indices.
        std::vector<int> nodeId(T.tris.size(), -1);
        auto useNode = [&](int tri) {
            if (nodeId[tri] < 0) {
                nodeId[tri] = (int)mesh.nodes.size();
                mesh.nodes.push_back(T.tris[tri].cc);
            }
            return nodeId[tri];
        };

        std::vector<std::vector<int>> rawCells; // triangle-id cells (pre-remap)
        std::vector<int> rawType;
        for (int v = 0; v < nReal; ++v) {
            auto tris = inc[v];
            if ((int)tris.size() < 3) continue;
            const Vec2 c = all[v];
            std::sort(tris.begin(), tris.end(), [&](int i, int k) {
                const Vec2 pi = T.tris[i].cc - c;
                const Vec2 pk = T.tris[k].cc - c;
                return std::atan2(pi.y, pi.x) < std::atan2(pk.y, pk.x);
            });
            rawCells.push_back(tris);
            rawType.push_back(real_[v].type);
        }

        // Boundary edges = cell edges used by only one real cell.
        std::unordered_map<std::int64_t, int> edgeCount;
        auto ekey = [](int u, int v) {
            if (u > v) std::swap(u, v);
            return ((std::int64_t)u << 32) ^ (std::uint32_t)v;
        };
        for (const auto& cell : rawCells) {
            const int m = (int)cell.size();
            for (int k = 0; k < m; ++k)
                edgeCount[ekey(cell[k], cell[(k + 1) % m])]++;
        }
        std::vector<char> boundaryNode(T.tris.size(), 0);
        for (const auto& cell : rawCells) {
            const int m = (int)cell.size();
            for (int k = 0; k < m; ++k) {
                const int u = cell[k], w = cell[(k + 1) % m];
                if (edgeCount[ekey(u, w)] == 1) {
                    boundaryNode[u] = 1;
                    boundaryNode[w] = 1;
                }
            }
        }

        // Build compact cells.
        for (size_t ci = 0; ci < rawCells.size(); ++ci) {
            std::vector<int> cell;
            cell.reserve(rawCells[ci].size());
            for (int tri : rawCells[ci]) cell.push_back(useNode(tri));
            mesh.cells.push_back(std::move(cell));
            mesh.cellType.push_back(rawType[ci]);
        }

        // Snap boundary nodes onto the domain boundary.
        for (int tri = 0; tri < (int)T.tris.size(); ++tri)
            if (nodeId[tri] >= 0 && boundaryNode[tri])
                mesh.nodes[nodeId[tri]] = dom_.snapToBoundary(mesh.nodes[nodeId[tri]]);

        cleanup(mesh);
        return mesh;
    }

    // Weld coincident nodes and drop degenerate edges/cells so the extruded
    // mesh is watertight (no non-manifold edges from snapped-together nodes).
    void cleanup(PolyMesh2D& mesh) const {
        const double tol = std::max(diag_ * 1e-8, 1e-12);
        const double inv = 1.0 / tol;
        std::unordered_map<std::int64_t, int> lut;
        std::vector<Vec2> nodes;
        std::vector<int> remap(mesh.nodes.size());
        auto key = [&](const Vec2& p) {
            auto qx = (std::int64_t)std::llround(p.x * inv);
            auto qy = (std::int64_t)std::llround(p.y * inv);
            return (qx << 32) ^ (qy & 0xffffffff);
        };
        for (int i = 0; i < (int)mesh.nodes.size(); ++i) {
            const auto k = key(mesh.nodes[i]);
            auto it = lut.find(k);
            if (it == lut.end()) {
                lut[k] = (int)nodes.size();
                remap[i] = (int)nodes.size();
                nodes.push_back(mesh.nodes[i]);
            } else {
                remap[i] = it->second;
            }
        }
        std::vector<std::vector<int>> cells;
        std::vector<int> types;
        for (std::size_t ci = 0; ci < mesh.cells.size(); ++ci) {
            std::vector<int> c;
            for (int id : mesh.cells[ci]) {
                const int r = remap[id];
                if (c.empty() || c.back() != r) c.push_back(r);
            }
            while (c.size() > 1 && c.front() == c.back()) c.pop_back();
            if (c.size() < 3) continue;                 // collapsed sliver
            if (std::abs(Domain::signedArea(polyOf(nodes, c))) < tol * tol) continue;
            cells.push_back(std::move(c));
            types.push_back(mesh.cellType[ci]);
        }
        mesh.nodes.swap(nodes);
        mesh.cells.swap(cells);
        mesh.cellType.swap(types);
    }

    static std::vector<Vec2> polyOf(const std::vector<Vec2>& n, const std::vector<int>& c) {
        std::vector<Vec2> p;
        p.reserve(c.size());
        for (int id : c) p.push_back(n[id]);
        return p;
    }
};

} // namespace poly2d
