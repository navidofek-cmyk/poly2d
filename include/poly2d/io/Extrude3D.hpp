#pragma once
// -----------------------------------------------------------------------------
// Extrude3D -- extrude a 2D polygonal mesh into a genuine 3D volume mesh
// (N cells along z) and write it as an OpenFOAM constant/polyMesh.
//
// Every 2D polygonal cell becomes a stack of N polyhedral prisms:
//   * each 2D edge  -> N vertical quad faces (internal or a named wall patch),
//   * between z-levels, each 2D cell polygon -> a horizontal internal face,
//   * the z=0 / z=Lz caps -> two end patches (e.g. inlet / outlet).
//
// Suitable for prismatic 3D geometries (e.g. a block pierced by a straight
// through-hole): mesh the cross-section in 2D, extrude here.
// -----------------------------------------------------------------------------
#include "../Domain.hpp"
#include "../Mesher.hpp"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace poly2d::io {

struct Extrude3DOptions {
    int nLayers = 20;          // cells along z
    double length = 1.0;       // total z length
    double scale = 1.0;        // xy scale
    std::string zmin = "inlet";
    std::string zmax = "outlet";
    std::function<std::string(const std::string&)> patchType =
        [](const std::string& n) -> std::string {
            if (n == "inlet" || n == "outlet" || n == "left" || n == "right") return "patch";
            return "wall";
        };
};

class Extrude3D {
public:
    Extrude3D(const Domain& dom, Extrude3DOptions opt) : dom_(dom), opt_(std::move(opt)) {}

    void write(const PolyMesh2D& m, const std::string& polyMeshDir) {
        build(m);
        writePoints(polyMeshDir + "/points");
        writeFaces(polyMeshDir + "/faces");
        writeLabelList(polyMeshDir + "/owner", owner_, "owner");
        writeLabelList(polyMeshDir + "/neighbour",
                       std::vector<int>(neigh_.begin(), neigh_.begin() + nInternal_), "neighbour");
        writeBoundary(polyMeshDir + "/boundary");
    }

    int nFaces() const { return (int)faces_.size(); }
    int nInternalFaces() const { return nInternal_; }
    int nCells() const { return nCells_; }
    int nPoints() const { return (int)P_.size(); }

private:
    struct Vec3 { double x, y, z; };
    struct Face { std::vector<int> pts; int owner; int neigh; int patch; };

    const Domain& dom_;
    Extrude3DOptions opt_;
    std::vector<Vec3> P_;
    std::vector<Face> faces_;
    std::vector<int> owner_, neigh_;
    int nInternal_ = 0, nCells_ = 0, nC2_ = 0, nN2_ = 0, N_ = 0;
    std::vector<std::string> patchNames_, patchTypes_;
    std::vector<int> patchStart_, patchCount_;

    int pid(const std::string& name) {
        for (int i = 0; i < (int)patchNames_.size(); ++i)
            if (patchNames_[i] == name) return i;
        patchNames_.push_back(name);
        patchTypes_.push_back(opt_.patchType(name));
        return (int)patchNames_.size() - 1;
    }
    int cell(int c2, int k) const { return k * nC2_ + c2; }
    int node(int i, int level) const { return level * nN2_ + i; }

    Vec3 newell(const std::vector<int>& f) const {
        Vec3 n{0, 0, 0};
        const int m = (int)f.size();
        for (int i = 0; i < m; ++i) {
            const Vec3& a = P_[f[i]];
            const Vec3& b = P_[f[(i + 1) % m]];
            n.x += (a.y - b.y) * (a.z + b.z);
            n.y += (a.z - b.z) * (a.x + b.x);
            n.z += (a.x - b.x) * (a.y + b.y);
        }
        return n;
    }
    void orient(std::vector<int>& p, const Vec3& ref) {
        const Vec3 n = newell(p);
        if (n.x * ref.x + n.y * ref.y + n.z * ref.z < 0.0)
            std::reverse(p.begin(), p.end());
    }

    void build(const PolyMesh2D& m) {
        nN2_ = (int)m.nodes.size();
        nC2_ = (int)m.cells.size();
        N_ = std::max(1, opt_.nLayers);
        nCells_ = nC2_ * N_;
        const double s = opt_.scale, dz = opt_.length / N_;

        P_.resize((N_ + 1) * nN2_);
        for (int L = 0; L <= N_; ++L)
            for (int i = 0; i < nN2_; ++i)
                P_[node(i, L)] = {m.nodes[i].x * s, m.nodes[i].y * s, L * dz};

        std::vector<Face> internal, boundary;

        // --- horizontal faces (cell polygons at each z-level) ---------------
        for (int c2 = 0; c2 < nC2_; ++c2) {
            const auto& poly = m.cells[c2];
            for (int L = 0; L <= N_; ++L) {
                std::vector<int> f;
                f.reserve(poly.size());
                for (int id : poly) f.push_back(node(id, L));
                if (L == 0) {
                    orient(f, {0, 0, -1});
                    boundary.push_back({f, cell(c2, 0), -1, pid(opt_.zmin)});
                } else if (L == N_) {
                    orient(f, {0, 0, 1});
                    boundary.push_back({f, cell(c2, N_ - 1), -1, pid(opt_.zmax)});
                } else {
                    // between layer L-1 (below) and L (above): owner below, +z
                    orient(f, {0, 0, 1});
                    internal.push_back({f, cell(c2, L - 1), cell(c2, L), -1});
                }
            }
        }

        // --- vertical side faces from 2D half-edges -------------------------
        struct HE { int cell, u, v; };
        std::map<std::int64_t, std::vector<HE>> he;
        auto ekey = [](int u, int v) {
            if (u > v) std::swap(u, v);
            return ((std::int64_t)u << 32) ^ (std::uint32_t)v;
        };
        for (int c2 = 0; c2 < nC2_; ++c2) {
            const auto& poly = m.cells[c2];
            const int mm = (int)poly.size();
            for (int k = 0; k < mm; ++k)
                he[ekey(poly[k], poly[(k + 1) % mm])].push_back({c2, poly[k], poly[(k + 1) % mm]});
        }
        auto edgeRef = [&](int u, int v) {
            const double dx = (m.nodes[v].x - m.nodes[u].x) * s;
            const double dy = (m.nodes[v].y - m.nodes[u].y) * s;
            return Vec3{dy, -dx, 0.0};
        };
        for (const auto& [key, list] : he) {
            for (int k = 0; k < N_; ++k) {
                if (list.size() == 2) {
                    const HE& a = list[0];
                    const HE& b = list[1];
                    const HE& own = (a.cell < b.cell) ? a : b;
                    const int oc = own.cell, nc = (a.cell < b.cell) ? b.cell : a.cell;
                    std::vector<int> q = {node(own.u, k), node(own.v, k),
                                          node(own.v, k + 1), node(own.u, k + 1)};
                    orient(q, edgeRef(own.u, own.v));
                    internal.push_back({q, cell(oc, k), cell(nc, k), -1});
                } else if (list.size() == 1) {
                    const HE& e = list[0];
                    std::vector<int> q = {node(e.u, k), node(e.v, k),
                                          node(e.v, k + 1), node(e.u, k + 1)};
                    orient(q, edgeRef(e.u, e.v));
                    const Vec2 mid{(m.nodes[e.u].x + m.nodes[e.v].x) * 0.5,
                                   (m.nodes[e.u].y + m.nodes[e.v].y) * 0.5};
                    boundary.push_back({q, cell(e.cell, k), -1, pid(dom_.patchAt(mid))});
                }
            }
        }

        std::sort(internal.begin(), internal.end(), [](const Face& a, const Face& b) {
            return a.owner != b.owner ? a.owner < b.owner : a.neigh < b.neigh;
        });
        nInternal_ = (int)internal.size();
        faces_ = internal;
        patchStart_.assign(patchNames_.size(), 0);
        patchCount_.assign(patchNames_.size(), 0);
        for (int p = 0; p < (int)patchNames_.size(); ++p) {
            patchStart_[p] = (int)faces_.size();
            for (const auto& f : boundary)
                if (f.patch == p) { faces_.push_back(f); patchCount_[p]++; }
        }
        owner_.resize(faces_.size());
        neigh_.resize(faces_.size());
        for (int i = 0; i < (int)faces_.size(); ++i) {
            owner_[i] = faces_[i].owner;
            neigh_[i] = faces_[i].neigh;
        }
    }

    static void header(std::ofstream& f, const std::string& cls, const std::string& obj) {
        f << "FoamFile\n{\n    version 2.0;\n    format ascii;\n    class " << cls
          << ";\n    location \"constant/polyMesh\";\n    object " << obj << ";\n}\n\n";
    }
    void writePoints(const std::string& path) {
        std::ofstream f(path);
        header(f, "vectorField", "points");
        f << P_.size() << "\n(\n";
        for (const auto& p : P_) f << "(" << p.x << ' ' << p.y << ' ' << p.z << ")\n";
        f << ")\n";
    }
    void writeFaces(const std::string& path) {
        std::ofstream f(path);
        header(f, "faceList", "faces");
        f << faces_.size() << "\n(\n";
        for (const auto& fc : faces_) {
            f << fc.pts.size() << "(";
            for (std::size_t i = 0; i < fc.pts.size(); ++i) f << (i ? " " : "") << fc.pts[i];
            f << ")\n";
        }
        f << ")\n";
    }
    void writeLabelList(const std::string& path, const std::vector<int>& v, const std::string& obj) {
        std::ofstream f(path);
        header(f, "labelList", obj);
        f << v.size() << "\n(\n";
        for (int x : v) f << x << '\n';
        f << ")\n";
    }
    void writeBoundary(const std::string& path) {
        std::ofstream f(path);
        header(f, "polyBoundaryMesh", "boundary");
        f << patchNames_.size() << "\n(\n";
        for (int i = 0; i < (int)patchNames_.size(); ++i) {
            f << "    " << patchNames_[i] << "\n    {\n";
            f << "        type            " << patchTypes_[i] << ";\n";
            f << "        nFaces          " << patchCount_[i] << ";\n";
            f << "        startFace       " << patchStart_[i] << ";\n    }\n";
        }
        f << ")\n";
    }
};

} // namespace poly2d::io
