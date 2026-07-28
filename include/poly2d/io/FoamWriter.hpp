#pragma once
// -----------------------------------------------------------------------------
// FoamWriter -- export a 2D polygonal mesh as an OpenFOAM constant/polyMesh.
//
// OpenFOAM has no native 2D mesh: a 2D case is a single-cell-thick 3D mesh whose
// front and back faces form an "empty" patch. We therefore extrude every 2D
// polygon into a prism (one z-layer):
//
//   * each 2D cell            -> one 3D polyhedral cell,
//   * each interior 2D edge   -> one internal quad face (owner < neighbour),
//   * each boundary 2D edge   -> one boundary quad face (named patch),
//   * each cell's z=0 / z=tz  -> two quad/polygon faces on the empty patch.
//
// Face point ordering is fixed up with a Newell normal so that every face normal
// points from owner to neighbour (internal) or out of the domain (boundary),
// exactly as OpenFOAM's checkMesh expects.
//
// Works with OpenFOAM 14 (classic ascii polyMesh format).
// -----------------------------------------------------------------------------
#include "../Domain.hpp"
#include "../Mesher.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace poly2d::io {

struct FoamOptions {
    double scale = 1.0;        // multiply x,y by this (e.g. 0.01 to go cm -> m)
    double thickness = 0.01;   // extrusion thickness in z, in *scaled* units
    std::string emptyPatch = "frontAndBack";
    // name -> OpenFOAM patch type; default: inlet/outlet/farfield => patch, else wall
    std::function<std::string(const std::string&)> patchType =
        [](const std::string& n) -> std::string {
            if (n == "inlet" || n == "outlet" || n == "farfield") return "patch";
            return "wall";
        };
};

class FoamWriter {
public:
    FoamWriter(const Domain& dom, FoamOptions opt) : dom_(dom), opt_(std::move(opt)) {}

    void write(const PolyMesh2D& m, const std::string& polyMeshDir) {
        build(m);
        writePoints(polyMeshDir + "/points");
        writeFaces(polyMeshDir + "/faces");
        writeLabelList(polyMeshDir + "/owner", owner_, "owner");
        writeLabelList(polyMeshDir + "/neighbour",
                       std::vector<int>(neigh_.begin(), neigh_.begin() + nInternal_),
                       "neighbour");
        writeBoundary(polyMeshDir + "/boundary");
    }

    int nInternalFaces() const { return nInternal_; }
    int nFaces() const { return (int)faces_.size(); }
    int nCells() const { return nCells_; }

private:
    struct Vec3 { double x, y, z; };
    struct Face { std::vector<int> pts; int owner; int neigh; int patch; }; // patch<0 => internal

    const Domain& dom_;
    FoamOptions opt_;

    std::vector<Vec3> P_;          // 3D points (scaled)
    std::vector<Face> faces_;      // final ordered faces
    std::vector<int> owner_, neigh_;
    int nInternal_ = 0;
    int nCells_ = 0;

    std::vector<std::string> patchNames_;
    std::vector<std::string> patchTypes_;
    std::vector<int> patchStart_, patchCount_;

    static std::int64_t ekey(int u, int v) {
        if (u > v) std::swap(u, v);
        return ((std::int64_t)u << 32) ^ (std::uint32_t)v;
    }

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
    void build(const PolyMesh2D& m) {
        const int N = (int)m.nodes.size();
        nCells_ = (int)m.cells.size();
        const double s = opt_.scale, tz = opt_.thickness;

        // 3D points: back copy (z=0) at i, front copy (z=tz) at i+N.
        P_.resize(2 * N);
        for (int i = 0; i < N; ++i) {
            P_[i]     = {m.nodes[i].x * s, m.nodes[i].y * s, 0.0};
            P_[i + N] = {m.nodes[i].x * s, m.nodes[i].y * s, tz};
        }

        // Directed half-edges: each 2D edge is stored with the winding of the
        // cell that owns it. The in-plane outward normal of a directed edge
        // u->v in a CCW cell is (dy,-dx) -- a robust orientation reference that
        // does not depend on cell convexity.
        struct HE { int cell, u, v; };
        std::map<std::int64_t, std::vector<HE>> he;
        for (int c = 0; c < nCells_; ++c) {
            const auto& cell = m.cells[c];
            const int mm = (int)cell.size();
            for (int k = 0; k < mm; ++k) {
                const int u = cell[k], v = cell[(k + 1) % mm];
                he[ekey(u, v)].push_back({c, u, v});
            }
        }

        auto patchId = [&](const std::string& name) {
            for (int i = 0; i < (int)patchNames_.size(); ++i)
                if (patchNames_[i] == name) return i;
            patchNames_.push_back(name);
            patchTypes_.push_back(opt_.patchType(name));
            return (int)patchNames_.size() - 1;
        };

        std::vector<Face> internal, boundary;

        auto orient = [&](std::vector<int>& pts, const Vec3& ref) {
            const Vec3 n = newell(pts);
            if (n.x * ref.x + n.y * ref.y + n.z * ref.z < 0.0)
                std::reverse(pts.begin(), pts.end());
        };
        auto sideQuad = [&](int u, int v) {
            return std::vector<int>{u, v, v + N, u + N};
        };
        auto edgeRef = [&](int u, int v) {   // in-plane outward normal (dy,-dx)
            const double dx = (m.nodes[v].x - m.nodes[u].x) * s;
            const double dy = (m.nodes[v].y - m.nodes[u].y) * s;
            return Vec3{dy, -dx, 0.0};
        };

        for (const auto& [key, list] : he) {
            if (list.size() == 2) {
                // owner = smaller cell index; use its winding for the normal.
                const HE& a = list[0];
                const HE& b = list[1];
                const HE& own = (a.cell < b.cell) ? a : b;
                const int o = own.cell, ne = (a.cell < b.cell) ? b.cell : a.cell;
                std::vector<int> quad = sideQuad(own.u, own.v);
                orient(quad, edgeRef(own.u, own.v));   // outward of owner -> toward neigh
                internal.push_back({quad, o, ne, -1});
            } else if (list.size() == 1) {
                const HE& e = list[0];
                std::vector<int> quad = sideQuad(e.u, e.v);
                orient(quad, edgeRef(e.u, e.v));       // outward of domain
                const Vec2 mid{(m.nodes[e.u].x + m.nodes[e.v].x) * 0.5,
                               (m.nodes[e.u].y + m.nodes[e.v].y) * 0.5};
                boundary.push_back({quad, e.cell, -1, patchId(dom_.patchAt(mid))});
            }
            // list.size()>2 would be non-manifold; skip defensively.
        }

        // Front/back (empty) faces, two per cell: back normal -z, front +z.
        const int emptyId = patchId(opt_.emptyPatch);
        patchTypes_[emptyId] = "empty";
        for (int c = 0; c < nCells_; ++c) {
            const auto& cell = m.cells[c];
            std::vector<int> back, front;
            back.reserve(cell.size());
            front.reserve(cell.size());
            for (int id : cell) back.push_back(id);
            for (int id : cell) front.push_back(id + N);
            orient(back, {0, 0, -1});
            orient(front, {0, 0, 1});
            boundary.push_back({back, c, -1, emptyId});
            boundary.push_back({front, c, -1, emptyId});
        }

        // Order faces: internal (sorted by owner,neigh) then boundary by patch.
        std::sort(internal.begin(), internal.end(), [](const Face& a, const Face& b) {
            return a.owner != b.owner ? a.owner < b.owner : a.neigh < b.neigh;
        });
        nInternal_ = (int)internal.size();

        faces_.clear();
        faces_.insert(faces_.end(), internal.begin(), internal.end());

        patchStart_.assign(patchNames_.size(), 0);
        patchCount_.assign(patchNames_.size(), 0);
        for (int pid = 0; pid < (int)patchNames_.size(); ++pid) {
            patchStart_[pid] = (int)faces_.size();
            for (const auto& f : boundary)
                if (f.patch == pid) { faces_.push_back(f); patchCount_[pid]++; }
        }

        owner_.resize(faces_.size());
        neigh_.resize(faces_.size());
        for (int i = 0; i < (int)faces_.size(); ++i) {
            owner_[i] = faces_[i].owner;
            neigh_[i] = faces_[i].neigh;
        }
    }

    // ---- file writers ------------------------------------------------------
    static void header(std::ofstream& f, const std::string& cls, const std::string& object) {
        f << "FoamFile\n{\n";
        f << "    version     2.0;\n";
        f << "    format      ascii;\n";
        f << "    class       " << cls << ";\n";
        f << "    location    \"constant/polyMesh\";\n";
        f << "    object      " << object << ";\n";
        f << "}\n\n";
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
            for (std::size_t i = 0; i < fc.pts.size(); ++i)
                f << (i ? " " : "") << fc.pts[i];
            f << ")\n";
        }
        f << ")\n";
    }
    void writeLabelList(const std::string& path, const std::vector<int>& v,
                        const std::string& object) {
        std::ofstream f(path);
        header(f, "labelList", object);
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
            if (patchTypes_[i] == "empty")
                f << "        inGroups        1(" << patchNames_[i] << ");\n";
            f << "        nFaces          " << patchCount_[i] << ";\n";
            f << "        startFace       " << patchStart_[i] << ";\n";
            f << "    }\n";
        }
        f << ")\n";
    }
};

} // namespace poly2d::io
