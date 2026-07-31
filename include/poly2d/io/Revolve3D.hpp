#pragma once
// -----------------------------------------------------------------------------
// Revolve3D -- revolve a 2D meridional (x, r>=0) mesh 360 deg about the x-axis
// into a genuine 3D volume polyMesh (M azimuthal sectors).
//
// This turns the axisymmetric body-of-revolution-in-a-sphere meridional mesh
// into the full 3D mesh. Handling:
//   * axis nodes (r=0) are welded to a single point (all sectors coincide),
//   * the revolution wraps (sector M-1 neighbours sector 0), so the meridional
//     faces are all internal -- no end caps,
//   * faces that collapse on the axis (zero area) are dropped,
//   * face orientation is fixed by owner->neighbour cell centroids.
//
// Note: cells near the axis become azimuthal slivers (inherent to a full
// revolution); a small-angle wedge avoids that if only a 2D-equivalent run is
// wanted.
// -----------------------------------------------------------------------------
#include "../Domain.hpp"
#include "../Mesher.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
#include <numbers>
#include <string>
#include <vector>

namespace poly2d::io {

struct Revolve3DOptions {
    int nSectors = 48;
    double scale = 1.0;
    double axisEps = 1e-7;   // r below this is treated as on-axis
};

class Revolve3D {
public:
    Revolve3D(const Domain& dom, Revolve3DOptions opt) : dom_(dom), opt_(std::move(opt)) {}

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
    Revolve3DOptions opt_;
    std::vector<Vec3> P_;
    std::vector<Face> faces_;
    std::vector<int> owner_, neigh_;
    std::vector<Vec3> cc_;              // 3D cell centroids
    int nInternal_ = 0, nCells_ = 0, nC2_ = 0, M_ = 0;
    std::vector<int> nonAxisLocal_, axisLocal_;   // per 2D node -> compact index (or -1)
    int cntNonAxis_ = 0, axisBase_ = 0;
    std::vector<std::string> patchNames_, patchTypes_;
    std::vector<int> patchStart_, patchCount_;

    int pid(const std::string& name) {
        for (int i = 0; i < (int)patchNames_.size(); ++i)
            if (patchNames_[i] == name) return i;
        patchNames_.push_back(name);
        patchTypes_.push_back(name == "farfield" ? "patch" : "wall");
        return (int)patchNames_.size() - 1;
    }
    int cell(int c2, int j) const { return j * nC2_ + c2; }
    int pt(int i, int j) const {   // 3D point index of 2D node i at sector j
        if (axisLocal_[i] >= 0) return axisBase_ + axisLocal_[i];
        return (j % M_) * cntNonAxis_ + nonAxisLocal_[i];
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
    Vec3 fcentroid(const std::vector<int>& f) const {
        Vec3 c{0, 0, 0};
        for (int id : f) { c.x += P_[id].x; c.y += P_[id].y; c.z += P_[id].z; }
        const double s = 1.0 / f.size();
        return {c.x * s, c.y * s, c.z * s};
    }
    static std::vector<int> dedup(std::vector<int> f) {
        std::vector<int> o;
        for (int id : f) if (o.empty() || o.back() != id) o.push_back(id);
        while (o.size() > 1 && o.front() == o.back()) o.pop_back();
        return o;
    }
    void orient(std::vector<int>& p, const Vec3& ref) {
        const Vec3 n = newell(p);
        if (n.x * ref.x + n.y * ref.y + n.z * ref.z < 0.0)
            std::reverse(p.begin(), p.end());
    }

    void build(const PolyMesh2D& m) {
        const int nN2 = (int)m.nodes.size();
        nC2_ = (int)m.cells.size();
        M_ = std::max(3, opt_.nSectors);
        nCells_ = nC2_ * M_;
        const double s = opt_.scale;

        // classify nodes on/off axis
        nonAxisLocal_.assign(nN2, -1);
        axisLocal_.assign(nN2, -1);
        for (int i = 0; i < nN2; ++i) {
            if (m.nodes[i].y <= opt_.axisEps) axisLocal_[i] = 0;   // mark; index later
            else nonAxisLocal_[i] = cntNonAxis_++;
        }
        int nAxis = 0;
        for (int i = 0; i < nN2; ++i) if (axisLocal_[i] == 0) axisLocal_[i] = nAxis++;
        axisBase_ = M_ * cntNonAxis_;

        // points
        P_.assign(axisBase_ + nAxis, {0, 0, 0});
        for (int j = 0; j < M_; ++j) {
            const double th = 2.0 * std::numbers::pi * j / M_;
            const double ct = std::cos(th), st = std::sin(th);
            for (int i = 0; i < nN2; ++i) {
                if (nonAxisLocal_[i] < 0) continue;
                const double x = m.nodes[i].x * s, r = m.nodes[i].y * s;
                P_[j * cntNonAxis_ + nonAxisLocal_[i]] = {x, r * ct, r * st};
            }
        }
        for (int i = 0; i < nN2; ++i)
            if (axisLocal_[i] >= 0)
                P_[axisBase_ + axisLocal_[i]] = {m.nodes[i].x * s, 0.0, 0.0};

        // 3D cell centroids
        cc_.assign(nCells_, {0, 0, 0});
        for (int c2 = 0; c2 < nC2_; ++c2) {
            const auto& poly = m.cells[c2];
            for (int j = 0; j < M_; ++j) {
                Vec3 g{0, 0, 0}; int cnt = 0;
                for (int id : poly)
                    for (int jj : {j, j + 1}) {
                        const Vec3& p = P_[pt(id, jj)];
                        g.x += p.x; g.y += p.y; g.z += p.z; ++cnt;
                    }
                cc_[cell(c2, j)] = {g.x / cnt, g.y / cnt, g.z / cnt};
            }
        }

        std::vector<Face> internal, boundary;

        // A) meridional faces: cell polygon at angle theta_j, between sector j-1 and j
        for (int c2 = 0; c2 < nC2_; ++c2) {
            const auto& poly = m.cells[c2];
            for (int j = 0; j < M_; ++j) {
                std::vector<int> f;
                f.reserve(poly.size());
                for (int id : poly) f.push_back(pt(id, j));
                f = dedup(f);
                if ((int)f.size() < 3) continue;
                const int o = cell(c2, (j - 1 + M_) % M_), ne = cell(c2, j);
                const int lo = std::min(o, ne), hi = std::max(o, ne);
                orient(f, {cc_[hi].x - cc_[lo].x, cc_[hi].y - cc_[lo].y, cc_[hi].z - cc_[lo].z});
                internal.push_back({f, lo, hi, -1});
            }
        }

        // B) revolved-edge faces (side walls / internal radial faces)
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
        for (const auto& [key, list] : he) {
            for (int j = 0; j < M_; ++j) {
                const HE& e0 = list[0];
                std::vector<int> q = dedup({pt(e0.u, j), pt(e0.v, j),
                                            pt(e0.v, j + 1), pt(e0.u, j + 1)});
                if ((int)q.size() < 3) continue;            // collapses on the axis
                if (list.size() == 2) {
                    const int o = cell(list[0].cell, j), ne = cell(list[1].cell, j);
                    const int lo = std::min(o, ne), hi = std::max(o, ne);
                    orient(q, {cc_[hi].x - cc_[lo].x, cc_[hi].y - cc_[lo].y, cc_[hi].z - cc_[lo].z});
                    internal.push_back({q, lo, hi, -1});
                } else if (list.size() == 1) {
                    const int oc = cell(e0.cell, j);
                    const Vec3 fc = fcentroid(q);
                    orient(q, {fc.x - cc_[oc].x, fc.y - cc_[oc].y, fc.z - cc_[oc].z});
                    const Vec2 mid{(m.nodes[e0.u].x + m.nodes[e0.v].x) * 0.5,
                                   (m.nodes[e0.u].y + m.nodes[e0.v].y) * 0.5};
                    boundary.push_back({q, oc, -1, pid(dom_.patchAt(mid))});
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
