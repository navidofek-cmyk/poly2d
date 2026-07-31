// -----------------------------------------------------------------------------
// poly2d -- general 2D polygonal (polyhedral) mesher with prism layers.
//
// Two demo cases (the mesher core is fully geometry-agnostic):
//   rect     : 50 x 15 cm channel with a centered circular hole (d = 4 cm)
//   airfoil  : NACA 0012 (chord 1 m) inside a circular farfield (d = 10 m),
//              15 prism layers on the airfoil, first cell 0.005 m
//
// Outputs per case in out/<case>/ :
//   mesh.vtk                  polygonal mesh for ParaView
//   mesh.svg                  standalone preview
//   constant/polyMesh/*       OpenFOAM 14 mesh (extruded, empty front/back)
// -----------------------------------------------------------------------------
#include "poly2d/Domain.hpp"
#include "poly2d/Mesher.hpp"
#include "poly2d/io/Extrude3D.hpp"
#include "poly2d/io/Revolve3D.hpp"
#include "poly2d/io/FoamWriter.hpp"
#include "poly2d/io/SvgWriter.hpp"
#include "poly2d/io/VtkWriter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <numbers>
#include <string>
#include <vector>

using namespace poly2d;
namespace fs = std::filesystem;

// NACA 0012 contour (chord along +x), cosine-clustered at LE/TE.
//   teCut = 1.0 : classic sharp (zero-thickness) trailing-edge point
//   teCut < 1.0 : truncate at x=teCut -> small blunt trailing edge of finite
//                 thickness 2*yt(teCut). The short vertical base is kept as a
//                 boundary edge, so the prism layer wraps continuously around
//                 the whole profile *including* that small edge.
static std::vector<Vec2> naca0012(double chord, Vec2 le, int nPerSide,
                                  double aoaDeg = 0.0, double teCut = 1.0) {
    const double t = 0.12;
    auto yt = [&](double x) {
        return 5.0 * t * (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x +
                          0.2843 * x * x * x - 0.1036 * x * x * x * x);
    };
    std::vector<double> xs(nPerSide + 1);
    for (int i = 0; i <= nPerSide; ++i) {
        const double th = std::numbers::pi * i / nPerSide;
        xs[i] = teCut * 0.5 * (1.0 - std::cos(th)); // 0 -> teCut, clustered at ends
    }
    const bool blunt = teCut < 1.0 - 1e-9;
    std::vector<Vec2> pts;
    for (int i = 0; i <= nPerSide; ++i) pts.push_back({xs[i], yt(xs[i])});          // upper LE->TE
    // keep the lower TE point for a blunt TE (adds the short base segment)
    const int loStart = blunt ? nPerSide : nPerSide - 1;
    for (int i = loStart; i >= 1; --i) pts.push_back({xs[i], -yt(xs[i])});          // lower TE->LE

    const double c = std::cos(-aoaDeg * std::numbers::pi / 180.0);
    const double s = std::sin(-aoaDeg * std::numbers::pi / 180.0);
    for (auto& p : pts) {
        double x = p.x * chord, y = p.y * chord;      // scale by chord
        double xr = x * c - y * s, yr = x * s + y * c; // angle of attack
        p = {le.x + xr, le.y + yr};
    }
    return pts;
}

static void runCase(const std::string& name, const Domain& dom,
                    const Mesher::Options& mopt, const io::FoamOptions& fopt,
                    double svgPx) {
    std::printf("=== case '%s' ===\n", name.c_str());
    Mesher mesher(dom, mopt);
    PolyMesh2D mesh = mesher.generate();

    std::size_t prism = 0;
    for (int t : mesh.cellType) prism += (t == 1);
    std::printf("  cells: %zu  (prism %zu / core %zu)   nodes: %zu\n",
                mesh.cells.size(), prism, mesh.cells.size() - prism, mesh.nodes.size());

    const std::string dir = "out/" + name;
    fs::create_directories(dir + "/constant/polyMesh");
    io::writeVtk(mesh, dir + "/mesh.vtk");
    io::writeSvg(mesh, dir + "/mesh.svg", svgPx);

    io::FoamWriter fw(dom, fopt);
    fw.write(mesh, dir + "/constant/polyMesh");
    std::printf("  OpenFOAM: %d faces (%d internal), %d cells -> %s/constant/polyMesh\n",
                fw.nFaces(), fw.nInternalFaces(), fw.nCells(), dir.c_str());
}

// ---- case 1: rectangular channel with a circular hole (centimetres) --------
static void caseRect() {
    Domain dom;
    const double W = 50.0, H = 15.0;
    PrismSpec wallPrism{5, 0.08, 1.3};
    PrismSpec cylPrism{6, 0.06, 1.3};
    dom.addRectangle(0, 0, W, H, "inlet", "outlet", "bottom", "top", 1.0, wallPrism);
    dom.addCircle({W / 2, H / 2}, 2.0, "cylinder", 0.35, /*hole*/ true, cylPrism);
    dom.build();

    Mesher::Options mo;
    const Vec2 cc{W / 2, H / 2};
    mo.sizeField = [cc](Vec2 p) {
        const double d = dist(p, cc) - 2.0;               // distance from cylinder
        return std::clamp(0.45 + 0.12 * d, 0.45, 1.1);    // fine near cylinder
    };
    mo.lloydIters = 4;    // regular, rounded polyhedral cells
    io::FoamOptions fo;
    fo.scale = 0.01;      // cm -> m
    fo.thickness = 0.01;  // 1 cm in z
    runCase("rect", dom, mo, fo, 18.0);
}

// ---- case 2: NACA 0012 in a circular farfield (metres) ---------------------
static void caseAirfoil() {
    Domain dom;
    const double chord = 1.0, R = 5.0;
    dom.addCircle({0, 0}, R, "farfield", 0.35, /*hole*/ false); // outer, no prism
    PrismSpec afPrism{15, 0.0012, 1.06};  // thin boundary layer, total ~0.028 (2.8% c)
    const double aoaDeg = 10.0;
    const int nPerSide = 160;
    // Small blunt trailing edge (truncate at 94% chord -> TE thickness ~1.7% c).
    // The prism layer grows continuously around the whole profile including the
    // small TE base, forming closed boundary-layer lanes that preserve the edge.
    // Sharp trailing edge: the prism generator rounds the convex TE vertex with
    // arc fans, so the boundary-layer lanes wrap smoothly around the trailing
    // edge (concentric lanes) instead of collapsing to a point.
    auto af = naca0012(chord, {-0.5, 0.0}, nPerSide, aoaDeg, /*teCut*/ 1.0);
    // Fine wall resolution; prism tangential spacing matches the near-airfoil
    // core size so the prism->polyhedral interface is continuous (1 prism : 1
    // polyhedron).
    dom.addPolyLoop(af, "airfoil", /*hole*/ true, 0.006, afPrism);
    dom.build();

    const int afLoop = 1; // farfield=0, airfoil=1
    const double band = afPrism.totalThickness();   // prism band thickness
    const double hWall = 0.006;                      // core size at the prism front (1:1)
    Mesher::Options mo;
    // Hold the core size equal to the prism tangential spacing until the prism
    // front, then grade outward. This makes the first polyhedral cells the same
    // size as the prism columns -> one prism connects to one polyhedron.
    mo.sizeField = [&dom, afLoop, band, hWall](Vec2 p) {
        const double d = dom.distanceToLoop(p, afLoop);
        const double dd = std::max(0.0, d - band);   // distance beyond the prism front
        return std::clamp(hWall + 0.13 * dd, hWall, 0.6);
    };
    mo.lloydIters = 5;          // regular, rounded polyhedral cells (ANSYS-like)
    mo.transitionRing = true;   // exact 1 prism : 1 polyhedron at the interface
    io::FoamOptions fo;
    fo.scale = 1.0;       // already metres
    fo.thickness = 0.05;
    runCase("airfoil", dom, mo, fo, 90.0);
}

// NACA 0012 half-thickness at chord fraction x in [0,1] (for a body of revolution).
static double naca0012_yt(double x) {
    const double t = 0.12;
    return 5.0 * t * (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x +
                      0.2843 * x * x * x - 0.1036 * x * x * x * x);
}

// ---- case 3: NACA 0012 body of revolution inside a sphere (axisymmetric) ----
// Meshed in the meridional half-plane (x, r>=0): the body sits on the axis
// (r=0) with radius r=yt(x); the sphere becomes an upper semicircle; the axis
// segments upstream/downstream of the body are the symmetry axis. Revolving
// this 2D mesh about the x-axis gives the 3D mesh (OpenFOAM: wedge patches).
static void caseSphereAxi() {
    const double R = 5.0, chord = 1.0, x0 = -0.5;
    const int nb = 180, na = 160;

    std::vector<Vec2> nodes;
    std::vector<std::string> patch;   // patch of segment nodes[i]->nodes[i+1]
    auto add = [&](Vec2 p, const std::string& seg) { nodes.push_back(p); patch.push_back(seg); };

    // Truncate the body at 95% chord -> small blunt tail (a flat disc on the
    // axis) so the trailing edge is not a mesh singularity, exactly like the
    // small blunt TE used on the 2D airfoil.
    const double teCut = 0.95;
    add({-R, 0.0}, "axis");                                   // (-R,0) -> nose : axis
    std::vector<Vec2> bodyPts;
    for (int i = 0; i <= nb; ++i) {                           // body surface nose -> tail
        const double xl = teCut * 0.5 * (1.0 - std::cos(std::numbers::pi * i / nb));
        const Vec2 p{x0 + chord * xl, chord * naca0012_yt(xl)};
        bodyPts.push_back(p);
        add(p, "body");                    // last body node's seg -> base (also "body")
    }
    const Vec2 teU = bodyPts.back();       // tail top corner (x_TE, r_TE)
    const Vec2 teL{teU.x, 0.0};            // tail on the axis
    add(teL, "axis");                      // base (teU->teL) is "body"; teL->(R,0) axis
    add({R, 0.0}, "farfield");             // (R,0) -> arc : farfield
    for (int j = 1; j < na; ++j) {                            // upper semicircle back to (-R,0)
        const double th = std::numbers::pi * j / na;
        add({R * std::cos(th), R * std::sin(th)}, "farfield");
    }

    Domain dom;
    Loop L;
    L.nodes = nodes;
    L.patch = patch;
    L.hole = false;
    L.hBnd = 0.006;                           // tangential step ~ core size (1:1)
    L.prism = PrismSpec{15, 0.0007, 1.05};    // thin layer, total ~0.015
    // prism only on the body surface (skip the axis, the farfield arc and the
    // small blunt-tail base)
    auto same = [](const Vec2& p, const Vec2& q) { return dist(p, q) < 1e-9; };
    L.prismSkip = [R, teU, teL, same](const Vec2& a, const Vec2& b) {
        const bool base = (same(a, teU) && same(b, teL)) || (same(a, teL) && same(b, teU));
        const Vec2 m{(a.x + b.x) * 0.5, (a.y + b.y) * 0.5};
        const bool body = (m.y > 1e-6) && (norm(m) < R - 1e-3) && !base;
        return !body;
    };
    dom.loops.push_back(std::move(L));
    dom.build();

    // Smooth radial grading from the body (incl. the blunt tail corner): no
    // thin on-axis refinement strip, which would create stretched wake cells.
    const double band = dom.loops[0].prism.totalThickness();
    const double hWall = 0.006;
    const double xTE = teU.x;
    Mesher::Options mo;
    mo.sizeField = [bodyPts, band, hWall, xTE](Vec2 p) {
        double d = 1e9;
        for (const auto& b : bodyPts) d = std::min(d, dist(p, b));
        // mild extra refinement in the near wake behind the tail (a region, not
        // a line) so the transition stays smooth and cells stay well-shaped.
        if (p.x > xTE && p.y < 0.25) d = std::min(d, 0.5 * (p.x - xTE) + p.y);
        return std::clamp(hWall + 0.12 * std::max(0.0, d - band), hWall, 0.5);
    };
    mo.lloydIters = 5;
    mo.transitionRing = true;

    std::printf("=== case 'sphere' (axisymmetric meridional plane) ===\n");
    PolyMesh2D mesh = Mesher(dom, mo).generate();
    std::size_t prism = 0;
    for (int t : mesh.cellType) prism += (t == 1);
    std::printf("  cells: %zu  (prism %zu / core %zu)   nodes: %zu\n",
                mesh.cells.size(), prism, mesh.cells.size() - prism, mesh.nodes.size());
    fs::create_directories("out/sphere");
    io::writeVtk(mesh, "out/sphere/mesh.vtk");
    io::writeSvg(mesh, "out/sphere/mesh.svg", 90.0);
    std::printf("  wrote out/sphere/mesh.{vtk,svg}\n");

    // Revolve the meridional mesh 360 deg -> full 3D mesh of the body of
    // revolution inside the sphere.
    io::Revolve3DOptions ro;
    ro.nSectors = 48;
    io::Revolve3D rev(dom, ro);
    fs::create_directories("out/sphere3d/constant/polyMesh");
    rev.write(mesh, "out/sphere3d/constant/polyMesh");
    std::printf("  3D revolve: %d cells, %d faces (%d internal), %d points -> out/sphere3d\n",
                rev.nCells(), rev.nFaces(), rev.nInternalFaces(), rev.nPoints());
}

// ---- case 4: square block with a square through-hole (duct cross-section) ---
// This 2D cross-section is prismatic along z, so extruding it gives the 3D mesh
// of a cube pierced by a square hole. Tests prism layers + arc fans on the
// right-angle hole corners.
static void caseDuct() {
    Domain dom;
    const double W = 2.0, H = 2.0;
    dom.addRectangle(0, 0, W, H, "left", "right", "bottom", "top", 0.05,
                     PrismSpec{8, 0.004, 1.2});
    // centered square hole
    const double a = 0.6;
    const Vec2 c{W / 2, H / 2};
    std::vector<Vec2> hole = {{c.x - a / 2, c.y - a / 2}, {c.x + a / 2, c.y - a / 2},
                              {c.x + a / 2, c.y + a / 2}, {c.x - a / 2, c.y + a / 2}};
    dom.addPolyLoop(hole, "hole", /*hole*/ true, 0.03, PrismSpec{8, 0.003, 1.2});
    dom.build();

    Mesher::Options mo;
    mo.sizeField = [c, a](Vec2 p) {
        const double d = std::max(0.0, std::max(std::abs(p.x - c.x), std::abs(p.y - c.y)) - a / 2);
        return std::clamp(0.03 + 0.12 * d, 0.03, 0.09);
    };
    mo.lloydIters = 5;
    mo.transitionRing = true;
    io::FoamOptions fo;
    fo.scale = 1.0;
    fo.thickness = 0.05;
    runCase("duct", dom, mo, fo, 260.0);

    // Extrude the cross-section into a genuine 3D volume mesh: a cube pierced by
    // a straight square through-hole (z = flow direction).
    Mesher::Options mo3 = mo;
    PolyMesh2D mesh = Mesher(dom, mo3).generate();
    io::Extrude3DOptions eo;
    eo.nLayers = 24;
    eo.length = 3.0;      // duct length along z
    io::Extrude3D ex(dom, eo);
    fs::create_directories("out/duct3d/constant/polyMesh");
    ex.write(mesh, "out/duct3d/constant/polyMesh");
    std::printf("  3D extrude: %d cells, %d faces (%d internal), %d points -> out/duct3d\n",
                ex.nCells(), ex.nFaces(), ex.nInternalFaces(), ex.nPoints());
}

int main(int argc, char** argv) {
    std::string which = (argc > 1) ? argv[1] : "both";
    if (which == "rect" || which == "both") caseRect();
    if (which == "airfoil" || which == "both") caseAirfoil();
    if (which == "sphere" || which == "both") caseSphereAxi();
    if (which == "duct" || which == "both") caseDuct();
    std::puts("done.");
    return 0;
}
