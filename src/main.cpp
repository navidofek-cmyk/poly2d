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

// NACA 0012 closed contour (chord along +x), cosine-clustered at LE/TE.
//   sharpTE = true : classic closed (zero-thickness) trailing edge point
//   sharpTE = false: open coefficient -> thin blunt TE (base ~0.0025 c) that
//                    the prism layers can land on instead of collapsing to a
//                    single point.
static std::vector<Vec2> naca0012(double chord, Vec2 le, int nPerSide,
                                  double aoaDeg = 0.0, bool sharpTE = true) {
    const double t = 0.12;
    const double c4 = sharpTE ? -0.1036 : -0.1015; // closed vs finite TE
    auto yt = [&](double x) {
        return 5.0 * t * (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x +
                          0.2843 * x * x * x + c4 * x * x * x * x);
    };
    std::vector<double> xs(nPerSide + 1);
    for (int i = 0; i <= nPerSide; ++i) {
        const double th = std::numbers::pi * i / nPerSide;
        xs[i] = 0.5 * (1.0 - std::cos(th)); // 0 -> 1, clustered at ends
    }
    std::vector<Vec2> pts;
    for (int i = 0; i <= nPerSide; ++i) pts.push_back({xs[i], yt(xs[i])});          // upper LE->TE
    // include the lower TE point for a blunt TE (adds the short base segment)
    const int loStart = sharpTE ? nPerSide - 1 : nPerSide;
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
    PrismSpec afPrism{15, 0.005, 1.05};  // thin boundary layer (total ~0.11 m)
    const double aoaDeg = 10.0;
    const int nPerSide = 160;
    // Sharp trailing edge -> the prism layer wraps continuously around the whole
    // airfoil (as in the reference), with no blunt base and no gap.
    auto af = naca0012(chord, {-0.5, 0.0}, nPerSide, aoaDeg, /*sharpTE*/ true);
    // Prism tangential spacing matches the near-airfoil core size so the
    // prism->polyhedral interface is continuous (no scalloping / kinks).
    dom.addPolyLoop(af, "airfoil", /*hole*/ true, 0.015, afPrism);
    dom.build();

    const int afLoop = 1; // farfield=0, airfoil=1
    const double band = afPrism.totalThickness();   // prism band thickness
    const double hWall = 0.015;                      // core size at the prism front
    Mesher::Options mo;
    // Hold the core size equal to the prism tangential spacing until the prism
    // front, then grade outward. This makes the first polyhedral cells the same
    // size as the prism columns -> one prism connects to one polyhedron.
    mo.sizeField = [&dom, afLoop, band, hWall](Vec2 p) {
        const double d = dom.distanceToLoop(p, afLoop);
        const double dd = std::max(0.0, d - band);   // distance beyond the prism front
        return std::clamp(hWall + 0.13 * dd, hWall, 0.6);
    };
    mo.lloydIters = 5;    // regular, rounded polyhedral cells (ANSYS-like)
    io::FoamOptions fo;
    fo.scale = 1.0;       // already metres
    fo.thickness = 0.05;
    runCase("airfoil", dom, mo, fo, 90.0);
}

int main(int argc, char** argv) {
    std::string which = (argc > 1) ? argv[1] : "both";
    if (which == "rect" || which == "both") caseRect();
    if (which == "airfoil" || which == "both") caseAirfoil();
    std::puts("done.");
    return 0;
}
