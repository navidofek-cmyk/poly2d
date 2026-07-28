#pragma once
// -----------------------------------------------------------------------------
// SvgWriter -- quick standalone preview (no ParaView needed).
//   prism cells: light orange   core cells: light blue
// -----------------------------------------------------------------------------
#include "../Mesher.hpp"

#include <algorithm>
#include <fstream>
#include <string>

namespace poly2d::io {

inline void writeSvg(const PolyMesh2D& m, const std::string& path,
                     double pxPerUnit = 0.0, double marginPx = 20.0) {
    Vec2 lo{1e300, 1e300}, hi{-1e300, -1e300};
    for (const auto& p : m.nodes) {
        lo.x = std::min(lo.x, p.x); lo.y = std::min(lo.y, p.y);
        hi.x = std::max(hi.x, p.x); hi.y = std::max(hi.y, p.y);
    }
    const double w = hi.x - lo.x, h = hi.y - lo.y;
    if (pxPerUnit <= 0.0) pxPerUnit = 1000.0 / std::max(w, h);
    const double W = w * pxPerUnit + 2 * marginPx;
    const double H = h * pxPerUnit + 2 * marginPx;

    auto X = [&](double x) { return marginPx + (x - lo.x) * pxPerUnit; };
    auto Y = [&](double y) { return H - marginPx - (y - lo.y) * pxPerUnit; }; // flip

    std::ofstream f(path);
    f << "<svg xmlns='http://www.w3.org/2000/svg' width='" << (int)W
      << "' height='" << (int)H << "' viewBox='0 0 " << W << ' ' << H << "'>\n";
    f << "<rect width='100%' height='100%' fill='white'/>\n";
    const double sw = std::max(0.15, 0.6 * pxPerUnit / 50.0);
    for (std::size_t i = 0; i < m.cells.size(); ++i) {
        const auto& c = m.cells[i];
        const char* fill = m.cellType[i] == 1 ? "#ffd8a8" : "#d0ebff";
        f << "<polygon points='";
        for (int id : c) f << X(m.nodes[id].x) << ',' << Y(m.nodes[id].y) << ' ';
        f << "' fill='" << fill << "' stroke='#333' stroke-width='" << sw << "'/>\n";
    }
    f << "</svg>\n";
}

} // namespace poly2d::io
