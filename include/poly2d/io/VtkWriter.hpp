#pragma once
// -----------------------------------------------------------------------------
// VtkWriter -- legacy VTK UnstructuredGrid (VTK_POLYGON) for ParaView.
// -----------------------------------------------------------------------------
#include "../Mesher.hpp"

#include <fstream>
#include <string>

namespace poly2d::io {

inline void writeVtk(const PolyMesh2D& m, const std::string& path) {
    std::ofstream f(path);
    f << "# vtk DataFile Version 3.0\n";
    f << "poly2d polygonal mesh\n";
    f << "ASCII\n";
    f << "DATASET UNSTRUCTURED_GRID\n";

    f << "POINTS " << m.nodes.size() << " double\n";
    for (const auto& p : m.nodes)
        f << p.x << ' ' << p.y << " 0\n";

    std::size_t total = 0;
    for (const auto& c : m.cells) total += c.size() + 1;
    f << "CELLS " << m.cells.size() << ' ' << total << '\n';
    for (const auto& c : m.cells) {
        f << c.size();
        for (int id : c) f << ' ' << id;
        f << '\n';
    }
    f << "CELL_TYPES " << m.cells.size() << '\n';
    for (std::size_t i = 0; i < m.cells.size(); ++i) f << "7\n"; // VTK_POLYGON

    f << "CELL_DATA " << m.cells.size() << '\n';
    f << "SCALARS region int 1\nLOOKUP_TABLE default\n";
    for (int t : m.cellType) f << t << '\n';
    f << "SCALARS nSides int 1\nLOOKUP_TABLE default\n";
    for (const auto& c : m.cells) f << c.size() << '\n';
}

} // namespace poly2d::io
