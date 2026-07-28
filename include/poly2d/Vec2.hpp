#pragma once
// -----------------------------------------------------------------------------
// Vec2 -- minimal 2D vector used across the mesher.
// -----------------------------------------------------------------------------
#include <cmath>

namespace poly2d {

struct Vec2 {
    double x{0.0};
    double y{0.0};

    constexpr Vec2() = default;
    constexpr Vec2(double x_, double y_) : x(x_), y(y_) {}

    constexpr Vec2 operator+(const Vec2& o) const { return {x + o.x, y + o.y}; }
    constexpr Vec2 operator-(const Vec2& o) const { return {x - o.x, y - o.y}; }
    constexpr Vec2 operator*(double s) const { return {x * s, y * s}; }
    constexpr Vec2 operator/(double s) const { return {x / s, y / s}; }
};

constexpr Vec2 operator*(double s, const Vec2& v) { return {v.x * s, v.y * s}; }

constexpr double dot(const Vec2& a, const Vec2& b) { return a.x * b.x + a.y * b.y; }
// 2D cross product (z-component); >0 means b is CCW from a.
constexpr double cross(const Vec2& a, const Vec2& b) { return a.x * b.y - a.y * b.x; }

inline double norm(const Vec2& v) { return std::sqrt(v.x * v.x + v.y * v.y); }
inline double norm2(const Vec2& v) { return v.x * v.x + v.y * v.y; }

inline Vec2 normalized(const Vec2& v) {
    const double n = norm(v);
    return (n > 0.0) ? Vec2{v.x / n, v.y / n} : Vec2{0.0, 0.0};
}

// Left normal of a direction (rotate +90 deg): points to the interior for a
// CCW-oriented boundary edge.
constexpr Vec2 leftNormal(const Vec2& dir) { return {-dir.y, dir.x}; }

inline double dist(const Vec2& a, const Vec2& b) { return norm(a - b); }

} // namespace poly2d
