// Small vector, quaternion and fixed-size matrix types. Enough for attitude filters, nothing more.
#pragma once
#include <array>
#include <cmath>
#include <cstddef>

namespace imufuse {

struct Vec3 {
  double x = 0, y = 0, z = 0;
  Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
  Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
  Vec3 operator-() const { return {-x, -y, -z}; }
  double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
  Vec3 cross(const Vec3& o) const { return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x}; }
  double norm() const { return std::sqrt(dot(*this)); }
  Vec3 normalized() const { const double n = norm(); return n > 0 ? *this * (1.0 / n) : *this; }
};

// Unit quaternion q = (w, x, y, z) rotating body vectors into the world frame: v_world = q * v_body * q^-1.
struct Quat {
  double w = 1, x = 0, y = 0, z = 0;

  static Quat identity() { return {}; }
  static Quat from_axis_angle(const Vec3& axis, double angle) {
    const Vec3 a = axis.normalized();
    const double s = std::sin(angle / 2);
    return {std::cos(angle / 2), a.x * s, a.y * s, a.z * s};
  }
  // small rotation vector -> quaternion (exact for any size, used for gyro integration)
  static Quat exp(const Vec3& rotvec) {
    const double n = rotvec.norm();
    if (n < 1e-12) return {1, rotvec.x / 2, rotvec.y / 2, rotvec.z / 2};
    return from_axis_angle(rotvec, n);
  }
  // ZYX Euler angles (yaw about z, then pitch about y, then roll about x), radians
  static Quat from_euler(double roll, double pitch, double yaw) {
    const double cr = std::cos(roll / 2), sr = std::sin(roll / 2), cp = std::cos(pitch / 2), sp = std::sin(pitch / 2);
    const double cy = std::cos(yaw / 2), sy = std::sin(yaw / 2);
    return {cr * cp * cy + sr * sp * sy, sr * cp * cy - cr * sp * sy, cr * sp * cy + sr * cp * sy, cr * cp * sy - sr * sp * cy};
  }

  Quat operator*(const Quat& o) const {
    return {w * o.w - x * o.x - y * o.y - z * o.z, w * o.x + x * o.w + y * o.z - z * o.y,
            w * o.y - x * o.z + y * o.w + z * o.x, w * o.z + x * o.y - y * o.x + z * o.w};
  }
  Quat conj() const { return {w, -x, -y, -z}; }
  double norm() const { return std::sqrt(w * w + x * x + y * y + z * z); }
  Quat normalized() const { const double n = norm(); return {w / n, x / n, y / n, z / n}; }
  double dot(const Quat& o) const { return w * o.w + x * o.x + y * o.y + z * o.z; }

  Vec3 rotate(const Vec3& v) const {  // body -> world
    const Vec3 u{x, y, z};
    const Vec3 t = u.cross(v) * 2.0;
    return v + t * w + u.cross(t);
  }
  Vec3 rotate_inverse(const Vec3& v) const { return conj().rotate(v); }  // world -> body

  // integrate a body-frame angular rate over dt
  Quat integrate(const Vec3& omega, double dt) const { return (*this * exp(omega * dt)).normalized(); }

  std::array<double, 3> to_euler() const {  // roll, pitch, yaw
    const double sinr = 2 * (w * x + y * z), cosr = 1 - 2 * (x * x + y * y);
    const double sinp = 2 * (w * y - z * x);
    const double siny = 2 * (w * z + x * y), cosy = 1 - 2 * (y * y + z * z);
    const double pitch = std::fabs(sinp) >= 1 ? std::copysign(3.14159265358979323846 / 2, sinp) : std::asin(sinp);
    return {std::atan2(sinr, cosr), pitch, std::atan2(siny, cosy)};
  }
  // angle between two attitudes, radians
  double angle_to(const Quat& o) const {
    const double d = std::fabs(dot(o));
    return 2 * std::acos(d > 1 ? 1 : d);
  }
};

// Fixed-size matrix with just the operations the EKF needs.
template <std::size_t R, std::size_t C>
struct Mat {
  std::array<double, R * C> a{};
  double& operator()(std::size_t r, std::size_t c) { return a[r * C + c]; }
  double operator()(std::size_t r, std::size_t c) const { return a[r * C + c]; }
  static Mat identity() { Mat m; for (std::size_t i = 0; i < R && i < C; ++i) m(i, i) = 1; return m; }
  Mat operator+(const Mat& o) const { Mat m; for (std::size_t i = 0; i < R * C; ++i) m.a[i] = a[i] + o.a[i]; return m; }
  Mat operator-(const Mat& o) const { Mat m; for (std::size_t i = 0; i < R * C; ++i) m.a[i] = a[i] - o.a[i]; return m; }
  Mat operator*(double s) const { Mat m; for (std::size_t i = 0; i < R * C; ++i) m.a[i] = a[i] * s; return m; }
  template <std::size_t K>
  Mat<R, K> operator*(const Mat<C, K>& o) const {
    Mat<R, K> m;
    for (std::size_t i = 0; i < R; ++i)
      for (std::size_t j = 0; j < K; ++j) {
        double s = 0;
        for (std::size_t k = 0; k < C; ++k) s += (*this)(i, k) * o(k, j);
        m(i, j) = s;
      }
    return m;
  }
  Mat<C, R> transposed() const { Mat<C, R> m; for (std::size_t i = 0; i < R; ++i) for (std::size_t j = 0; j < C; ++j) m(j, i) = (*this)(i, j); return m; }
  // Gauss-Jordan with partial pivoting; only used on small symmetric positive matrices
  Mat inverse() const {
    static_assert(R == C, "square only");
    Mat m = *this, inv = identity();
    for (std::size_t c = 0; c < R; ++c) {
      std::size_t p = c;
      for (std::size_t r = c + 1; r < R; ++r) if (std::fabs(m(r, c)) > std::fabs(m(p, c))) p = r;
      if (p != c) for (std::size_t k = 0; k < R; ++k) { std::swap(m(c, k), m(p, k)); std::swap(inv(c, k), inv(p, k)); }
      const double d = m(c, c);
      for (std::size_t k = 0; k < R; ++k) { m(c, k) /= d; inv(c, k) /= d; }
      for (std::size_t r = 0; r < R; ++r) {
        if (r == c) continue;
        const double f = m(r, c);
        if (f == 0) continue;
        for (std::size_t k = 0; k < R; ++k) { m(r, k) -= f * m(c, k); inv(r, k) -= f * inv(c, k); }
      }
    }
    return inv;
  }
};

inline Mat<3, 3> skew(const Vec3& v) {
  Mat<3, 3> m;
  m(0, 1) = -v.z; m(0, 2) = v.y; m(1, 0) = v.z; m(1, 2) = -v.x; m(2, 0) = -v.y; m(2, 1) = v.x;
  return m;
}

}  // namespace imufuse
