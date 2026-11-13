// Three attitude estimators with the same interface:
//   void update(const Vec3& gyro, const Vec3& accel, const Vec3& mag, double dt);
//   Quat attitude() const;  Vec3 gyro_bias() const;
// Frames: world is z-up with magnetic north along +x. The accelerometer measures specific force, so at rest it
// reads +g along the body's "up". Rates are rad/s, accelerations m/s^2 (only the direction is used), mag any unit.
#pragma once
#include "imufuse/math.hpp"

namespace imufuse {

constexpr Vec3 kUp{0, 0, 1};
constexpr Vec3 kNorth{1, 0, 0};

// Tilt-compensated yaw from a body-frame magnetometer reading and the current roll/pitch.
inline double yaw_from_mag(const Vec3& mag, double roll, double pitch) {
  const double cr = std::cos(roll), sr = std::sin(roll), cp = std::cos(pitch), sp = std::sin(pitch);
  const double mx = mag.x * cp + mag.y * sr * sp + mag.z * cr * sp;
  const double my = mag.y * cr - mag.z * sr;
  return std::atan2(-my, mx);
}

// Roll and pitch from the accelerometer alone (valid only when not accelerating).
inline std::array<double, 2> tilt_from_accel(const Vec3& a) {
  return {std::atan2(a.y, a.z), std::atan2(-a.x, std::sqrt(a.y * a.y + a.z * a.z))};
}

// 1. Complementary filter on Euler angles: integrate the gyro, pull slowly towards accel/mag angles.
class Complementary {
 public:
  explicit Complementary(double alpha = 0.98) : alpha_(alpha) {}
  void update(const Vec3& gyro, const Vec3& accel, const Vec3& mag, double dt) {
    const auto e = q_.integrate(gyro, dt).to_euler();
    const auto [ra, pa] = tilt_from_accel(accel);
    const double ya = yaw_from_mag(mag, ra, pa);
    const double r = blend(e[0], ra), p = blend(e[1], pa), y = blend(e[2], ya);
    q_ = Quat::from_euler(r, p, y);
  }
  Quat attitude() const { return q_; }
  Vec3 gyro_bias() const { return {}; }  // it does not estimate one
  void set_attitude(const Quat& q) { q_ = q; }

 private:
  double blend(double gyro_angle, double ref_angle) const {  // shortest-way blend so yaw can wrap
    double d = ref_angle - gyro_angle;
    while (d > 3.14159265358979323846) d -= 2 * 3.14159265358979323846;
    while (d < -3.14159265358979323846) d += 2 * 3.14159265358979323846;
    return gyro_angle + (1 - alpha_) * d;
  }
  double alpha_;
  Quat q_;
};

// 2. Mahony: feed the error between measured and predicted reference directions back into the rate, PI style.
class Mahony {
 public:
  Mahony(double kp = 2.0, double ki = 0.05) : kp_(kp), ki_(ki) {}
  void update(const Vec3& gyro, const Vec3& accel, const Vec3& mag, double dt) {
    const Vec3 a = accel.normalized();
    const Vec3 up_b = q_.rotate_inverse(kUp);
    Vec3 err = a.cross(up_b);
    if (mag.norm() > 0) {
      // only the horizontal part of the field is informative, so compare against the projected reference
      const Vec3 m_w = q_.rotate(mag.normalized());
      const Vec3 ref_w = Vec3{std::sqrt(m_w.x * m_w.x + m_w.y * m_w.y), 0, m_w.z};
      err = err + mag.normalized().cross(q_.rotate_inverse(ref_w));
    }
    bias_ = bias_ - err * (ki_ * dt);
    const Vec3 omega = gyro - bias_ + err * kp_;
    q_ = q_.integrate(omega, dt);
  }
  Quat attitude() const { return q_; }
  Vec3 gyro_bias() const { return bias_; }
  void set_attitude(const Quat& q) { q_ = q; }

 private:
  double kp_, ki_;
  Quat q_;
  Vec3 bias_;
};

// 3. Multiplicative EKF. Nominal state is the quaternion and the gyro bias; the error state is
// [small rotation (3), bias error (3)]. Accelerometer and magnetometer each give a direction measurement.
class Ekf {
 public:
  struct Noise {
    double gyro = 0.01;      // rad/s/sqrt(Hz)-ish, per step
    double bias_walk = 1e-4; // how fast the bias may drift
    double accel = 0.1;      // direction noise, unit vector units
    double mag = 0.2;
  };
  explicit Ekf(Noise n = {}) : n_(n) {
    P_ = Mat<6, 6>::identity() * 0.5;
    for (int i = 3; i < 6; ++i) P_(static_cast<std::size_t>(i), static_cast<std::size_t>(i)) = 1e-3;
  }

  void update(const Vec3& gyro, const Vec3& accel, const Vec3& mag, double dt) {
    predict(gyro, dt);
    correct(accel.normalized(), q_.rotate_inverse(kUp), n_.accel);
    if (mag.norm() > 0) {
      const Vec3 m_w = q_.rotate(mag.normalized());
      const Vec3 ref_w = Vec3{std::sqrt(m_w.x * m_w.x + m_w.y * m_w.y), 0, m_w.z};
      correct(mag.normalized(), q_.rotate_inverse(ref_w), n_.mag);
    }
  }
  Quat attitude() const { return q_; }
  Vec3 gyro_bias() const { return bias_; }
  const Mat<6, 6>& covariance() const { return P_; }
  void set_attitude(const Quat& q) { q_ = q; }

 private:
  void predict(const Vec3& gyro, double dt) {
    const Vec3 omega = gyro - bias_;
    q_ = q_.integrate(omega, dt);
    // F = [ I - [w x] dt   -I dt ; 0  I ]
    Mat<6, 6> F = Mat<6, 6>::identity();
    const Mat<3, 3> S = skew(omega) * dt;
    for (std::size_t i = 0; i < 3; ++i)
      for (std::size_t j = 0; j < 3; ++j) { F(i, j) -= S(i, j); F(i, j + 3) = i == j ? -dt : 0; }
    Mat<6, 6> Q;
    for (std::size_t i = 0; i < 3; ++i) { Q(i, i) = n_.gyro * n_.gyro * dt; Q(i + 3, i + 3) = n_.bias_walk * n_.bias_walk * dt; }
    P_ = F * P_ * F.transposed() + Q;
  }

  // z = measured unit direction in body frame, h = predicted one; H = [ [h x]  0 ]
  void correct(const Vec3& z, const Vec3& h, double sigma) {
    Mat<3, 6> H;
    const Mat<3, 3> hx = skew(h);
    for (std::size_t i = 0; i < 3; ++i) for (std::size_t j = 0; j < 3; ++j) H(i, j) = hx(i, j);
    Mat<3, 3> R = Mat<3, 3>::identity() * (sigma * sigma);
    const Mat<3, 3> S = H * P_ * H.transposed() + R;
    const Mat<6, 3> K = P_ * H.transposed() * S.inverse();
    const Vec3 r = z - h;
    Mat<3, 1> rv; rv(0, 0) = r.x; rv(1, 0) = r.y; rv(2, 0) = r.z;
    const Mat<6, 1> dx = K * rv;
    const Vec3 dtheta{dx(0, 0), dx(1, 0), dx(2, 0)};
    q_ = (q_ * Quat{1, dtheta.x / 2, dtheta.y / 2, dtheta.z / 2}).normalized();
    bias_ = bias_ + Vec3{dx(3, 0), dx(4, 0), dx(5, 0)};
    const Mat<6, 6> IKH = Mat<6, 6>::identity() - K * H;
    P_ = IKH * P_ * IKH.transposed() + K * R * K.transposed();  // Joseph form keeps P symmetric
  }

  Noise n_;
  Quat q_;
  Vec3 bias_;
  Mat<6, 6> P_;
};

}  // namespace imufuse
