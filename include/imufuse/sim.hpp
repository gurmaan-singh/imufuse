// A synthetic flight: smooth angular rates, true attitude by integration, and an IMU with noise and a gyro bias.
#pragma once
#include <random>
#include <vector>

#include "imufuse/math.hpp"
#include "imufuse/filters.hpp"

namespace imufuse {

struct Sample {
  double t;
  Quat truth;
  Vec3 gyro, accel, mag;
};

struct ImuSpec {
  double gyro_noise = 0.005;   // rad/s, white, per sample
  double accel_noise = 0.15;   // m/s^2
  double mag_noise = 0.02;     // unit-vector units
  Vec3 gyro_bias{0.02, -0.015, 0.01};
  double g = 9.81;
};

// Rates are sums of sines so the motion is smooth and never stops; wobble sets the amplitude.
inline Vec3 true_rate(double t, double wobble = 1.0) {
  return Vec3{0.8 * std::sin(0.5 * t), 0.6 * std::sin(0.3 * t + 1.0), 0.4 * std::sin(0.2 * t + 2.0)} * wobble;
}

inline std::vector<Sample> make_flight(double seconds = 60.0, double dt = 0.005, const ImuSpec& spec = {}, unsigned seed = 1,
                                       double wobble = 1.0) {
  std::mt19937 rng(seed);
  std::normal_distribution<double> n(0.0, 1.0);
  std::vector<Sample> out;
  Quat q = Quat::from_euler(0.2, -0.1, 0.5);
  const auto steps = static_cast<std::size_t>(seconds / dt);
  out.reserve(steps);
  for (std::size_t i = 0; i < steps; ++i) {
    const double t = static_cast<double>(i) * dt;
    const Vec3 w = true_rate(t, wobble);
    q = q.integrate(w, dt);
    Sample s;
    s.t = t;
    s.truth = q;
    s.gyro = w + spec.gyro_bias + Vec3{n(rng), n(rng), n(rng)} * spec.gyro_noise;
    s.accel = q.rotate_inverse(kUp * spec.g) + Vec3{n(rng), n(rng), n(rng)} * spec.accel_noise;  // no linear acceleration
    const Vec3 field = Vec3{std::cos(1.1), 0, -std::sin(1.1)};  // 63 degrees of inclination, like northern Europe
    s.mag = q.rotate_inverse(field) + Vec3{n(rng), n(rng), n(rng)} * spec.mag_noise;
    out.push_back(s);
  }
  return out;
}

// Attitude error split the way it matters: tilt = angle between the two "up" directions in the body frame,
// heading = angle between the two body x axes projected on the horizontal plane. Euler differences wrap and
// blow up near 90 degrees of pitch; these do not.
inline double tilt_error(const Quat& est, const Quat& truth) {
  const Vec3 a = est.rotate_inverse(kUp), b = truth.rotate_inverse(kUp);
  const double d = a.dot(b);
  return std::acos(d > 1 ? 1 : d < -1 ? -1 : d);
}
inline double heading_error(const Quat& est, const Quat& truth) {
  const Vec3 a = est.rotate({1, 0, 0}), b = truth.rotate({1, 0, 0});
  return std::atan2(a.y, a.x) - std::atan2(b.y, b.x) > 3.14159265358979323846
             ? std::atan2(a.y, a.x) - std::atan2(b.y, b.x) - 2 * 3.14159265358979323846
         : std::atan2(a.y, a.x) - std::atan2(b.y, b.x) < -3.14159265358979323846
             ? std::atan2(a.y, a.x) - std::atan2(b.y, b.x) + 2 * 3.14159265358979323846
             : std::atan2(a.y, a.x) - std::atan2(b.y, b.x);
}

struct Result {
  double rms_deg = 0;      // attitude error over the scoring window
  double max_deg = 0;
  double bias_err = 0;     // |estimated - true| gyro bias, rad/s, at the end
  std::vector<double> err_deg;  // per sample, for plotting
};

template <class Filter>
Result run(Filter& f, const std::vector<Sample>& flight, const ImuSpec& spec = {}, double score_from = 10.0) {
  Result r;
  r.err_deg.reserve(flight.size());
  double acc = 0; std::size_t n = 0;
  double prev = flight.empty() ? 0 : flight.front().t;
  for (const Sample& s : flight) {
    const double dt = s.t - prev > 0 ? s.t - prev : 0.005;
    prev = s.t;
    f.update(s.gyro, s.accel, s.mag, dt);
    const double e = f.attitude().angle_to(s.truth) * 180.0 / 3.14159265358979323846;
    r.err_deg.push_back(e);
    if (s.t >= score_from) { acc += e * e; ++n; if (e > r.max_deg) r.max_deg = e; }
  }
  r.rms_deg = n ? std::sqrt(acc / static_cast<double>(n)) : 0;
  r.bias_err = (f.gyro_bias() - spec.gyro_bias).norm();
  return r;
}

}  // namespace imufuse
