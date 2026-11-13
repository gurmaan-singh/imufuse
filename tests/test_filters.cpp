#include <gtest/gtest.h>

#include "imufuse/sim.hpp"

using namespace imufuse;
constexpr double kPi = 3.14159265358979323846;

namespace {
const ImuSpec kSpec;
const std::vector<Sample>& flight() { static const auto f = make_flight(40.0, 0.005, kSpec); return f; }
Quat sixty_off(const Quat& truth) { return truth * Quat::from_axis_angle({1, 1, 0}, 60 * kPi / 180); }
}  // namespace

TEST(Static, TiltFromAccelMatchesTheAttitude) {
  const Quat q = Quat::from_euler(0.3, -0.2, 0.0);
  const Vec3 a = q.rotate_inverse(kUp * 9.81);
  const auto [roll, pitch] = tilt_from_accel(a);
  EXPECT_NEAR(roll, 0.3, 1e-9); EXPECT_NEAR(pitch, -0.2, 1e-9);
}

TEST(Static, YawFromMagIsTiltCompensated) {
  const Quat q = Quat::from_euler(0.3, -0.2, 1.1);
  const Vec3 field{std::cos(1.1), 0, -std::sin(1.1)};
  const Vec3 m = q.rotate_inverse(field);
  EXPECT_NEAR(yaw_from_mag(m, 0.3, -0.2), 1.1, 1e-9);
}

TEST(Complementary, ConvergesFromSixtyDegreesOff) {
  Complementary f(0.98);
  f.set_attitude(sixty_off(flight().front().truth));
  const Result r = run(f, flight(), kSpec);
  EXPECT_LT(r.rms_deg, 3.0);
}

TEST(Mahony, ConvergesAndFindsTheBias) {
  Mahony f(2.0, 0.5);
  f.set_attitude(sixty_off(flight().front().truth));
  const Result r = run(f, flight(), kSpec);
  EXPECT_LT(r.rms_deg, 1.5);
  EXPECT_LT(r.bias_err, 0.3 * kSpec.gyro_bias.norm());
}

TEST(Ekf, ConvergesAndFindsTheBias) {
  Ekf f;
  f.set_attitude(sixty_off(flight().front().truth));
  const Result r = run(f, flight(), kSpec);
  EXPECT_LT(r.rms_deg, 1.0);
  EXPECT_LT(r.bias_err, 0.2 * kSpec.gyro_bias.norm());
}

TEST(Ekf, CovarianceShrinksAndStaysSymmetric) {
  Ekf f;
  const double before = f.covariance()(0, 0);
  run(f, flight(), kSpec);
  const auto& P = f.covariance();
  EXPECT_LT(P(0, 0), before * 0.1);
  for (std::size_t i = 0; i < 6; ++i) for (std::size_t j = 0; j < 6; ++j) EXPECT_NEAR(P(i, j), P(j, i), 1e-9);
}

TEST(Ekf, WithoutAMagnetometerYawDriftsButTiltHolds) {
  auto fl = flight();
  for (auto& s : fl) s.mag = Vec3{};
  Ekf g;
  g.set_attitude(sixty_off(fl.front().truth));
  double tilt = 0, heading = 0; std::size_t n = 0;
  for (const Sample& s : fl) {
    g.update(s.gyro, s.accel, s.mag, 0.005);
    if (s.t < 10) continue;
    const double te = tilt_error(g.attitude(), s.truth), he = heading_error(g.attitude(), s.truth);
    tilt += te * te; heading += he * he; ++n;
  }
  EXPECT_LT(std::sqrt(tilt / static_cast<double>(n)) * 180 / kPi, 1.0);
  EXPECT_GT(std::sqrt(heading / static_cast<double>(n)) * 180 / kPi, 3.0);  // heading is unobservable without the field
}

TEST(AllFilters, GyroOnlyDriftsBecauseOfTheBias) {
  Quat q = flight().front().truth;
  for (std::size_t i = 1; i < flight().size(); ++i) q = q.integrate(flight()[i].gyro, 0.005);
  // 0.027 rad/s for 40 s would be 62 degrees about a fixed axis; the motion mixes the axes, so it ends 10-20 off
  EXPECT_GT(q.angle_to(flight().back().truth) * 180 / kPi, 10.0);
}
