#include <gtest/gtest.h>

#include "imufuse/math.hpp"

using namespace imufuse;
constexpr double kPi = 3.14159265358979323846;

TEST(Quat, NinetyDegreesAboutZTakesXToY) {
  const Quat q = Quat::from_axis_angle({0, 0, 1}, kPi / 2);
  const Vec3 v = q.rotate({1, 0, 0});
  EXPECT_NEAR(v.x, 0, 1e-12); EXPECT_NEAR(v.y, 1, 1e-12); EXPECT_NEAR(v.z, 0, 1e-12);
  const Vec3 back = q.rotate_inverse(v);
  EXPECT_NEAR(back.x, 1, 1e-12); EXPECT_NEAR(back.y, 0, 1e-12);
}

TEST(Quat, ProductMatchesComposition) {
  const Quat a = Quat::from_axis_angle({1, 0, 0}, 0.7), b = Quat::from_axis_angle({0, 1, 0}, -0.4);
  const Vec3 v{0.3, -0.2, 0.9};
  const Vec3 lhs = (a * b).rotate(v), rhs = a.rotate(b.rotate(v));
  EXPECT_NEAR(lhs.x, rhs.x, 1e-12); EXPECT_NEAR(lhs.y, rhs.y, 1e-12); EXPECT_NEAR(lhs.z, rhs.z, 1e-12);
}

TEST(Quat, EulerRoundTrip) {
  const double r = 0.3, p = -0.5, y = 2.0;
  const auto e = Quat::from_euler(r, p, y).to_euler();
  EXPECT_NEAR(e[0], r, 1e-12); EXPECT_NEAR(e[1], p, 1e-12); EXPECT_NEAR(e[2], y, 1e-12);
}

TEST(Quat, IntegratingAConstantRateGivesTheExpectedAngle) {
  Quat q;
  for (int i = 0; i < 1000; ++i) q = q.integrate({0, 0, 0.5}, 0.001);  // 0.5 rad/s for 1 s
  EXPECT_NEAR(q.to_euler()[2], 0.5, 1e-9);
  EXPECT_NEAR(q.norm(), 1.0, 1e-12);
}

TEST(Quat, AngleBetweenAttitudes) {
  const Quat a = Quat::from_axis_angle({0, 1, 0}, 0.2), b = Quat::from_axis_angle({0, 1, 0}, 0.2 + 0.6);
  EXPECT_NEAR(a.angle_to(b), 0.6, 1e-12);
  EXPECT_NEAR(a.angle_to(a), 0.0, 1e-9);
}

TEST(Mat, InverseTimesItselfIsIdentity) {
  Mat<3, 3> m;
  m(0, 0) = 4; m(0, 1) = 1; m(0, 2) = 0.5; m(1, 0) = 1; m(1, 1) = 3; m(1, 2) = -0.2; m(2, 0) = 0.5; m(2, 1) = -0.2; m(2, 2) = 2;
  const Mat<3, 3> id = m * m.inverse();
  for (std::size_t i = 0; i < 3; ++i) for (std::size_t j = 0; j < 3; ++j) EXPECT_NEAR(id(i, j), i == j ? 1 : 0, 1e-12);
}

TEST(Mat, SkewMatrixReproducesTheCrossProduct) {
  const Vec3 a{1, 2, 3}, b{-2, 0.5, 4};
  const Mat<3, 3> S = skew(a);
  Mat<3, 1> bv; bv(0, 0) = b.x; bv(1, 0) = b.y; bv(2, 0) = b.z;
  const Mat<3, 1> c = S * bv;
  const Vec3 ref = a.cross(b);
  EXPECT_NEAR(c(0, 0), ref.x, 1e-12); EXPECT_NEAR(c(1, 0), ref.y, 1e-12); EXPECT_NEAR(c(2, 0), ref.z, 1e-12);
}
