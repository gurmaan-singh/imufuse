// Runs the three filters on the same synthetic flight and prints a table; with --csv writes per-sample errors.
#include <cstdio>
#include <cstring>
#include <string>

#include "imufuse/sim.hpp"

using namespace imufuse;

int main(int argc, char** argv) {
  const char* csv = nullptr;
  double seconds = 60.0;
  bool no_mag = false;
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--csv") == 0 && i + 1 < argc) csv = argv[++i];
    else if (std::strcmp(argv[i], "--seconds") == 0 && i + 1 < argc) seconds = std::stod(argv[++i]);
    else if (std::strcmp(argv[i], "--no-mag") == 0) no_mag = true;
    else { std::printf("usage: imufuse_sim [--seconds N] [--no-mag] [--csv out.csv]\n"); return 2; }
  }
  const ImuSpec spec;
  auto flight = make_flight(seconds, 0.005, spec);
  if (no_mag) for (auto& s : flight) s.mag = Vec3{};

  Complementary comp(0.98);
  Mahony mahony(2.0, 0.5);
  Ekf ekf;
  // all three start 60 degrees off so the plot shows convergence
  const Quat wrong = flight.front().truth * Quat::from_axis_angle({1, 1, 0}, 60 * 3.14159265358979323846 / 180);
  comp.set_attitude(wrong); mahony.set_attitude(wrong); ekf.set_attitude(wrong);

  const Result rc = run(comp, flight, spec), rm = run(mahony, flight, spec), re = run(ekf, flight, spec);
  std::printf("%-14s %10s %10s %14s\n", "filter", "rms (deg)", "max (deg)", "bias err (rad/s)");
  std::printf("%-14s %10.2f %10.2f %14.4f\n", "complementary", rc.rms_deg, rc.max_deg, rc.bias_err);
  std::printf("%-14s %10.2f %10.2f %14.4f\n", "mahony", rm.rms_deg, rm.max_deg, rm.bias_err);
  std::printf("%-14s %10.2f %10.2f %14.4f\n", "ekf", re.rms_deg, re.max_deg, re.bias_err);
  std::printf("(scored from t = 10 s; %zu samples at 200 Hz; true gyro bias %.3f rad/s)\n", flight.size(), spec.gyro_bias.norm());

  if (csv) {
    std::FILE* f = std::fopen(csv, "w");
    if (!f) { std::perror(csv); return 1; }
    // per-sample total error of each filter, plus the EKF's error split into tilt and yaw, plus the true angles
    Ekf again; again.set_attitude(wrong);
    std::fprintf(f, "t,complementary,mahony,ekf,ekf_tilt,ekf_yaw,roll,pitch,yaw\n");
    double prev = flight.front().t;
    for (std::size_t i = 0; i < flight.size(); ++i) {
      const double dt = i ? flight[i].t - prev : 0.005;
      prev = flight[i].t;
      again.update(flight[i].gyro, flight[i].accel, flight[i].mag, dt);
      if (i % 4) continue;
      const auto e = flight[i].truth.to_euler();
      const double deg = 180 / 3.14159265358979323846;
      std::fprintf(f, "%.3f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", flight[i].t, rc.err_deg[i], rm.err_deg[i], re.err_deg[i],
                   tilt_error(again.attitude(), flight[i].truth) * deg, heading_error(again.attitude(), flight[i].truth) * deg, e[0], e[1], e[2]);
    }
    std::fclose(f);
  }
  return 0;
}
