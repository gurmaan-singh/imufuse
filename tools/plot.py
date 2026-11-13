"""Draw the error curves from the simulator's CSVs and write docs/run.json for the portfolio replay.

    build/Release/imufuse_sim --csv out/errors.csv
    build/Release/imufuse_sim --no-mag --csv out/errors-nomag.csv
    python tools/plot.py
"""
import csv, json, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "img"; OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": True, "grid.alpha": 0.3})

def load(name):
    rows = list(csv.DictReader(open(ROOT / "out" / name, encoding="utf-8")))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}

d = load("errors.csv")
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
for k, label in [("complementary", "complementary"), ("mahony", "Mahony"), ("ekf", "multiplicative EKF")]:
    ax[0].plot(d["t"], d[k], label=label, lw=1); ax[1].plot(d["t"], d[k], label=label, lw=1)
ax[0].set(xlim=(0, 8), ylim=(0, 65), xlabel="time (s)", ylabel="attitude error (deg)", title="Convergence from 60 degrees off")
ax[1].set(xlim=(10, 60), ylim=(0, 5), xlabel="time (s)", title="Steady state (gyro bias 0.027 rad/s)")
ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "errors.png"); plt.close(fig)

n = load("errors-nomag.csv")
fig, ax = plt.subplots(figsize=(9, 3.4))
ax.plot(n["t"], n["ekf_tilt"], label="EKF tilt error", lw=1)
ax.plot(n["t"], [abs(v) for v in n["ekf_yaw"]], label="EKF heading error", lw=1)
ax.plot(n["t"], n["mahony"], label="Mahony total error", lw=1)
ax.set(xlim=(0, 60), ylim=(0, 30), xlabel="time (s)", ylabel="error (deg)", title="No magnetometer: tilt is observable, heading is not")
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(OUT / "no-mag.png"); plt.close(fig)

step = max(1, len(d["t"]) // 320)
run = {k: [round(v, 3) for v in d[k][::step]] for k in ("t", "complementary", "mahony", "ekf")}
(ROOT / "docs" / "run.json").write_text(json.dumps(run, separators=(",", ":")), encoding="utf-8")
print("wrote", OUT / "errors.png", OUT / "no-mag.png", ROOT / "docs" / "run.json")
