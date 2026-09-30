"""Tabulate key metrics of all validation runs against the experiments.

Usage: python summarize.py  -> prints Markdown tables, writes results/summary.md"""

import os
import json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
from common import RESULTS, run_dir  # noqa: E402
EXP = json.load(open(os.path.join(HERE, "data", "experiments.json")))
MODELS = ["SA", "SA-Edwards", "SA-comp"]
out = []


def load(path):
    return json.load(open(path)) if os.path.exists(path) else None


# --- 2-D flat plate ----------------------------------------------------------
out.append("### Mach 8.2 cold-wall flat plate, x = -11 cm (Kussoy & Horstman inflow)\n")
out.append("| Model | theta [cm] | Cf x1e3 | q_w [W/cm^2] |")
out.append("|---|---|---|---|")
for m in MODELS:
    d = load(os.path.join(run_dir("adflow", "plate", "medium", m), "bl.json"))
    if d:
        r = d["x=-0.11"]
        out.append(f"| {m} | {r['theta_cm']:.4f} | {r['Cf'] * 1e3:.3f} | {r['qw_W_cm2']:.3f} |")
c = EXP["kussoy_fin15"]["conditions"]
out.append(f"| **Experiment** | {c['inflow_theta_cm']} | {c['inflow_Cf'] * 1e3:.3f} | {c['inflow_qw_W_cm2']} |\n")

# --- Sharp fin -----------------------------------------------------------------
e = EXP["kussoy_fin15"]
out.append("### Mach 8.2 sharp fin (peak / plateau values)\n")
out.append("| Model | plate P/Pinf peak (X=18.2) | plate P plateau (Z=4-6 cm) | plate Q/Qinf peak (X=16.45) | plate Cf peak x1e3 (X=15.5) | fin P/Pinf (Y=6-10 cm) | fin Q/Qinf (Y=6-10 cm) |")
out.append("|---|---|---|---|---|---|---|")


def band(d, key, lo, hi, xk):
    x = np.array(d[key][xk])
    v = np.array(d[key]["value"])
    return float(np.mean(v[(x >= lo) & (x <= hi)]))


for m in MODELS:
    d = load(os.path.join(run_dir("adflow", "fin", "medium", m), "lines.json"))
    if d:
        out.append(
            f"| {m} | {max(d['plate_p_X18.2']['value']):.2f} | {band(d, 'plate_p_X18.2', 4, 6, 'Z_cm'):.2f} | "
            f"{max(d['plate_q_X16.45']['value']):.2f} | {max(d['plate_cf_X15.5']['value']) * 1e3:.2f} | "
            f"{band(d, 'fin_p_X17.72', 6, 10, 'Y_cm'):.2f} | {band(d, 'fin_q_X16.13', 6, 10, 'Y_cm'):.2f} |"
        )


def eband(key, xk, yk, lo, hi):
    x = np.array(e[key][xk])
    v = np.array(e[key][yk])
    return float(np.mean(v[(x >= lo) & (x <= hi)]))


out.append(
    f"| **Experiment** | {max(e['plate_p_X18.2']['P_over_Pinf']):.2f} | {eband('plate_p_X18.2', 'Z_cm', 'P_over_Pinf', 4, 6):.2f} | "
    f"{max(e['plate_q_X16.45']['Q_over_Qinf']):.2f} | {max(e['plate_cf_X15.5']['Cf']) * 1e3:.2f} | "
    f"{eband('fin_p_X17.72', 'Y_cm', 'P_over_Pinf', 6, 10):.2f} | {eband('fin_q_X16.13', 'Y_cm', 'Q_over_Qinf', 6, 10):.2f} |"
)
p = e["paper_SA_curves_approx"]
out.append(f"| Edwards & Chandra SA (approx.) | {p['plate_p_peak']} | {p['plate_p_plateau']} | {p['plate_q_peak']} | {p['plate_cf_peak'] * 1e3:.1f} | {p['fin_p_plateau']} | {p['fin_q_plateau']} |\n")

# --- Cylinder / flare -------------------------------------------------------
w = EXP["wideman_flare"]
t = w["text_metrics"]
out.append("### Mach 2.89 cylinder / offset flare\n")
out.append("| Model | delta99 x=-2 cm [cm] | Cf upstream x1e3 | phi=0 P rise | phi=180 x_sep [cm] | phi=180 x_att [cm] | phi=180 Cf min x1e3 | phi=90 Cf max x1e3 |")
out.append("|---|---|---|---|---|---|---|---|")
for m in MODELS:
    d = load(os.path.join(run_dir("adflow", "flare", "medium", m), "lines.json"))
    if d:
        p0, p90, p180 = d["phi0"], d["phi90"], d["phi180"]
        x0 = np.array(p0["x_cm"])
        prise = float(np.max(np.array(p0["P_over_Pinf"])[(x0 > 0) & (x0 < 10)]) / p0["P_upstream_x-10cm"])
        zc = p180["Cfx_zero_crossings_cm"]
        x180 = np.array(p180["x_cm"])
        c180 = np.array(p180["Cfx"])
        x90 = np.array(p90["x_cm"])
        c90 = np.array(p90["Cfx"])
        out.append(
            f"| {m} | {d['incoming_BL']['delta99_cm']:.2f} | {p0['Cfx_upstream_x-10cm'] * 1e3:.3f} | {prise:.2f} | "
            f"{zc[0] if zc else float('nan'):.2f} | {zc[1] if len(zc) > 1 else float('nan'):.2f} | "
            f"{c180[(x180 > 3) & (x180 < 9)].min() * 1e3:.2f} | {c90[(x90 > 5) & (x90 < 14)].max() * 1e3:.2f} |"
        )
out.append(
    f"| **Experiment** | {w['conditions']['delta_cm']} | {w['conditions']['Cf_upstream_LISF'] * 1e3:.2f} (LISF), {w['conditions']['Cf_upstream_profile'] * 1e3:.2f} (profile) | "
    f"{t['phi0']['p_rise_ratio']} | {t['phi180']['x_sep_cf_cm']} | {t['phi180']['x_att_cm']} | {t['phi180']['cf_min'] * 1e3:.2f} | {t['phi90']['cf_max'] * 1e3:.2f} |"
)

txt = "\n".join(out)
print(txt)
open(os.path.join(RESULTS, "summary.md"), "w").write(txt + "\n")
