# V&V of the SA modifications in this ADflow branch

Verification and validation of the two Spalart-Allmaras modifications
implemented in `src/turbulence/saCorrections.F90` (branch `sa-comp-edwards`):

* **SA-Edwards** (`turbulenceModel = "SA-Edwards"`) - Edwards & Chandra,
  AIAA 94-2275 (1994), Sec. 4.3, Eqs. 20, 24, 25.
* **SA-comp** (`useCompressibilitySA = True`, constant `SAc5`) - Spalart,
  AIAA 2000-2306: -C5 nuTilde^2/a^2 (du_i/dx_j)^2, C5 = 3.5.

Everything needed to reproduce the evidence below is in this folder
(`run_all.sh`), except the unit/regression tests, which live with ADflow's
tests in `tests/reg_tests/test_sa_corrections.py`.

## Summary of the evidence

| Check | Result |
|---|---|
| Equations vs. the papers, the NASA Turbulence Modeling Resource and the SU2 source | Match (SU2 differs from the paper in one strain term; this code follows the paper) |
| Independent code review | No errors; findings addressed |
| Unit/regression tests (`test_sa_corrections.py`) | 11 pass, 1 skip (N/A): blockette vs block residual 1e-12, forward AD vs FD (1st-order FD convergence), forward/reverse/reverse-fast dot products 1e-11, C5 linearity |
| Upstream ADflow regression tests | 103 / 103 pass |
| **Code-to-code (ADflow vs SU2), SA-Edwards increment** | Same effect in both codes to 0.1-0.6 percentage points on all three grids |
| Grid refinement, 2-D plate (3 levels) | Grid independent: Cf, q_w <= 0.2 % change, GCI <= 0.2 % |
| Grid refinement, fin / flare (2 levels) | <= 1.5 % change (except the flare separation-bubble Cf minimum, ~7 %) |
| Validation vs. Kussoy & Horstman (M 8.2 sharp fin) | Pressures within a few %; SA reproduces the paper's SA results; heating overpredicted by all models, as in the paper |
| Validation vs. Wideman et al. (M 2.89 cylinder / flare) | Separation / pressure plateau locations captured; same model deficiencies as in the paper |
| SA-comp in hypersonic wall-bounded flow | Cf and q_w 42 % too low on the Mach 8 plate - do not use for such flows (consistent with Xue et al., AIAA J. 2024) |

All documents are in `results/` (`summary.md`, `grid_study.md`,
`su2_comparison.md`) and all figures in `results/figures/`:

* `plate_ADflow_vs_SU2.png` - wall Cf and q_w along the plate, ADflow vs SU2, three grids
* `plate_profiles_ADflow_vs_SU2.png` - boundary-layer profiles at x = -11 cm
* `fin_vs_experiment.png` - Kussoy & Horstman sharp fin
* `flare_vs_experiment.png` - Wideman et al. cylinder / flare

## 1. Verification

### 1.1 Implementation checks

* **Equations** read from Edwards & Chandra (paper in `../research_docs/`),
  Spalart (2000) via the NASA Turbulence Modeling Resource, and compared term
  by term with SU2's implementation (`SA_OPTIONS= EDWARDS` / `COMPRESSIBILITY`).
* **Regression tests** (`../tests/reg_tests/test_sa_corrections.py`, MDO tutorial
  RANS wing, SA-Edwards / SA+comp / SA-Edwards+comp): consistency of the two
  residual paths, forward AD vs finite differences, forward vs reverse vs
  reverse-fast AD (dot-product tests), and that the SA-comp term changes only
  the turbulence residual, linearly in C5. The adjoint code for the new
  routines was generated with Tapenade 3.16 (`src/adjoint/output*/saCorrections_*.f90`).
* **Log-layer check** of the Edwards form: Stilde and r match standard SA where
  standard SA is well behaved; Stilde stays positive where standard SA's
  Stilde goes negative (chi ~ 5).

### 1.2 Code-to-code: ADflow vs SU2 (2-D Mach 8.2 cold-wall plate)

`su2/` is the same flat plate for SU2 v8.5.0 (built in `~/packages/SU2`), on
exactly the same grids, with the physics matched to ADflow (SA without ft2 /
Edwards option, Sutherland 1.716e-5 / 273.15 K / 110.55 K, Pr 0.72, Pr_t 0.9,
nuTilde_inf = 1.342 nu_inf, T_w = 300 K; freestream state agrees to 0.001 %).
SU2 is run first order, then second order (Roe, MUSCL + van Albada). Its global
residual stalls because of the corner node at the plate leading edge; away
from the leading edge the density residual is 1e-5 (coarse) to 1e-7
(medium/fine) and the wall quantities change by <= 0.6 % over the last 2000
iterations.

**SA-Edwards relative to SA at x = -11 cm** (the direct test of the Edwards
implementation, independent of each code's baseline discretization):

| grid | Cf ADflow / SU2 | q_w ADflow / SU2 | theta ADflow / SU2 |
|---|---|---|---|
| coarse (7k cells) | -7.49 / -7.61 % | -7.66 / -7.90 % | -11.6 / -12.1 % |
| medium (27k) | -7.52 / -7.62 % | -7.69 / -7.87 % | -11.4 / -12.1 % |
| fine (109k) | -7.49 / -7.83 % | -7.66 / -8.24 % | -11.4 / -12.2 % |

**Baseline (SA) difference between the codes.** SU2 gives lower Cf/q_w than
ADflow on the same grid (x = -11 cm: Cf -6.9 / -5.4 / -3.9 %, q_w -4.8 / -4.2 /
-2.8 % on coarse / medium / fine). The difference is largest at the plate
leading edge (~-30 % 4 cm downstream of it), decays downstream, and decreases
with grid refinement, while ADflow is already grid converged: SU2 converges
more slowly towards the same solution on these grids. It is not the freestream
turbulence level (nuTilde_inf = 3 nu_inf changes it by < 0.3 %) nor ADflow's
strain-based production (vorticity changes ADflow's result by 0.02 %).
Profiles agree near the wall (nuTilde and T within ~1 % for y < 1 mm); the
differences are in the outer layer. Each code's reported heat flux matches
its own k dT/dy at the wall (`su2/wallgrad.py`). Full tables:
`results/su2_comparison.md`.

### 1.3 Grid refinement (`grid_study.py` -> `results/grid_study.md`)

Grid families from `meshes.py` (`ref` scales the number of intervals and
divides all spacings, first wall spacing included); the medium grids are the
ones used for all other results.

* **2-D plate** (r = 2; 7k / 27k / 109k cells; SA and SA-Edwards): Cf and q_w
  change by <= 0.2 % between levels (GCI_fine <= 0.2 %); theta converges
  monotonically with observed order 3-4 (GCI_fine <= 0.03 %).
* **Sharp fin** (r = sqrt 2; 0.34 M and 0.97 M cells): all extracted quantities
  change by <= 1.5 % from coarse to medium.
* **Cylinder / flare** (r = sqrt 2; 0.22 M and 0.63 M cells): <= 1.4 %, except the
  phi = 180 deg separation-bubble Cf minimum (-7 % SA, -6 % SA-Edwards) and the
  phi = 90 deg Cf maximum (-3 %).

## 2. Validation (cases of Edwards & Chandra, AIAA 94-2275)

| Case | Conditions | Medium grid | Data |
|---|---|---|---|
| 2-D plate | M 8.2, Re 5.32e6/m, T 76.89 K, T_w 300 K | 171 x 80 x 2 cells | Kussoy & Horstman inflow survey, x = -11 cm |
| Sharp fin (15 deg), Kussoy & Horstman | same, T_w/T 3.9 | 0.97 M cells | Edwards & Chandra Figs. 2-6 |
| Cylinder / offset 20 deg flare, Wideman et al. | M 2.89, p 5.54 kPa, T 105 K, adiabatic | 0.63 M cells | NASA CR-202617 Figs. 11-13 + text |

Full tables: `results/summary.md` (medium grids).

* **Incoming boundary layers** (SA): Kussoy x = -11 cm: theta 0.094 cm (exp. 0.094),
  Cf 0.95e-3 (1.0e-3), q_w 1.10 W/cm^2 (1.04); Wideman x = -2 cm: delta99 1.14 cm
  (1.10), upstream Cf 1.67e-3 (1.44-1.56e-3).
* **Mach 8.2 plate:** SA-Edwards lowers Cf and q_w by ~8 % relative to SA;
  **SA-comp lowers them by ~42 %**, far outside the data.
* **Sharp fin:** pressures are nearly model independent and agree with the data
  (plate plateau 2.4 vs 2.46, fin 10.4 vs 10.1). SA reproduces the SA curves of
  Edwards & Chandra (peak plate heating 6.6 vs 6.6, peak Cf 5.2e-3 vs 5.2e-3);
  SA-Edwards gives higher peak plate heating (7.7; measured 5.0). All models
  overpredict the fin heating (~7.8 vs 3.1), as Edwards & Chandra report for all
  of theirs (Boussinesq relation / laminar fin boundary layer in the experiment).
* **Cylinder / flare:** separation and pressure-plateau locations are captured;
  as in Edwards & Chandra, the phi = 180 deg separation is stronger than measured
  (Cf min -1.0e-3 vs -0.43e-3) and the phi = 90 deg recovery is overpredicted. At
  Mach 2.89 with an adiabatic wall the modifications change little.

## Folder layout

| Path | Contents |
|---|---|
| `README.md` | this report |
| `results/summary.md` | validation tables (medium grids) |
| `results/grid_study.md` | grid-refinement tables |
| `results/su2_comparison.md`, `.json` | ADflow vs SU2 tables |
| `results/figures/` | all plots |
| `results/<code>/<case>/<level>/<model>/` | solutions (`code` = adflow \| su2, `case` = plate \| fin \| flare) |
| `run_all.sh` | regenerates everything (`./run_all.sh grids plate fin flare su2 post`) |
| `run.sh` | `./run.sh NPROC script args`: runs with this ADflow build |
| `meshes.py` | grid generation (flat plate, sharp fin, cylinder/flare; `ref` for grid families) |
| `common.py` | conditions, solver options, `MODELS`, results paths, boundary-layer initial condition, RK/ANK strategy, readers |
| `run_flatplate.py`, `run_case.py`, `postprocess.py`, `dump_restart.py` | running and post-processing ADflow cases |
| `summarize.py`, `grid_study.py`, `plot_results.py` | tables and plots in `results/` |
| `data/experiments.json` | digitized experimental data (sources and uncertainty inside) |
| `su2/` | SU2 case generator and runner, `compare.py`, `profiles.py`, `wallgrad.py` |
| `grids/`, `logs/` | generated grids and run logs |

Solution folders are not versions of the same result: each is one model on
one grid. `level` is

* `medium` - the reported result (all validation tables and figures);
* `coarse`, `fine` - used only by the grid-refinement study and the SU2
  comparison (fine exists for the 2-D plate only);
* `sensitivity` - checks cited in Sec. 1.2 (`SA_nu3`: nuTilde_inf = 3 nu_inf,
  `SA_vort`: vorticity-based production).

`grids/`, `logs/`, `results/adflow/` and `results/su2/` are large and not
tracked by git.

## Solver notes

* From a uniform freestream the Mach 8.2 cases do not get through the start-up
  transient with ADflow's ANK solver. The flow is initialized with an
  approximate turbulent boundary layer (`init_boundary_layer`) and, at Mach 8.2,
  1500 explicit Runge-Kutta iterations are run before ANK (`STAGE1_RK`).
* On the fin the ANK coupled turbulence solve is switched on early
  (`ANKCoupledSwitchTol` 1e-3); on the coarse flare grid it must stay decoupled
  (1e-16). Restarting a different model from a converged solution sometimes
  stalls in ANK; those runs were repeated from the boundary-layer initial condition.
* The fin case with SA-comp could not be converged and is not reported; the
  effect of SA-comp is shown by the 2-D plate.

## Known limitations

* Experimental data are digitized from published figures (NASA TM-103838 for the
  fin is not available online); expect +-0.3 cm / +-3-5 %.
* Three-level grid study for the 2-D plate (including the SU2 comparison); two
  levels for the fin and the flare (their fine grids were not run to completion).
* Wind-tunnel walls, the Wideman sting/nose shape and the fin height are not
  modelled.
