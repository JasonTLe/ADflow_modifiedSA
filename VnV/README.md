# V&V of the SA modifications in this ADflow branch

Verification and validation of the Spalart-Allmaras modifications
implemented in `src/turbulence/saCorrections.F90` (branch `sa-comp-edwards`):

* **SA-Edwards** (`turbulenceModel = "SA-Edwards"`) - Edwards & Chandra,
  AIAA 94-2275 (1994), Sec. 4.3, Eqs. 20, 24, 25 (SA-noft2-Edwards on the
  NASA Turbulence Modeling Resource, TMR).
* **SA-comp** (`useCompressibilitySA = True`, constant `SAc5`) - Spalart,
  AIAA 2000-2306: -C5 nuTilde^2/a^2 (du_i/dx_j)^2, C5 = 3.5.
* **Density corrections of the diffusion term** (`densityCorrectionSA`):
  * `"Catris-Aupoix"` - Catris & Aupoix, Aerosp. Sci. Technol. 4(1), 2000
    (TMR: SA-noft2-Catris): diffused quantity sqrt(rho) nuTilde, laminar
    diffusion of mu, and the original non-conservative left hand side
    D(rho nuTilde)/Dt (which adds nuTilde div(u) to ADflow's equation for nuTilde);
  * `"Catris-Aupoix-cons"` - the same diffusion term with the conservative
    left hand side (TMR: SA-noft2-CatrisCons, Waligura et al. AIAA 2022-0587);
  * `"conservative"`, `"Catris-1999"`, `"Seror"` - the density-gradient cross
    terms that turn ADflow's non-conservative diffusion term into the
    conservative form of the TMR / SU2, the form of Catris (1999) and the form
    of Seror et al. (J. Aircraft 42(1), 2005); see Xue et al., AIAA J. 62(1),
    2024, Eqs. 5-6 and 10.

  All of them can be combined with `"SA"` or `"SA-Edwards"` and with SA-comp.
  Model keys used in this folder (`common.MODELS`): `SA`, `SA-Edwards`,
  `SA-comp`, `SA-Catris`, `SA-Catris-cons`, `SA-Catris1999`, `SA-Seror`,
  `SA-cons`, `SA-Edwards-Catris` (all `SA*` variants without ft2).

Everything needed to reproduce the evidence below is in this folder
(`run_all.sh`), except the unit/regression tests, which live with ADflow's
tests in `tests/reg_tests/test_sa_corrections.py`.

## Summary of the evidence

| Check | Result |
|---|---|
| Equations vs. the papers, the NASA Turbulence Modeling Resource and the SU2 source | Match (SU2 differs from the paper in one strain term; this code follows the paper). Catris-Aupoix transcribed from the TMR equation and checked against the expanded form of Waligura et al. |
| Independent code review | No errors; findings addressed |
| Unit/regression tests (`test_sa_corrections.py`, 9 model combinations) | 39 pass, 15 skip (N/A) under testflo's process-per-test mode (a single-process runner gives 44 pass / 10 skip because `setUp` leaks options between the parameterized classes; see `AUDIT_2026-10-08.md`): blockette vs block residual 1e-12, forward AD vs FD (1st-order FD convergence), forward/reverse/reverse-fast dot products 1e-11, C5 linearity, density corrections change only the turbulence residual, density-corrected diffusion = standard SA at constant density (1e-10) |
| Full ADflow regression test suite (`tests/reg_tests`, with the new build) | 181 pass, 0 fail, 15 skip (the N/A combinations above) |
| **Code-to-code (ADflow vs SU2), SA-Edwards increment** | Same effect in both codes on all three grids: Cf within 0.1 pt, q_w within 0.2 pt, theta within 0.9 pt (SU2 converged to rms[rho] 1e-11; rerun 2026-10-09) |
| Grid refinement, 2-D plate (3 levels) | Grid independent: Cf, q_w <= 0.2 % change, GCI <= 0.2 % |
| Grid refinement, fin / flare (2 levels) | <= 1.5 % change (except the flare separation-bubble Cf minimum, ~7 %) |
| Validation vs. Kussoy & Horstman (M 8.2 sharp fin) | Pressures within a few %; heating overpredicted by all models, as in the paper. NB the paper's 'Spalart-Allmaras' is its Edwards-modified SA (its Sec. 4.3): like for like, ADflow SA-Edwards gives peak plate heating 7.7 vs 6.6 and peak Cf 5.6e-3 vs 5.2e-3, with an incoming theta 11 % thinner than the paper's matched inflow (see `AUDIT_2026-10-08.md`) |
| Validation vs. Wideman et al. (M 2.89 cylinder / flare) | Separation / pressure plateau locations captured; same model deficiencies as in the paper |
| SA-comp in hypersonic wall-bounded flow | Cf and q_w 42 % too low on the Mach 8 plate - do not use for such flows (consistent with Xue et al., AIAA J. 2024) |
| Density corrections on the Mach 8.2 plate | Catris-Aupoix (original form): Cf, q_w +6 % vs SA, theta +7 %; conservative-form Catris-Aupoix, Catris-1999, Seror and conservative: within 1 % of SA. The incoming boundary-layer profiles of all of them match the survey; the wall heat flux stays ~10 % above the measurement for all models except SA-comp |
| Density corrections on the 3-D cases (Catris-Aupoix, both forms) | Fin: original form lowers peak plate heating 6.6 -> 5.9 (measured 5.0) and peak Cf 5.2 -> 4.9e-3 (5.9 +- 1.3e-3), fin heating 7.8 -> 7.6 (3.0); conservative form within 2 % of SA. Flare (adiabatic, M 2.89): conservative form within 1 % of SA; original form moves separation 0.4 cm upstream and raises Cf behind the flare expansion by ~25 % (away from the data) because its nuTilde div(u) term is a source in expansions |

All documents are in `results/` (`summary.md`, `grid_study.md`,
`su2_comparison.md`) and all figures in `results/figures/`:

* `plate_vs_experiment.png` - boundary-layer profiles (u, T, M) at x = -11 cm and theta, Cf, q_w vs the Kussoy & Horstman survey, all models
* `plate_ADflow_vs_SU2.png` - wall Cf and q_w along the plate, ADflow vs SU2, three grids
* `plate_profiles_ADflow_vs_SU2.png` - boundary-layer profiles at x = -11 cm, ADflow vs SU2 (superseded by the SU2 overlays in `by_model/plate_SA-Edwards.png`; kept because `su2/profiles.py` also prints the near-wall numbers quoted in Sec. 1.2)
* `fin_vs_experiment.png` - Kussoy & Horstman sharp fin
* `flare_vs_experiment.png` - Wideman et al. cylinder / flare
* `by_model/<case>_<model>.png` - one figure per modification (plate, fin, flare): the modification (red, dashed) vs SA (green) vs experiment (black circles); on the plate SA-Edwards figure the converged SU2 solutions are overlaid dashed (SU2 SA purple, SU2 SA-Edwards blue); `plate_summary_bars.png` - theta, Cf, q_w of all models as ratios to the experiment, with SU2's SA and SA-Edwards values as hatched bars; `VnV_by_model_all.pdf` - all of them plus the along-the-plate ADflow-vs-SU2 figure in one file (17 pages)

## 1. Verification

### 1.1 Implementation checks

* **Equations** read from Edwards & Chandra (paper in `../research_docs/`),
  Spalart (2000) via the NASA Turbulence Modeling Resource, and compared term
  by term with SU2's implementation (`SA_OPTIONS= EDWARDS` / `COMPRESSIBILITY`).
* **Regression tests** (`../tests/reg_tests/test_sa_corrections.py`, MDO tutorial
  RANS wing; SA-Edwards, SA+comp, SA-Edwards+comp, the five density
  corrections and SA-Edwards+Catris+comp): consistency of the two
  residual paths, forward AD vs finite differences, forward vs reverse vs
  reverse-fast AD (dot-product tests), and that the SA-comp term changes only
  the turbulence residual, linearly in C5; that the density corrections
  change only the turbulence residual; and that the diffusion terms of
  `Catris-Aupoix-cons`, `Catris-1999` and `conservative` reduce to the
  standard SA term when the density is set to a constant (relative
  difference < 1e-10). The last test caught a difference in the clipping of
  the face diffusion coefficients at wall cells during development. The
  adjoint code for the new routines was generated with Tapenade 3.16
  (`src/adjoint/output*/saCorrections_*.f90`).
* **Form of the equation.** ADflow solves the SA equation for nuTilde in
  non-conservative form with the local kinematic viscosity. Multiplied by rho,
  this is exactly the conservative form of Allmaras, Johnson & Spalart (2012)
  including its cross term -(1/sigma)(nu + nuTilde) grad(rho).grad(nuTilde),
  which the TMR considers the correct compressible form (its Note 4 warns
  against moving rho inside the derivatives without that term). The
  `conservative` option therefore *adds* the cross term to recover the
  rho-inside-the-derivative form used by the TMR tables and SU2; on the Mach
  8.2 plate the two differ by 0.3 % in Cf.
* **Log-layer check** of the Edwards form: Stilde and r match standard SA where
  standard SA is well behaved; Stilde stays positive where standard SA's
  Stilde goes negative (chi ~ 5).

### 1.2 Code-to-code: ADflow vs SU2 (2-D Mach 8.2 cold-wall plate)

`su2/` is the same flat plate for SU2 v8.5.0 (built in `~/packages/SU2`), on
exactly the same grids, with the physics matched to ADflow (SA without ft2 /
Edwards option, Sutherland 1.716e-5 / 273.15 K / 110.55 K, Pr 0.72, Pr_t 0.9,
nuTilde_inf = 1.342 nu_inf, T_w = 300 K; freestream state agrees to 0.001 %).
SU2 is run first order (CFL 5, 8000 iterations, initial guess only), then
second order (Roe, MUSCL + Venkatakrishnan limiter frozen after 500 iterations,
adaptive CFL 2 -> 50) restarted from it, to rms[rho] = 1e-11 on every grid
(stage 2: 640 / 760 / 2400-2700 iterations on coarse / medium / fine). This
recipe dates from 2026-10-09: the earlier runs (MUSCL + van Albada edge
limiter, fixed CFL 2) stalled at rms[rho] 1e-4 (coarse) to 1e-5 (medium/fine)
because of a persistent residual in the few nodes at the plate leading edge
(slip -> no-slip, 77 K -> 300 K at Mach 8.2); freezing that limiter or
tightening the linear solver did not help, the Venkatakrishnan limiter (first
order at extrema) does. Those earlier results are kept in
`results/su2_v1_vanalbada/`; converging the fine-grid solutions moved its
Cf / q_w at x = -11 cm by -0.4 to -1.0 % and the Edwards increments by
0.2-0.4 pt. On the fine grid the wall quantities change by < 0.02 % between
the last two snapshots (720 / 414 iterations apart).

**SA-Edwards relative to SA at x = -11 cm** (the direct test of the Edwards
implementation, independent of each code's baseline discretization):

| grid | Cf ADflow / SU2 | q_w ADflow / SU2 | theta ADflow / SU2 |
|---|---|---|---|
| coarse (7k cells) | -7.49 / -7.59 % | -7.66 / -7.88 % | -11.6 / -12.2 % |
| medium (27k) | -7.52 / -7.62 % | -7.69 / -7.87 % | -11.4 / -12.1 % |
| fine (109k) | -7.49 / -7.61 % | -7.66 / -7.83 % | -11.4 / -12.3 % |

Code-to-code difference of the increment (SU2 - ADflow): Cf 0.10 / 0.10 / 0.12
pt, q_w 0.22 / 0.18 / 0.17 pt, theta 0.6 / 0.7 / 0.9 pt on coarse / medium /
fine. (theta is compared at the nearest ADflow cell column, 2-5 mm apart
between grid levels; theta grows ~0.5 %/cm here.)

**Baseline (SA) difference between the codes.** SU2 gives lower Cf/q_w than
ADflow on the same grid (x = -11 cm: Cf -7.2 / -5.4 / -4.6 %, q_w -5.0 / -4.1 /
-3.8 % on coarse / medium / fine). The difference is largest near the plate
leading edge, decays downstream, and decreases with grid refinement while
ADflow is already grid converged (Cf 0.954 / 0.952 / 0.954e-3). SU2's Cf
sequence 0.885 / 0.901 / 0.910e-3 is monotonic (observed order ~0.8) but
extrapolates to ~0.92e-3, ~3 % below ADflow, so the two baselines do not
converge to the same solution: the codes differ in the SA diffusion form,
limiter and node- vs cell-centred discretization, which is why the Edwards
implementation is compared on the increment and not on the absolute values.
It is not the freestream
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
| 2-D plate | M 8.2, Re 5.32e6/m, T 76.89 K, T_w 300 K | 171 x 80 x 2 cells | Kussoy & Horstman inflow survey, x = -11 cm (KUSSOY3.DAT) |
| Sharp fin (15 deg), Kussoy & Horstman | same, T_w/T 3.9 | 0.97 M cells | KUSSOY3.DAT (NASA TM-103838) as listed in NASA CR-177638 |
| Cylinder / offset 20 deg flare, Wideman et al. | M 2.89, p 5.54 kPa, T 105 K, adiabatic | 0.63 M cells | NASA TM-108824, Tables 1-37 |

Full tables: `results/summary.md` (medium grids). The experimental data are
the tabulated NASA data (`data/experiments_nasa.json`, see `data/README.md`);
the hand-digitized data of the first round (`data/experiments.json`) are
superseded. Note that the free-stream conditions of the fin data file (M 8.18,
T 81 K, Re 4.9e6/m, Cf 0.98e-3, q_w 0.99 W/cm^2) differ slightly from the
conditions of Edwards & Chandra used for the runs (table above); the plate
leading edge is placed 1.64 m upstream of the fin apex (fully turbulent plate;
physical distance 1.76 m, with laminar/transitional flow near the leading
edge in the experiment) so that theta at x = -11 cm matches the survey. The
Reynolds number of Wideman's Table 1 (1.5e7/m) is used; Edwards & Chandra
list 1.7e7/m.

* **Incoming boundary layers** (SA): Kussoy x = -11 cm: theta 0.094 cm (exp. 0.094),
  Cf 0.95e-3 (0.98e-3), q_w 1.10 W/cm^2 (0.99); the u, T and M profiles follow
  the survey (`plate_vs_experiment.png`). Wideman x = -2 cm: delta99 1.14 cm
  (1.10), upstream Cf 1.67e-3 (1.44-1.56e-3).
* **Mach 8.2 plate:** SA-Edwards lowers Cf and q_w by ~8 % relative to SA;
  **SA-comp lowers them by ~42 %**, far outside the data. The original
  Catris-Aupoix model raises Cf and q_w by 6 % and theta by 7 % (the nuTilde
  div(u) term of its non-conservative left hand side acts as a source in the
  boundary layer); its conservative form, Catris-1999, Seror and the
  conservative form of standard SA are all within 1 % of SA. None of the
  density corrections removes the ~10 % overprediction of the wall heat flux.
* **Sharp fin:** pressures are nearly model independent and agree with the data
  (plate plateau 2.4 vs 2.5, fin 10.4 vs 10.2). The paper's only
  Spalart-Allmaras model is its Edwards-modified form (its Sec. 4.3); its
  curves (peak plate heating 6.6, peak Cf 5.2e-3, digitised from Figs. 2-6,
  +-0.1-0.2) sit between ADflow's SA (6.6, 5.2e-3) and SA-Edwards (7.7,
  5.6e-3). The paper matched the incoming theta for each model, whereas one
  virtual origin is used for all models here, so SA-Edwards arrives with theta
  11 % thinner (see `AUDIT_2026-10-08.md`). SA-Edwards gives higher peak plate
  heating (7.7; measured 5.0). The measured
  peak Cf (5.9e-3 +- 1.3e-3) is matched within its uncertainty. All models
  overpredict the fin heating (~7.8 vs 3.0), as Edwards & Chandra report for all
  of theirs (Boussinesq relation / laminar fin boundary layer in the experiment).
  The original Catris-Aupoix model reduces the peak plate heating and skin
  friction (5.9 and 4.9e-3, vs 6.6 and 5.2e-3 for SA) and the fin heating
  slightly (7.6), i.e. it moves towards the measured heating; its conservative
  form stays within 2 % of SA.
* **Cylinder / flare:** separation and pressure-plateau locations are captured
  (plateau 2.75 vs 2.62-2.68 measured); as in Edwards & Chandra, the phi = 180
  deg separation is stronger than measured (Cf min -1.0e-3 vs -0.43e-3) and the
  phi = 90 deg recovery is overpredicted. The latter is not an ADflow problem:
  the SA curve of Edwards & Chandra (their Fig. 10) reaches the same 3.5-4e-3
  at x = 10-14 cm, against 2.2e-3 measured, and they attribute it to the SA
  destruction term, which only sees wall-normal variations. At Mach 2.89 with
  an adiabatic wall the modifications change little.
  The conservative Catris-Aupoix form is within 1 % of SA everywhere. The
  original form separates 0.4 cm earlier at phi = 180 deg (x_sep 4.5 vs 4.9 cm,
  measured 5.2) and gives 25 % higher Cf on the afterbody behind the flare
  expansion (x > 10 cm), where the data lie below SA already: the nuTilde div(u)
  term of its non-conservative left hand side is a source of nuTilde in
  expansions (and a sink in shocks). For adiabatic supersonic flow the
  conservative form is the safer choice; on the cold-wall Mach 8.2 cases the
  original form is closer to the measured heating.

## Folder layout

| Path | Contents |
|---|---|
| `README.md` | this report |
| `results/summary.md` | validation tables (medium grids) |
| `results/grid_study.md` | grid-refinement tables |
| `results/su2_comparison.md`, `.json` | ADflow vs SU2 tables |
| `results/figures/` | all plots |
| `results/<code>/<case>/<level>/<model>/` | solutions (`code` = adflow \| su2, `case` = plate \| fin \| flare) |
| `run_all.sh` | regenerates everything (`./run_all.sh grids plate fin flare plate_dens fin_dens flare_dens su2 post`) |
| `run.sh` | `./run.sh NPROC script args`: runs with this ADflow build |
| `meshes.py` | grid generation (flat plate, sharp fin, cylinder/flare; `ref` for grid families) |
| `common.py` | conditions, solver options, `MODELS`, results paths, boundary-layer initial condition, RK/ANK strategy, readers |
| `run_flatplate.py`, `run_case.py`, `postprocess.py`, `dump_restart.py` | running and post-processing ADflow cases |
| `summarize.py`, `grid_study.py`, `plot_results.py` | tables and plots in `results/` |
| `plot_by_model.py` | per-modification comparison figures and the all-model bar chart (`results/figures/by_model/`, `plate_summary_bars.png`, `VnV_by_model_all.pdf`) |
| `data/experiments_nasa.json` | tabulated NASA experimental data (`data/README.md`); `experiments.json`: superseded digitized data |
| `su2/` | SU2 case generator and runner, `compare.py`, `profiles.py`, `wallgrad.py` |
| `grids/`, `logs/` | generated grids and run logs (not tracked) |
| `archive/` | superseded runs, not tracked (`su2_v1_vanalbada/`: the SU2 runs with the van Albada limiter that stalled at rms[rho] 1e-4..1e-5, see Sec. 1.2) |
| `AUDIT_2026-10-08.md` | independent audit of this V&V (equations, tests, data, SU2 comparison, convergence) with the 2026-10-09 SU2-convergence addendum |

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

* The fin and flare computations use the free-stream conditions of Edwards &
  Chandra, which differ by a few % from those listed with the data (see above);
  the plate leading edge is a fully turbulent virtual origin.
* Fin heat transfer is compared at the thermocouple station S = 16.70 cm;
  the Schmidt-Boelter data at S = 19.62 cm are shown for reference only.
* Three-level grid study for the 2-D plate (including the SU2 comparison); two
  levels for the fin and the flare (their fine grids were not run to completion).
* Wind-tunnel walls, the Wideman sting/nose shape and the fin height are not
  modelled.
* One virtual origin (plate leading edge at -1.64 m, tuned so that standard SA
  gives theta = 0.094 cm at x = -11 cm) is used for every model, whereas Edwards
  & Chandra matched the inflow per model; SA-Edwards therefore reaches the fin
  with theta 11 % thinner (Catris-Aupoix 7 % thicker) than the experiment.
* With two grid levels the fin and flare tables carry an unquantified grid
  uncertainty: assuming second order, the medium-grid GCI is ~1.3-1.8 % for the
  fin peak heating / Cf and ~9 % for the flare Cf minimum, larger than some of
  the model-to-model differences in the tables.
