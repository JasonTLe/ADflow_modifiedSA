# Experimental data

## `experiments_nasa.json` (used by all scripts)

Tabulated data transcribed from the NASA reports in `../../research_docs/`
(every data set carries a `_source` entry with report, table and page):

* **Kussoy & Horstman, Mach 8.2 flat plate / 15 deg sharp fin** (NASA
  TM-103838, 1991). The TM itself is not on NTRS; its data file KUSSOY3.DAT
  (including corrected skin-friction data "supplied by the experimenters") is
  listed in Settles & Dodson, *Hypersonic Shock/Boundary-Layer Interaction
  Database: New and Corrected Data*, NASA CR-177638 (1994), pp. 25-28
  (rotated page images in
  `../../research_docs/Settles_Dodson_1994_KUSSOY3_fin_data_pages/`). Contents:
  free-stream conditions, upstream boundary-layer profile, plate pressure and
  heat transfer at X = 18.19 / 16.45 cm (two sets each), plate skin friction at
  X = 15.5 cm with uncertainties, fin pressure (S = 18.34 cm) and heat transfer
  (thermocouples S = 16.70 cm, Schmidt-Boelter gauges S = 19.62 cm), surface
  streamline angles; the same for the 5, 7.5, 10 and 12.5 deg fins. Z is the
  distance from the fin surface in the plate plane; S the distance along the
  fin from its leading edge (X = S cos 15 deg is used in the figure titles, as
  in Edwards & Chandra).

  Note that the free-stream conditions listed in the data file (M 8.18, T 81 K,
  p 430 Pa, Re 4.9e6/m, Cf 0.98e-3, q_w 0.99 W/cm^2) differ slightly from the
  conditions tabulated by Edwards & Chandra (AIAA 94-2275: M 8.2, T 76.89 K,
  Re 5.32e6/m, Cf 1.0e-3, q_w 1.04 W/cm^2), which are the ones run here so that
  the computations can also be compared with their results.
* **Wideman, Brown, Miles & Ozcan, Mach 2.89 cylinder / offset 20 deg flare**:
  NASA TM-108824, *Surface Documentation of a 3-D Supersonic
  Shock-Wave/Boundary-Layer Interaction* (1994), Tables 1-37 (the complete
  data report; the AIAA 94-0314 / NASA CR-202617 paper only has figures).
  Contents: tunnel conditions, junction locations, oil-flow separation /
  attachment lines, LDV profile of the undisturbed boundary layer, LISF skin
  friction (streamwise along phi = 0, 90, 180 deg; transverse at 90 deg;
  undisturbed boundary layer) with the uncertainty analysis, surface pressure
  along phi = 0 ... 180 deg on the flare and along the cylinder side. P/Pinf
  was computed from the tabulated P/Pt with the measured Pinf/Pt = 0.0332 of
  the report. x = 0 is the innermost (phi = 0) cylinder/flare junction.

## `experiments.json` (superseded)

The data digitized by hand from the figures of Edwards & Chandra (AIAA
94-2275) and Wideman et al. (AIAA 94-0314) that were used for the first round
of results. Kept only for the approximate readings of Edwards & Chandra's own
SA curves (`paper_SA_curves_approx`), which `summarize.py` still reports.
Compared with the tabulated data, the digitized values agree within ~1-4 %
where the curves are flat, but show position shifts of 0.1-0.5 cm in steep
gradients and were 6-10 % low in the undisturbed plate heat transfer and
pressure.
