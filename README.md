# Hydrate Formation Conditions Web App

© 2026 [Qazi Nasir](https://qazinasir.com) — https://qazinasir.com

Interactive, self-contained implementation of:

> Moradi & Khosravani, *Modeling of hydrate formation conditions for CH4, C2H6, C3H8, N2, CO2
> and their mixtures using the PRSV2 equation of state and obtaining the Kihara potential
> parameters*, Fluid Phase Equilibria 338 (2013) 179–187.

## Quick start

Open **`index.html`** in any modern browser — no server, no network, no dependencies required.
(Everything — CSS, the JavaScript thermodynamic engine, the app logic and the ~70 experimental
anchor points — is embedded in that single file.)

## What the app does

- **P(T) dissociation curves** for pure gases and arbitrary 5-gas mixtures (H–I–V, H–LW–V,
  H–LW–LC), with automatic sI/sII selection (lowest equilibrium pressure wins).
- **All-gas overview** mode: the five pure-component curves on one chart with anchors.
- **P(x) mode**: structure-transition scan along a binary composition axis (paper Fig. 5).
- **State details**: Langmuir constants, cage occupancies, water/gas fugacities, hydration
  number at a chosen (T, P).
- **Paper tab**: reproductions of the paper's Fig. 1 (gas solubility in water), Fig. 5
  (CH4+C3H8 sI→sII transition), Fig. 7 (natural-gas mixture), Tables 6–8 (occupancies).
- **Compare**: model vs. experimental anchors with AADP, deviation-vs-T and log–log parity
  panels (±10/±25% bands).
- **Params / Validation tabs**: every calibrated constant, and an honest list of the paper's
  internal inconsistencies and this implementation's known limitations.

## Model

- van der Waals–Platteeuw hydrate model with McKoy–Sinanoğlu cell potential and
  Kihara (spherical-core) guest parameters.
- PRSV2 equation of state (Stryjek–Vera κ1κ2κ3 from the DWSIM tables) with the paper's
  van Laar-type two-parameter kij mixing rule and temperature-dependent kij (Eq. 10, Table 1).
- Water reference fugacity: Martín–Peters fβ(T) form, A/B calibrated here (see below).
- **Stage-2 aqueous correction (this work):** k(g,w) → k(g,w) + dK0 + dK1·(T−273.15 K),
  applied only in the water-rich phase. Needed because the paper's printed kij give
  CO2 solubility ~1e-17 (negative-a artifact) while the paper's own Fig. 1 shows ~1e-2.

## Repository layout

```
index.html          the complete app (single file, open in browser)
dev/
  engine.py         Python research prototype (reference implementation)
  engine_port.js    line-by-line JavaScript port (verified to ~1e-10 vs Python)
  app.js / app.css  UI logic and styling
  build_html.py     inlines engine+css+js+anchors into ../index.html
  calibrate.py      experimental anchor set (~70 points, Sloan & Koh / Zenodo compilations)
  calib4.py         cached fugacity machinery used by all fits
  fitv13.py         stage-2 objective: root-space residuals + solubility anchors (dK-only)
  gnfit13.py        damped Gauss–Newton driver that produced the final v8 parameters
  scan_v8.py        guard-first parameter scan (CH4/C2H6 dK0 selection)
  val8.py           full validation harness (per-gas AADP, structures, occupancies, Fig.5, NG)
  calib10.npy / params_v7.json   frozen stage-1 parameters (Kihara + fβ)
  calib18.npy / params_v8.json   FINAL parameters (v7 + aqueous dK corrections)
  patch_v8.py       injects params_v8 into engine_port.js and regenerates JS test refs
  test_port.js      Node parity tests (Python vs JS engine) — run: node test_port.js
  domtest.js        headless DOM smoke test of index.html — run: node domtest.js
  CALIBRATION.md    full calibration history: method, rejected approaches, final numbers
  fitv3–12, gnfit, calib2–9, params_v1–13 …  earlier (superseded) iterations, kept for provenance
```

## Rebuilding / testing

```
cd dev
python3 build_html.py     # regenerate ../index.html from sources
node test_port.js         # JS↔Python parity (needs /tmp/jsref.json; patch_v8.py regenerates it)
node domtest.js           # headless UI smoke test
python3 val8.py calib18.npy   # validation report for the final parameters
```

## Accuracy (v8, vs. ~70 experimental anchors)

| gas | AADP this work | AADP paper (Table 5) |
|---|---|---|
| CH4 | 7.6% | 1.37% |
| C2H6 | 72% (3 LC points dominate; V/ice avg ≈7%) | 2.47% |
| C3H8 | 4.7% | — |
| N2 | 39% (281–293 K dominate; ≤277 K ≤19%) | 6.84% |
| CO2 | 30% | 2.89% |

Gas solubilities (paper Fig. 1 anchors): model/exp = 0.83–1.02 for C2H6, C3H8, N2, CO2;
0.32–0.45 for CH4 (deliberate trade — see Validation tab). Known limitations (C2H6 LC,
N2 high-T) are documented honestly in the app's Validation tab and in `dev/CALIBRATION.md`.
