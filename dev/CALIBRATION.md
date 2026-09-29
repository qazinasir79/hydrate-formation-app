# Hydrate Model Calibration Report
## Moradi & Khosravani (FPE 338 (2013) 179–187) — PRSV2 + vdW–P + Kihara

## 1. Model (implemented per paper)
- PRSV2 EOS (Stryjek–Vera κ1κ2κ3 from DWSIM's genuine prsv2.dat table), van Laar-type
  two-binary-parameter mixing rule: a_ij = √(a_i a_j)(1 − k_ij k_ji/(x_i k_ij + x_j k_ji)),
  k_ij(T) = p_ij T + q_ij with the paper's water–gas parameters (assigned water_first:
  in water-rich liquid the (gas,water) direction dominates; in gas-rich phase the
  (water,gas) direction dominates).
- van der Waals–Platteeuw: f_w^H = f_w^β · exp(Σ ν_m ln(1−θ_m)), θ = Cf/(1+ΣCf),
  McKoy–Sinanoğlu cell potential w(r) with Kihara (ε, σ, a), 200-node Gauss–Legendre.
- f_w^β = exp(A − B/T) · exp(V_β(T,P)(P − f_β)/RT); V_β per paper Eq. (26)/(27) incl. P terms.
- H–I–V below 273.15 K via exact ice bridge: f_β,I = f_β,LW · f_ice/f_liq,
  psat_ice = exp(28.89 − 6134.9/T) Pa, V_ice = 19.65 cm³/mol.
- Structure selection: lowest equilibrium pressure among sI/sII.
- Flash: successive substitution w/ Rachford–Rice; zero-feed components kept exactly zero
  (an earlier 1e-16 composition floor amplified trace lnφ_L of absent gases — e.g.
  lnφ_L(CO2)≈+39 in water — into phantom vapor fractions; fixed 2026-09).

## 2. FINAL CALIBRATED PARAMETERS (calib10.npy / params_v7.json)
Kihara (ε/k [K], σ [Å], a [Å]):
  CH4   161.05  3.1595  0.3317
  C2H6  180.46  3.1638  0.4799
  C3H8  210.41  3.4554  0.6502
  N2    151.38  2.9596  0.3410
  CO2   172.37  2.9973  0.6786
Empty-lattice fugacity f_β = exp(A − B/T) [Pa]:
  sI:  A = 28.7933, B = 5953.7   → f_β(273.15) = 1093.1 Pa
  sII: A = 28.6056, B = 5922.6   → f_β(273.15) = 1015.5 Pa
Fit: ~70 literature anchors (Sloan & Koh / Zenodo CO2-ice data), slope-scaled residuals
(r/|dr/dlnP| ≈ ln P-error), per-anchor structure hinges, occupancy anchors (Tables 6–8),
weak P–P Langmuir regularization, Fig-5 transition state. Objective ≈ paper Eq. (27) but
scaled so flat-slope (saturated-cage) regions carry their true Peq sensitivity.

## 3. Paper internal inconsistencies found (documented per user instruction)
1. **Table 2 Kihara parameters do not reproduce the paper's own results.** Used verbatim
   in a standard McKoy–Sinanoğlu implementation with Eqs. (16)-(17) f_β they give:
   θ_S(CH4) ≈ 0.01 vs Table 6's 0.935; Peq(CH4) off by ~×100; C3H8/CO2 Langmuir constants
   vanish (σ too large vs cage radius). Their printed values cannot be what their code used.
2. **Table 6 occupancies + Martín–Peters f_β (Eqs. 16-17) are mutually inconsistent**:
   θ = (0.969, 0.935) at 273.2 K with f_β,sI(273.15) = 1295 Pa gives f_w^H = 748 Pa vs
   f_w^L ≈ 600 Pa (25% off-equilibrium). f_β,sI(273.15) ≈ 1050–1100 Pa (our fitted value)
   reconciles Table 6 with H–LW–V data. Data-driven fits land at (28.4–28.8, 5800–5960) —
   close to but distinguishable from Eqs. (16)-(17).
3. Table 7 caption duplicates Table 6's (says CH4; the 276.7 K data are indeed CH4).

## 4. Validation (this work vs anchors; paper AADP in parentheses)
Per-region AADP of dissociation pressure:
- CH4  (1.37%): H-LW-V 273–293: 0.2–12% (mean ~5%); ends: 263.15 ice +18.8%, 298.15 +19.3%.
- C2H6 (2.47%; LC 11.15%): H-I/LW-V 263–287.75: 3–23% (mean ~7%); H-LW-LC 293–303:
  +180…+408% — NOT REPRODUCIBLE with shared-B f_β: in the LC region |dr/dlnP| ≈ 0.02,
  and the residual slope dr/dT ≈ +0.01/K (from B/T² − dlnf_LW/dT) cannot be flattened
  without tilting CH4/CO2 V-region curves (tested: B1 shifts, D/T² curvature, per-gas
  a-core, weight boosts — all trade one region for another; the tension is structural).
- C3H8 (3.83%): 4.75% overall (ice 268/271 ±15%, V 273–278 ±5%).
- CO2  (2.89%; LC 10.2%): ice 217–253 (Zenodo): ±6%; H-LW-V 263–283: ±4.5% (paper-level);
  Q2 282.65–285.15: −15…−24%; H-LW-LC 287–291: +131…+271% (same structural tension).
- N2   (6.84%): 263–277: −6…−19%; 281–293: +44…+77% (same tension; sII structure ✓).
- Structure selection: CH4 sI ✓, C2H6 sI ✓, C3H8 sII ✓, CO2 sI ✓ (incl. ice branch),
  N2 sII ✓, natural gas sII ✓, CH4–CO2 sI ✓.
- Occupancies (Tables 6–8): CH4@273.2: 0.9737/0.9232 (paper 0.969/0.935);
  CH4@276.7: 0.9774/0.9330 (0.981/0.917); CO2@273.2: 0.9849/0.5860 (0.986/0.581) — ≤1.3%.
- Fig. 5 (CH4+C3H8, 277.6 K): sI→sII onset at y_C3H8 ≈ 3e-4, P ≈ 4.0 MPa
  (paper: y_CH4 = 0.9995 at 4.0 MPa) ✓.
- Fig. 6 (CH4+C2H6, 277.6 K): our model gives sI over the whole composition range
  (smooth P–x curve); the paper reports transitions at y_CH4 = 0.753/0.724 (an sII lens).
  Not reproduced — with physically calibrated C(T) the sII branch never dips below sI here.
- Fig. 7 (natural gas): sII curve, e.g. 1.09 MPa @277.7 K → 9.85 MPa @296.7 K,
  qualitatively matching the experimental trend shown in the paper.

## 5. Known limitations to surface in the app
- H-LW-LC equilibria (T above the Q2 point of C2H6/CO2) carry large errors (documented above);
  the app should flag results in that region.
- CH4 above ~295 K and ice-branch endpoints (~±19%) are the weakest CH4 points.
- Gas solubility in water is underpredicted (~×10 for CH4 at 2.6 MPa: x ≈ 1e-4 vs 1.1e-3)
  because the paper's water–gas van Laar k's give very large infinite-dilution lnφ; this
  barely affects f_w^L (x_w ≈ 0.9999) and hence equilibrium, but absolute solubilities
  shown in the app are indicative only.

## v8 — aqueous-phase binary correction (current, adopted 2026-09-28)

**Motivation.** v7 (frozen Kihara/fβ above) reproduced the P(T) curves but gas solubilities
in water were 2–100× too low (CO₂ x≈1e-17 — a negative-a artifact of the paper's printed
van Laar k_ij for CO₂–water; see Validation tab item 4), making the paper's own Fig. 1
irreproducible. LC-region and N₂ high-T errors also persisted.

**Method (stage 2).** Freeze v7 Kihara + fβ entirely. Fit only 10 parameters
DK_AQ[g] = (dK0, dK1): aqueous-phase-only shift k(g,w) → k(g,w) + dK0 + dK1·(T−273.15),
active only when the flash finds a water-rich phase (x[IW] > 0.5). Vapor phase untouched
⇒ ice branch, Langmuir constants, occupancies, Fig. 5 and structure selection cannot move.
Objective (fitv13.py/gnfit13.py): root-space anchor residuals ln(P_root/P_exp) with
multi-root-robust bracketing ([P/1.8, 4P], 9-pt sign-change scan, brentq), w=1 on all
anchors above 273.2 K; 10 solubility anchors w=0.6 clip ±3; NO mixture anchors (invented
NG Pexp values corrupted an earlier run — never anchor on made-up data). Custom damped GN,
~70 s/iter, converged cost 6.59 → manual guard-first override for CH₄ (scan_v8.py:
dK0 ∈ [−0.04,−0.10] keeps CH₄ AADP at v7 level; −0.164 buys solubility 0.9 but AADP 11.7).

**Final dK:** CH₄ (−0.100, 0) · C₂H₆ (−0.1128, +0.002866) · C₃H₈ (−0.0277, −0.001056)
· N₂ (−0.4968, +0.007746) · CO₂ (−1.9906, −0.008089). Files: calib18.npy, params_v8.json.

**Results (val8):** AADP CH₄ 7.63 (v7 7.71) · C₂H₆ 72.2 (72.1) · C₃H₈ 4.74 (4.75) ·
N₂ 39.4 (38.7) · CO₂ 30.1 (42.1). Solubility ratios model/exp: CH₄ 0.32/0.45,
C₂H₆ 0.96/1.00, C₃H₈ 1.00, N₂ 0.98/0.97, CO₂ 0.83/1.02 (v7: 0.09–0.45, N₂/CO₂ ≈ 0.00).
Structures, Fig5, NG, occupancies unchanged from v7. CO₂ LC improved (+129/+207% at
289/291 K, was +270/+536%); C₂H₆ LC and N₂ 281–293 K remain the documented structural
limits (single fβ(T) + frozen paper kij cannot flatten dr/dT there without wrecking V).

**Rejected on the way:** 23-param root-space GN with free Kihara/fβ/dK (cost 56: Kihara σ
drift to bounds + fβ B→5987 destroyed LC/N₂ and flipped Fig5 to sII); solubility-clip
gradient death; invented NG anchors.

**Ported:** engine_port.js (patch_v8.py), node parity ALL PASS incl. 3 solubility checks
(/tmp/jsref_v8.json); app gained Overview mode, Fig. 1 solubility mode, deviation+parity
panels on Compare; index.html params tables now auto-generated from params_v8.json.
