# Ewald reformulation — derivation and calibration note (REVIEW-CORRECTED)

**Status: REVIEW-CORRECTED derivation plus exploratory calibration.** This closes the on-paper math
needed to write Step 3. The appended MATLAB results were produced by scratch code that was not
retained; they support parameter selection and the Γ convention but do **not** count as prospective
Gate-A/B/C passes. The frozen preregistration must require retained test/oracle code, exact inputs
and outputs, and fresh formal reruns. §F (integration surfaces) lives in
`docs/invzp_ewald_integration_map.md`.

Pinned convention (from design §3, re-derived and confirmed by the controller):
- `T_ab(x) = (3 x_a x_b − |x|²δ_ab)/|x|⁵ = +∂_a∂_b(1/|x|)`  (the `+` sign is load-bearing).
- Primitive returns `MF_dipole`'s quantity and sign: `dip_ab,nm(q) = −Σ'_R T_ab(R+d_nm)
  exp[−i q·(R+d_nm)]`, `d_nm = τ_m − τ_n`, **total-displacement phase**, Å⁻³.
- Reciprocal phase is `exp(+i G·d_nm)` (Poisson on the total-displacement sum: shifting `u=x+d_nm`
  yields `e^{+iG·d_nm}·ĝ(q+G)`), NOT `exp[−i(q+G)·d]`.
- Gauge covariance: `dip_nm(q+K) = exp(−iK·d_nm)·dip_nm(q)` ⇒ raw sublattice matrices are
  gauge-covariant; sorted eigenvalue multisets are periodic (matches Phase-1 item-2's sorted-branch
  test).
- The initial boundary is `conducting_k0_omitted`: the exact reciprocal `k=0` term and any Ewald
  surface term are absent from the primitive. Lorentz reconciliation and caller-owned ellipsoidal
  demagnetization are separate.
- For a finite reusable reciprocal candidate set, reduce `q` to
  `qbar=q-floor(q+1/2) in [-1/2,1/2)^3`, evaluate at `qbar`, and restore
  `exp[-i2πK·(tau_m-tau_n)]`. The candidate union must be complete for the entire canonical cube,
  not only the first q used to build a cache.

---

## §A. The three convergent parts (from design §3)

```
dip^(r)_ab,nm(q)  = −Σ'_{R:|x|≤r_cut} g_ab(x) exp(−i q·x),          x = R + d_nm
dip^(G)_ab,nm(q)  = +(4π/Vc) Σ_{G:k≠0,|k|≤g_cut} (k_a k_b/|k|²) e^{−|k|²/4α²} e^{+i G·d_nm},  k=q+G
dip^(self)_ab,nm  = −δ_nm δ_ab · 4α³/(3√π)
dip = dip^(r) + dip^(G) + dip^(self),  the exact k=0 reciprocal vector omitted and no surface
term added ("conducting_k0_omitted").
```

## §B. Blocker 1 — the screened real-space tensor g_ab (DERIVED)

For a radial `φ(r)`: `∂_a∂_b φ = (φ″/r² − φ′/r³) x_a x_b + (φ′/r) δ_ab`. With `φ(r)=erfc(αr)/r`:

```
φ′(r)  = −(2α/√π) e^{−α²r²}/r − erfc(αr)/r²
φ″(r)  = (4α/√π) e^{−α²r²}(α² + 1/r²) + 2 erfc(αr)/r³
```

Assembling and regrouping onto the bare-tensor basis gives the standard damped dipolar tensor:

```
g_ab(x) = ∂_a∂_b[erfc(αr)/r]
        = P(r)·(3 x_a x_b − r² δ_ab)/r⁵  +  Q(r)· x_a x_b
P(r) = erfc(αr) + (2αr/√π) e^{−α²r²}
Q(r) = (4α³/√π) e^{−α²r²} / r²
```

Equivalently `g_ab = A(r)δ_ab + B(r)x_a x_b`, `A=−erfc/r³ − (2α/√π)e^{−α²r²}/r²`,
`B = 3erfc/r⁵ + (2α/√π)e^{−α²r²}(2α² + 3/r²)/r²`.

Checks: (i) α→0 ⇒ `P→1, Q→0` ⇒ `g_ab → (3x_ax_b − r²δ)/r⁵ = T_ab` (bare tensor). (ii) Controller
spot-check at α=1, r=1: `φ′=−0.572406, φ″=1.975037, A=−0.572406, B=2.547442, P=0.572406,
Q=0.830224` vs the closed forms `P=0.572406, Q=0.830219` (agree). **[Full finite-difference
calibration across (r,α), and the sign of every term against a small explicit enumerated sum, is in
the appended calibration section — the crux the design §3 flags as "a named paper alone is not
enough."]**

## §C. Blocker 2 — k=0 / Γ / Lorentz mapping (STRUCTURE derived; one constant is numerical)

The reciprocal G=0 term (only k=0 is omitted, so at nonzero grid q the G=0 term with k=q is present)
is `+(4π/Vc)(q_a q_b/|q|²)e^{−q²/4α²} → +(4π/Vc) q̂_a q̂_b` as q→0. Hence, splitting off the regular
(direction-independent) part `dip_reg = dip^(r) + Σ_{G≠0} dip^(G) + dip^(self)`,

```
dip_ab(q→0, q̂) = dip_reg,ab + (4π/Vc) q̂_a q̂_b .
```

In coupling units (consumer forms `J_cc = −gfac·dip_cc + …`):

```
J_cc(q→0, q̂) [Ewald conducting_k0_omitted] = J_reg,cc − gfac (4π/Vc) q̂_c²
J_cc(q→0, q̂) [existing convention] = J_reg′,cc + gfac (4π/Vc) (1/3 − q̂_c²)
```

The **directional parts are identical** (`−gfac(4π/Vc)q̂_c²`). The two forms differ **only** by an
isotropic constant `gfac·4π/(3Vc) = lorz` between the two regular parts:
`J_reg,cc[Ewald] = J_reg′,cc + lorz`.

To prevent that relation from being hidden in one overloaded metadata name, the integration exports
two matrices:

```text
Jgamma_cc[Ewald]    = −gfac*dip_reg,cc + Jex,cc
Jpath_base_cc       = Jgamma_cc − lorz*ones(4)
J(q→0,qhat)         = Jpath_base_cc
                      + gfac*(4π/Vc)*(1/3−qhat_c²)*ones(4).
```

For brute force, `Jpath_base_cc` is its existing no-Lorentz `Greg` and
`Jgamma_cc=Greg+lorz*ones(4)`. Thus the same q-path formula is used without mistaking the raw Ewald
regular matrix for the legacy-normalized path base.

**Pre-calibration decision test for the integration (design §4.2/4.3):** whether the Ewald branch
adds `0`, `lorz`, or another constant at exactly Γ reduces to a single, well-posed numerical
question — **does
`−gfac·dip_reg,cc(0) + 4·J12` already equal the physical `Jcc0 = 0.006424435656` meV, or is it short
by exactly `lorz`?** If equal ⇒ the Ewald regular part already contains the isotropic term ⇒ the Ewald
branch adds **0** at Γ (retaining the caller's `+lorz` would double-count). If short by `lorz` ⇒ the
Ewald branch keeps the caller's `+lorz` unchanged. At exactly Γ the uniform mode uses the isotropic
`(1/3)δ_ab` average (the sphere/Lorentz convention), and demag (`dm_cc`/`dm_aa` via `ellipsoid_demagn`)
rides on top exactly as today — the four-bullet semantics (`Jcc0`/`Tc` demag-invariant, demag exported)
are unaffected either way. **[The numerical determination of this one constant, plus the §4.3 multi-ray
reconciliation (a*/b*/c* + two generic rays, full Cartesian `δ_ab/3 − q̂_aq̂_b` identity, both `Jcc0`
and `Jaa0`), is in the appended calibration section; Calibration D resolves the decision to
“adds 0.”]**

## §D. Blocker 3 — independent numerical reference (design Gate B) [calibration section]

Plan: implement a scalar-Coulomb Ewald potential `V(x)=Σ' erfc(α|R+x|)/|R+x| + reciprocal + self`,
and obtain an independent dipolar-tensor oracle by numerically differentiating `−∂_a∂_b V` w.r.t. the
sublattice displacement, compared against the direct dipolar-Ewald `dip` above with an explicit
convention/sign map. (A brute-force spherical extrapolation at generic q bounded away from Γ is a
secondary cross-check only; it cannot validate the Γ boundary term.) **[Built and run in the appended
calibration section.]**

## §E. Blocker 4 — parameter ladders + tolerances (design §5.1/§5.2) [calibration section]

To be calibrated by the convergence study (append): a deterministic default `alpha` (Å⁻¹) rule, an α
bracket, ≥3 `r_cut` and ≥3 `g_cut` values with **separate one-axis** tests (grow `r_cut` at fixed
generous `g_cut` and vice-versa, so cancelling truncation errors cannot fake convergence), a joint
refinement, minimum accuracy guards, separate hard enumeration/resource bounds, and raw-tensor +
derived-coupling tolerances (expressed via the frozen `J_ref`). **[Recommended values from the study
appear in Calibration E; exact norms, samples, and resource caps remain Step-3 deliverables.]**

## §F. Blockers 5–6 — integration surfaces

Wrapper API (physics anchors) and the `invz_jq_path` Ewald exact-Γ/metadata contract are mapped in
`docs/invzp_ewald_integration_map.md` (separate note).

---

**Historical pre-append state:** the §C Γ-correction constant, finite-difference confirmation,
independent-reference agreement, and ladder calibration were open before the appended work. They
are now closed as design calibration, but remain prospective formal gates. Nothing here authorizes
production code, default changes, or Phase-3/3A/3B/4 work.

---
---

# APPENDED CALIBRATION (scratch-only MATLAB; no production code touched)

**Status: exploratory design calibration, not formal gate evidence.** The scratch work addressed the
§B finite-difference coefficients (A), a candidate primitive (B), an independent scalar reference
(C), the §C Γ decision (D), and parameter ladders (E). It used throwaway scripts under MATLAB R2025a
against the live `invz_ion()`/`MF_dipole`/`invz_jq_modes`; **no reusable repo test/helper or complete
machine-readable output was retained**. Every number below must therefore be regenerated by the
post-preregistration implementation. The calibration found no inconsistency in the pinned
convention (`T=+∂∂(1/r)`; `dip=−Σ'Te^{−iq·(R+d)}`; reciprocal phase `e^{+iG·d}`; self
`−δ_nm δ_ab·4α³/(3√π)`; conducting `k=0` omission).

Reference constants (live): `gfac=0.08388`, `Vc=287.8917187 Å³` (`Vc^{1/3}=6.6030`),
`4π/Vc=0.0436496425424 Å⁻³`, `lorz=4π/(3Vc)·gfac=1.22044400549e-3 meV` (`4·lorz=4.88177602194e-3`),
`J12=−1e-4` (`4·J12=−4e-4`). Live `invz_jq_modes(invz_ion(),[0 0 0])` (dpRng=30):
`info.Jcc0 = 0.00642443565570149 meV`, `info.Jcc0_dipole = 0.00682443565570149`,
`info.Jaa0 = 0.0035104462050649`, `info.Jaa0_dipole = 0.0039104462050649` (target `Jcc0` confirmed).

## Calibration A — §B screened Hessian `g_ab` (provisional pass)

Central-difference `∂_a∂_b[erfc(αr)/r]` (h=1e-4) vs the closed `g_ab = P(r)(3x_ax_b−r²δ)/r⁵ +
Q(r)x_ax_b`, over 200 random `x` per α (|x|∈[0.4,6] Å) at α∈{1e-6,0.05,0.1,0.2,0.3,0.5,1,2}:

- **Global max relative error = 1.34e-6** (per-component max ≈1.0–1.3e-6, uniform across all 9
  components). **Sign mismatches on non-negligible components = 0** (every component and sign confirmed).
- The 1.3e-6 is finite-difference truncation, not a formula defect: **Richardson extrapolation**
  (leading O(h²) removed) at the hardest small-r points matches the closed form to **4.88e-9** (abs).
- **α→0 limit:** at α=1e-8, `g_ab` reproduces the bare tensor `T_ab=(3x_ax_b−r²δ)/r⁵` to **0.0e0**
  (exact; `P→1, Q→0`) over 2000 points. Controller spot-check reproduced at α=1,r=1:
  `P=0.572406`, `Q=0.830215`, `g(1,1)=A+B=1.975028`.

## Calibration B — scratch `dip_ewald(q,ion,alpha,r_cut,g_cut)` (provisional pass)

`dip = dip^{(r)}[g_ab, |R+d|≤r_cut] + dip^{(G)}[+(4π/Vc)Σ_{k≠0,|k|≤g_cut}(k_ak_b/k²)e^{−k²/4α²}e^{+iG·d}]
+ dip^{(self)}`, `[3,3,4,4]`, Å⁻³, `k=q+G`. Sanity vs `MF_dipole` (N=40) at a generic q:
`max|Δ|=8.49e-6` (abs; scale 5.4e-2 — brute-force conditional-truncation level).

| Check | Result |
|---|---|
| (ii) conjugation `dip(:,:,m,n)=conj(dip(:,:,n,m))`, 5 generic q | **2.9e-17** |
| (ii) Γ: block Hermiticity `‖B−Bᴴ‖`, and `‖Im‖` | **2.4e-17 / 1.8e-26** |
| (i) gauge covariance `dip(q+K)=e^{−iK·d}dip(q)` (6 K × 5 q) | **8.9e-17** |
| (i) sorted-eig periodicity `sort(eig[cc])(q+K)=…(q)` | **5.2e-17** |
| (iii) α-independence, bracket [0.16,0.40], α-matched fixed cutoffs | **2.70e-15** |
| (iv) r_cut convergence (separate axis): residual ≈ `erfc(α·r_cut)` | rc=20→1.2e-10, rc=24→1.1e-14 (α=0.24) |
| (iv) g_cut convergence (separate axis): residual ≈ `exp(−g_cut²/4α²)` | gc=2.5→3.6e-12, gc=3→1.4e-16 (α=0.24) |

α- and cutoff-independence hold to machine precision once cutoffs are matched to α (the earlier
"residuals" were pure truncation tails: e.g. α=0.08 at rc=34 leaves `erfc(2.72)≈1.3e-4`). The
separate one-axis growth (grow r_cut at fixed generous g_cut, and vice-versa) drives each residual
monotonically to machine zero along its own tail law, so cancelling truncation cannot fake convergence.

## Calibration C — candidate Gate-B independent reference (provisional pass)

Independent scalar-Coulomb Ewald potential `φ_lat(x;q)=Σ_R e^{−iq·R} erfc(α|x+R|)/|x+R| +
(4π/Vc)Σ_{k≠0} e^{−k²/4α²}/k² e^{ik·x}` (phase on **R only**; separate code path, no reuse of the
tensor formulas). Identity `dip_ab,nm(q) = −e^{−iq·d_nm}·∂_a∂_b φ_lat|_{x=d_nm}` reproduces
`dip^{(r)}+dip^{(G)}` exactly; the scalar k=0 term is a constant and drops from the Hessian.

- **Off-diagonal (n≠m) blocks**, 3 q-values × 12 blocks, FD-Hessian (h=1e-4) vs `dip_ewald`:
  **worst |Δ| = 9.17e-9** (FD-limited). Independently confirms the reciprocal phase `e^{+iG·d}`, the
  `+(4π/Vc)` projector, and the boundary structure.
- **Self term isolated** on n=m (R=0 excluded by index): `dip_nn − ref_nn` matches the analytic
  `−δ_ab·4α³/(3√π)` to **4.14e-9** (FD-limited; magnitude 1.175e-2).
- **Secondary** brute-force `MF_dipole` extrapolation at benign q=[0.137 0.291 0.453] (away from Γ):
  max|MF−dip_ewald| = 2.50e-5 (N=20) → 1.53e-5 → 8.49e-6 → 6.37e-6 → **3.39e-6 (N=60)**, monotone
  toward `dip_ewald`.

## Calibration D — §C Γ-constant decision + multi-ray reconciliation (provisional)

Regular part `dip_reg = dip^{(r)}(0)+Σ_{G≠0}dip^{(G)}(0)+dip^{(self)}` = `dip_ewald([0 0 0])` (k=0
omitted; α=0.25, rc=30, gc=5 ⇒ converged to machine precision). Uniform projection exactly as
`invz_jq_modes` (`v=[1 1 1 1]/2`):

```
P_ewald,cc = v'·(−gfac·dip_reg_cc)·v = 0.00682166181275765 meV
   vs "adds 0"    target  Jcc0_dipole            = 0.00682443565570149   (Δ = −2.77e-6)
   vs "adds lorz" target  Jcc0_dipole − 4·lorz   = 0.00194265963375772   (Δ = +4.88e-3)
P_ewald,cc + 4·J12 = 0.00642166181  vs physical Jcc0 = 0.00642443566   (Δ = −2.77e-6)
```

### DESIGN DECISION: the Ewald `conducting_k0_omitted` branch adds **0** at Γ (NOT `lorz`).

`v'·(−gfac·dip_reg_cc)·v + 4·J12` reproduces the physical `Jcc0` to **2.77e-6** while the "adds-lorz"
alternative is off by **4.88e-3** (three orders larger) — unambiguous. The residual 2.77e-6 is the
**dpRng=30 target's own conditional-convergence error**, not a missing constant: the brute-force
sphere projection `+4·lorz` converges to `P_ewald,cc` as N grows (N=40→60: 0.0068218 → 0.0068213 →
**0.0068213**, i.e. `|P_ewald − converged brute| ≈ 3.2e-7`). **The Ewald regular part already contains
the isotropic Lorentz term** that the brute-force branch obtains by adding `lorz`; retaining the
caller's `+lorz` on the Ewald branch would **double-count**. Structurally, per element,
`dip_reg_ab(Ewald) = dip_sphere_ab(brute) − (4π/Vc)(1/3)δ_ab` (the isotropic average of the directional
term); verified per element over all 16 sublattice pairs to **6.2e-6 (cc) / 7.8e-6 (aa)** (brute
level). Consequently the physical Γ couplings are `Jcc0 = v'·(−gfac·dip_reg_cc)·v + 4·J12` and
`Jaa0 = v'·(−gfac·dip_reg_aa)·v + 4·J12` with **no extra constant**; `Jaa0` channel:
`P_ewald,aa = 0.00391183` vs `Jaa0_dipole 0.00391045` (Δ=1.39e-6). Demag (`dm_cc`/`dm_aa`) and the
four-bullet semantics ride on top unchanged.

**Multi-ray test (§4.3), uniform (macroscopic) projection** — approach Γ along a*, b*, c*, and two
generic rays [1 1 1], [2 1 −1]; the macroscopic block `dip^{(G=0)}=(4π/Vc)q̂q̂` is identical for all 16
pairs (`e^{+iG·d}=1` at G=0), i.e. the all-sublattice uniform block `ones(4,4)` (v-projection factor
`(Σv)²=4`):

- Recovered `(4π/Vc)q̂_aq̂_b` (full Cartesian) matches the analytic projector with **q² scaling**:
  fullCartErr ≤ 6.3e-6 at s=1e-3, ≤ **6.3e-8 at s=1e-4** (100× → confirms it is the analytic term in
  every Cartesian component, i.e. the full `δ_ab/3 − q̂_aq̂_b` structure, not just cc).
- Uniform projection is real (`‖Im‖≈1e-21`, even in q).
- Same regular matrix from **every ray**: `Jcc0(ray)−Jcc0` and `Jaa0(ray)−Jaa0` are identical to
  −2.774e-6 / +1.387e-6 across all 5 rays (s=1e-4); recovered `dip_reg` spread across rays = **1.1e-8**.
- Macroscopic term **couples only to the uniform mode**: `v'·Δ·v = 1.813e-2` (= predicted
  `4·(4π/Vc)q̂_c² = 1.813e-2`), while the non-uniform block `V⊥'·Δ·V⊥ = 1.2e-6` carries no
  `4π/Vc`-scale term (nonuniform sublattice modes unaffected, as theory requires).

## Calibration E — parameter ladders + candidate tolerances (≤1e-8)

**Deterministic default rule:** `alpha0 = √π / Vc^{1/3} = 0.268431 Å⁻¹`;
`r_cut = C_r/alpha0`, `g_cut = C_g·alpha0`. Tail laws (empirical, matching theory):
raw-tensor REL residual `≈ O(30)·erfc(alpha·r_cut)` (real axis) and `≈ O(3)·exp(−(g_cut/2·alpha)²)`
(reciprocal axis). Coupling REL (vs frozen `J_ref = 0.006424435656`) tracks the raw-tensor REL to
within a factor ≈3. **LiHoF₄ defaults** `C_r=5.5, C_g=11` ⇒ `r_cut0=20.49 Å` (`erfc(5.5)=7.4e-15`),
`g_cut0=2.953 Å⁻¹` (`exp(−30.25)=7.3e-14`):

- **Default residual vs 1.7×-refined:** raw-tensor **3.62e-13**, coupling **2.07e-13** — ~5 orders
  below the 1e-8 target.

**α bracket** `[0.6, 1.5]·alpha0 = [0.161, 0.403] Å⁻¹`; with α-matched generous cutoffs the
α-independence residual is **3.1e-15**.

**r_cut ladder** (≥3, *separate axis*: grow r_cut at fixed generous g_cut):

| C_r | r_cut (Å) | erfc(C_r) | raw REL | coupling REL |
|---|---|---|---|---|
| 4.5 | 16.76 | 2.0e-10 | 3.45e-9 | 1.87e-9 |
| 5.0 | 18.63 | 1.5e-12 | 5.38e-11 | 5.94e-11 |
| 5.5 | 20.49 | 7.4e-15 | 2.04e-13 | 2.17e-13 |

**g_cut ladder** (≥3, *separate axis*: grow g_cut at fixed generous r_cut):

| C_g | g_cut (Å⁻¹) | exp(−(C_g/2)²) | raw REL | coupling REL |
|---|---|---|---|---|
| 9  | 2.416 | 1.6e-9  | 5.47e-9 | 2.67e-9 |
| 10 | 2.684 | 1.4e-11 | 7.79e-11 | 2.21e-11 |
| 11 | 2.953 | 7.3e-14 | 2.53e-13 | 1.45e-13 |

**Candidate tolerances:** raw-tensor REL ≤ **1e-8**; derived-coupling
≤ **1e-8·J_ref ≈ 6.4e-11 meV** (equivalently ≤1e-8 REL). Step 3 must replace “raw-tensor REL” with
an explicit absolute-plus-relative norm that is defined at zero components. Both candidate
tolerances are met at the finest rung and default with ≥3 orders of margin; the coarsest ladder
rung (`alpha*r_cut=4.5`, `g_cut/(2*alpha)=4.5`) sits near the target and therefore defines the
minimum permitted accuracy guard. The earlier suggested 4.4/4.2 guards are rejected: their own
tail model gives approximately `1.47e-8` and `6.55e-8`, above the target. Hard maximum vector,
memory, and elapsed-time caps are a distinct Step-3 item and were not calibrated here.

The on-paper items are closed only as inputs to preregistration. Nothing in this append constitutes
a formal pass or authorizes production code, default changes, cache-schema changes, or
Phase-3/3A/3B/4 work.
