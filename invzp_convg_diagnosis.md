# invzp — why the 1/z susceptibility does not converge in the ordered state

**Independent diagnosis, 2026-07-27.**
Method: the `invz_projected` source itself plus `invz_projected/README.html` for architecture, then direct
numerical experiments in MATLAB R2025a. No document in `docs/`, and no other pre-existing analysis
(`Task-17_failure.md`, the independent QCP analysis, or the Codex reviews) was read before or during the
investigation, so nothing below is inherited from an earlier conclusion. Every number quoted is measured
in this session; the scripts are listed in §8.

**Revised and consolidated on the same day.** First after reading the prior record (§9) and then an
independently written QCP analysis in full (§9.5): that pass withdrew two recommendations, corrected an
outright error about `Dq` versus `D_uni`, and imported its source-backed findings. The standalone QCP
file was subsequently removed after its unique measurements, architectural conclusions, and regression
requirements were preserved here. A second pass followed an external review of the
result, which was substantially right and is the reason §10.1 no longer claims that the strict
truncation is controlled *iff* the resummed scheme is (my own §3.1 refutes it), no longer conflates the
coupling-support edge with the phase boundary, and no longer describes a 16384-pole rational function as
having a branch cut; R3 is withdrawn entirely and D4 demoted from defect to documentation. The final
consolidation also corrects the NaN-site count, defers the proposed h-grid regrading, scopes `omit_max`
to the closure that defines it, and points the continuation experiment at the explicit defactored
square residual in `invzp_convg_fix.md`. Withdrawn historical claims are identified rather than
silently reused. Passages marked as belonging to §1–§8 were written before any of that reading.

**Current-status update, 2026-07-28.** The historical 0.31 K diagnosis
below remains valid for the field cut on which it was measured, but it is
not the current global availability statement. Four later results now
govern:

1. The response evaluator is still healthy once an accepted state exists.
   On the research-priority 0.10 K, `q=[0 0 0]`, 4.60--4.90 T window, the
   unchanged production equations returned 61/61 finite susceptibility
   columns, including both sides of the 1/z transition.
2. The low-field ordered failure is substantially numerical but not a
   one-tolerance problem. A branch-free simultaneous audit removed one
   false rejection, while continuation exposed reproducible folds and
   multiple real roots. Equation closure and physical/Jensen branch
   selection are therefore separate problems.
3. A smooth QCP-connected auxiliary-`h` branch at 4.05 T can be followed
   accurately but has no nonzero Jensen endpoint of its own on the complete
   allowed interval. Its state cannot borrow the legacy branch's free-energy
   crossing. This is the cheapest decisive gate on the proposed smooth-`r`
   backup and prevents an off-shell diagnostic from becoming a production
   spectrum.
4. The accepted near-QCP ordered band is a finite-grid computability
   sliver. Across `12^3--24^3` its contiguous width scales approximately
   as `N^-1.076`, and `Bc_1z` moves by `0.02067 T`. The `16^3` spectrum is
   a visual regression, not a grid-converged critical result.

Section 11 records the measurements and revised net diagnosis. Sections
1--10 preserve the causal derivation and historical decision trail.

All runs: `ion = invz_ion()`, `T = 0.31 K`, transverse field along `a`, 16×16×16 BZ grid (16384 (q,ν)
modes, Γ dropped), `dpRng = 30`, `transverse_mf = 'legacy_x'`, `J(0) = J0eff = 6.42444 μeV`,
`max_q J_ν(q) = 5.98514 μeV`.

---

## 1. The symptom, stated precisely

`invz_spectra_map` masks the ordered columns: `S.phase_1z = 0`, `S.chiz(:,k) = NaN`. The mask comes from
`invz_solve_point_ordered(..., 'ordered_mode','jensen')` returning

```
is_ordered = false, converged = false, hmf_status = 'node_failed'
```

i.e. `invz_hmf_ordered` never finds the H_MF root, so `hstar` is `NaN` and the solver takes its
paramagnetic early return (`invz_solve_point_ordered.m:202-220`).

Measured over field at T = 0.31 K (default `static_medium = 'resummed'`):

| B (T) | ordered (jensen) leg | paramagnetic leg (`invz_solve_point`, `resummed` medium) |
|---|---|---|
| 1.0 | `node_failed` | **not converged**, 200 iters, Σ(0) = 2.823 |
| 2.0 | `node_failed` | **not converged**, 200 iters, Σ(0) = −0.675 |
| 3.0 | `node_failed` | **not converged**, 200 iters, Σ(0) = 0.127 |
| 3.6 | `node_failed` | converged, 86 iters, Σ(0) = 0.1267, crit = −0.1355 |
| 3.8 | ok, m₀ = 1.687 | — |
| 3.9 | ok, m₀ = 1.243 | — |
| 4.0 | ok, m₀ = 0.492, D_uni = 0.0111 | converged, 13 iters, crit = −0.00488 |
| 4.2 | paramagnetic (m₀ = 0) | converged, 13 iters, crit = +0.0453 |

At the driver level (`invz_spectra_map`, T = 0.31 K, 35-point ω grid) that is what the user sees — the
ordered columns carry **no spectrum at all**:

```
   B   phase phase_1z  phase_1z_reason        ordered_diag_reason   finite chi'' points
 1.00     1       0    solver_failed          solver_failed              0 / 35
 2.00     1       0    solver_failed          solver_failed              0 / 35
 3.00     1       0    solver_failed          solver_failed              0 / 35
 3.60     1       0    solver_failed          solver_failed              0 / 35
 3.80     1       1    ordered                accepted                  35 / 35
 4.00     1       1    ordered                accepted                  35 / 35
 4.20     1       2    pm                     not_attempted             35 / 35
 4.60     2       2    pm                     not_attempted             35 / 35
```

Note `phase = 1` on every ordered column: the bare mean-field dispatcher *does* find a moment at 1–3.6 T.
Only the 1/z leg fails.

**The name of the symptom is misleading, and it matters for where to look.** Nothing fails inside the
response function. `invz_chi_realaxis` reads `opts.Jfull` at `invz_chi_realaxis.m:57` and **no caller in
the repository ever sets it** — `invz_spectra_map`'s `copts` (`:517`) and `invz_spectra_qpath`'s do not —
so the `npass` loop breaks on its first pass. The production real-axis evaluation is a single shot with
no nonlinear medium loop that *could* fail to converge. Handed a valid state it is healthy; measured on a
60-point ω grid at T = 0.31 K:

| B (T) | state | Σ(0) | finite χ''(ω) | finite Σ(ω) | `npass` 3 → 12 |
|---|---|---|---|---|---|
| 1.0 | `bare` ordered | −0.0607314 | 60 / 60 | 60 / 60 | 0.000e+00 |
| 3.8 | `jensen` ordered | +0.0604428 | 60 / 60 | 60 / 60 | 0.000e+00 |
| 4.0 | `jensen` ordered | +0.0965479 | 60 / 60 | 60 / 60 | 0.000e+00 |
| 4.6 | paramagnet | +0.0733638 | 60 / 60 | 60 / 60 | 0.000e+00 |

So "χ does not converge in the ordered state" has to be read as **"the thermodynamic state that χ needs
is not available from the current solver"**. Everything below is about state construction, strictly
upstream of the response.
(This restatement was also the executive conclusion of the independent QCP analysis consolidated in
§9.5; the table is my independent re-measurement and reproduces its Σ(0) values at both anchors to every
digit.)

Two facts that frame the historical 0.31 K cut below:

1. The ordered leg works **only in a ~0.2 T sliver just below B_c**, and fails everywhere deeper in the
   ordered phase. It is not "slow to converge deep in the FM phase" — it is the opposite of the usual
   expectation, in which the deeply ordered state is the easy one.
2. **The paramagnetic solver fails in the same region.** This is not an ordered-branch bug. It is a
   property of the effective-medium closure that both branches share.

Neither statement means the present production solver fails throughout the
ordered phase at every temperature. Section 11.4 gives the later 0.10 K
QCP-window result, where all sampled ordered and paramagnetic columns are
available.

---

## 2. Root cause

### 2.1 The effective medium is a Brillouin-zone average over the RPA denominator

Both medium solvers average `1/(RPA denominator of mode q)` over the grid:

* dynamic / paramagnetic slot — `invz_emt_scalar.m:58`
  ```matlab
  A(w) = mean_q [ J_q / (1 + Sigma(w) + J_q*G0(w)) ]
  K    = A.*(1+Sigma)./(1 - A.*G0)
  ```
  with `G0 = -chi0^cc`, so the denominator is `1 + Sigma(0) - J_q*chi0^cc(0)`.

* ordered static sector — `invz_emt_static_ordered.m:62-65`
  ```matlab
  Gq     = Gstat ./ (1 + (J_q - K0).*Gstat)
  K0_new = mean(J_q.*Gq) / mean(Gq)
  ```
  denominator `Dq = 1 + (J_q - K0)*Gstat`.

In both cases the denominator **is** the Landau/RPA stability function of mode `q`. It is positive in the
paramagnet and changes sign, mode by mode, as the state goes unstable — which is the definition of being
inside the ordered phase. (Strictly, the grid modes go negative slightly *after* the uniform mode does,
because Γ is excluded from the sum; that lag is the sliver of §3.3, quantified in §5/D1 and §9.5.)

### 2.2 In the ordered region a large fraction of the grid sits past that pole

Measured at T = 0.31 K. The last column probes the mode count at a representative Σ(0) = 0.1; the
self-consistent values are given underneath.

| B (T) | χ₀^cc(0) (meV⁻¹) | Δ (meV) | 1 − J(0)·χ₀ | Σ(0) needed to make even the *worst grid mode* marginal | # of 16384 modes with negative denominator at Σ(0) = 0.1 |
|---|---|---|---|---|---|
| 1.0 | 951.7 | 0.0486 | −5.114 | **+4.696** | 5976 |
| 3.0 | 251.4 | 0.2659 | −0.615 | +0.505 | 480 |
| 3.6 | 196.5 | 0.3294 | −0.262 | +0.176 | 80 |
| 4.0 | 172.3 | 0.3692 | −0.107 | +0.031 | 0 |
| 4.4 | 154.1 | 0.4070 | +0.010 | −0.078 | 0 |

At 1 T you would need Σ(0) ≈ +4.7 to push the *most dangerous* grid mode back onto the physical branch.
The self-consistent value (taken from the strict scheme of §3.1, which is the only one that produces a
value at all here) is Σ(0) = 1.652 at 1 T, 0.441 at 2 T, 0.165 at 3 T, 0.111 at 3.6 T, 0.090 at 4.0 T.
At every one of those values thousands of modes remain past the pole below ~3.8 T — at 1 T, ~2500 of
16384 even at Σ(0) = 1.65 (interpolating the scan of §2.3), and ~6000 at the O(0.1) values that the
iteration actually visits on its way there.

Comparing the last two columns gives the criterion directly: the medium is well posed only where the
self-consistent Σ(0) exceeds the marginality threshold — true at 4.0 and 4.4 T, false at 3.6 T and
below.

### 2.3 Consequence: the outer map is meromorphic, pole-sensitive, and not a contraction

With sign-straddling denominators the BZ mean is a difference of large numbers and `A(0)`, hence `K(0)`,
becomes a meromorphic function of Σ(0), with a dense set of finite-grid poles and steep sign-changing
intervals. One application of the outer map at B = 1 T, holding Σ(ω_n≠0) = 0 and scanning Σ(0):

| Σ(0) in | # negative denominators | A(0) | K(0) (meV) | **F(Σ)(0) out** |
|---|---|---|---|---|
| −0.5 | 7320 | −1.29e−3 | +2.813e−3 | **+2.226** |
| −0.2 | 6584 | +3.78e−4 | +2.222e−4 | **+0.198** |
| 0.0 | 6376 | −1.23e−3 | +7.066e−3 | **+5.554** |
| 0.05 | 6264 | +1.43e−3 | +6.369e−4 | **+0.523** |
| 0.10 | 5976 | −5.61e−4 | −1.324e−3 | **−1.012** |
| 0.20 | 5912 | −1.86e−3 | +2.906e−3 | **+2.299** |
| 1.00 | 4176 | −1.11e−3 | +4.171e−2 | **+32.66** |
| 2.00 | 1456 | −4.29e−4 | −2.175e−3 | **−1.677** |

There is no fixed point *in this scan* and no monotonic trend; the sign of `A(0)` flips between adjacent
samples. Damping cannot repair this — a damped Picard iteration converges only on a map that is at least
continuous and locally contracting.

Two limits on how far that last sentence may be pushed. First, the scan is one-dimensional (Σ(ω_n≠0)
held at 0), so it bounds the behaviour of the ω = 0 slot, not the existence of a fixed point of the full
map. Second, the independent run preserved in §9.5 measured exactly that: at this field,
`mix_outer = 0.02` with a
1000-iteration budget **does** converge, after 879 iterations, to `crit = −3.669`. So a real algebraic
fixed point of the continued paramagnetic problem exists here — it is simply a deeply unstable one, and
it is reachable only by damping the iteration into a quasi-continuous limit. It does not rescue the
ordered path: with the same `mix_outer = 0.02`, `max_outer = 2000` and only nine nonzero H_MF nodes, the
1 T profile still returned `node_failed` (8/10 nodes accepted, 10,354 outer iterations). The claim that
survives is the operational one — **no mixing constant makes this map usable** — not the stronger
algebraic one.

### 2.4 The two-level self-energy multiplies the damage by ~10³

To linear order in `K(0)` (only the ω_n = 0 Matsubara slot carries `K(0)`):

```
dSigma(0)/dK(0) = (M^2/n01^2) * [ g0^2/beta - (g0 + beta(1-n01^2))*g0/(2*beta) + (g0/beta - (1-n01^2))*g0 ]
```

At (0.31 K, 1 T) — M² = 29.273, n01 = 0.721249, g0 = 29.661 meV⁻¹, 1/β = k_BT = 0.02671 meV — this is
**≈ 7.8 × 10² meV**. That prediction reproduces the scan above almost exactly:

| Σ(0) in | K(0) | 782·K(0) | measured F(Σ)(0) |
|---|---|---|---|
| 0.00 | 7.0658e−3 | 5.53 | 5.554 |
| 0.10 | −1.3240e−3 | −1.04 | −1.012 |
| 1.00 | 4.1709e−2 | 32.6 | 32.66 |

So the outer map is, to excellent accuracy, **Σ(0) ← C·K(0;Σ(0)) with C ≈ 782 meV**, where `K(0)` is the
singular BZ average of §2.1. The prefactor scales as `M²g₀²/(n01²β)`, i.e. it grows as the doublet
splitting Δ shrinks: Δ = 48.6 μeV at 1 T versus 407 μeV at 4.4 T. Deep in the ordered phase the
*paramagnetic* doublet is nearly degenerate, so the stiffness and the pole density peak together.

The resulting iteration is chaotic, not slow. Measured PM loop at B = 1 T (mix 0.7):

```
 iter    dS        Sigma(0)     K(0)        #neg denominators
   1    5.554        0        +7.07e-3       6376
   2   19.41     +3.888      -1.99e-2        184
   3    9.275    -9.702      -5.54e-4      16384
   4    4.312    -3.210      +1.44e-3      13968
   5    1.448    -0.191      -2.13e-3       6584
  ...
```

`dS` must reach 1e−8 for the loop to stop. Over the 60 iterations logged here it never falls below
0.393; in the ordered node solver's own 200-iteration trace at the same point it never falls below
0.123. Note `#neg` swinging between 0 and all 16384 modes from one iteration to the next — the iterate
is jumping across the pole surface, not approaching anything.

### 2.5 Why the *ordered* leg inherits this even though it replaces slot 1

`invz_ordered_node_solve` does discard `invz_emt_scalar`'s ω=0 slot (`K(1) = K0s`,
`invz_ordered_node_solve.m:230`) — but the replacement, `invz_emt_static_ordered` in its default
`'resummed'` branch, has the *same* pole in its own q-average (§2.1). Worse, the Jensen H_MF construction
(`invz_hmf_ordered`, J 2.31–2.33) evaluates that medium at **every molecular field h from ≈ 0 to
1.25·h_bare**, starting from an explicit h = 0 predictor node (`invz_hmf_ordered.m:211`). The h → 0 end of
that path *is* the unstable paramagnet.

Node-by-node at B = 1 T (h_bare ≈ 0.0348 meV, h* ≈ h_bare):

| profile node | h (meV) | outcome | D_uni | Dq_min | Σ(0) | r |
|---|---|---|---|---|---|---|
| predictor (h=0) | 0 | `max_iter` (200) | −1.908 | −1.804 | — | — |
| 1 | 4.35e−5 | `max_iter` | −0.812 | −0.755 | 0.936 | 3.33 |
| 13 | 5.80e−4 | `max_iter` | −1.927 | −0.739 | 0.571 | 3.99 |
| 21 | 3.26e−3 | `max_iter` | −0.035 | +0.057 | 2.583 | 4.26 |
| 25 | 7.74e−3 | `max_iter` | +0.570 | +0.602 | **135.0** | 2.99 |
| 27 | 1.19e−2 | `max_iter` | +0.754 | −0.213 | −1.242 | −0.54 |
| **28** | 1.48e−2 | **converged** (20 iters) | +0.779 | +0.794 | −0.255 | 0.764 |
| 30 | 2.28e−2 | converged (15) | +0.954 | +0.957 | −0.102 | 0.896 |
| 33 | 4.35e−2 | converged (13) | +0.989 | +0.990 | −0.048 | 0.951 |

**27 of 33 profile nodes plus the predictor fail; only the 6 nodes nearest the root converge — and those
converge easily** (13–20 iterations, D_uni and Dq_min ≈ 0.8–0.99, r a smooth 0.76 → 0.95).

Two aggravating design choices make this fatal rather than merely lossy:

* `invz_hmf_grid` uses a **geometric grid clustered at h → 0** (`hmin_frac = 1e-3`), so ~80 % of the
  quadrature nodes land below h*/3, i.e. exactly in the region where the medium is pole-crossed. That
  clustering is right near B_c (where h* → 0) and exactly wrong deep in the ordered phase (where
  h* ≈ 0.8·h_max).
* `invz_hmf_status.m:33` maps **any** non-accepted node to `node_failed`, and `path_from_nodes`
  (`invz_hmf_ordered.m:396-401`) NaNs the whole path if only the h = 0 predictor fails. A single bad
  quadrature node destroys the column.

### 2.6 The medium's inner iteration "converges" to unphysical roots

A separate, important observation from the same trace: `resid_static` at the failing predictor node is
frequently ~6e−11, i.e. `invz_emt_static_ordered` *did* satisfy its closure `|mean_q Gq − Gstat| < 1e-10`
— while reporting `Dq_min = −155` and `D_uni = −164`. Once the denominator changes sign the closure
equation `mean_q [1/(1+(J_q−K0)Gstat)] = 1` acquires additional roots, and the damped iteration has no
branch constraint, so it lands on whichever one it happens to reach. **A small `resid_static` is
therefore not evidence of a physical medium.** The `medium_status` machinery reports domain events only
in the strict branch; the resummed branch has no such test.

---

## 3. Three single-variable experiments that confirm the cause

### 3.1 Remove the pole → the ordered leg converges

`static_medium = 'strict_1z_dyson_ref'` replaces the ω = 0 slot with the one-shot moment closure
`K0 = J̄ − μ₂·G_ref` (`invz_medium_moment_closure.m:72`), which has **no q-denominator at all**. Nothing
else changes. Same fields, same grid, same T:

**Ordered (jensen) leg:**

| B (T) | `resummed` | `strict_1z_dyson_ref` |
|---|---|---|
| 1.0 | `node_failed` | `medium_out_of_domain` |
| 2.0 | `node_failed` | **ok**, m₀ = 4.948, D_uni = 0.853, 33/33 nodes, path_omit_max = 0.217 |
| 3.0 | `node_failed` | **ok**, m₀ = 3.784, D_uni = 0.510, 33/33, omit = 0.0995 |
| 3.6 | `node_failed` | **ok**, m₀ = 2.501, D_uni = 0.239, 33/33, omit = 0.0668 |
| 3.8 | ok, m₀ = 1.687 | ok, m₀ = 1.864, omit = 0.0595 |
| 4.0 | ok, m₀ = 0.492 | ok, m₀ = 0.880, omit = 0.0534 |

**Paramagnetic leg** (`invz_solve_point`) **under the same switch:**

| B (T) | `resummed` | `strict_1z_dyson_ref` |
|---|---|---|
| 1.0 | 200 iters, not converged | **converged, 15 iters**, Σ(0) = 1.652 |
| 2.0 | 200 iters, not converged | **converged, 14 iters**, Σ(0) = 0.441 |
| 3.0 | 200 iters, not converged | **converged, 13 iters**, Σ(0) = 0.165 |
| 3.6 | 86 iters | converged, 13 iters, Σ(0) = 0.111 |
| 4.0 | 13 iters | converged, 13 iters, Σ(0) = 0.0896 |

Changing exactly one thing — whether the ω = 0 medium carries a q-denominator — turns a permanently
non-convergent iteration into a 13–15 iteration solve at every field. That is the cause, isolated.

At the driver level the same switch restores the spectra outright (`invz_spectra_map`, same 35-point ω
grid as §1):

```
   B   phase_1z  phase_1z_reason        finite chi'' points   m_1z     D_ord
 1.00     0      medium_out_of_domain        0 / 35            NaN       NaN
 2.00     1      ordered                    35 / 35          4.948    0.8527
 3.00     1      ordered                    35 / 35          3.784    0.5103
 3.60     1      ordered                    35 / 35          2.501    0.2394
 3.80     1      ordered                    35 / 35          1.864    0.1416
 4.00     1      ordered                    35 / 35          0.880    0.0351
```

`Bc_1z = 4.1 T` under both schemes — but that is a **0.2 T field bracket**, not a demonstration that the
boundaries coincide, and it should not be read as one. The properly resolved comparison is §9.4(iii):
at T = 0.10 K the two PM boundaries agree to 0.014 T, while at T = 1.00 K they differ by ~0.1–0.2 T
(§9.4 viii). So the switch does not *grossly* move the boundary, which is all this table supports; the
claim that matters is the narrower one, that the ordered side becomes computable at fields where it
previously was not.

### 3.2 An ordered endpoint is perfectly well conditioned; it is the *path* that is not

`ordered_mode = 'bare'` solves **a** deeply ordered state — at the self-consistent mean-field molecular
field, with no H_MF path integral — still under `'resummed'`:

| B (T) | converged | iters | m₀ | Σ(0) | crit |
|---|---|---|---|---|---|
| 1.0 | yes | 15 | 5.420 | −0.0607 | +0.923 |
| 2.0 | yes | 15 | 4.976 | −0.0739 | +0.801 |
| 3.0 | yes | 16 | 4.001 | −0.0583 | +0.546 |
| 3.6 | yes | 14 | 3.053 | −0.0128 | +0.355 |
| 4.0 | yes | 13 | 2.122 | +0.0325 | +0.215 |

At B = 1 T the *ordered* state has χ₀^cc(0) = 2.50 meV⁻¹ (back-computed from `crit`), against
**951.7 meV⁻¹** for the paramagnetic state at the same field — a factor of 380. The molecular field
splits the doublet, the single-ion response collapses, every RPA denominator is comfortably positive,
and the medium is trivially solvable. The jensen leg fails not at its answer but on the paramagnetic end
of its own quadrature path.

This is *not* an argument for reverting to `'bare'`: that mode's moment is the bare mean-field order
parameter, so it onsets at the MF boundary rather than the 1/z one — the documented Option-A caveat, and
the reason the H_MF construction exists.

**And the bare state is not the state the jensen leg is looking for.** The two moments differ
substantially — at 4.0 T, bare `m₀ = 2.122` against jensen `m₀ = 0.492`; at 3.6 T, 3.053 against 2.501
(strict) — so this experiment does not show "the jensen endpoint is well conditioned". What it does show
is the narrower thing that the argument actually needs: **a deeply ordered state at the same (T, B) is
trivially solvable, so being deep in the ordered phase is not per se what breaks the medium.** The
distinguishing feature of the jensen leg is that it reaches its endpoint through a quadrature path whose
lower end is the unstable paramagnet — and that is where it fails.

### 3.3 The surviving working window is a q-grid artifact

The window exists only because Γ is dropped from the medium sum, so `max_q J_ν(q) < J(0)`. At strict
q = 0 the denominator would be `1 + Σ(0) − J(0)χ₀ = crit`, negative everywhere below B_c. So the grace
region has width set by `(J(0) − max_q J_ν(q))·χ₀` — a discretisation quantity that vanishes on
refinement:

| grid | modes | max_q J_ν (μeV) | J(0) − max_q J_ν (μeV) | resummed jensen at 3.4 / 3.6 / 3.8 / 4.0 T |
|---|---|---|---|---|
| 8³ | 2048 | 5.46319 | 0.9612 | ok / ok / ok / (PM) |
| 12³ | 6912 | 5.82319 | 0.6012 | fail / fail / ok / ok |
| 16³ | 16384 | 5.98514 | 0.4393 | fail / fail / ok / ok |
| 20³ | 32000 | 6.08578 | 0.3387 | fail / fail / ok / ok |

The predicted threshold is where the *worst* grid mode goes marginal,
`1 + Σ(0) − max_q J_ν · χ₀^cc(0) = 0`. On the 16³ grid, using the strict-scheme Σ(0) and the measured
χ₀^cc(0), that expression is −0.0652 at 3.6 T and +0.0586 at 4.0 T, i.e. a zero at **B ≈ 3.81 T** —
which is exactly where the observed window edge sits (fails at 3.6, works at 3.8).

The honest reading of the table: the coarsest grid (8³) has a visibly wider working region, while the
12³/16³/20³ edges all sit between 3.6 and 3.8 T — a separation my 0.2 T field resolution cannot resolve
further. The numerics therefore *support* the mechanism rather than sharply demonstrate the limit. The
analytic statement is what is load-bearing, and it is grid-independent: any apparent success of the
default scheme below B_c exists only because Γ is excluded from the medium sum, and it cannot be made
better by converging the lattice sum.

---

## 4. Summary of the causal chain

1. The effective-medium closure is defined as a resummed BZ average of `1/(RPA denominator)`.
2. That denominator changes sign across the grid once the *worst grid mode* goes marginal — which lags
   the uniform instability by the Γ-exclusion gap `(J(0) − max_q J_ν)·|Gstat|` (§5/D1, §9.5), a
   discretisation quantity that shrinks on refinement. So: throughout the ordered phase, bar a sliver
   below B_c that is a grid artifact (§3.3).
3. The BZ average therefore becomes a sign-indefinite, pole-sensitive meromorphic function of Σ(0);
   the closure has multiple pole-separated roots and the code has no branch constraint.
4. The two-level self-energy amplifies `K(0)` by ≈ 7.8 × 10² meV at (0.31 K, 1 T), so the outer damped
   Picard map has an O(10²–10³), sign-indefinite effective gain.
5. The outer loop therefore never converges; it hits `max_outer = 200` with a chaotic residual.
6. `invz_hmf_ordered` needs that medium at every h along a quadrature path whose lower end is precisely
   the unstable paramagnet, and whose grid deliberately clusters there.
7. `invz_hmf_status` converts any single failed node into `node_failed`; `hmf_star = NaN`;
   `invz_solve_point_ordered` returns its paramagnetic early return; `invz_spectra_map` masks the column.
8. The chain ends there. `invz_chi_realaxis` is never reached on a masked column, and when it *is*
   reached it does not fail (§1) — it is a single-shot evaluation with no medium loop, because no caller
   passes `Jfull`. The reported symptom is one step removed from the defect.

---

## 5. Secondary defects found on the way

Each of these is independent of the main cause. Some are implementation defects; others are
instrumentation or candidate-design constraints and must not be turned into default-path gates.

**D1 — no branch diagnostic on the resummed static closure.** `invz_emt_static_ordered`'s resummed
branch accepts any root of its closure equation. Measured: `resid_static = 6e−11` ("converged") with
`Dq_min = −155`. Record `min_q Dq`, `min_q |Dq|`, negative-mode count, and branch interval per
iteration. Do not turn `min_q Dq > 0` into a universal acceptance gate: intermediate H_MF nodes are
off-shell, and the plan has no physical prescription that would justify rejecting every negative
denominator. Today a medium evaluated past its own pole is indistinguishable in the trace from a good
one; that is the defect.

Note the test is `min_q Dq`, **not** `D_uni`, and the two are not the same condition. Since `Gstat < 0`
and `J_ν(q) < J(0)` for every grid mode (Γ excluded),

```
Dq - D_uni = (J(0) - J_nu(q))*|Gstat|  >  0        =>   D_uni < min_q Dq
```

so `Dq ≤ 0 ⇒ D_uni ≤ 0`, but **not** conversely. The gap between them is exactly the Γ-exclusion sliver
of §3.3 and the `155.7 / 167.1` ladder of §9.1 — the same object three times over — and it is the entire
region in which the default scheme still works below B_c. See §9.3 for what this test can and cannot be
expected to do.

**D2 — the node gate destroys information (but is the right gate).** `invz_hmf_status.m:33`
(`any(~[allnodes.accepted])` → `node_failed`) plus `path_from_nodes`'s NaN-everything on a failed
predictor. A run with 33/34 good nodes and a resolved bracket is indistinguishable at the caller from a
run where nothing worked. That is a **reporting** defect and worth fixing.

It is *not* an argument for relaxing the gate, and an earlier draft of R3 wrongly made it one. See
§10.2: the failed nodes sit exactly where the branch is pole-sensitive, so accepting a path with holes
in it — or interpolating `r(h)` across them — fabricates the Jensen integral rather than evaluating it.
The all-node acceptance rule is severe and, as things stand, scientifically correct.

**D3 — failure evidence was discarded; RESOLVED as instrumentation.**
The original jensen early return exported `pt.hmf_status` but not
`pt.hmf_prof`, losing the per-node evidence on exactly the masked columns.
Commit `261b12b` now preserves the complete profile and copies the
deterministic binding node's `medium_status`, reference denominator, and
margin to the point-level result. This changes reporting, not acceptance,
and does not make an incomplete path admissible.

**D4 — `mix_outer` does not damp the whole outer state** (documentation, *not* a defect). The outer state
is (Σ, λ, K0s), but only Σ is damped (`invz_ordered_node_solve.m:239`). `lam` is refreshed undamped every
pass (`:232`) and re-enters the next static closure through ξ (`invz_gstat_ordered.m:45`), and `K0s` is
warm-carried with its own inner mixing. In the paramagnet λ is a pure function of Σ, so damping Σ damps
everything; in the ordered sector it is not.

I originally filed this as a defect. It is not one: with `λ` derived from the current medium and `K0s`
solved by an inner closure, damping only Σ is a legitimate **block-Picard** scheme, and nothing measured
here shows that damping `λ` or `K0s` separately would change the branch behaviour. It is recorded
because it is a live surprise for anyone tuning `mix_outer` expecting it to damp the whole state, and
because it is the reason a *coupled* full-residual solve (§10.2) is the informative experiment rather
than a mixing sweep.

**D5 — NaN-ignoring convergence test; RESOLVED.** MATLAB's bare `max`
skips NaN entries, so the former `dS=max(abs(Sigma_new-Sigma))` at four
projected-path sites could in principle accept a partially non-finite
self-energy. Commit `c6a5fc4` replaced all four with the shared
NaN-propagating `invz_finite_max_abs`. This fail-closed hardening did not
alter the measured branch topology.

**D6 — candidate-only arithmetic note: `stable_form` is enabled only under a strict scheme.**
`invz_gstat_ordered`'s default branch computes
`Gtil0 = Gstat/(1-K0*Gstat); r = G0bare/Gtil0`, which evaluates `Inf/(-Inf) = NaN` exactly where
`gstat_local_denom = 1 + Σ(0) + K0·G0inel0` crosses zero — an event the H_MF path is *designed* to
traverse (measured negative at B = 1 T profile nodes 15, 18, 19, 24, 27). The `stable_form`
reassociation is documented as exact, not a regulariser. Nevertheless, making it unconditional would
violate the frozen default-path bit-identity gate and does not supply a branch prescription. Preserve
it as a requirement to assess for a successor candidate, not as a Phase-0 edit.

---

## 6. Recommendations, in order of leverage

**R1 — default promotion of the pole-free static medium: WITHDRAWN.** Switching the ordered/Jensen leg
and its PM probe to `static_medium = 'strict_1z_dyson_ref'` is a strong causal experiment: §3.1 shows
that it converts "no ordered χ below ~3.8 T" into "ordered χ from ~2 T upward" and fixes the PM
iteration at the same time. It is not a production recommendation because that exact candidate failed
its preregistered gate (§9.2).

The resummation `⟨1/(1+Σ+J_q G0)⟩` generates terms beyond retained order, and the one-shot truncation
removes the denominator that carries the numerical pole. That establishes mechanism, not correctness:
the candidate's reference leaves its domain at large moment, and its own omitted-order ratios fail on
the registered fixture.

**R2 — immediate h-grid regrading: DEFERRED.** The geometric grid does put ~80% of its nodes near
`h → 0` deep in the ordered phase, but the proposed "coarse bracket, then refine inside" remedy is not
valid: `F(h)` contains the integral from zero to `h`, so a bracket cannot be computed without the
low-`h` path, and refining only inside it leaves the integral error uncontrolled. Revisit the grid only
after a branch prescription exists, as an adaptive quadrature over the full `[0,h*]` interval.

**R3 — ~~downgrade the node gate from all-or-nothing~~. WITHDRAWN.** As first written, R3 proposed
allowing a bounded number of failed *interior* nodes, replacing them by interpolation of `r(h)`,
propagating a quadrature uncertainty into a reported error on `h*`, and making the `h = 0` predictor
non-fatal via the m → 0 identity `r(0) = 1 + Σ(0)`. Both halves are wrong and are withdrawn — the failed
nodes are exactly where the branch may cross poles or change identity, so interpolation manufactures the
integral, and the predictor identity needs the very `Σ(0)` that failed to converge. The full argument
is in §10.2; the independent analysis consolidated in §9.5 reached the same conclusion first.
**What survives is reporting only:** a `node_failed` column should say *which* nodes failed, where they
sat relative to the bracket, and why — without changing what is accepted.

**R4 — instrument D1; do not gate on it or make D6 unconditional.** Pole proximity and branch identity
belong in the trace. A default-path `stable_form` change is excluded by the bit-identity invariant and
would not select the physical branch.

**R5 — NaN-propagating residuals (D5); and solve the coupled (Σ, K0) state as an *experiment* (D4).**
The NaN hardening is straightforward and worth doing, though it is not the cause of anything measured
here. The coupled solve is a diagnostic, not a maintenance fix: as §5/D4 now records, damping Σ alone is
a legitimate block-Picard scheme, and no evidence here says that damping λ or K0s separately would help.
The converged nodes take 13–20 iterations, so a well-posed loop is cheap; the 200-iteration cap is only
ever reached by genuinely non-convergent nodes, and raising it would not help.

**R6 — preserve the failure evidence (D3).** Attach `hmf_prof` (and the per-node `medium_status`,
`Dq_min`, `ref_denom` arrays) to the ordered solver's early return, so a masked column carries its own
diagnosis instead of a bare status string.

**R7 — infer a strict-scheme validity floor from this cut: WITHDRAWN.** The strict PM leg returns
Σ(0) = 1.652 at 1 T and 0.441 at 2 T, and the current closure's `path_omit_max` reaches 0.217 at 2 T
and 0.0995 at 3 T. Those are useful scheme-specific diagnostics, but they do not amend a frozen
promotion gate or establish a universal floor. Report them; do not convert this one cut into a new
acceptance rule.

A structural note for the longer term: the fact that `ordered_mode = 'bare'` converges effortlessly at
every field while the Jensen path does not is worth taking seriously. The H_MF self-consistency
(J 2.31–2.33) is what makes the moment onset at the 1/z boundary rather than the mean-field one — a real
improvement near B_c — but its integral representation forces an excursion through states that the
theory's own medium cannot describe. A formulation that integrates the *stable* branch (e.g. downward
from the ordered endpoint, or in a variable in which the unstable interval is not traversed) would remove
the failure mode at its source rather than making the unstable region merely evaluable.

---

## 7. What this diagnosis does **not** claim

* It does not claim the strict scheme is *physically* correct deep in the ordered phase — §3.1's own
  `path_omit_max` numbers show why convergence is insufficient, and R1/R7 are withdrawn as production
  proposals.
* It is measured at T = 0.31 K on the 16³ grid. The mechanism (§2.1–2.3) is T-independent in form, and
  §3.3 shows the grid dependence goes the wrong way for the default scheme, but the specific field
  thresholds quoted are for this cut only.
* The B = 1 T strict failure is reported as `medium_out_of_domain`; because of D3 the failing node could
  not be attributed from the returned struct. Fixing D3 is a prerequisite for diagnosing that one
  further.
* §1's response check establishes that `invz_chi_realaxis` returns **finite** values on four anchors and
  a 60-point ω grid. Finite is not correct: nothing here tests whether the ordered-form Σ(ω) (J 2.26,
  `realaxis_sigma`) is right, and the two ordered anchors sit in the grid-artifact sliver of §3.3. The
  claim is only that the response evaluator is not the failing component.
* Three measurements came from the independent run rather than being reproduced in the original
  source-first investigation: the Ewald/brute-force comparison, the `mix_outer = 0.02` converged run,
  and the 3.6 T node's residual block. Their full values and provenance are now preserved in §9.5.
* **It does not claim the ordered state fails to exist.** Everywhere the text says the medium is
  "pole-crossed", "not evaluable on the physical branch", or "has no physical-branch root", the claim is
  about the closure equation as the code poses it at a given iterate — not about the physics. A real
  algebraic fixed point of the continued paramagnetic problem was found at 1 T (§2.3); convergence did
  not make it physical, but it does rule out literal non-existence.
* **§10.1's `x` is not a pre-solve screen.** It requires a converged `Σ(0)`; two of the seven fields in
  its table did not converge and are marked. It brackets the window edge; it does not predict it in
  advance.
* **Nothing here establishes where any scheme's domain *ends*.** §10.1 gives the form of the strict
  truncation's expansion parameter, not a bound on its radius at a state that has not been computed.
  Domain limits in this report are measured (§3.1, §9.4), not derived.

---

## 8. Historical reproduction record

The original scripts and logs were removed from the working tree on
2026-07-28 after their results were consolidated here. They hard-code the
author's workstation path and several depend on scratch `.mat` inputs that
were never committed, so they are provenance rather than a current
reproduction suite. They remain recoverable under
`31a7fd0:docs/diagnostics/claude_convg_2026-07-27/`. If restored, edit
`SD`/`LOG` before invoking `matlab -batch "addpath('<dir>'); diagN"`:

| script | log | what it establishes |
|---|---|---|
| `diag1.m` | `diag1.log` | §1 field map: ordered leg vs PM leg, T = 0.31 K |
| `diag2.m` | `diag2_B1.log` | §2.5 node-by-node profile + per-outer-iteration trace at B = 1 T (`invz_hmf_ordered` with `opts.trace = true`). Also writes a 46 MB `diag2_B1.mat` holding the full `prof`/`trc` structs — not committed; re-run to regenerate. |
| `diag3.m` | `diag3.log` | §3.1 resummed vs `strict_1z_dyson_ref`, ordered leg |
| `diag4.m` | `diag4.log` | §2.2–2.4 pole census, map scan, PM iteration trace |
| `diag5.m` | `diag5.log` | §3.2 `ordered_mode='bare'`, strict PM leg, §3.3 grid refinement |
| `diag6.m` | `diag6.log` | driver-level symptom via `invz_spectra_map` (per-column `phase_1z_reason`) |
| `diag7.m` | `diag7.log` | §9.4(i)/(iii) head-to-head at the **Gate-0 fixture** (T = 0.10 K, `nH = 33`) + the strict candidate's own PM boundary |
| `diag8.m` | `diag8.log` | §9.4(vii) B_c^{1/z}(T) over 0.10–1.50 K — locates the temperature the Gate-0 field set fits |
| `diag9.m` | `diag9.log` | §9.4(viii) the Gate-0 ordered field set + PM controls re-run at T = 1.00 K |
| `diag10.m` | `diag10.log` | §10.1 verification of the two Stieltjes identities and the `x` vs `supp(ρ)` window-edge criterion |
| `diag11.m` | `diag11.log` | §1 the response evaluator given a valid state (60/60 finite, `npass` a no-op) and §9.5's exact Γ-gap `min_q Dq − D_uni = (J(0) − max_q J_ν)·\|Gstat\|` |

Minimal reproduction of the core failure:

```matlab
ion = invz_ion();
[Jnu, info] = invz_bz_couplings(ion, struct('grid',[16 16 16],'dpRng',30,'cache',true));
o = struct('J0eff', info.Jcc0, 'Jxx0', info.Jaa0, 'ordered_mode','jensen');
pt = invz_solve_point_ordered(ion, 0.31, [1 0 0], Jnu(:), o);      % -> node_failed
o.static_medium = 'strict_1z_dyson_ref';
pt = invz_solve_point_ordered(ion, 0.31, [3 0 0], Jnu(:), o);      % -> ok, m0 = 3.784
```

---

## 9. Relation to the prior record (added after reading `docs/`)

§1–§8 above were written without reading any prior analysis, deliberately. They were then compared
against the since-retired convergence-fix attempt record (rev. 2, 2026-07-27),
`docs/invzp_strict_medium_gate0_report.md` (the Gate-0 verdict), the
since-retired ordered-state map
(`31a7fd0:docs/invzp_ordered_1z_state.md`),
`docs/INVZ-DEVELOPMENT-RECORD.md`, a since-retired external review, and an
independently written QCP analysis (the last read in full, now consolidated
into §9.5).
This section records what survives that comparison, what is new, and — most importantly — **what in
§5–§7 above had to be withdrawn.**

### 9.1 The root cause was already established. §2 is a re-derivation, not a discovery.

The retired attempt record's §2 states the same mechanism, in a sharper form than mine — as a ladder
of thresholds in the bare local static weight `|G0|` (its table, with `D ≈ 1`, `K0 ≈ 0`):

| `\|G0\|` | event |
|---|---|
| 155.7 | `D_uni ≤ 0` — physical uniform instability |
| 167.1 | first q-average pole, `D + J_max·G0 = 0` — "the resummed closure dies" |
| 416.2 | one-shot `mu2` pole |

with the remark that "the uniform instability precedes the resummed pole by ~7 %".

Those two numbers are exactly `1/J(0) = 1/6.424436e−3 = 155.66` and
`1/max_q J_ν = 1/5.985139e−3 = 167.08`, so their ratio is precisely `J(0)/max_q J_ν = 1.0734` — the
**same Γ-exclusion gap** my §3.3 measured in field space and attributed to dropping Γ from the medium
sum. Their analysis and mine are the same statement in two coordinate systems. Theirs came first
(spec 2026-07-25).

Three further items in my §2 are also already on record:

* The retired ordered-state map
  (`31a7fd0:docs/invzp_ordered_1z_state.md`), "Durable gotchas":
  *"Since `Gstat = −chi < 0` and
  `max(Jnu) < J0eff`, `Dq = D_uni + (J0eff−J(q))·chi ≥ D_uni`, so any `Dq ≤ 0` ⇒ `D_uni ≤ 0`… A
  Picard/EMT loop reporting `converged` proves nothing physical without checking `D_uni`/`min(Dq)`."*
  That is my §2.6 and the inequality I used in §5/D1, already stated.
* `docs/INVZ-DEVELOPMENT-RECORD.md` (Stage 2c, 2026-07-23): *"Stable accepted nodes exist at all four
  fields (6–10 of 34)."* That is my §2.5's "only the 6 nodes nearest the root converge", measured
  earlier at a different temperature.
* The strict one-shot closure tested in §3.1 and withdrawn from promotion in §6/R1 **is their
  hypothesis**, specified, preregistered before any strict run, built as 15 additive files, and
  measured to conclusion.

### 9.2 R1 is withdrawn as a promotion proposal.

This is the substantive correction. `strict_1z_dyson_ref` was put through a preregistered Gate 0 on
2026-07-27 and **FAILED**, on clauses (a), (c) and (e):

* **(a)** `ref_denom_nonpositive` at `B ∈ {0.05, 0.25, 0.5, 1} T`, all three `nH`;
* **(c)** `max(omit_max) > omit_promote = 0.10` on **all 24 rows**, smallest anywhere `0.17782`
  (3 T, `nH = 129`, a row that solves cleanly with a finite root);
* **(e)** both PM controls (3.1, 3.5 T) return converged, in-domain, **negative** mass
  (`crit = −0.884482`, `−0.528219`), violating the frozen positive-mass control.

The frozen preregistration §3 states that on failure the run **stops at diagnosis**, and that changing
`Gref`, carrying another moment, or truncating other Matsubara sectors is a *new theory candidate
requiring a fresh preregistration — never an in-run fallback*. My R1 asked to promote that exact
rejected candidate to the production default. That is not a recommendation I am entitled to make from
convergence evidence, and **R1 should be read as evidence about the mechanism, not as a proposal.**

My R7 ("report `path_omit_max`, flag columns that exceed the gate") is likewise **not** an answer to
clause (c). A promotion gate is a gate, not a label. The prior record's own lesson — *"converged is not
physical"*, stated twice in its §8 — applies to my recommendation exactly as it applied to theirs: I
measured that strict makes the iteration converge, which says nothing about whether the state it
converges to is right.

### 9.3 Two further corrections to §5.

**D1's remedy was over-sold** (and, in its first draft, misstated — see §9.5). The record's durable
gotcha gives `Dq ≤ 0 ⇒ D_uni ≤ 0`, one-way, not an equivalence; the converse fails precisely inside the
Γ-exclusion sliver. So a `min_q Dq > 0` test would *not* flag every intermediate H_MF node — the
existing two-tier verdict (`res.stability` computed for every node, gating none) deliberately does not
reject those nodes, and this test would not either. But that is also the limit of what it buys, in two
directions:

* it is **not sufficient**. The independent run preserved in §9.5 measures a node at 3.6 T that failed
  the outer map (residual 4.0e−2) with a *terminal* `Dq_min = +0.002475` and a static residual of
  4.9e−11. Earlier iterates crossed poles even though the last one did not; a terminal-iterate test says
  nothing about the path the iteration took.
* where it does fire, it reports something no test can repair: **on the branch the H_MF path must
  traverse deep in the ordered phase, the resummed closure has no root on the physical branch — the
  problem is not conditioning.**

So keep D1, but as *reporting* — don't let `resid_static ≈ 1e-11` read as "physical", which is already
the record's durable gotcha — not as a fix, and score it per-iteration rather than at the endpoint.

**D6 conflicts with a frozen gate.** Making `stable_form` unconditional would break gate G9 (the
default `resummed` path is bit-identical, and is asserted so). Given the paragraph above it is also
moot for that path. Keep D6 as a note for a successor candidate, not as an action.

**§7's third bullet overstated the 1 T strict failure.** The prior record's §5.1 — itself corrected
after external review — establishes only that the warm-started *and* cold-retry Picard attempts
*encounter* `1 + Σ(0) ≤ 0` in a bounded interior `h` window (their nodes 21–24 at 1 T, with solved `ok`
nodes on both sides). It does **not** establish that `1 + Σ(0)` passes through zero as a property of a
converged solution, nor that no positive-denominator fixed point exists there. My "the Dyson reference
itself goes out of domain" should carry the same qualification.

### 9.4 What these measurements do add

**(i) The head-to-head the prior record explicitly lists as missing.** Its §9 states: *"Not
established: that the strict scheme extends the ordered domain relative to resummed at 2–3 T.
`resummed` was measured on the real multiset only at 1 T (`node_failed`) and at the two PM controls
(non-convergent). Its status at 2, 2.5, 2.9 and 3.0 T was never measured head-to-head, so the
comparison at those fields is inference, not data."*

Run at the Gate-0 fixture (**T = 0.10 K**, 16³, dpRng 30, `nH = 33`, `Ecut = 40`, `legacy_x`,
`demag = 0`):

| B (T) | `resummed` | `strict_1z_dyson_ref` | Gate-0's own `max(omit_max)` |
|---|---|---|---|
| 0.25 | `node_failed` | `medium_out_of_domain` | (n/a) |
| 0.50 | `node_failed` | `medium_out_of_domain` | (n/a) |
| 1.00 | `node_failed` | `medium_out_of_domain` | (n/a) |
| 2.00 | **`node_failed`** | **`ok`, m₀ = 4.988, omit = 0.660** | 0.66027 |
| 2.50 | **`node_failed`** | **`ok`, m₀ = 4.612, omit = 0.317** | 0.3165 |
| 2.90 | **`node_failed`** | **`ok`, m₀ = 4.213, omit = 0.198** | 0.19766 |
| 3.00 | **`node_failed`** | **`ok`, m₀ = 4.100, omit = 0.178** | 0.17782 |

Those four rows are the ones that were inference. They are now data: **on the Gate-0 ordered field set
at T = 0.10 K the resummed scheme's ordered domain is empty, and the strict scheme's is 2.0–3.0 T.**
The strict scheme does extend the ordered domain there. My `omit` column reproduces the Gate-0 table
to every quoted digit, which cross-validates the two runs against each other and confirms the two
`max(omit_max)` definitions (profile-only here, ledger-wide there) agree to the precision quoted.

This strengthens rather than weakens the Gate-0 verdict's standing: the same measurement that shows
strict extends the domain also reproduces the `omit` ratios that failed clause (c) at every one of
those fields.

**(ii) `max(omit_max)` is strongly temperature-dependent, and Gate 0 froze the hardest temperature.**
Same quantity, same grid, same coupling multiset:

| B (T) | `max(omit_max)`, **T = 0.10 K** (Gate-0) | `path_omit_max`, **T = 0.31 K** (this work) |
|---|---|---|
| 2.0 | 0.660 | 0.217 |
| 3.0 | 0.178 | **0.0995** |
| 3.6 | — | 0.0668 |
| 3.8 | — | 0.0595 |
| 4.0 | — | 0.0534 |

Because `omit_cubic ∝ mu2·Gref²` and `|Gref|` falls as `χ₀` falls, the ratio drops by 2–4× between
0.10 K and 0.31 K, and at 0.31 K it goes **below** the frozen `omit_promote = 0.10` from ~3 T upward.
This does **not** overturn clause (c): the preregistration is frozen at T = 0.10 K, and the verdict
there stands. What it does establish is that clause (c)'s failure is a statement about the *registered
fixture temperature*, not a universal property of the truncation — the region where the moment
truncation sits inside its own gate is non-empty. A successor preregistration should register a
**(T, B) domain of validity** rather than a single cut, and should pick its fixture deliberately: the
frozen one is the worst case.

**(iii) The strict candidate's own PM boundary — and a re-reading of clause (e).** The record's open
question 3 reads: *"A necessary first step is to determine the strict candidate's own PM boundary,
which was never established."* Measured (`invz_solve_point`, T = 0.10 K, same fixture):

| B (T) | `resummed` crit (iters) | `strict_1z_dyson_ref` crit (iters) |
|---|---|---|
| 2.50 | −1.92075 (200, **not converged**) | −1.88732 (15) |
| 3.00 | −1.01463 (200, **not converged**) | −1.00207 (13) |
| **3.10** | −0.859869 (200, **not converged**) | **−0.884562** (13) |
| **3.50** | −0.515711 (200, **not converged**) | **−0.528272** (13) |
| 4.00 | −0.238144 (13) | −0.244039 (13) |
| 4.50 | −0.0528097 (13) | −0.0580326 (13) |
| 5.00 | **+0.0734666** (13) | **+0.071907** (13) |
| 6.00 | +0.24066 (14) | +0.240161 (14) |
| 7.00 | +0.344126 (14) | +0.343889 (14) |

Three things follow.

*First*, the strict candidate's PM boundary is **≈ 4.72 T** at T = 0.10 K (linear interpolation on this
0.5 T grid), and the resummed scheme's is **≈ 4.71 T**. The two agree to 0.014 T, and agree on crit
itself to 2–10 % at every field where the resummed leg converges at all. **The strict truncation does
not move the PM boundary.**

*Second*, my strict values reproduce the record's §5.3 to five digits (−0.884562 vs their −0.884482 at
3.1 T; −0.528272 vs −0.528219 at 3.5 T), so this is the same measurement, not a different one.

*Third* — and this is the consequence — at T = 0.10 K the frozen PM controls at **3.1 and 3.5 T lie
~1.2 T inside the ordered phase**, not in the paramagnet. A negative PM mass there is the physically
correct answer, and the resummed scheme returns the same sign at the same fields (it simply fails to
converge, for the §2 reason, which is why it could not serve as the cross-check the record wanted). So
the PM half of clause (e) fired on **controls placed on the ordered side of the boundary** — a
control-placement defect in the frozen preregistration, in the same class as the three mis-specified
gates the record itself identifies in its §7, and not evidence against the candidate.

**This does not overturn the Gate-0 FAIL.** Clause (a) still fires (`ref_denom_nonpositive` at
0.05–1 T), clause (c) still fires at the registered temperature, and the *ordered* half of clause (e)
still fires because those same fields do not return `status = 'ok'`. What changes is that one of the
three adverse legs — the one the record calls *"the cleanest single anomaly in the run"* — is most
likely an artifact of where the controls were placed. The record's own §9 was right to refuse to call
it "wrong sign above `Bc_1z`": the sign is right, and `Bc_1z ≈ 3.025 T` appears not to be the boundary
at the registered temperature (see §9.4(vii)).

**(iv) The grid dependence of the resummed *window edge*.** The record's probe (iii) reports *"grid
sweep at 1 T: `node_failed` at N = 8, 12, 16, 24 — grid-independent, not a quadrature-density
artifact"*. That is correct and I reproduce it — but 1 T is far inside the failure region, where no
grid could help. Measured at the *edge* (§3.3, 3.4–4.0 T at T = 0.31 K) the resummed window **is**
grid-dependent: 8³ solves at 3.4 and 3.6 T where 12³/16³/20³ fail. Not a contradiction — a different
field regime — but "grid-independent" should be qualified to "grid-independent deep in the ordered
phase".

**(v) Per-column driver evidence, and a wording discrepancy inside the prior record.**
The retired ordered-state map's TL;DR said the panel masks *"at every
field"*; the retired attempt record's §0 said *"masks below `Bc_1z ≈ 3.025
T`"*. At T = 0.31 K the second is right and the first is too strong: the
resummed leg returns complete ordered spectra at 3.8 and 4.0 T (35/35 finite
ω points, `phase_1z = 1`, `phase_1z_reason = 'ordered'`; §1 and §3.1
tables). This stale contradiction is one reason the ordered-state map was
removed from the current documentation set.

**(vii) `Bc_1z ≈ 3.025 T` is not the boundary at the temperature Gate 0 registered.** Measured 1/z PM
mass zero-crossing (`invz_solve_point`, resummed, 16³, 0.5 T sampling, linear interpolation):

| T (K) | 0.10 | 0.31 | 0.60 | 0.90 | 1.20 | 1.50 |
|---|---|---|---|---|---|---|
| **B_c^{1/z} (T)** | **4.709** | 4.021 | 3.562 | 3.192 | 2.717 | 2.001 |

`B_c = 3.025 T` — the Stage-2 closure figure (`INVZ-DEVELOPMENT-RECORD.md`, 2026-07-22: *"`Bc_1z =
3.025 T` versus the PM `crit`-zero at `3.033 T`"*) and the value the Gate-0 control placement assumes —
corresponds to **T ≈ 1.0 K**, not to the registered `T = 0.10 K`, where the boundary is 4.71 T. The
trend is monotone and physical (B_c falls as T rises), so there is no ambiguity about which cut 3.025 T
belongs to.

The Gate-0 field set is therefore **the field set of a T ≈ 1 K cut, evaluated at T = 0.10 K**: 8
"ordered fields" spanning 0.05–3.0 T against a boundary at 4.71 T rather than 3.03 T, and two "PM
controls" at 3.1 / 3.5 T that are 1.2–1.6 T *inside* the ordered phase. Every field is ~1.7 T deeper
below the boundary than the set's construction implies.

That bears directly on two of the three failing clauses, because both scale with `|Gref| ∝ χ₀`, which
is precisely what grows as a field moves deeper below B_c:

* **clause (e)** — the PM half fired on ordered-side controls (§9.4 iii);
* **clause (c)** — the `omit_max` ratios were evaluated much further below the boundary than intended,
  and §9.4(ii) already shows those ratios fall steeply as the state approaches the boundary.

**(viii) The Gate-0 field set re-run at the temperature it fits.** Identical field set, identical
grid, `nH = 33`, `T = 1.00 K`:

| B (T) | `resummed` | `strict_1z_dyson_ref` | strict `path_omit_max` | same at T = 0.10 K |
|---|---|---|---|---|
| 0.05 | `node_failed` | `medium_out_of_domain` | — | — |
| 0.25 | `node_failed` | `medium_out_of_domain` | — | — |
| 0.50 | `node_failed` | `medium_out_of_domain` | — | — |
| 1.00 | `node_failed` | `medium_out_of_domain` | — | — |
| 2.00 | `node_failed` | `medium_out_of_domain` | — | 0.660 |
| 2.50 | `node_failed` | **ok**, m₀ = 3.998 | **0.0749** ✓ | 0.317 ✗ |
| 2.90 | ok, m₀ = 2.258 | **ok**, m₀ = 2.834 | **0.0606** ✓ | 0.198 ✗ |
| 3.00 | ok, m₀ = 1.375 | **ok**, m₀ = 2.400 | **0.0579** ✓ | 0.178 ✗ |

PM controls at the same temperature (prereg requires converged, finite, **positive** mass):

| B (T) | `resummed` crit | `strict` crit |
|---|---|---|
| 3.10 | +0.00690 | **−0.0398** |
| 3.50 | +0.0840 | **+0.0587** ✓ |

Reading the three failing clauses against this:

* **Clause (c) is largely a fixture-temperature artifact.** At the temperature the field set fits,
  `max(omit_max)` is **0.058–0.075 at every field that solves** — comfortably inside the frozen
  `omit_promote = 0.10`. The same three fields give 0.178–0.317 at T = 0.10 K. The record's finding
  that the ratio "exceeds 0.10 on ALL 24 rows" is a statement about evaluating a T ≈ 1 K field set
  1.7 T deeper below the boundary than its construction implies.
* **Clause (e)'s PM half shrinks to a small, genuine scheme difference.** At T = 1.00 K, 3.5 T passes
  (+0.0587) and only 3.1 T fails (−0.0398) — because the strict scheme's own B_c sits ~0.1–0.2 T above
  the resummed one there, so a control placed 0.07 T above the resummed boundary lands just below the
  strict one. That is a measurable difference between two schemes' boundaries, not "the candidate
  predicts the PM state is unstable through 3.5 T".
* **Clause (a) is real and temperature-robust.** `medium_out_of_domain` fires at 0.05–1 T at 0.10 K and
  at 0.05–**2.0** T at 1.00 K. Changing the temperature does not remove it; the reference
  `Gref = G0bare/(1+Σ(0))` leaves its domain deep in the ordered phase either way. This is the
  substantive adverse finding, and it is exactly where the record's §11 item 1 already points:
  *"Where does the low-field breakdown actually live? … The `1+Sigma0` reference was chosen for its
  `m → 0` cross-leg identity; that choice is untested at large `m`."*

**The Gate-0 FAIL verdict stands** — clause (a) fires at both temperatures, and the *ordered* half of
clause (e) fires with it. What these measurements change is the **weight**: two of the three adverse
legs are substantially artifacts of the fixture temperature, and a successor's effort belongs on the
reference construction at large `m`, not on the moment truncation's omitted order.

**This is not a re-verdict.** I ran `nH = 33` at one temperature through the production solver, not the
frozen Gate-0 protocol (all three `nH`, the aggregator, the digest checks, the B = 0 control). The
preregistration is frozen at T = 0.10 K and stays frozen; changing the fixture temperature requires a
dated amendment, and this section is input to that decision, not a substitute for it.

**(vi) Four solver-mechanics findings I do not find anywhere in the record** — D2 (the all-or-nothing
node gate), D3 (the early return discarding `hmf_prof`), D4 (damping applied to Σ but not to λ/`K0s`),
D5 (NaN-ignoring `max` at four projected-path sites) — plus R2's observation that
`invz_hmf_grid`'s geometric clustering at `h → 0` puts ~80 % of the quadrature nodes in the region
where the medium is pole-crossed. The observation stands; the proposed immediate regrading does not
(§6/R2).

### 9.5 Independent replication and consolidated QCP findings

A separate QCP analysis was written from the same code, at the same temperature, in the same week,
without either investigation seeing the other's working. The standalone file has now been retired; this
section preserves its unique evidence and conclusions. The two investigations agree on the mechanism
and cross-validate every shared quoted number:

| quantity | independent QCP run | source-first run in this report |
|---|---|---|
| accepted HMF nodes, 1 T, 0.31 K | 6 / 34 | 6 of 33 profile records + failed predictor (§2.5) |
| accepted HMF nodes, 3.6 T | 33 / 34 | `node_failed` (§1) |
| strict `omit_max`, 0.31 K, 2/3/3.6/4 T | 0.217 / 0.0995 / 0.0668 / 0.0534 | identical (§3.1) |
| Σ(0), `bare` ordered at 1 T | −0.0607314 | −0.0607314 (§1, §3.2) |
| Σ(0), `jensen` at 4 T | +0.0965479 | +0.0965479 (§1) |
| finite response samples given a state | 60 / 60 | 60 / 60 (§1) |
| `\|G0\|` ladder 155.7 / 167.1 | via the earlier attempt record | 155.656 / 167.081 (§9.1, diag11) |

Measurements unique to that run, preserved here rather than left behind in the retired file:

| configuration | result |
|---|---|
| then-current driver anchor, T = 0.10 K, B = 2.85 T | `node_failed`; 9/34 accepted; 5,159 outer iterations, 5,011 with at least one `Dq <= 0`; predictor exhausted 200 iterations with terminal static residual 9.11e−11 |
| T = 0.31 K, B = 1 T | 5,691 outer iterations; 3,869 encountered at least one `Dq <= 0` mode |
| T = 0.31 K, B = 3.6 T failed node | outer residual 4.007e−2, static residual 4.928e−11, `D_uni = −0.1482`, `Dq_min = +0.002475` |
| isolated PM continuation, 1 T | `mix_outer = 0.02` converged after 879 iterations to `crit = −3.669`; 0.05–1.0 did not converge within 1,000 |
| shortened ordered profile, 1 T | the same 0.02 damping still gave `node_failed`, 8/10 accepted after 10,354 outer iterations |
| Ewald versus brute-force, 16³ | relative sorted-multiset difference 1.168e−3 and `Jcc0` difference 4.318e−4; both failed at 1 T with 1/5 nodes accepted |

These measurements impose three important bounds on the diagnosis:

1. The symptom is a missing thermodynamic state, not a nonlinear failure in `invz_chi_realaxis`.
2. A real algebraic continued-PM fixed point can exist even where ordinary Picard is unusable;
   convergence does not make that deeply unstable root the physical ordered continuation.
3. A positive terminal `Dq_min` and a closed static leaf do not imply that the coupled outer node
   succeeds or that its branch was selected consistently.

The third item also corrected an earlier algebra error. For `Gstat < 0`,

```
min_q Dq - D_uni = (J0eff - max_q J_nu)*|Gstat| > 0.
```

At the 4.0 T Jensen endpoint the measured gap is 0.0854417 and the right-hand side is 0.0854417.
Therefore `Dq <= 0` implies `D_uni <= 0`, but not conversely. The gap is the same Γ-exclusion sliver
seen in §3.3, the `155.7 / 167.1` ladder in §9.1, and the support-edge construction in §10.1.

The surviving architectural conclusions from the original 2026-07-21 QCP diagnosis are also retained:

1. RPA and 1/z have distinct pole conditions, `1-J*chi0 = 0` and
   `1+Sigma-J*chi0 = 0`; they require separate phase/state provenance and must not be forced to close
   at one field.
2. A boundary-linearized modified field is not a controlled finite-moment ordered theory.
3. Pole or minimum-singular-value continuation is a more reliable critical-mode diagnostic than the
   brightest pixel of a broadened finite-frequency map.
4. The later implementation of J 2.28–2.29, J 2.31–2.33, and the J 2.34 audit superseded the original
   implementation-gap observations, but not conclusions 1–3.

The present source-first analysis adds the Stieltjes reduction in §10.1, the post-hoc
`x=(1+Sigma(0))/chi0` support diagnostic, and the temperature/fixture analysis in §9.4. It does **not**
make the current `omit_max` generic: that quantity is the `mu3`/`mu4` remainder ratio of the existing
one-shot moment closure and its chosen `Gref`. A successor closure or functional must derive its own
omitted-order diagnostic.

Two reconciliation decisions are binding:

* `ordered_mode = 'bare'` is acceptable as an explicitly labelled diagnostic product, never as a
  silent replacement for a requested `jensen` column or as closure of this defect.
* Failed HMF nodes may not be skipped or interpolated. They sit where poles or branch changes occur;
  filling those holes manufactures the Jensen integral. The all-node gate stays, while its reporting is
  improved.

Consolidated regression requirements for any claimed fix:

1. Use real-coupling anchors: deep ordered, failure edge, near-boundary accepted, and PM.
2. Require an accepted thermodynamic state before response evaluation, with state-construction failure
   distinct from `response_failed`.
3. For a continuation claim, require the explicit defactored square residual; for every accepted state,
   require the independent A–D audit. A small static-leaf residual is insufficient.
4. Require forward/reverse and cold/warm branch agreement.
5. Separate coupling-grid/backend stability from nonlinear-solver stability.
6. Gate every strict/moment candidate with an error diagnostic derived from its own expansion on every
   required path node.
7. Specify and test any analytic, principal-value, or complex continuation before using it in
   thermodynamics.
8. Preserve separate RPA/1z boundaries; `Bc_auto` must not be presented as `Bc_rpa`.
9. Refine the real-axis mesh and broadening only after the Matsubara state is accepted.

### 9.6 Net position after reconciliation

The **diagnosis** is unchanged and now **triply** sourced: the attempt record's `|G0|` ladder (§9.1),
the independent pole-crossing run preserved in §9.5, and §2 here, each reached without the others.

The **recommendations** change. R1 and R7 are withdrawn as proposals (§9.2);
D1 and D6 are corrected (§9.3). The small implementation hardening is now
delivered: R6/D3 preserves `hmf_prof` and binding-node provenance, D2 reports
the failed-node census without relaxing the gate, and R5/D5 makes all four
self-energy convergence norms NaN-propagating. These changes touch neither
the frozen preregistration, the production default, nor the theory. D4's
coupled solve has also been executed as a diagnostic; its result is in
§11.2 and is not a maintenance fallback.

R2 is no longer on this list. The current grid may be inefficient, but a coarse root bracket still
requires the low-`h` integral. Grid redesign is deferred until a branch prescription permits
error-controlled adaptive quadrature on the whole interval.

Two contributions to the successor question the record leaves open:

1. **If pursuing a strict successor, investigate the reference before adding moments.** Section
   9.4(viii) separates the three Gate-0 clauses by temperature robustness: clause (c) (omitted order)
   and the PM half of clause (e) are largely fixture artifacts, while clause (a) —
   `Gref = G0bare/(1+Σ(0))` leaving its domain deep in the ordered phase — fires at every temperature
   tested. This motivates deriving the large-`m` reference before assuming that carrying another moment
   addresses the limiting failure.
2. **Re-examine the frozen fixture before the next preregistration.** `B_c(T)` is measured in
   §9.4(vii); the field set and controls should be placed against `B_c` at whatever temperature is
   registered, and the registration should cover a **(T, B) domain** rather than a single cut.

And the structural point stands: the current H_MF path must traverse the pole-sensitive region in which
the resummed branch is unresolved (§9.3). A promising successor is therefore a **reformulation that
does not require evaluating the unstable path**, rather than merely another numerical iteration.
The independent analysis in §9.5 reaches the same conclusion from the other side and names a concrete
route: one truncated free-energy/effective-action functional, with χ from the stationary Hessian.
Section 10.2 adopts it as Step 3.

---

## 10. Recommendation

### 10.1 The whole medium is one scalar function, and the problem is where its argument sits

Both medium sectors reduce to the **Stieltjes transform of the coupling density**

$$S(x) \;=\; \Big\langle \frac{1}{J_\nu(\mathbf q) - x} \Big\rangle_{\mathbf q,\nu} \;=\; \int \frac{\rho(J)}{J-x}\,\mathrm dJ ,
\qquad \rho(J) = \frac1N\sum_{\mathbf q\nu}\delta\big(J-J_\nu(\mathbf q)\big).$$

*Dynamic / PM sector* (`invz_emt_scalar`). With `D = 1 + Σ`, `G0 = −χ₀`,

```
D + J*G0 = G0*(J - x),      x = -D/G0 = (1 + Sigma) / chi0^cc
A = <J/(D + J*G0)>          = (1/G0) * [ 1 + x*S(x) ]
```

so `A`, hence `K` and `G`, are explicit functions of `S(x)` alone.

*Ordered static sector* (`invz_emt_static_ordered`). With `y = K0 − 1/Gstat`,

```
1 + (J - K0)*Gstat = Gstat*(J - y)   =>   Gq = 1/(J - y)
mean_q Gq = S(y)     and    K0 = <J/(J-y)>/S(y) = 1/S(y) + y
```

so the entire resummed closure is the **single scalar equation `S(y) = Gstat(K0)`**.

**What `S` actually is on this problem.** For the real 16384-value multiset, `ρ` is an empirical measure
and

$$S_N(z) \;=\; \frac1N\sum_{i=1}^{N}\frac{1}{J_i-z}$$

is a **rational function with 16384 simple poles on the real axis** — meromorphic, not cut. A branch cut
appears only after `N → ∞` or after replacing the empirical measure by a smooth density. That
distinction is not pedantic and §10.2 turns on it: any density-based treatment is a *new approximation*,
not a stabler evaluation of the present scheme.

**Where the argument sits, precisely.** Two thresholds must be kept apart, and my first draft of this
section merged them:

| threshold | condition | what it marks |
|---|---|---|
| **physical** | `x = J(0) = J0eff` | the 1/z uniform instability — the phase boundary |
| **computational** | `x = max_q J_ν` | the edge of the finite-grid pole cloud — where the solver breaks |

The first is an identity, not a discovery: `crit = 1 + Σ(0) − J(0)·χ₀^cc(0)`, so `crit = 0 ⟺ x = J(0)`,
and §9.4(vii)'s `B_c^{1/z}` *is* the `crit` zero-crossing. The Stieltjes variable therefore carries the
physical boundary exactly — but at `J(0)`, which the grid excludes, **not** at `max_q J_ν`. The two
differ by `J(0) − max_q J_ν = 0.439 μeV` — the Γ-exclusion gap of §3.3 and §9.5 once again.

So the correct sentence is *not* "the ordered phase is where `x` enters the support". Sweep the field
down from the paramagnet at T = 0.31 K: `x` falls monotonically, and it crosses the two thresholds **in
this order**,

* `x = J(0) = 6.4244 μeV` at **B ≈ 4.02 T** — the state orders (and §9.4 vii's `B_c^{1/z}(0.31 K)` is
  4.021 T, necessarily, since these are the same equation);
* `x = max_q J_ν = 5.9851 μeV` at **B ≈ 3.73 T** — `x` enters the computed pole cloud and the solver
  breaks.

Between them lies the grace sliver: 3.73–4.02 T is **ordered and still computable**, which is exactly
where the 3.8 and 4.0 T resummed jensen columns live. So:

> **The state becomes ordered when `x` passes `J(0)`. The finite-grid medium becomes pole-crossed — and
> the solver fails — only later, when `x` reaches `max_q J_ν`. The window between the two is the
> Γ-exclusion sliver, and it shrinks on grid refinement (§3.3).**

**Both identities are verified numerically** (`diag10`, on the real 16384-point multiset, at converged
production states):

```
(1) A = (1/G0)[1 + x S(x)]  vs the direct BZ mean, all Matsubara slots, B = 4.4 T:
      max abs diff 8.24e-14 ; max rel diff 6.89e-10
      K rebuilt from S alone vs med.K:  max rel diff 6.89e-10
(2) ordered static at the converged B = 4.0 T jensen state:
      mean_q Gq  = -194.496445881
      S(y)       = -194.496445881     rel diff  0.000e+00   (bit-identical)
      closure |S(y) - Gstat|/|Gstat|  = 3.05e-13
      1/S(y) + y = 0.00134004468185   vs K0 = 0.00134004468185   rel diff 1.17e-12
```

*(The 6.9e-10 in (1) is cancellation in the rewritten form when `x` is far outside the support, not an
identity failure. Identity (2) is exact.)*

And the support crossing brackets the observed window edge, at T = 0.31 K against
`max_q J_ν = 5.98514 μeV`:

| B (T) | 3.00 | 3.40 | 3.60 | **3.80** | 4.00 | 4.20 | 4.40 |
|---|---|---|---|---|---|---|---|
| `x` (μeV) | *4.482* | *5.416* | 5.735 | **6.118** | 6.396 | 6.703 | 7.013 |
| `x − max_q J_ν` (μeV) | *−1.503* | *−0.569* | −0.251 | **+0.133** | +0.411 | +0.718 | +1.028 |
| outside supp(ρ)? | no | no | no | **YES** | yes | yes | yes |
| resummed jensen leg | fail | fail | **fail** | **solves** | solves | solves | (PM) |
| PM leg converged? | **no** | **no** | yes (86 it) | yes | yes (13) | yes | yes |

The sign of `x − max_q J_ν` flips between 3.6 and 3.8 T; the resummed ordered leg's window edge is
between 3.6 and 3.8 T. §2's "grace sliver", §3.3's Γ-exclusion argument and the record's `|G0|` ladder
are all this one crossing, seen three ways.

**Two caveats on this table, both real.** *(i)* `x` needs `Σ(0)`, which comes from `invz_solve_point` —
so it is **not** a pre-solve predictor, and at 3.00 and 3.40 T (italicised) the PM leg did **not**
converge, so those two `x` values are read off a terminal iterate and are indicative only. The rows that
locate the edge (3.60, 3.80) are converged, so the bracket itself stands, but the quantity is a
*diagnostic evaluated on a solved state*, not a cheap screen. *(ii)* The genuinely pre-solve version is
the `Σ = 0` one, `x_RPA = 1/χ₀` — a different and cruder threshold, and the distinction is precisely
the independent QCP conclusion preserved in §9.5: the RPA and 1/z boundaries must not be conflated.

**What this does and does not say about the strict truncation.** The strict closure is the leading term
of the geometric expansion of the *same* closure equation. Writing `u_q = (J_q − K0)·Gstat = Dq − 1`,

```
<1/(1+u)> = 1   =>   <u> = <u^2> - <u^3> + ...   =>   K0 = Jbar - mu2*Gref + O(Gref^2)
```

The geometric series is uniformly justified when `max_q |Dq − 1| < 1`, i.e. `0 < Dq < 2` for every
mode. The implemented `omit_max` is narrower: it compares the `mu3` and `mu4` contributions with the
retained `mu2` term for the current one-shot closure and its chosen `Gref`
(`invz_medium_moment_closure.m:56-58`). It is **not** a direct measurement or bound on
`max_q |Dq−1|`.

My first draft concluded from this that strict converges *iff* `x` is outside the support, hence "no
truncation order can be controlled where the resummed scheme fails". **That is wrong, and my own §3.1
refutes it**: at 3.6 T, `x` is inside the computed support and the resummed jensen path fails, yet the
strict path solves the whole profile with `omit_max = 0.0668`. The error was to evaluate a
state-dependent condition as if the two schemes shared a state. They do not — removing the
q-denominator changes the fixed point, so each scheme's convergence condition is evaluated at *its own*
solution, and the resummed scheme's failure at a field says nothing directly about the strict scheme's
expansion parameter there.

What survives is weaker and more useful: `omit_max` is a local, scheme-specific omitted-order
diagnostic evaluated at the strict solution. It should be kept for that closure, not inherited by a
different reference, density prescription, or functional without a new derivation. Its temperature and
field dependence in §9.4(ii)/(viii) is still informative, and the record's question whether some
truncation has a window where both domain and error are acceptable remains open.

### 10.2 What I recommend, in order

**Step 0 — diagnostic hardening: COMPLETE.** D3/R6 now attaches
`hmf_prof` and binding-node provenance to the early return; D2 reports the
failed-node census without changing acceptance; and D5 is NaN-propagating
at all four sites. This converts "masked column, no information" into
"this node, this reason". It does not fix or reinterpret the physics.

**What Step 0 no longer contains, and why.** My earlier R3 proposed relaxing the all-node gate — accept
a bounded number of failed *interior* nodes, interpolate `r(h)` across them, propagate a quadrature
uncertainty — and making the `h = 0` predictor non-fatal via the identity `r(0) = 1 + Σ(0)`.
**Withdrawn, both halves.** The failed nodes are not randomly placed: they sit exactly where the branch
is pole-sensitive and may change identity, so interpolating across them *manufactures* the Jensen
integral instead of evaluating it — the integrand is smooth only where it is computable, which is
circular. And the predictor identity does not help, because the `Σ(0)` it needs is the h = 0
paramagnetic self-energy, i.e. precisely the quantity that fails to converge at these fields (§1). The
all-node gate is severe, but as things stand it is **scientifically correct**, and relaxing it would
convert a visible failure into an invisible fabrication. The independent analysis in §9.5 reached this
conclusion first.

R2 is also absent. Because `F(h)` contains the full integral from zero, neither a coarse bracket nor
inside-bracket refinement avoids the low-`h` path. Grid redesign waits for a branch prescription and an
error-controlled quadrature over the complete interval.

D4 also leaves this list. Damping Σ while `λ` is recomputed from the current medium and `K0s` is solved
by an inner closure is a legitimate **block-Picard** iteration, not a defect; nothing measured here
shows that damping `λ` or `K0s` separately would help. It stays in §5 as documentation — `mix_outer`
does not damp the whole state, which is a live surprise for anyone tuning it — and the coupled solve
belongs in the diagnostic experiment below, not in a maintenance list.

Two additions from the independent analysis in §9.5 remain central. First,
the delivered failure classification exports the failed `h` nodes and
their binding cause. Second, the safeguarded coupled solve has now been run
as a diagnostic, continuing from certified states and retaining
forward/reverse/cold evidence where available.

The implementation must use the explicit defactored square residual, not
`invz_ordered_residual`, as the vector equation for Newton:
`invz_ordered_residual` returns scalar audit blocks. Its newer
`formulation='coupled'` does, however, audit the simultaneous Sigma and
defactored-static equations without re-entering the nested static solver.
The vector solve uses `RΣ=Σmap-Σ` and the defactored K closure, with the
coupled A--D formulation as an independent post-corrector audit. Section
11.2 records the false nested-audit veto this removed and the multiple roots
the corrected experiment exposed. This is not a production patch: a Newton
corrector can identify algebraic branches, but it cannot supply the missing
physical off-shell prescription.

**Step 1 — adopt `S(x)` as the reporting diagnostic now; treat the density continuation as a
preregistered research candidate, not as a numerical fix.** My earlier draft made this the primary
recommendation and described it as "evaluating the medium as `S(x)` on the coupling density rather than
as a 16384-term sum", as if it were an implementation change. It is not, and §10.1 now says why: for
the real multiset `S_N` is **meromorphic**, so replacing the empirical measure by a smooth `ρ` is a new
approximation and a new continuation prescription in one move. Split accordingly:

*1a — free, and worth doing with Step 0.* Report `x`, `x − max_q J_ν`, and the per-iteration count of
modes past the pole. This is pure instrumentation of quantities the code already forms, it changes no
result, and it turns "the outer loop did not converge" into "the argument was inside the pole cloud from
iteration 3 onward". §10.1 shows this bracket is genuine; §10.1's caveat *(i)* — it needs a solved
`Σ(0)` — limits it to post-hoc diagnosis, which is exactly what the instrumentation phase is for.

*1b — a candidate requiring its own spec and preregistration.* Replace `S_N` by `∫ρ(J)/(J−x)dJ` for a
smoothed `ρ`, evaluated inside the support by singularity subtraction,
`PV∫ρ(J)/(J−x)dJ = ∫[ρ(J)−ρ(x)]/(J−x)dJ + ρ(x)·ln|(J_max−x)/(x−J_min)|`, with `S(x ∓ i0) = PV ∓ iπρ(x)`
as the one-sided alternative. What recommends it: outside the support its difference from the empirical
transform can be measured and must vanish under a declared density/grid limit; it replaces a difference
of ~16384 near-singular terms, whose sign pattern changes steeply across a dense pole set (§2.3), by a
smooth transform; and it makes the continuation an explicit, named choice with testable limits — the
`m → 0` cross-leg identity, the bare/MF limit as `Σ → 0`, and continuity of `crit(h)` through
`crit = 0`.

What must **not** be claimed for it, and what I claimed too freely before: that it is merely a stabler
evaluation of the existing equations (it is a different medium, by an amount set by the smoothing); that
it settles which continuation is physical (it does not — it only makes the choice explicit); or that the
imaginary part `−πρ(x)` is "the decay rate of the mode that has gone unstable". That last was wrong as
stated: `ρ` here is the density of *coupling eigenvalues over q*, so `−πρ(x)` is the imaginary part of a
Hilbert transform of a static distribution. Reading it as a damping rate on this off-shell thermodynamic
path needs a dynamical derivation that nothing here supplies.

**Step 2 — specify, then test, a successor strict reference candidate.** Section 9.4(viii) isolates the
temperature-robust failure: `Gref = G0bare/(1+Σ(0))` leaves its domain at large `m`. "Build it from the
ordered local propagator" is motivation, not yet an implementable formula. Before code, derive the
reference, retained order, domain, cross-leg `m → 0` identity, and its own omitted-order estimator.
Dependence on `K0`, `lambda`, or `xi` may reintroduce self-consistent denominator feedback. Any such
formula is a new candidate with a fresh specification and preregistration.

**Step 3 — the durable fix.** The independent analysis's recommendation 4, preserved in §9.5: derive
the ordered thermodynamics *and* the response from one truncated free-energy/effective-action
functional, so that the stationary ordered solution, its field derivative and χ all come from the same
diagram set, with χ from the stationary Hessian. This dominates the other candidates on conceptual
correctness for a precise reason that §10.1 cannot fix:

> Steps 1 and 2 both keep the architecture in which the order parameter is built by **integrating a
> separately resummed susceptibility across the unstable branch** (J 2.31–2.33). Step 1b makes that
> integral computable and makes the continuation choice explicit; Step 2 changes what is truncated.
> Neither *removes the need to choose a continuation through the pole-crossed region* — and that choice
> cannot be settled by convergence data. A functional formulation obtains the ordered state from
> stationarity instead, so the unstable interval is never traversed and the question does not arise.

The steps are complements, not alternatives, and the ordering is by cost, not by merit. Step 0 unblocks
diagnosis of everything else. Step 1a is instrumentation and rides along with it.
Steps 1b and 2 are each a theory candidate needing a preregistration, and on present evidence **Step 2
has the more localized motivation**: its failure mode is isolated to one identified construct
(§9.4 viii), whereas 1b changes the medium itself by an amount nobody has bounded. Neither is ready
without a derivation and specification. Step 3 is a research programme, with the consolidated checklist
in §9.5 — strict retained-order vacuum/one-point/
two-point diagrams, validation against a small exact finite-cluster oracle, analytic handling of the
elastic zero-frequency sector, resummation only where it follows from stationarity of the same
functional, PM/FM boundary closure on the real coupling density. `rigorous_z^1_extension_Codex.md`
is the existing entry point.

One constraint applies to every candidate: derive its error diagnostic from its own approximation.
The current `omit_max` belongs to the present `mu2` closure, its `mu3`/`mu4` remainder, and its chosen
reference. It cannot be inherited by a new reference, a density continuation, or a common functional
merely because each formula contains powers of a propagator.

### 10.3 What I recommend against

* **Do not raise `max_outer`, soften `mix_outer`, or widen a tolerance.** Measured (§2.3–2.4): the map
  is meromorphic and pole-sensitive with an O(10²–10³) sign-indefinite gain, not slowly converging. The
  strongest available test of the contrary is the independent run in §9.5: `mix_outer = 0.02` *does*
  reach a fixed point of the isolated 1 T PM problem after 879 iterations — a deeply unstable one,
  `crit = −3.669`. With `max_outer = 2000` and only nine nonzero H_MF nodes, the 1 T ordered profile
  still returns `node_failed` after 10,354 outer iterations. So the defensible claim is the measured
  one — **no tested mixing constant produced a usable ordered path** — and what heavy damping buys is
  an unstable continued-PM root, not an ordered solution. It is not established that no constant could.
* **Do not promote `strict_1z_dyson_ref` on convergence evidence.** It failed a preregistered gate;
  §9.4 refines *why*, but refining a verdict is not overturning it.
* **Do not use `ordered_mode = 'bare'` as the answer to this problem** — though it is fine as a labelled
  product. It converges everywhere (§3.2) precisely because it never evaluates the medium on the
  unstable branch, but its moment onsets at the mean-field boundary, which is the defect the H_MF
  machinery exists to remove. Shipping it explicitly labelled as a bare-HMF moment-form 1/z response,
  never silently substituted into a requested `jensen` column, is a reasonable interim plotting
  decision. It does not close this item.
* **Do not switch dipolar backends.** The independent measurement in §9.5 found Ewald versus
  brute-force on the same 16³ grid differ by `1.168e-3` in `norm(sort J)` and `4.318e-4` in `Jcc0`,
  and both return `node_failed` at 1 T with 1/5 nodes accepted. Ewald is the correct cure for the
  conditionally
  convergent real-space sum and should be pursued on its own merits; it is not related to this failure.
* **Do not regularise the pole** (broadening, an `ε` in the denominator, a floor on `Dq`). The pole is
  evidence of an unresolved off-shell continuation, not a numerical nuisance. Smoothing it manufactures
  a number under an undeclared theory change rather than computing the existing equations.

---

## 11. Post-diagnosis numerical findings and current interpretation

### 11.1 Simple numerical remedies: useful discriminators, not a global fix

The later work deliberately tested inexpensive remedies before extending the
theory:

- Raising the iteration budget and using heavier mixing did not produce a
  complete usable ordered path. It either left failed nodes or converged to
  an unstable continued-paramagnetic root.
- Reusing a nearby field's self-energy/profile did not remove the
  nonuniformity. A complete 4.05 T profile still left 4.04 and 4.00 T
  incomplete when used as their warm seed; in one measured sequence,
  field-to-field seeding reduced the accepted-node count.
- A defactored fixed-node Newton corrector repaired the three missing 3.6 T
  nodes in 4/4/8 iterations and completed that profile, but at 1.5 T it
  repaired only the predictor and two nodes. This establishes a local
  Picard-stability defect without establishing a global branch.

These results explain why the failure looks stochastic even though no
random-number path is involved. Small changes in field, seed, ordering, or
roundoff change which basin a noncontractive pole-sensitive iteration
reaches. A larger iteration cap cannot remove folds or decide between
multiple algebraic roots.

### 11.2 The simultaneous audit removed a false veto and exposed real branches

The simultaneous unknown is `[Sigma(:);K0]`. The earlier Block-A audit
replayed the nested static Picard solver; in a multi-root region that replay
could select a different `K0` branch and report the difference as failure of
the supplied root. The optional coupled formulation now independently
repeats the simultaneous Sigma and defactored-static equations, rebuilds and
checks `K`/`lambda`, and invokes no nested static solve unless an explicitly
report-only legacy diagnostic is requested.

At seven representative roots, coupled A/C/D were between zero and
`7.3e-15`, and defactored B was below `7.72e-10`; the old nested replay on
the same states falsely reported `3.40e-3` to `3.12e-2`. With the corrected
audit, a 4.05 T continuation retained 126 consecutive accepted states from
`h=0.0041142866` to `0.0117280521 meV`, with maximum coupled residual
`9.76e-10`.

This correction did not restore uniqueness. At the identical
`h=0.01171732294388717 meV`, two A--D-accepted roots have
`r=0.768127507` and `r=0.822169537`, despite a scaled full-state separation
of only `0.00103702`. Continuation and natural-`h` Newton selected different
ones. The audit now answers “does this state solve the declared equations?”;
it does not answer “which solution supplies the Jensen path?”

The retained adaptive fixed-`h` handoff follows the same distinction. It
uses two certified states for its predictor, shrinks and retries after a
state- or `r`-tube rejection, retains every attempt, and has bounded point
and attempt budgets. It is a guarded numerical primitive, not a bridge
across a fold or a physical selector.

### 11.3 Low-field topology and endpoint evidence

At 1.5 T, a 25-seed zero-field census returned 16 accepted attempts grouped
into seven distinct `[Sigma;K0]` roots. Pseudo-arclength traces established
multiple regular folds, including root 4 near `5.56814e-6 meV` and root 6
near `9.54114e-6 meV`; the quoted turn coordinates are diagnostics, not
certified enclosures.

Root 6 provides genuine endpoint evidence. Its returning leg crossed zero
and matched the zero-field census root within `1.42e-14` in the frozen
cluster metric after a predictor correction of `7.43e-12`. The opposite leg
was replayed on 33 decreasing fixed-`h` nodes: all passed A--D with no event
failure, the largest predictor distance was `0.0247702`, and its endpoint
matched root 6 within `7.92e-14` with residual `7.52e-14`.

This closes the zero-field endpoint identity for the two observed root-6
legs. It does not construct the Jensen path: finite accepted samples do not
prove continuous existence, uniqueness, nonsingularity, or event clearance
between nodes. Other legs remain untraced, the root-7 boundary layer remains
unresolved, no complete single-valued section covers a candidate's full
Jensen interval, and no production continuation dispatch is authorized.

### 11.4 The QCP-connected smooth branch fails its own Jensen-endpoint gate

The smooth low-`r` component at 4.05 T was continued from zero field in the
auxiliary `h` coordinate to the common HMF ceiling
`h=0.02313633386565907 meV` (1369 retained states). At that ceiling it had

```text
F_path = -0.00293612268
r      = 0.429100409
m      = 3.27662221
```

and no nonzero increasing `F_path=0` crossing anywhere on its complete
interval. By contrast, the legacy high-`r` component has
`F=+0.00342692213` and `r=0.996459592` at the same ceiling and crosses at
its own `hstar=0.01473374100101341 meV`.

The latter crossing cannot be borrowed: the Jensen integral and therefore
`F_path` are branch-specific. Continuity of the physical state across a
continuous transition as physical field `B` passes the QCP does not imply
that every smooth auxiliary-`h` component is thermodynamically admissible.
The smooth component is an accurate simultaneous-equation branch, but it
returns `no_admissible_endpoint` and remains off-shell for spectra.

### 11.5 The research-priority QCP susceptibility is available as a finite-grid preview

A separate direct production run used the accepted default equations at
`T=0.10 K`, `q=[0 0 0]`, `B=4.60--4.90 T`, and `0--6 GHz`, with 61 field
points, `0.01 GHz` frequency spacing, and `5e-5 meV` real-axis HWHM. It
returned:

- 19 Jensen-ordered and 42 stable-paramagnetic columns;
- zero masked, suspect, or non-finite-peak columns;
- `Bc_1z=4.6925 T`, bracketed by `[4.690,4.695] T`;
- a continuous V-shaped mode: about `1.00 GHz` at 4.60 T,
  `0.174 GHz` at 4.690 T, `0.110 GHz` at 4.695 T, and about
  `1.02 GHz` at 4.90 T.

The finite sampled minimum is set by the field/frequency mesh and
broadening and is not evidence for a finite critical gap. A closer
asymptotic comparison requires local field/frequency and HWHM refinement,
not a new low-field state solver.

This is a deliberately narrow, finite-`16^3` availability claim. The
verified interval is 4.60--4.90 T; a user-broadened 3--6 T driver range has
not thereby acquired the same all-column verification, and the accepted
ordered width must not be interpreted as grid-converged.

### 11.6 The QCP coverage is a finite-grid computability sliver

A prediction-driven grid gate then separated physical criticality from
finite-grid state availability. It retained the exact legacy grid route and
evaluated no real-axis response. On one common `4.000:0.025:4.700 T` mesh:

| grid | solver-grade `Bc_1z` (T) | total accepted | contiguous accepted width below `Bc` (T) |
|---:|---:|---:|---:|
| 12³ | 4.682284546 | 23 | 0.382284546 |
| 16³ | 4.692758179 | 14 | 0.267758179 |
| 20³ | 4.699093628 | 11 | 0.224093628 |
| 24³ | 4.702957153 | 11 | 0.177957153 |

The contiguous widths scale approximately as `N^-1.076`; `N*width` stays
between `4.27` and `4.59 T`. The lower edge is quantized on a `0.025 T`
mesh, so the exponent is descriptive, but the trend strongly supports that
the near-QCP ordered availability band is tied to the excluded-Γ finite-grid
edge and shrinks toward the continuum. Accepted low-field islands remain,
so this is not a single monotone iteration boundary.

The preregistered expectation that `Bc` would move by at most `0.01 T`
between `12^3` and `24^3` was falsified: it moves `0.02067 T`. Linear and
quadratic four-grid `1/N` fits give `4.72386` and `4.72224 T`, but their
agreement is model sensitivity, not a rigorous continuum error bar.

The coupling-only quantities explain both trends. With
`W=Jmax-Jmin`, the excluded-Γ gap falls approximately as `N^-1.103`;
`S_N(J0)` moves from `-196.7302` to `-199.3844 meV^-1`. At the `16^3`
mass root the PM static propagator is `-198.0070 meV^-1`, matching
`S_16(J0)=-198.0061 meV^-1`. The four-grid extrapolation of the latter is
still model-sensitive (`-201.81` linear versus `-202.84 meV^-1`
quadratic).

A phase-aligned susceptibility gate then evaluated ten common offsets
`B-Bc(N)` from `-0.080` to `+0.080 T`. All 40 response columns were finite
and correctly phased. Across `12^3--24^3`, the parabolically refined
soft-mode peaks differ by only `0.38--0.53%`; a `0.001 GHz` refinement on
the extreme grids reproduced the peak spreads. Thus the mode shape after
alignment is already sub-percent stable, while the dominant unresolved grid
effect is the horizontal `Bc(N)` shift. Spectral weight and an independent
analytic pole solve were not graded by this gate.

A traced `16^3` edge pair also localizes the numerical failure. At both
4.400 T (rejected) and 4.425 T (accepted), every iterate stays in the same
rightmost static interval `y>Jmax`; no pole-index switch occurs. At
4.400 T the inner static residual is about `5e-11`, while the outer
Σ-map predictor oscillates with asymptotic residual ratio about `0.916`
and stops at `8.41e-6` after 200 iterations. At 4.425 T the same predictor
reaches `6.37e-9` in 13 iterations. A 1000-iteration cap accepts 4.400 T
but still fails at 4.300 T with Block-A residual `0.00833`: it moves the
edge rather than curing it.

This also rejects a proposed simple inner replacement as currently stated.
`invz_emt_static_ordered` recomputes `Gstat(K0)` inside each K0 iteration;
the implemented scalar equation is not `S_N(y)=constant`, and the measured
failed block already has a converged inner closure. A rightmost fixed-G
bisection would therefore solve a different equation and target the wrong
block unless the coupled scalar reduction and monotonicity are first
derived.

Finally, the exact identity
`F=h0-J0*m=integral(crit dh)` is useful as a cancellation/quadrature
oracle, not as a missing-state fix. At 4.05 T, increasing the HMF grid from
33 to 65 nodes reduces the route difference by about fourfold, but the
65-node profile hits one newly sampled `node_failed` point despite an inner
static residual `4.24e-11`. The complete evidence and reproducible fixtures
are in `docs/diagnostics/invzp_qcp_grid_2026-07-28/`.

### 11.7 Revised net diagnosis

The apparent ordered-state “susceptibility nonconvergence” consists of two
separable issues:

1. **Observable availability near the QCP:** available for visual inspection
   on the verified finite-`16^3` window. The response code is not the failing
   component, and the phase-aligned peak curve is sub-percent stable over the
   tested grids. The absolute field alignment, ordered coverage, and `Bc`
   are not yet grid-converged.
2. **Complete low-field Jensen coverage:** open and secondary. The nested
   Picard map is noncontractive and badly conditioned in parts of the
   auxiliary path; simultaneous equations have multiple roots and folds;
   coordinate changes improve numerical reach but do not select the
   thermodynamic branch; and a smooth root component may fail its own
   branch-specific Jensen endpoint.

Thus the failure is nonuniform, deterministic numerical sensitivity of the
outer Σ closure, amplified by a finite-grid pole edge and combined with an
unresolved off-shell branch prescription—not universal critical slowing,
not a failure of real-axis susceptibility evaluation, and not something an
iteration limit or a fixed-G inner bisection alone can cure.
