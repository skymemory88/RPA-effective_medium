# invzp ordered-solver reformulation: a strict-order static medium

**Status:** design-approved basis. **Not executable** until the numerical preregistration in §6.0 is
filled and frozen, and until the scope decisions in §7 are settled. This is **not** approval to change
a production default.
**Date:** 2026-07-25. **Branch** `invzp-stage2c-diagnostic`, **HEAD** `2ee310b`. No production file
touched yet.
**Conventions:** `G = -chi` (meV^-1), ferromagnetic positive `J`.

**Task.** Un-mask the 1/z ordered spectra panel below `Bc_1z`, keep the window
`Bc_1z < B < Bc_bare` reported as PM, and do it with a defensible 1/z stability criterion rather than
node-convergence luck.

**Scope.** Production phase dispatch is **PM-mass-first** (`invz_spectra_map.m:288-333`). A *successful*
PM probe with `crit_pm > 0` owns the window; a *failed* PM probe is **unknown**, not evidence for order.
The defect to be fixed is the unexplained `phase_1z = 0` region below `Bc_1z` — without converting
solver availability into a phase criterion.

### Record of revisions

This document combines and supersedes two working layouts (`invzp_solver_section1-3_layout_Claude.md`,
`invzp_solver_section4-6_layout_Claude.md`), both now retired. The corrections from two review rounds
are folded into the body rather than kept as separate verdict lists. For provenance, the substantive
reversals were:

- `r` is not generally `1 + Sigma0` at finite moment; the ill-defined path content is `integral (r-1) dh`
  (§A).
- `crit` is a dimensionless mass, not an inverse susceptibility; the dimensional object is
  `crit/chi_path` (§A, §1).
- `D_q^-1` carries the leading collective/RPA stability physics; only its *feedback into `K0`* exceeds
  retained order (§0).
- The reference propagator must be named and constructed by its own primitive; a closure leaf handed
  only the quotient `Gref` cannot diagnose the denominator behind it (§4.1).
- A `mu2`-only truncation omits `mu3*Gref^2` **first**, before the cubic term; both ratios are reported,
  and neither is excused by the near-zero skewness of one fixture (§B, §4.1).
- Block B is revised **in place**, not duplicated; only load-bearing residuals enter
  `res.finite`/`res.accepted` (§4.4).
- A single shared error classifier must govern **every** catch on the strict path, not just the node
  solver (§5.1).
- A failed/unknown PM probe may run the ordered solver diagnostically but may not emit a production
  phase label (§4.4).
- `invz_deltaF_ordered` is a partial, non-gating diagnostic and cannot be promoted into a J 2.34 gate
  (§A, G6d).
- `|K(iw_1) - K(0)|` mixes physical dispersion with the scheme change (G7).
- Every blocking gate needs a Boolean predicate; "reported" is an output. An adverse result fails
  promotion and starts a new preregistration cycle — never an in-run fallback (§6.0).
- The 1/z leg must not stay nested under the auto/overlay `phase==1` branch (§4.4, §7).
- The `Gstat` singularity is **removable in the integrand**, which changes both the risk assessment and
  the promotion rule (§3, §6.0).
- Gate 0's negative-outcome rule gains a fifth condition, (e): a gate that passes while nothing orders
  cannot decide whether stage 4 is worth planning (§6.0).

---

## 0. Decisions of record

1. **Keep Jensen's path-integral self-consistency; reformulate only the static medium.** The earlier
   proposal to replace the full Jensen condition with a new local order-parameter fixed point is
   superseded. This change removes a beyond-order static resummation; it does not invent a new
   ordered-state stationarity condition.
2. **The pole-free static closure applies to BOTH legs (PM and ordered), consistently** — not to the
   ordered leg alone:
   - a formally `O(1/z^2)` PM/FM mismatch is dangerous at a continuous boundary, where the target mass
     is exactly zero; it can become the leading residual gap or critical-field mismatch;
   - the same static approximation on both sides preserves the exact `m -> 0` PM/ordered identity
     within the chosen computational scheme; the ordered-leg-only variant would deliberately abandon it;
   - the ordered-leg-only variant does not in fact have a one-function blast radius —
     `invz_ordered_residual.m:109` independently re-demands the full resummed static closure, and Gate
     C1 requires exact equality with the PM medium, so both contracts need revision either way;
   - preserving the old PM numbers is useful for comparison but is not sufficient reason to join two
     different resummation schemes across one transition. The resummed PM calculation is retained as a
     **clearly labelled legacy comparator**, not as the production path.
3. **One-shot, no `K0`-denominator feedback:** `K0 = Jbar - mu2*Gref`, evaluated once on a
   `K0`/`lambda`-independent `O(1)` reference propagator, in both legs, at `omega_n = 0` only. The
   frozen convention is the **Dyson reference**, `Gref = G0bare/(1+Sigma0)` (§4.1).
   *Documented artifact of the `omega_n = 0` restriction:* `K(1)` (strict) and `K(2)` (exact resummed)
   then differ while `omega_1 -> 0` as `T` falls, and both feed the same `invz_lambdas` sums — a scheme
   discontinuity across adjacent grid points, not a physical one. It must be measured (G7), not assumed
   negligible.
4. **Domain-switched prescriptions are rejected** — the result would depend on where the switch occurs.
   A separately derived local `Phi_1z(m)` is unnecessary unless the consistent strict-order construction
   fails; §A already supplies the Landau-like functional.

### The `D_q` correction, absorbed

`D_q^-1` is **not** an `O(1/z^2)` artifact — it carries the leading collective/RPA stability physics.
What begins beyond retained order is the **feedback of the fully resummed denominator into `K0`**,
beyond `K0^(1/z) = Jbar - mu2*Gref`. `D_q` is therefore kept and reported as the
susceptibility/stability observable; only its use *as the definition of the medium* is dropped.

---

## A. Theory basis

Jensen's spontaneous condition (framework §9.3, J 2.31–2.33) is
`h0(hmf) = integral_0^hmf r dh'` with `r = G0bare/Gtil0`, and `F = h0 - J0eff*m = 0`.
Using `F(0) = 0` and J 2.31 (`dm/dh = -G0bare`):

```
F'(h) = r(h) + J0eff*G0bare(h) =: crit(h)          [ dimensionless mass ]

F(h)  = integral_0^h crit(h') dh'
```

`crit(0)` is exactly `invz_hmf_ordered.m:138`'s `slope_pred = r0n + J0eff*Gb0`, i.e. the existing
`prof.slope0`, generalised to finite `h`.

**Consequence 1 — Jensen's root supplies a one-dimensional Landau-like stationarity condition.**
Setting `Phi_path(m) = integral_0^m F dm'` gives `dPhi_path/dm = F` and `dF/dm = crit/chi_path`. Since
`chi_path = -G0bare > 0`, the ordered root is a local minimum **iff `crit(h*) > 0`**. `crit` is
dimensionless; the dimensional inverse susceptibility / curvature is `crit/chi_path`, so
`chi_1z,uni = chi_path/crit`. The existing bisection (`invz_hmf_ordered.m:216`, sign pattern
`F: - -> +`) is *intended* to select such a locally stable crossing, but the endpoint derivative must be
evaluated, not inferred from a finite-grid sign change. If multiple stable roots exist, compare
`Phi_path`.

Before calling `Phi_path` a thermodynamic free energy, its normalization must be cross-checked against
Jensen's independent relation (J 2.34). **That cross-check is not available from existing code.**
`invz_deltaF_ordered.m` is by its own docstring a "PARTIAL HYBRID DIAGNOSTIC", explicitly "NOT
dF(m=0)", meaningful only at a common cutoff, carrying a documented ~13.7% same-retained-order
static-elastic residual, and it designates the closed 2x2 model as the J 2.34 comparison route — "an
order-consistency check plus a non-gating approximation-fingerprint regression, never 'validated'". It
is validation-only with no production dependant. It therefore becomes *runnable* again once the profile
solves, but it cannot be promoted into a "two routes agree" gate: a genuine J 2.34 normalization check
is separate derivation work, out of scope here.

**Consequence 2 — the obstruction localises to the static elastic correction, not to `Sigma0` alone.**

```
crit = [ 1 - J0eff*chi_path(h') ] + [ r(h') - 1 ]
```

The first bracket is a direct diagonalization, well defined in the unstable window
(`invz_hmf_ordered.m:130-135`). At `m = 0` the exact hybrid identity gives `r = 1 + Sigma0`; at finite
ordered moment it does **not** — from `invz_gstat_ordered`, `r` also contains the elastic factor `xi`
and so depends on `K0` and `lambda(1:2)`. The ill-defined content of the current construction is
accordingly `integral [r(h)-1] dh` over the unstable window. `integral Sigma0 dh` remains a useful
component diagnostic but is not the whole correction except at onset.

**Binding caution.** The ~0.3% PM boundary shift (`Bc_1z 3.025` vs bare `3.033`) does **not** bound
either `integral Sigma0(h) dh` or the full `integral [r(h)-1] dh` deep in the ordered phase. Both must be
measured (G5), not inferred from the boundary shift.

**Consequence 3 — re-anchoring the integral above the unstable window is impossible.** `h0(h*)`
genuinely requires `r` on all of `[0, h*]`. There is no re-anchoring of the ODE `dH0/dHmf = r` that
avoids it: `H0 = H + J0eff*m` is a state function, so the integration constant is fixed at `Hmf = 0`,
and continuation in `B` still needs `dr/dB` over the same window. What needs a pole-free definition is
the finite-order static correction `r - 1` throughout that path — not a replacement for Jensen's
condition.

---

## B. The moment expansion, and why both legs truncate identically

The ordered static closure (`invz_emt_static_ordered.m`) is
`Gq = Gs/(1 + (J-K0)*Gs)`, `mean_q Gq = Gs`, `K0 = mean_q(J*Gq)/mean_q(Gq)`.
Expanding both closure conditions independently in `Gs` (`J = Jbar + d`, `mu_n = E[d^n]`, `mu_1 = 0`) —
they agree identically:

```
K0 = Jbar - mu2*Gs + mu3*Gs^2 + (2*mu2^2 - mu4)*Gs^3 + (mu5 - 5*mu2*mu3)*Gs^4 + O(Gs^5)
```

Every coefficient is polynomial in the moments: **no q-denominator remains after formal expansion**.
Under the high-density power counting `mu2 ~ 1/z`, so the strict `O(1/z)` medium is
`K0 = Jbar - mu2*Gref`, with `Gref` evaluated without feeding the resulting `K0` back into itself. The
fully resummed `1/D_q` feedback enters the *medium definition* beyond retained order. This is a formal
order statement, **not** a proof that the numerical moment series converges near the critical path.

**The PM leg has the same pole, and the same expansion.** `invz_emt_scalar.m` eliminates `K` in closed
form, but its `A = mean_q J/(D + J*G0)` (`:37-50`, `D = 1+Sigma`) carries the identical q-denominator —
`D + J*G0` is proportional to `D_q` — and the file already flags `1 - A*G0 -> 0` as a "vanishing RPA
denominator" (`:79-83`). Re-expressing that closed-form `K` in terms of its **own local dressed `G`**
gives, verified term-by-term to `O(G^4)` (difference exactly 0):

```
K_PM(G_local) = Jbar - mu2*G + mu3*G^2 + (2*mu2^2 - mu4)*G^3 + (mu5 - 5*mu2*mu3)*G^4
              = K_ordered(G_local)                                   <-- SAME FUNCTION
```

So the two legs' exact closures are the same function of `(G_local; mu_n)`. Truncating both at the same
order **and using the same explicitly defined `Gref` at `m -> 0`** makes the cross-phase static-medium
identity exact within the selected numerical scheme. Merely sharing the polynomial is insufficient if
the two callers feed it differently: the shared primitive must receive an already-constructed `Gref`,
and Gate C1 must test the complete caller wiring (§4.2, G2).

### Measured multiset

**Provenance (load-bearing — these numbers are not universal):** production default grid `[16 16 16]`,
`dpRng 30`, `bruteforce` dipole backend, no grid-policy fields. `n = 16384` entries = 4096 q x 4
branches. The symmetric even-count grid contains no exact Gamma point, so `J0eff` never enters the
q-average and only reaches `D_uni`/`crit`.

| quantity | value |
|---|---|
| `J0eff` (`info.Jcc0`) | `6.42444e-3` meV |
| `Jbar = mean_q J` | `1.20766e-4` |
| `sqrt(mu2)` (rms spread) | `2.3415e-3` |
| `mu2` | `5.48264e-6` |
| `mu3` | `-3.42228e-11` (skewness `-0.0027`) |
| `mu4` | `2.3894 * mu2^2` |
| `J_max`, `J_min` | `5.98514e-3`, `-6.7631e-3` |
| `(J_max - Jbar)/sqrt(mu2)` | `2.5045` sigma |

On this multiset `Jmax < J0eff`; together with the physical static sign `Gstat < 0` this gives
`D_q - D_uni = (J(q) - J0eff)*Gstat > 0`, hence `D_q >= D_uni`. The inequality follows from those two
facts, **not** from the omission of Gamma alone.

**These moments support truncation diagnostics, not a convergence claim.** The `mu3` contribution is
negligible *on this multiset* — that is a measured property of one grid/cutoff/backend, and
generalising it is precisely the inference error that produced the synthetic-`Jnu` defect being
repaired here. Both omitted ratios are therefore always reported (§4.1). The first omitted cubic
coefficient is `(2*mu2^2 - mu4) = -0.3894*mu2^2`; relative to the retained `mu2*Gref` term its
magnitude is approximately 5.2%, 6.0%, 9.3% and 37.1% at `|G| = 155.7, 167.1, 208.2, 416.2`, and
becomes equal only near `|G| ~ 684`. "Less than one" is not a small-error guarantee, and still higher
moments can matter.

### Breakdown thresholds, `D = 1+Sigma ~ 1`, `|G0|` = bare local static weight

| `\|G0\|` | event | kind |
|---|---|---|
| **155.7** | `D_uni = 1 + (J0eff - K0)*G <= 0` — physical uniform instability (Gamma) | **real physics; keep as observable** |
| **167.1** | first q-average pole, `D + J_max*G0 = 0` — **resummed closure dies** | arithmetic |
| **208.2** | `mu2` closure, *self-consistent PM/Dyson proxy*: branch point `\|D + Jbar*G0\| = 2*sqrt(mu2)*\|G0\|` | arithmetic |
| **416.2** | `mu2` closure, *one-shot PM/Dyson coordinate*: pole `1 + Jbar*G0 - mu2*G0^2 = 0` at `D=1` | arithmetic |

This is consistent with the recorded diagnostic ordering: the uniform instability precedes the first
resummed q-average pole by about 7%, which is why "below `Bc` => node dies". The hierarchy is monotone
in **how much beyond-`O(1/z)` denominator feedback is retained** — the direct justification for
decision 0.3. The last two numbers are PM ordinary-Dyson estimates only; they are **not** universal
ordered-elastic thresholds, because `invz_gstat_ordered` has a different local function. They motivate
removing feedback but do not replace an actual ordered-path domain map (G0).

**The rejected self-consistent variant, for the record.** Merging `K0 = Jbar - mu2*G` with
`G = G0/(D + K0*G0)` gives a quadratic with physical branch

```
G = [ (D + Jbar*G0) - sqrt( (D + Jbar*G0)^2 - 4*mu2*G0^2 ) ] / (2*mu2*G0),   K0 = Jbar - mu2*G
```

(`-> G0/(D + Jbar*G0)` as `mu2 -> 0`), real only while `|D + Jbar*G0| > 2*sqrt(mu2)*|G0|`. It extends the
domain over the resummed closure iff `J_max - Jbar > 2*sqrt(mu2)`; measured gain `1.2523`, versus
`2.49` for the one-shot form. It is **not** a normal `static_medium` scheme — it needs different inputs
(`G0`, `D`, and a branch choice) — so it belongs in a separately named diagnostic comparator, keeping
it unselectable from the production resolver and keeping the one-shot leaf's API honest.

---

## 1. The construction

Jensen's condition is kept **verbatim**. Only the `omega_n = 0` medium changes.

```
chi_path = -G0bare = dm/dh > 0
crit(h)  = r(h) + J0eff*G0bare(h)                                [ = F'(h), dimensionless ]
chi_1z,uni(h) = chi_path(h)/crit(h)                              [ where stable ]
F(h)     = integral_0^h crit dh'   = h0(h) - J0eff*m(h)          [ Jensen, unchanged ]
h*       : F(h*) = 0,  crit(h*) > crit_tol      [ Landau minimum; tolerances per §1 two-tier verdict ]

Gref     = G0bare/(1 + Sigma0)                                   [ selected Dyson-reference convention ]
K0       = Jbar - mu2 * Gref                                     [ ONE-SHOT, closed form ]
r        = G0bare * (1/Gstat - K0)                               [ see the arrangement note below ]
```

`Jbar` and `mu2` are taken over the same multiset `invz_emt_scalar` already averages. No q-denominator
survives in the medium definition. `Gref` is deliberately **not** the physical ordered `Gstat`: the
latter continues to use Jensen's hybrid elastic expression with the solved `Sigma`, `lambda` and
one-shot `K0`. `D_q = 1 + (J(q) - K0)*Gstat` and `D_uni = 1 + (J0eff - K0)*Gstat` are still built from
that physical `Gstat` and reported in full as the collective/RPA observables.

**Two-tier verdict (do not collapse these):**

- **Path-node consistency:** the complete outer `Sigma`/`lambda` residual contract passes and every
  exported quantity required by `r` is finite. Intermediate nodes may have `crit`, `D_uni` or `Dq_min`
  non-positive, because they are explicitly the analytic continuation through the unstable Landau
  interval.
- **Physical endpoint stability:** a consistent root additionally has `crit > crit_tol`,
  `D_uni > D_tol` and `Dq_min > Dq_tol`, with boundary-scaled tolerances frozen **before the first
  strict-mode run** (§6.0). Tolerances chosen after seeing strict-mode output would be tolerance-tuning
  through the back door, which the durable gotcha forbids. Node convergence remains mandatory; only the
  obsolete inner resummed-closure convergence flag disappears.

**Invariants that survive by construction, not by tuning:**

- `r = 1 + Sigma0` at `m = 0`, for **any** `K0`:
  `1/Gtil0 = 1/Gstat - K0 = (1+Sigma0)/G0inel0 + K0 - K0`. The `K0` dependence cancels identically, so
  the pinned `invz_gstat_ordered` identity cannot be broken by the truncation.
- The static `K0` **stops being an inner iteration.** There is no fixed point and no branch margin in
  the one-shot default. The outer `Sigma <-> lambda <-> medium` Picard loop and its complete residual
  acceptance remain.

**Arrangement note — the `Gstat` singularity is removable in the integrand.** Because
`1/Gtil0 = 1/Gstat - K0`, as `Gstat -> ±Infinity`

```
r    -> -G0bare*K0             (same limit from both sides; r is continuous through the pole)
Gtil0 -> -1/K0                 crit -> G0bare*(J0eff - K0)
```

so the divergence **cancels in the quantity Jensen's condition integrates**. It survives only in
`D_uni` and `D_q`, which are diagnostics. But `invz_gstat_ordered.m:30-31` currently computes
`Gtil0 = Gstat/(1 - K0*Gstat)` and `r = G0bare/Gtil0`, which at `Gstat = Inf` evaluates as
`Inf/(-Inf) = NaN` and loses precision approaching it. The algebraically identical
`r = G0bare*(1/Gstat - K0)` returns the exact limit (verified: the two agree to `|Gstat| ~ 1e18`, then
the as-written form goes indeterminate while the rearranged one is exact). This is a **reassociation,
not a regulariser** — identical value, no broadening, no tolerance — so it is admissible under the
durable gotcha, and it is required for the removability above to hold in floating point.

---

## 2. Components

| unit | change | why this is the boundary |
|---|---|---|
| `invz_common/invz_coupling_moments.m` | **new**, pure | `(Jbar, mu2, mu3, mu4, n)` from `Jnu_flat`; per-column for the `[nJ,nw]` retardation interface. Derived at call time — **no coupling-cache schema change**. |
| `invz_common/invz_static_medium_reference.m` | **new**, pure | Constructs `Gref` **and owns its denominator metadata**. A closure leaf handed only the quotient cannot reconstruct `1+Sigma0` or its margin, so reference construction cannot live inside the closure. |
| `invz_common/invz_medium_moment_closure.m` | **new**, pure | `[K0, info]`: both legs call the same algebraic map with an explicit caller-supplied `Gref`. Carries scheme/order, the strict comparison value, and **both** omitted-order ratios. |
| `invz_common/invz_check_static_medium.m` | **new** | Sole public authority for the scheme id; stamps it into both leg option structs so they cannot disagree. |
| `invz_common/invz_is_recoverable_solver_error.m` | **new**, pure | One shared error classifier for every catch on the strict path (§5.1). |
| `invz_emt_scalar.m` | `omega = 0` slot only, mode-gated | Owns the **PM** static medium. `K(2:end)`/`G(2:end)` and both the vector and `[nJ,nw]` branches untouched; legacy path stays bitwise-reproducible. |
| `invz_emt_static_ordered.m` | scheme-gated: `resummed` \| strict | Owns the **ordered** static medium. Strict mode removes the damped `K0` iteration entirely. |
| `invz_gstat_ordered.m` | `r`/`Gtil0` arrangement only | `r = G0bare*(1/Gstat - K0)`, exact reassociation, so the removable singularity is removable in floating point (§1). No physics change; the pinned identities are unaffected. |
| `invz_ordered_residual.m` | block B revised **in place** | `:109` currently re-demands the discarded resummation. Under strict mode block B's load-bearing residual becomes `\|K0 - Kstrict(Gref)\|`; the resummed residual survives only as an opt-in diagnostic. |
| `docs/invz_ordered_residual_contract.md` | revised block-B section | The contract is binding; it must state the new acceptance rule. |
| `invz_hmf_ordered.m` / `invz_ordered_node_solve.m` | thread the scheme; separate node consistency from endpoint stability | Every predictor/profile/bisection node must use the same static convention, while unstable intermediate nodes stay admissible as analytic-continuation data. |
| `invz_solve_point.m` / `invz_solve_point_ordered.m` | thread the scheme + provenance | PM and ordered point outputs must identify static scheme, reference convention and moments. |
| `invz_spectra_map.m` | PM-mass-first three-way dispatch; expose provenance and two-tier failures | The phase label must come from a converged PM mass where available, never from ordered solver availability. |
| `invz_projected/invz_static_domain_scan.m` | **new**, diagnostic | Gate-0 prospective scan (§7 item 2 settles its grid ownership). |
| tests and saved-result metadata | new strict-mode gates; legacy fixtures retained | Gate C1, PM anchors, ordered residual tests, q-path/spectra wiring and any serialized result schema must distinguish legacy-resummed from strict-order output. |

The ordered leg already replaces `invz_emt_scalar`'s `K(1)`/`G(1)` (`invz_solve_point_ordered.m:217`,
`invz_ordered_node_solve.m:190`), so the strict PM slot must not leak into ordered lambdas before that
replacement. Pinned by a behavioural test (G13), not assumed from statement order.

---

## 3. What this fixes, and what it does not

**Production phase dispatch already has the right hierarchy.** In `invz_spectra_map.m:288-333` the bare
ordered solve first supplies the overlay state; the 1/z leg then performs a strict PM probe. A converged
finite PM point with `crit_pm > 0` is labelled PM (`phase_1z = 2`), and only when that fails does the
transverse Jensen ordered solve run. Therefore:

- `Bc_1z < B < Bc_bare`: the **PM solver and `crit_pm`** own the phase label. The ordered HMF routine
  should normally not be used to manufacture this PM result.
- `B < Bc_1z`: the PM mass is non-positive or no longer a valid stable PM endpoint, so the Jensen
  ordered leg must return a residual-consistent root whose endpoint passes the stability gates.
- `B > Bc_bare`: the existing automatic solve is already PM; the HMF routine's bare-order ceiling guard
  (`:95`) is only a direct-call safeguard.

**`pm_valid` conflates two different causes, and the ordered leg cannot tell them apart.**
`pm_valid = ~isempty(ptp) && ptp.converged && isfinite(ptp.Sigma0) && isfinite(crit_pm) && crit_pm > 0`
is false both when the PM state is *genuinely unstable* (`crit_pm <= 0`) and when the PM probe merely
*failed or threw*. Both fall through to the ordered leg, and today a PM numerical failure plus a
converging jensen yields `phase_1z = 1` — solver availability acting as a phase criterion. §4.4 replaces
this with a three-way dispatcher.

**The `B = 0` guard is a splitting guard, and narrower than "B = 0".**
`invz_twolevel_ordered.m:18-21` raises `invz:degenerateDoublet` when `tl.Delta < 1e-4` meV, with `Delta`
evaluated at **the node's own molecular field** `hp` (`:16-17`, via `hz_fixed`). Because the `h`-grid is
geometric and clustered at 0 (`:144-145`), the predictor and lowest nodes are at risk whenever `Bx` is
small — not only at exactly `Bx = 0`. And `:281` calls the constructor **outside any try/catch**, so the
throw escapes `eval_node`, escapes `invz_hmf_ordered`, and is absorbed by `invz_spectra_map`'s `invz:*`
catch (`:325-328`) leaving `ptj = []` and `phase_1z = 0` — indistinguishable from a solver failure.
Production sweeps start at `B = 0` (`linspace(0,9,101)`), so this path is hit on every run. `Delta` is a
direct diagonalization, so the column can be labelled "outside the two-level domain" rather than
"solver failed" (§5.3).

**Inside `invz_hmf_ordered`,** `slope_pred < 0` still predicts a nonzero root and the negative-to-positive
`F` bracket still selects the local minimum. But `idx` empty with `status = 'ok'` returns `hstar = NaN`,
and `invz_solve_point_ordered` wraps that in `early_return`, whose `converged` field is **false**
(`:272`) — it is not a genuine PM point. The production PM-mass-first dispatcher must be preserved and
the onset identity `slope_pred = crit_pm` gated separately (G2).

**Why the strict construction is the right hypothesis.** The backend-independent diagnostics localize
the observed masking to the resummed static-medium failure: even the `h = 0` predictor fails, and the
outer solve limit-cycles as the inner closure meets widespread non-positive `Dq`. It is **not** a
demonstrated fix until the full residual and path gates pass (G16).

**What it does not guarantee.** The outer `Sigma <-> medium` Picard loop remains an iteration and can
still fail. This removes the inner resummed-closure pole from that map; it does **not** guarantee
smoothness across Jensen's remaining local denominators, nor prove the outer loop converges.

**Residual risk, revised.** The one-shot PM ordinary-Dyson coordinate has an arithmetic pole at
`D + Jbar*G0 - mu2*G0^2/D = 0` (`|G0| ~ 416` in the `D = 1` estimate only). The ordered elastic `Gstat`
instead retains the local denominators of J 2.28–2.29 including the `xi` denominator, and has no
universal `416` threshold. Critically, **the `Gstat` singularity is removable in the integrand** (§1):
`r` and `crit` are continuous through it, so a `Gstat` denominator crossing is not by itself fatal to
Jensen's condition — provided the arrangement fix is in place. The **reference** denominator `1+Sigma0`
has no such cancellation and keeps a hard rule. Gate 0 must therefore map the *actual solved* path and
distinguish the two denominator classes, rather than condemning any node whose local response diverges.

**Durable gotcha, still in force.** Do not broaden or regularize the static response, add a pole
regularizer, flip a sign, or widen a tolerance as a convergence patch. Nothing above does: the change is
a truncation to the order the theory already retains (§B), plus an exact reassociation (§1), and its own
breakdown is reported rather than smoothed.

---

## 4. Data flow

### 4.1 The new pure leaves

All pure, all in `invz_common/`. Reference construction is separated from moment closure because the
latter cannot reconstruct a denominator or its margin once it receives only the quotient `Gref`.

```
mom = invz_coupling_moments(Jnu_flat)
    -> struct('Jbar','mu2','mu3','mu4','n')
       Vector Jnu_flat  -> scalar fields.
       [nJ,nw] matrix   -> length-nw field vectors (per-frequency, T2.1 retardation interface).
                           Every static caller uses column/index 1 only; it must never flatten all
                           frequency columns into one static multiset.

[Gref, ref] = invz_static_medium_reference(G0bare, Sigma0, scheme)
    -> strict_1z_dyson_ref:  denom = 1 + Sigma0;  Gref = G0bare/denom
       strict_1z_bare_ref:   denom = 1;           Gref = G0bare
       ref.denom, ref.margin, ref.status, ref.scheme
       ref.status = 'ok' | 'ref_denom_nonpositive' | 'ref_denom_small' | 'nonfinite'

[K0, info] = invz_medium_moment_closure(Gref, mom_static, scheme)
    -> K0    = mom.Jbar - mom.mu2*Gref                      (one-shot schemes)
       info.scheme        the canonical scheme id
       info.retained      'mu2'
       info.Kstrict       the strict value K0 is checked against (identical for one-shot)
       info.omit_mu3      |mu3*Gref^2| / |mu2*Gref|          <-- the FIRST omitted term
       info.omit_cubic    |(2*mu2^2-mu4)*Gref^3| / |mu2*Gref|
       info.omit_max      max(info.omit_mu3, info.omit_cubic)
       info.status        'ok' | 'nonfinite'
```

All moments use the **population** normalization implied by the BZ average, `mu_n = mean((J-Jbar).^n)`.
MATLAB's default sample-normalized `var(J)` is **not** interchangeable with `mu2`: the `N-1` difference
is 6e-5 relative at `N = 16384` but **4% at the `N = 24` synthetic fixtures**, i.e. largest exactly
where it would go unnoticed. The helper rejects non-finite/non-real inputs and preserves the exact
q/branch weighting of its input.

Ratios use an explicit zero convention: if `mu2*Gref == 0`, report zero only when the corresponding
omitted numerator is also zero, otherwise `Inf`. They are evaluated at the actual state, never inferred
globally. They do **not** reject a path node — the formal polynomial stays defined — but
`max(info.omit_max)` is a **promotion gate** against the preregistered control bound. A large ratio
never triggers a scheme switch. The `'resummed'` path bypasses both new strict leaves entirely.

### 4.2 Scheme resolution: one authority, stamped into both legs

`invz_check_static_medium.m` validates and canonicalises, following the established shared-validator
idiom (`invz_check_coupling_opts.m` from Step-5 Task 7, `invz_check_solve_opts.m`,
`invz_common/invz_check_transverse_mf.m`).

```
scheme ids:  'resummed'                  (DEFAULT — legacy, bit-identical, no flip in this phase)
             'strict_1z_dyson_ref'       (the selected production candidate, §0.3)
             'strict_1z_bare_ref'        (systematic comparator: Gref without the (1+Sigma0))
unknown id -> error('invz:staticMedium', ...)
```

`opts.static_medium` is the sole public authority. At a spectra/phase driver it is a top-level driver
option and is forbidden inside `opts.solve_opts`; direct point solvers accept the same field directly.
Public solver entries validate it and stamp the resolved id into both internal option structs — `eopts`
for `invz_emt_scalar`, `eso` for `invz_emt_static_ordered`. Direct user values in
`opts.emt.static_medium` or `opts.emt_static.static_medium` are rejected as conflicts rather than
allowed to split the two legs. Neither leaf defaults to anything but `'resummed'`, so an unstamped
legacy call stays numerically bit-identical. Gate C1 tests the public caller wiring end-to-end, not only
the shared polynomial.

### 4.3 Moment lifetime

`mean` plus central moments over the production 16384-entry multiset are cheap once but must not be
recomputed inside an outer iteration. Moments are computed **after** the point's actual coupling
spectrum has been resolved (important for field-dependent ODD modes), then threaded unchanged.

- `invz_solve_point` and `invz_solve_point_ordered` compute `mom` once per resolved point coupling. The
  latter passes it into `invz_hmf_ordered` and reuses it for the final root solve; HMF does not
  recompute when supplied.
- A direct `invz_hmf_ordered` call derives `mom` once at entry as a compatibility fallback and places it
  in every `node` struct.
- `node` gains a 14th field, `Jmom`. `invz_ordered_node_solve.m:119-125` and
  `invz_ordered_residual.m:58-63` each carry their own `req_node` list and are already kept in lockstep;
  both extend together. Missing `Jmom` is a wiring error **only after** scheme resolution selected a
  strict mode; legacy/resummed direct-node fixtures may derive or omit it without changing their
  numerical path.
- `invz_emt_scalar` accepts `opts.Jmom` **optionally** and derives it from `Jnu_flat` when absent, so
  every existing caller keeps working unchanged.

For `[nJ,nw]` retardation input, `mom` is per column, while the ordered static closure and all static
`Dq` diagnostics consume **only** `Jnu_flat(:,1)`.

> **Pre-existing latent defect, not introduced here.** `invz_emt_static_ordered.m:43` does
> `Jf = Jnu_flat(:)`, which for a `[nJ,nw]` input flattens *all* frequency columns into the static
> q-average today, and `invz_ordered_node_solve.m:164` deliberately preserves a matrix. Existing
> `Jnu_flat(:)` flattening is not an allowed strict-mode interpretation. If ordered-retarded support is
> not validated in this phase, the resolver must reject that combination explicitly rather than silently
> average frequencies together. Labelling it pre-existing matters so it receives its own fix decision
> instead of being read as strict-mode fallout.

### 4.4 New state carried out — the two-tier split, in data terms

The two tiers of §1 must be *separate fields*, never one flag. Collapsing them re-masks the ordered
phase, because intermediate path nodes are the unstable Landau interval by construction.

**`invz_ordered_residual` (`res`):**

| field | tier | role |
|---|---|---|
| `res.blockB` | consistency | Block B is revised **in place**. Under a strict scheme its load-bearing residual is `\|K0 - Kstrict(Gref)\|`, independently recomputed from the exported state, with fields `resid`, `scale_abs`, `scale_rel`, `pass`, `status`, `scheme`, `ref_denom`, `omit_mu3`, `omit_cubic`. Gate on the preregistered `K_atol + K_rtol*max(\|K0\|,\|Kstrict\|,Jscale)` — **not** on `eps(Jbar)` or on a correction term that can vanish as `Gref -> 0`. `Jscale` already exists at `:88`. Under `'resummed'`, block B keeps its current meaning. |
| `res.blockB.resid_resummed` | diagnostic | under a strict scheme the old resummed-closure residual is **nullable and computed only under an explicit debug/trace flag**. It never contributes to `.finite` or `.accepted`: the analytic-continuation path deliberately crosses its pole. The normal strict checker must not run the discarded inner solve — doing so would restore the very iteration and pole exposure this design removes. |
| `res.blockA/C/D` | consistency | unchanged in form and role. |
| `res.finite` / `res.accepted` | consistency | `.finite` covers the exported state and load-bearing A/B/C/D residuals only. `.accepted = finite && A && B && C && D`. Diagnostic `resid_resummed`, omitted-order ratios and stability signs are excluded. |
| `res.stability` | **stability** | raw `crit`, `D_uni`, `Dq_min`, plus a stable/unstable/boundary-band classification against the frozen margins. Computed for every node, folded into `res.accepted` for **none**. |

`crit` is available at node level from `info.so.r` and `node.G0bare0`: `crit = r + J0eff*G0bare0`.

**`invz_ordered_node_solve` (`info`):** gains `info.medium` (reference + closure metadata) and passes
`info.res.stability` through unchanged. If **reference** construction is out of domain, the attempt
stops before `invz_lambdas`/`invz_sigma_ordered` rather than propagating a non-finite `K0` through the
outer map. `info.term_reason` gains `'medium_out_of_domain'`.

**`invz_hmf_ordered` (`prof`):** gains per-node `crit`, `r_minus_1`, `Delta`, `medium_status`; and at the
root `G0bare_star`, `crit_star`, plus the two path integrals.

- `prof.crit_star = prof.r_star + J0eff*prof.G0bare_star`. This needs the root call at `:243` to stop
  discarding `Gbk` (currently `~` in output position 6) — a one-token change.
- `prof.int_Sigma0` and `prof.int_r_minus_1` use the **same first-panel seeding as `h0`** (`:152`,
  `cumtrapz([0 hgrid],[r0n rv])`): `int_Sigma0` seeded with `Sigma0_pm0`, `int_r_minus_1` with `r0n-1`.
  Both are free — no extra solve — and must be recomputed after extension or redensification exactly
  when `h0` is.
- `prof.status` gains `'medium_out_of_domain'` and `'degenerate_doublet'`, both distinct from
  `'node_failed'` and `'unresolved'`.

> **Status precedence must be specified, or the new statuses are unreachable.** An out-of-domain node is
> also a non-accepted node, so the existing `if any(~cnv)` check (`:211-214`) would claim it as
> `'node_failed'` first. Resolution order, most specific first: `'degenerate_doublet'` ->
> `'medium_out_of_domain'` -> `'node_failed'` -> `'unresolved'` -> `'ok'`. The status is resolved from the
> recorded per-node reasons, not from the order the checks appear in the file, and a test pins each
> status as individually reachable.

**`invz_solve_point_ordered` (`pt`):** gains `pt.stable_1z`, `pt.crit_1z` (`= prof.crit_star`),
`pt.static_medium`, `pt.medium_margin`, `pt.omit_mu3`, `pt.omit_cubic`, `pt.omit_max`. `pt.converged`
keeps its existing meaning (`info.accepted`, the consistency tier), so no existing consumer shifts
meaning under it. The production ordered label requires **both** `pt.converged` and `pt.stable_1z`.

**`invz_solve_point` (PM leg, `pt`):** gains `pt.static_medium`, `pt.Jmom`, `pt.medium_margin`.

**`invz_spectra_map` (`S`):** gains `S.static_medium` (scalar provenance) and per-column `S.crit_pm`,
`S.pm_probe_status`, `S.stability_1z`, `S.phase_1z_reason`. `invz_spectra_qpath` gains scheme provenance
only, unless it is separately upgraded to the PM-mass-first dispatcher; its current auto/bare-ordered
flow cannot honestly populate the ordered-stability fields.

#### The three-way dispatcher

1. converged finite PM with `crit_pm > crit_tol` -> **PM**;
2. converged finite PM with `crit_pm < -crit_tol` -> ordered solve is **phase-eligible**;
3. PM non-convergence / non-finite / recoverable error, or `|crit_pm| <= crit_tol` ->
   **unknown / boundary-indeterminate**.

Case 3 may run the ordered solver for diagnostics but **cannot emit `phase_1z = 1`** without a
separately validated free-energy/branch-selection rule. This is what prevents solver availability from
becoming the phase criterion. Preserve scalar `S.Bc_1z` only for a valid ordered/stable-PM bracket; add
`S.Bc_1z_interval` and `S.Bc_1z_status` so unknown points widen or invalidate the bracket rather than
being interpolated through.

`S.phase_1z_reason` is categorical, one value per column:
`'ordered'` | `'pm'` | `'unstable_endpoint'` | `'medium_out_of_domain'` | `'degenerate_doublet'` |
`'solver_failed'` | `'pm_probe_unknown'` | `'boundary_indeterminate'` |
`'not_attempted_longitudinal'` | `'bare_not_ordered'` | `'bare_escape_hatch'`.

`S.pm_probe_status` uses the narrower vocabulary `stable` | `unstable` | `boundary_band` |
`nonconverged` | `nonfinite` | `recoverable_error` | `not_attempted`, with a parallel error-id field.
Fatal/unclassified errors never become a status, because they rethrow.

> **The `ordered_1z = 'bare'` escape hatch needs an explicit rule.** `invz_spectra_map.m:340` sets
> `phase_1z = 1` from the *bare* solve with no PM verdict at all (the documented Stage-1 diagnostic
> escape hatch, also taken under a longitudinal tilt). Left unstated, the stopgap silently violates the
> new labelling contract. Rule: the escape hatch keeps `phase_1z = 1` — it is explicitly a bare
> diagnostic, not a 1/z claim — but must carry `S.phase_1z_reason = 'bare_escape_hatch'` and
> `S.static_medium = 'n/a'`, so no consumer can mistake it for a strict-order ordered label.

#### Sequencing constraint on the dispatcher and on un-nesting

The three-way dispatcher, and the `one_field` reorganisation that lifts the 1/z leg out of the
auto/overlay `phase == 1` branch, both change `phase_1z` **on the default path**. Two consequences:

1. They conflict with G9's "every pre-existing numerical field bitwise under `'resummed'`".
2. Worse than formal: under `'resummed'` the PM probe below `Bc_1z` fails for the same pole reason, so
   the three-way rule would classify it *unknown* -> no ordered label -> **more** masking than today.

**Therefore both changes are gated on a strict scheme being active.** Under `'resummed'` the historical
dispatch and nesting are preserved exactly. This is coherent because the three-way rule is only
*meaningful* once `crit_pm` is resolvable below `Bc_1z` — which is precisely what decision 0.2 buys.
The `one_field` separation itself is a scope decision (§7 item 4).

When active, `one_field` splits into two genuinely independent legs:

- the existing `invz_solve_auto` result owns only the auto/RPA overlay and `phase`;
- for every transverse field the 1/z leg runs its PM probe regardless of auto/overlay success, then
  applies the three-way verdict and calls Jensen only when phase-eligible (or diagnostically when
  unknown);
- longitudinal fields retain their documented separate route.

This is a control-flow separation, not a fallback: the ordered 1/z result never substitutes for a failed
overlay, and an overlay result never votes on `phase_1z`.

> **Compatibility decision, deliberate:** `phase_1z` keeps its existing `{0,1,2}` enum. Plotting code
> consumes it, and widening it to carry "converged but unstable" would silently change every existing
> figure and saved result. New information goes in the parallel arrays above. A column that is
> consistent but fails the stability tier stays `phase_1z = 0` (masked) with
> `S.phase_1z_reason = 'unstable_endpoint'` — masked, but no longer *unexplained*.

### 4.5 What does not change

- **The v5 coupling cache and its `cacheMeta`.** Moments are derived at call time from `Jnu_flat`, so no
  cache key, prefix or schema change, and Step-5's exact-`isequaln(cacheMeta)` validation, the
  `jq5_<backend>_` separation and the collision/genuine-read gates are all untouched.
- **`invz_jq_modes` / `invz_bz_couplings` / `invz_jq_path` and the Ewald backend.** Out of scope.
  Production dipole default stays `bruteforce` (that flip is a separate step).
- **`K(2:end)` in both legs**, `invz_emt_scalar`'s vector and `[nJ,nw]` branches, and its `opts.debug`
  closure diagnostic.
- **`pt.crit`'s historical ordinary-Dyson definition** (`invz_solve_point_ordered.m:60-62`, contract
  P2-G). It stays a legacy diagnostic; it is not the ordered pole mass.

**Result-level provenance does change.** Every returned point/map carries the canonical scheme and
reference convention, and any project-owned save/cache key includes them. Do **not** introduce
`invz:staticMediumConflict` until an API actually accepts a precomputed solver result for reuse — there
is no current ingestion path on which it could truthfully fire (the spectra drivers ingest precomputed
*couplings*, not precomputed solver results). If such a path is added later, scheme mismatch is a hard
conflict.

### 4.6 Call graph

```
drivers (run_spectra / spectra_map / spectra_qpath / run_phase_diagram)
  |  opts.static_medium -> invz_check_static_medium -> stamps eopts + eso
  |
  +-- invz_solve_auto                                  [independent auto/RPA overlay]
  |
  +-- invz_solve_point                                 [1/z PM probe, every transverse field]
  |     |  mom (once/point) via opts.Jmom
  |     +-- invz_emt_scalar   (omega = 0 slot only)
  |           +-- invz_static_medium_reference
  |           +-- invz_medium_moment_closure
  |
  +-- invz_solve_point_ordered                         [1/z ordered leg, Jensen when eligible]
        |  mom (once/point)
        +-- invz_hmf_ordered      mom -> node.Jmom for every node
              |
              +-- eval_node
              |     +-- invz_twolevel_ordered          [Delta domain mode, §5.3]
              |
              +-- invz_ordered_node_solve
                    +-- invz_emt_scalar                (K(1) discarded by the ordered leg)
                    +-- invz_emt_static_ordered
                    |     +-- invz_static_medium_reference
                    |     +-- invz_medium_moment_closure
                    |     +-- invz_gstat_ordered       (r arrangement fix, §1)
                    +-- invz_ordered_residual
                          -> res.accepted (tier 1) + res.stability (tier 2)
```

The `K(1)` replacement ordering is **already correct** at all three sites —
`invz_ordered_node_solve.m:190-192` (`K(1) = K0s` before `invz_lambdas`),
`invz_ordered_residual.m:203-204` (same), and block D excluding `K(1)` by design. G13 pins that existing
behaviour; it is not fixing a defect.

---

## 5. Error handling

### 5.1 Three-way taxonomy, and one shared classifier

| class | examples | behaviour |
|---|---|---|
| **wiring / programming** | unknown or conflicting scheme (`invz:staticMedium`), missing strict-mode `node.Jmom`, `Gref` not supplied, `Jnu_flat` shape mismatch (`invz:emtJnu`), `hz_fixed` not held (`invz:hzFixed`) | **throw, loudly, and do not absorb** |
| **physics non-convergence** | outer Picard `max_iter`, block A/C/D residual fail | return non-accepted; never throw (existing policy) |
| **out of domain** (new class) | reference denominator non-positive or inside its frozen zero margin; the sampled path brackets a *reference* denominator crossing; `Delta < 1e-4` meV | return a **distinct status**; neither an error nor a convergence failure |

**The narrowing, and why one site is not enough.** `invz_ordered_node_solve.m:213-216` absorbs *every*
`invz:*` identifier. Its docstring establishes why that was safe: across its whole chain
(`invz_emt_scalar` / `invz_emt_static_ordered` / `invz_gstat_ordered` / `invz_lambdas` /
`invz_sigma_ordered` / `invz_sigma`) the only `error()` site is `invz:emtJnu`. That premise fails once
the strict scheme adds throw sites in the same chain — a wiring error would be silently downgraded to
"node not accepted", i.e. a masked column. **That is exactly the failure mode that let the original
defect hide for a whole stage.** And there are at least four such absorbers on the path:
`invz_ordered_residual.safe_eval:178-192`, `invz_spectra_map:297-300` and `:325-328`, and
`invz_solve_auto`. Narrowing one merely relocates the swallow.

Use one shared predicate, `invz_is_recoverable_solver_error(id)`, with an **explicit whitelist**;
unclassified ids rethrow by default. Apply it at *every* catch reached by the strict path. Unit tests
inject one fatal and one recoverable id at every layer (G15).

The initial whitelist is limited to the existing branch/domain signals `invz:orderedPhase` and
`invz:degenerateDoublet`; strict-medium domain outcomes return statuses instead of throwing. The inner
node map and residual checker therefore have **no expected recoverable throw**: remove their broad catch
or rethrow every unclassified id. Adding an identifier to the whitelist is a reviewed contract change,
not a convenience response to a failing run.

> **Pre-execution check.** Removing the broad catch changes existing behaviour for `invz:emtJnu`
> (currently absorbed into a non-accepted node; afterwards it throws). Per `invz_emt_scalar`'s docstring
> no vector-`Jnu` fixture exercises that path, so exposure should be nil — but this must be *verified
> against the suite*, not assumed.

### 5.2 The out-of-domain path, and what must not be done to it

- `invz_static_medium_reference` never throws on a physical/domain denominator event. It returns
  `Gref = NaN` plus denominator metadata and a non-`ok` status. The medium leaf is not called on that
  node.
- `invz_emt_static_ordered` propagates it as `out.medium_status`; `invz_emt_scalar` as a field on `med`.
- `invz_ordered_node_solve` and the PM loop stop the attempt before lambdas or `Sigma` consume an
  invalid reference.
- `invz_ordered_residual` sets `res.blockB.pass = false` when the status is not `'ok'`, so the node is
  not accepted.
- `invz_hmf_ordered` sets `prof.status = 'medium_out_of_domain'` — **distinct from `'node_failed'`**,
  because the two demand different responses: one says the iteration did not settle, the other says the
  truncation has no value there. `prof.medium_status` records which nodes and the margin.
- **Two denominator classes, deliberately treated differently.** The **reference** denominator
  `1 + Sigma0` has no cancellation available and is a hard out-of-domain event. A **local** `Gstat`
  denominator crossing is *removable in the integrand* (§1): `r` and `crit` stay continuous through it,
  so it is recorded and its regularity verified (G17), not treated as invalidating the node.
- **Forbidden responses**, per the durable gotcha: broadening or regularising the static response,
  adding a pole regulariser, flipping a sign, or widening a tolerance to cross a boundary. The result
  fails promotion; any revised theory starts a new design/preregistration cycle (§6.0).

### 5.3 The `invz:degenerateDoublet` domain contract

The jensen leg's domain requires `Delta(Bx, h) >= 1e-4` meV at every node it evaluates (mechanism and
current failure mode in §3).

1. **Single-evaluation domain mode.** Extend `invz_twolevel_ordered` with an internal
   `domain_policy = 'return'` mode; the default stays the current throwing behaviour for compatibility.
   In return mode it performs the existing electronic diagonalization **once**, records `Delta` and
   `valid = false`, and returns before evaluating `g0` or any division that assumes a resolved doublet.
   `eval_node` records `prof.Delta` and maps invalidity to `degenerate_doublet`. Do **not** pre-screen by
   duplicating the same diagonalization and then calling the constructor again.
2. **No broad catch.** The default constructor throw remains a fatal signal if it unexpectedly escapes a
   return-mode call — that is a wiring defect, not a second recoverable path.
3. **`B = 0` is a documented hard exclusion, labelled, not silently proxied.** The driver marks the
   column `S.phase_1z_reason = 'degenerate_doublet'` and moves on. It does **not** substitute a
   small-field proxy: `README.html:208` already establishes that at exactly `B = 0` the solve throws
   *deliberately*, that 0.05 T is the practical proxy for a different (Tier-1/ODD) route, and that
   near-degenerate results are "directional only, not quantitative". Substituting a proxy inside the
   spectra driver would bury that caveat inside a production figure. Choosing the field grid is the
   caller's decision.
4. **The floor is not a new number.** `1e-4` meV is the constructor's own existing threshold, reused.

### 5.4 Warning policy

Unchanged where it exists: `eso.warn = false` inside sweep loops (`invz_hmf_ordered.m:57`,
`invz_solve_point_ordered.m:188`) so the per-node console never floods; `invz:hmfUnresolved` retained for
bracket and refinement failures. Strict mode adds **no** per-node warning — there is no closure iteration
left to warn about. One driver-level summary warning per sweep (`invz:staticMediumDomain`) reports how
many columns were out-of-domain or degenerate, with counts.

### 5.5 No silent caps

If Gate 0's scan skips nodes, or the profile drops any, the count is logged. A truncated scan reported as
a clean scan reads as "covered everything" when it did not — the discipline Step 5 applied to coverage.

---

## 6. Gate matrix

### 6.0 Preregistration — frozen before the first strict-mode run

1. `crit_tol`, `D_tol`, `Dq_tol` — the endpoint stability tolerances, boundary-scaled.
2. `K_atol`, `K_rtol`, the reference-denominator margin, and the rule for classifying
   `|crit_pm| <= crit_tol` as boundary-indeterminate.
3. **The Gate-0 negative-outcome rule.** Promotion fails if any required solved-path node has an invalid
   **reference** denominator, any skipped node is unaccounted for, or `max(omit_max)` exceeds the control
   bound. A **local** `Gstat` denominator crossing does not by itself fail promotion, provided G17 shows
   `r` and `crit` remain finite and continuous through it; a crossing where they do not is a failure.
   **(e)** Promotion also fails if any required ordered field does not return `status = 'ok'`, a finite
   nonzero root `hstar > 0`, and an endpoint that passes the frozen `crit_tol`/`D_tol`/`Dq_tol` margins; or
   if either designated PM control field fails to return a converged, finite, positive-mass PM state.
   Without (e), conditions (a)–(d) could all pass while every field remained masked, and a Gate 0 that
   reports PASS without the ordered leg producing a stable root cannot serve its one stated purpose —
   deciding whether stage 4 is worth planning. Exact `B = 0` is a labelled hard-domain control, not a
   required ordered-path node, so it does not fire (e). On failure the run stops at diagnosis. Carrying
   another moment, changing `Gref`, or truncating other Matsubara sectors is a **new theory candidate
   requiring a new spec and fresh preregistration** — never an in-run fallback. Regularisation,
   broadening and tolerance widening remain forbidden.
4. The `omit_mu3`/`omit_cubic` reporting threshold and the stricter `omit_max` promotion bound. Neither
   affects path-node algebraic acceptance nor switches schemes, but the promotion bound is a real gate.
5. The strict-vs-resummed control tolerances for G8, quantity by quantity: one norm cannot cover a
   critical field, a dimensionless mass, a self-energy and a spectrum.
6. Numerical convergence tolerances for G1/G5/G14 and the complete field/HMF/T/Matsubara/q grids.
7. The `Delta` floor: `1e-4` meV, inherited from `invz_twolevel_ordered`, not re-derived.

Until these entries carry numbers or formulas rather than adjectives such as "boundary-scaled" or
"negligible", the design is reviewable but **not executable**.

### 6.1 The gates

| id | question | method | pass | blocking for |
|---|---|---|---|---|
| **G0** | does the solved path stay in a controlled strict-order domain? | (a) bare proxy on the identical proposed grid; (b) actual solved-path **reference** and **local** denominator margins, kept as separate classes, plus coverage counters, `omit_mu3`, `omit_cubic` | every required node accounted for; no invalid *reference* denominator; local crossings all satisfy G17; `max(omit_max)` below the frozen bound; every required ordered field returns `status = 'ok'` with a stable finite root; both designated PM controls return a converged, finite, positive-mass PM state | default flip |
| **G1** | are the differential identities implemented? | panelwise `dm/dh = -G0bare`; matched-panel `Delta F/Delta h` against the trapezoidal average of `crit`; `dF/dm = crit/chi_path`; repeated under `nH` refinement | all three errors converge at the preregistered quadrature order and fall below frozen tolerances | build |
| **G2** | do the legs coincide at `m -> 0`? | `slope_pred = crit_pm`; `Gref_ord = Gref_PM`; `K(1)_ordered = K(1)_PM` at a non-degenerate zero-moment fixture under the same public scheme, through the full caller wiring | within frozen `K`/mass tolerances; bitwise only where the arithmetic path is intentionally identical | build |
| **G3** | does the pinned identity survive? | existing `r = 1 + Sigma0` at `m = 0` test, under strict mode | stays green (§1: holds for any `K0`) | build |
| **G4** | are the two tiers separate? | endpoint requires `crit`/`D_uni`/`Dq_min` above frozen margins; **negative** fixtures assert residual-consistent intermediate nodes with negative stability masses remain admitted | endpoint rejection and intermediate-node admission both pass | build |
| **G5** | are the path corrections measured reliably? | `prof.int_Sigma0` and `prof.int_r_minus_1` at several fields, repeated under profile refinement | both finite and quadrature-converged; magnitude reported, and explicitly **not** bounded by the 0.3% boundary shift | default flip |
| **G6** | are multiple ordered roots selected consistently? | if multiple locally stable roots occur, compare `Phi_path` on the same path and record all roots | selected root has the lowest converged `Phi_path`; vacuous if only one root | default flip |
| **G6d** | what does the existing thermodynamic diagnostic say? | run `invz_deltaF_ordered` at a common cutoff; retain the closed-2x2 fingerprint | labelled **partial/non-gating** exactly as the helper contract requires; **no** "two routes agree" production claim | diagnostic only |
| **G7** | does changing only `wn = 0` create a material scheme sensitivity? | `omega_n != 0` is unchanged by construction, so `K(2)_strict == K(2)_resummed` and the scheme jump at `omega = 0` is **exactly** `K(1)_strict - K(1)_resummed`. Compare that shift against the physical dispersion `K(2) - K(1)`, vs `T`; optionally add an explicitly labelled all-frequency strict diagnostic | converged measurement produced, with artifact and dispersion separated; no assumption that either is negligible | diagnostic only |
| **G8** | is strict order controlled where resummed is also valid? | common stable states only: compare `K0`, `Sigma0`, mass, integrated path quantities, `Bc_1z`, spectra, with quantity-specific metrics | every promotion metric within its frozen bound; otherwise the candidate is reported uncontrolled | default flip |
| **G9** | is the legacy path preserved? | `'resummed'` default/off path against frozen references, **including the historical dispatch and nesting**, which the §4.4 sequencing constraint keeps inactive under default | every pre-existing numerical field bitwise; new provenance/status fields checked separately | build |
| **G10** | are recomputed outputs labelled? | recompute `Bc_1z`, PM anchors, spectra; every point/map/save key carries canonical scheme and reference convention | no unlabeled result, no project-owned result-cache collision across schemes | default flip |
| **G11** | do fixtures exercise the real thing? | retain synthetic algebra fixtures; add real-coupling anchors from one fully explicit provenance tuple, pinning both its coupling hash and recomputed moments; use the static column and full multiset correctly | synthetic + real anchors pass; hard-coded moments accepted **only** with matching full provenance (the §B table is valid solely for grid `[16 16 16]`, `dpRng 30`, `bruteforce`, no grid-policy fields) | build |
| **G12** | are phase decisions and masks honest? | per-column `phase_1z_reason`, `pm_probe_status`, `crit_pm`; inject stable, unstable, boundary-band, PM-failure, auto-overlay-failure, domain, longitudinal and bare-escape-hatch cases | no `phase_1z = 0` without a reason; no unknown PM probe votes ordered; auto failure does not suppress an otherwise valid 1/z leg; boundary estimate exposes unknown intervals | default flip |
| **G13** | can the PM slot leak into ordered lambdas? | behavioural sentinel: make PM slot 1 distinguishable from ordered `K0s` and verify recomputed lambdas use the latter. **Not** a source-text-order assertion | sentinel never leaks at the live loop, the independent residual map, or the final refresh | build |
| **G14** | is it numerically converged? | forward/reverse continuation in `B`; refine field grid, `nH`, root tolerance, `Ecut`/Matsubara, q-grid; count all skipped/invalid nodes | every preregistered observable stable within its frozen tolerance, coverage complete | default flip |
| **G15** | does the error taxonomy hold end-to-end? | inject fatal and recoverable ids through node solver, residual checker, auto solver and spectra dispatcher; separately inject domain statuses; verify the `invz:emtJnu` behaviour change against the suite | fatal ids always escape; recoverable/domain outcomes retain their exact category and never become generic `node_failed` | build |
| **G16** | is the original spectra defect actually fixed? | preregistered field sweep across both sides of `Bc_1z`, excluding only declared hard-domain/boundary-band points; inspect phase coverage, finite spectra, `D_ord -> 0+`, `crit_pm -> 0+`, tracking the same denominator pole/critical branch rather than the brightest broadened pixel | no unexplained or solver-failed ordered columns in scope; ordered and PM masses close at one boundary; the pole-tracked critical mode softens continuously within resolution | default flip |
| **G17** | is the removable singularity removable in code? | synthetic states driving `Gstat` through and past its denominator zero (including exactly `Inf`): assert `r`, `Gtil0` and `crit` finite and continuous, matching the limits `-G0bare*K0`, `-1/K0`, `G0bare*(J0eff-K0)`; assert the as-written arrangement is gone | finite and continuous through the crossing at every tested magnitude | build |

"Blocking for build" = must pass before the phase is called complete. "Blocking for default flip" = must
pass before `'resummed'` stops being the default — which is **not** in this phase's scope.

### 6.2 Sequencing

Mirrors the Ewald Step-4/5/7 discipline this project has already run twice.

1. **Additive primitives + their own gates.** `invz_coupling_moments`,
   `invz_static_medium_reference`, `invz_medium_moment_closure`, `invz_check_static_medium`,
   `invz_static_domain_scan`, `invz_is_recoverable_solver_error`. Nothing wired. G9 holds trivially.
2. **Opt-in wiring, default `'resummed'`.** Both leg leaves, the `invz_gstat_ordered` arrangement fix,
   the `node.Jmom` schema bump, the residual two-tier split, the contract-doc revision. The three-way
   dispatcher and any `one_field` separation land here **gated on a strict scheme being active** (§4.4),
   so the default path is untouched. Gates G1–G4, G9, G11, G13, G15, G17.
3. **Gate 0 and the measurements.** G0, G5, G7 — diagnose-first, before any claim that the panel is
   fixed.
4. **End-to-end + provenance.** G6/G6d, G8, G10, G12, G14, G16.
5. **Default flip** — a separate, later decision, out of this phase's scope.

Stage 3 can falsify the candidate but cannot demonstrate the masking fix. The fix is demonstrated only
when stage 4 clears G16 end-to-end; until then the strict construction is a well-founded *hypothesis*.

### 6.3 Cost

The static sector loses an iteration (up to `maxit = 200` `invz_gstat_ordered` calls per outer iteration
collapse to one), so strict mode should be *cheaper* than resummed per node, and the moment computation
is amortised per resolved point / HMF call. Gate 0's proxy scan is bare diagonalizations only. The
expensive new items are G8 (both schemes over a field sweep) and G14 (convergence sweeps), both one-off.

---

## 7. Open items — settle these when approving for execution

1. **§6.0 still lacks numerical values and formulas.** They must be chosen and frozen as part of
   approving this spec for execution, not during it. Until then the design is approved but
   execution-blocked.
2. **`invz_static_domain_scan`'s grid ownership.** A standalone scanner must accept an explicit `hgrid`;
   it must not independently recreate HMF's adaptive extension/redensification logic. Factor the initial
   geometric grid into a small shared helper for prospective scans, and evaluate solved-path margins
   directly on `prof.hgrid` after adaptation. "Grid identity" then means shared initial input plus
   complete accounting of every actual solved node — not two duplicated algorithms that happen to agree
   in one test.
3. **G2 floating-point equality.** The algebra is exact at `m = 0`, but the two callers need not reach it
   through bitwise-identical arithmetic (`G0_PM(0)/(1+Sigma0)` vs `G0bare0/(1+Sigma0)`, equal only
   because `G0el0 -> 0`). Use the frozen `K_atol`/`K_rtol` unless a dedicated fixture proves intentional
   bitwise identity. Do not leave this as an execution-time choice.
4. **Whether the `one_field` un-nesting is in scope for this phase.** Lifting the 1/z leg out of the
   auto/overlay `phase == 1` branch is the largest structural change in this spec, and it is **not** on
   the critical path for the masking fix: below `Bc_1z` the bare set orders, so `phase == 1` holds, and
   there is already a separate transverse-PM branch for `phase == 2` (`invz_spectra_map.m:372+`). It is a
   genuine robustness improvement — an overlay failure should not pre-mask the theory being repaired —
   but it could be deferred to its own change without weakening G16. Recommend: **defer**, and keep only
   the three-way dispatcher (which *is* required, since it is what stops solver availability from
   labelling a phase).
5. **Whether the pre-existing `[nJ,nw]` static-flattening defect (§4.3) is fixed here or split out.**
   Recommend: reject the ordered-retarded combination in the resolver for this phase, and split the fix.
