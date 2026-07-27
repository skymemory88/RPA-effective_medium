# DRAFT preregistration amendment — proposed §11

**Status: DRAFT. NOT APPLIED. Requires user approval before anything is changed.**

Prereg §9 amendment rule: *"Changing one afterwards requires a new dated re-registration section,
never an edit in place."* Prereg §8: *"Changing these grids after seeing strict output requires a
dated preregistration amendment."*

This document is written to be **appended verbatim** to
`docs/invzp_strict_medium_prereg.md` as a new §11 if approved. The preregistration has not been
edited and must not be until then.

---

## What triggered this

Task 17 (`a923691`) ran the G1/G2/G13 gates for the first time. Three gate *definitions* proved
defective — one is algebraically redundant, two assert invariants the code does not promise. The
defects were found by measurement, and every relevant constant is frozen, so correcting them
requires this amendment rather than a test edit.

Evidence: `Task-17_failure.md` §4, §5, §6 and its verification log §8.

## What this amendment does NOT change

Stated first, because it is the part that matters most.

| frozen item | status |
|---|---|
| §3 Gate-0 negative-outcome predicate (a)–(e) | **UNCHANGED** |
| §4 `omit_report = 0.05`, `omit_promote = 0.10` | **UNCHANGED** |
| §2 `K_atol = 1e-14`, `K_rtol = 1e-12` (values) | **UNCHANGED** |
| §8 coupling fixture, digest, required field set, `nH = [33 65 129]` | **UNCHANGED** |
| §3 prohibition on regularisation, broadening, tolerance widening | **UNCHANGED** |

**The Gate-0 verdict is not reopened.** `strict_1z_dyson_ref` fails clauses (a), (c) and (e) on the
1 T real-coupling anchor. Nothing below bears on that result, and nothing below may be cited in
support of promotion. Every item here concerns *identity and contract gates*, none of which is an
input to the §3 predicate.

No numeric threshold proposed below was selected by inspecting a failing value. The one new number
(A1.2) is an a-priori midpoint discriminator, justified independently of the data and stated as such.

---

## A1 — §5 G1 identity gates: ALL THREE forms are one identity, and the threshold tests the grid

> **Revision 2026-07-27, after external review (`invzp_fix_attmpt_review_Codex.md` §1).** A1.1 below
> originally claimed only that form 2 is dependent, and proposed retaining forms 1 and 3 as "the two
> independent G1 identity checks". **That proposal is invalid — form 3 is dependent too.** With
> `a = Δm/Δh`, `g = midpoint-avg(G0bare)`, `ρ = midpoint-avg(r)`, `J = J0eff`:
>
> ```
> form 1  =  a + g
> form 2  =  ΔF/Δh − crit_avg           =  −J·(a + g)          =  −J · form 1
> form 3  =  ΔF/Δm − crit_avg/chi_path  =  ρ·(a + g)/(a·g)     =  ρ/(a·g) · form 1
> ```
>
> Since `chi_path = −G0bare > 0` gives `g ≠ 0`, and `a > 0` along the profile, forms 2 and 3 vanish
> exactly when form 1 does. **Form 1 is the sole load-bearing identity; forms 2 and 3 are derived
> consistency restatements.** Everywhere below, read "form 2 is DERIVED" as "forms 2 **and 3** are
> DERIVED". Registering §5 as three independent checks overstates its coverage by a factor of three,
> not two. A genuinely independent G1 check would have to be a local or thermodynamic test, not
> another rearrangement of the same secant identity; **whether to add one is an open decision for the
> owner**, not something this amendment resolves.

### A1.1 The second form is algebraically the first

Current §5 lists three forms as though independent:

```
|dm/dh + G0bare|                        <= 1e-6 * max(1, |G0bare|)
|Delta F/Delta h - trapz_avg(crit)|     <= 1e-6 * max(1, |crit|)
|dF/dm - crit/chi_path|                 <= 1e-6 * max(1, |crit/chi_path|)
```

`invz_hmf_ordered.m:397-401` builds `h0 = cumtrapz(hgrid, r)` and `F = h0 - J0eff*m`; `:539` sets
`crit = r + J0eff*G0bare`. Because `h0` is a cumulative trapezoid of `r`, the secant
`diff(h0)./diff(h)` is *exactly* the midpoint average `r_avg`. Therefore

```
dF/dh    = r_avg - J0eff * dm/dh
crit_avg = r_avg + J0eff * G0bare_avg
---------------------------------------------------
dF/dh - crit_avg  ==  -J0eff * (dm/dh + G0bare_avg)      exact, to reassociation
```

Confirmed numerically at `5.7e-15`. **Form 2 is `J0eff` times form 1.** It cannot fail unless form 1
fails and carries no independent information.

**Proposed:** retain form 1 and form 3 as the two independent G1 identity checks. Retain form 2
explicitly labelled **DERIVED** — asserted, but recorded as a consistency restatement of form 1
rather than a third gate. Registering it as independent overstates the coverage of §5.

### A1.2 The absolute threshold tests the grid, not the identity

All three forms compare **secant/trapezoidal** quantities against **differential** identities, so
their residual is O(h²) panel error. On the frozen `nH = [33 65 129]` ladder:

| nH | form 1 | form 2 | form 3 |
|---:|---:|---:|---:|
| 33 | 1.8216829803e-02 | 5.8431222186e-03 | 2.2609487299e-04 |
| 65 | 4.8713648324e-03 | 1.4726351883e-03 | 6.4840455639e-05 |
| 129 | 1.2528393091e-03 | 3.6828822520e-04 | 1.7383399960e-05 |

Reduction ratios: form 1 `3.7396, 3.8882`; form 2 `3.9678, 3.9986`; form 3 `3.4869, 3.7301`.

The identity holds — the residual falls at second order. The `1e-6` threshold is not mathematically
impossible; it is **incompatible with the frozen `nH = 129` grid**, requiring roughly
`sqrt(368) × 129 ≈ 2500` nodes at second order. Since both the threshold (§5) and the grid (§8) are
frozen, no dial can be turned without this amendment.

**Proposed replacement for the profile identities.** Assert the *convergence order*, not an absolute
floor:

```
G1 profile identities (forms 1 and 3), evaluated on the frozen nH ladder [33 65 129]:
   the max-norm residual must be MONOTONICALLY DECREASING across the ladder, and
   each successive reduction ratio must satisfy  ratio >= 3.
   Form 2 is DERIVED from form 1 (identity above) and is asserted as a restatement.
   The existing synthetic smooth nested-grid quadrature fixture is RETAINED unchanged
   as the control that the quadrature machinery itself is second order.
```

**Justification of `ratio >= 3`, with its status stated plainly:** first-order convergence gives a
ratio of 2, second-order gives 4. Three is the midpoint and discriminates between them. The
*reasoning* is structural, not fitted. **But the criterion was selected AFTER the failing run was
observed, and must not be described as blind or a priori.** It is being introduced through this dated,
transparent amendment precisely because it is post-hoc. (For the record, the measured minimum is
`3.4869`; the choice is not tight against it, and had the data been fitted to, a value nearer `3.4`
would have been picked. That is evidence of good faith, not a substitute for the disclosure.)

### Why a ratio test is admissible here, and the guard that keeps it so

Prereg §5 currently warns: *"do not require a [3,5] error ratio from non-nested adaptive grids."*
That caveat is **honoured, not overridden** — it is conditional on non-nesting, and the condition
must be checked rather than assumed.

*Structural argument.* `invz_hmf_grid` builds `hgrid = hmax * ratio.^((nH-1):-1:0)` with
`ratio = hmin_frac^(1/(nH-1))`. With `hmax` and `hmin_frac` fixed, `ratio_65 = ratio_33^(1/2)`, so
the `33` nodes are a subset of the `65`, which are a subset of the `129`; the largest panel scales
as `1/(nH-1)`, and `32 → 64 → 128` is an exact halving.

*But that describes the INITIAL grid.* G1 runs against `prof.hgrid`, which is the grid **after**
`invz_hmf_ordered`'s adaptive extension and redensification. If adaptation fires, nesting is broken
and the prereg's caveat applies. Measured on the G1 fixture (`T = 0.31 K`, `Bx = [2.85 0 0]`):

```
nH= 33  status=ok  n_extend=0  redensified=0  numel(hgrid)=33
nH= 65  status=ok  n_extend=0  redensified=0  numel(hgrid)=65
nH=129  status=ok  n_extend=0  redensified=0  numel(hgrid)=129

nesting of the ACTUAL adapted grids (tolerance 1e-14 relative):
  33 nodes found in  65 grid:  33/33   true
  33 nodes found in 129 grid:  33/33   true
  65 nodes found in 129 grid:  65/65   true
```

Adaptation did not fire, and the grids the tests actually use are exactly nested. This is measured
at one fixture only; at other fields adaptation may fire.

**Therefore the ratio assertion is registered with an explicit guard:**

```
The ratio >= 3 assertion applies only when adaptation did NOT fire on any rung of the
ladder, i.e. prof.n_extend == 0 and prof.redensified == false for all three nH. That
condition must be ASSERTED and RECORDED, not assumed. If adaptation fired on any rung,
the ladder is not nested, prereg SS5's existing caveat governs, and only the monotonic-
decrease requirement applies -- with the adaptation status reported alongside.
```

**Alternative considered and not proposed:** a Richardson-extrapolated continuum test retaining the
`1e-6` target. It is applicable (the ladder is nested), but the extrapolated values are not uniform
across the three forms, and selecting a threshold against them would be fitting to the failing run.
If the continuum target is judged worth keeping, the extrapolation rule and its threshold should be
registered *before* the amended test is run.

**Note for implementation:** relabelling form 1 as a "local derivative test" while retaining `1e-6`
does not solve the problem. Any finite difference on the profile inherits panel error; reaching
`1e-6` would need Richardson extrapolation or a much finer local stencil.

---

## A2 — §2/§7.3 G2: bound cross-solver comparisons by the solver tolerance, not the block-B gate

### The defect

§2 defines `gate = K_atol + K_rtol*max(|K0|, |Kstrict|, Jscale)` and states its purpose: *"A one-shot
formula recomputed from the exported state should agree to floating-point reassociation only; this
gate catches mis-wiring, not physics."* It is a **block-B recomputation** gate, scoped to one
exported state.

G2 as implemented applies it to `|p.K0_pm0 - ptp.K(1)|` — two **independently terminated outer
iterations** (the ordered predictor's node solve and the PM leg's converged `invz_emt_scalar`).
Measured `1.021112e-12` against a gate of `1.600000e-14`. That is `3.4e-10` relative: cross-solver
iteration residue, which this gate was never designed to bound.

The in-scope comparisons all pass, two of them better than the design promises. §7.3 states that
bitwise identity *"is not a property the design guarantees"*; it held anyway:

```
|p.slope0 - ptp.crit|   = 0.000000e+00
|GrefO - GrefP|/|GrefP| = 0.000000e+00
|KO - KP|               = 0.000000e+00     (bitwise)
```

### Proposed

```
G2 is split into three registered assertions:
  (i)   CLOSURE IDENTITY -- the strict closure recomputed from each leg's exported state
        (KO vs KP) under K_atol/K_rtol. This is the block-B gate applied in its designed
        scope and is UNCHANGED.
  (ii)  CALLER PROVENANCE -- both legs report the strict scheme and a converged state.
  (iii) CROSS-SOLVER EXPORTED K -- retained as coverage, bounded by a value DERIVED from the
        declared outer-solver tolerance (tol_outer) of the two callers, not by K_atol/K_rtol.
```

**Open item requiring a decision before this clause is registered:** the exact functional form of
the (iii) bound. It must be derived from `tol_outer` by a stated rule, registered before the amended
test runs. It must not be set to the observed `1.02e-12`.

Full-caller coverage is deliberately **not** dropped; only its bound is moved to the tolerance that
actually governs it.

---

## A3 — G13: replace a stale-versus-refreshed equality with a behavioural discriminator

### The defect

G13 asserts `state.lam == invz_lambdas(state.K, …)` at `RelTol 1e-10`; measured `4.705009e-10`.

Verified in source: `state.lam` is computed **inside** the outer loop at
`invz_ordered_node_solve.m:232`, whereas `state.K(1)` is overwritten **after** the loop by the
post-loop static refresh at `:257` (`K(1) = K0s`). The assertion compares a pre-refresh quantity
against a post-refresh one. **Exact agreement is not an invariant the code promises**, so the
residual is structural rather than numerical noise.

The sentinel's actual purpose is satisfied: the ordered `K(1)` differs materially from the PM slot
(`4.711546e-06`, relative `1.734262e-03`), and the counterfactual leaked lambdas sit
`2.974382e-05` away in absolute terms — about `1.7e6` times the observed absolute residual
(`1.771636e-11`).

### Proposed

```
G13 asserts, in place of the RelTol equality:
  (i)   the ordered and PM static slots are distinguishable (RETAINED unchanged);
  (ii)  the observed lambdas are CLOSER to the ordered-slot reconstruction than to the
        PM-slot-leaked reconstruction  -- a scale-free comparative test with no threshold
        to register;
  (iii) info captured from invz_ordered_node_solve reports info.accepted;
  (iv)  the residual/final-refresh route is exercised explicitly rather than discarded.
```

Clause (ii) carries no numeric constant, so nothing here can be fitted to the failing run.

---

## Approval

Registered only when approved and dated by the project owner, and appended to
`docs/invzp_strict_medium_prereg.md` as §11 without editing any existing section.

Until then:

- the two test files stay as committed and the branch stays red (`P=544 F=6 I=29 T=579`);
- Task 18 may proceed with the Gate-0 **FAIL** report, which depends on none of the above;
- the A2(iii) bound and the A1.2 Richardson alternative remain open decisions.

**Signature line (owner):** `Amendment §11 approved and frozen: ____________  date: __________`
