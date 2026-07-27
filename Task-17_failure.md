# Task 17 — gate-test failure record

**Date:** 2026-07-27 (revised 2026-07-27 after external review)
**Branch:** `invzp-stage2c-diagnostic`
**Commit under test:** `a923691` — *test(invzp): add strict-medium identity and contract gates (G1-G3, G11, G13, G15)*
**Previous green commit:** `c525d42` (Task 16), suite `P=540 F=0 I=28 T=568`
**Review incorporated:** Codex review of this record, 2026-07-27 (source file since deleted at the
owner's request). Every correction in that review was verified against source before adoption; §8
below records what was checked, so the audit trail survives the file's removal.

**Later correction (2026-07-27, rev. 2 of the successor record).** §4 below states that G1a and G1b
are the same test and implies G1c is independent. **G1c is dependent too** —
`G1c = ρ·(a+g)/(a·g) = ρ/(a·g) · G1a`, so all three forms vanish together. See
`docs/invzp_convg_fix_attmpt_claude.md` §7 for the derivation and §12 for the verification log.

Task 17 added two test-only files (`test_invz_strict_identities.m`, `test_invz_strict_contracts.m`,
254 insertions, no production file touched). **Seven of its eleven gate tests fail.** The branch is
red and Task 18 has not been dispatched.

---

## 0. Verdict

**Gate 0 FAILS for the preregistered candidate `strict_1z_dyson_ref`, on three independent clauses
of the frozen predicate (prereg §3).**

| clause | prereg text | measured |
|---|---|---|
| **(a)** | any required solved-path node has a non-`ok` REFERENCE denominator status | 4 of 33 nodes report `ref_denom_nonpositive` |
| **(c)** | `max(omit_max)` over the solved path exceeds `omit_promote` | `1.20077` vs frozen `0.10` (attempted path) |
| **(e)** | any required ordered field does not return `status='ok'`, a finite nonzero root, and a stable endpoint | `status = medium_out_of_domain`, `hstar = NaN`, no stable endpoint |

Each is independently load-bearing. The no-root outcome (e) is not a consequence of (c) — the
candidate fails both **operationally** and **asymptotically**.

Prereg §3 continues: *"On failure the run STOPS AT DIAGNOSIS. Carrying another moment, changing
Gref, or truncating other Matsubara sectors is a NEW theory candidate requiring a new spec and fresh
preregistration — never an in-run fallback. Regularisation, broadening and tolerance widening remain
forbidden."*

Accordingly Task 18 may produce a formal **FAIL** diagnostic report. It must not be used to reopen
promotion. Investigating the reference-denominator window is legitimate diagnosis toward a *future*
theory candidate, not an in-run rescue of this one.

---

## 1. Suite state

**Controller-run authoritative full suite at `a923691` (default mode):**

```
SUITE P=544 F=6 I=29 T=579
FAILEDNAME test_invz_strict_identities/test_G1a_dm_dh_equals_minus_G0bare
FAILEDNAME test_invz_strict_identities/test_G1b_F_prime_equals_crit
FAILEDNAME test_invz_strict_identities/test_G1c_dF_dm_equals_crit_over_chi_path
FAILEDNAME test_invz_strict_identities/test_G1d_second_order_convergence
FAILEDNAME test_invz_strict_identities/test_G2_onset_coincidence_at_m_zero
FAILEDNAME test_invz_strict_contracts/test_G13_pm_slot_does_not_leak_into_ordered_lambdas
INCOMPLETECOUNT 29
```

`F=6` in default mode because **G11 is `INVZ_SLOW`-gated and registers as *incomplete*, not failed**.
Under `INVZ_SLOW=1` the count is **7 failures**.

Nothing outside the two new files broke. The incomplete set is the 28 pre-existing entries — all
present by name, none lost — plus G11. `P` rose by exactly the four passing new tests.

**`INVZ_SLOW=1` run of the two files (166.4097 s wall):**

| test | result | duration |
|---|---|---|
| `test_G1a_dm_dh_equals_minus_G0bare` | **FAIL** | 17.564 s |
| `test_G1b_F_prime_equals_crit` | **FAIL** | 16.595 s |
| `test_G1c_dF_dm_equals_crit_over_chi_path` | **FAIL** | 16.572 s |
| `test_G1d_second_order_convergence` | **FAIL** | 30.794 s |
| `test_G1d_nested_trapezoid_is_second_order` | pass | 0.0082 s |
| `test_G2_onset_coincidence_at_m_zero` | **FAIL** | 5.5699 s |
| `test_G3_r_equals_one_plus_sigma_under_strict` | pass | 0.0265 s |
| `test_G11_real_coupling_ordered_anchor` | **FAIL** | 72.438 s |
| `test_G13_pm_slot_does_not_leak_into_ordered_lambdas` | **FAIL** | 0.6014 s |
| `test_G15_fatal_ids_escape_every_layer` | pass | 0.2069 s |
| `test_G15b_domain_outcome_keeps_its_category` | pass | 4.2564 s |

Per file: identities `n=7 P=2 F=5` (87.1292 s); contracts `n=4 P=2 F=2` (77.5025 s).

---

## 2. G11 — the Gate 0 result

Real `bz_couplings` (`grid [16 16 16]`, `dpRng 30`, `dipole 'bruteforce'`, `cache false`),
`T = 0.1 K`, `B = [1 0 0] T`, `static_medium = 'strict_1z_dyson_ref'`.

**Provenance assertions all passed** — the fixture is the frozen one:

```
digest = ddb9532d11326458554b93b2ce09c80a3299cca9baa69202dc311f138b4fae17
digest matches frozen prereg = 1
info.dipole.backend = bruteforce ; isfield(info,'grid') = 0
mom.n = 16384 ; mom.Jbar = 1.2076644328e-04 ; mom.mu2 = 5.4826376529e-06
info.Jcc0 = 6.4244356557e-03
```

The digest reproducing at `a923691` — 35 commits after it was recorded at `2ee310b` — confirms this
branch did **not** perturb the coupling computation.

**The failing outcome:**

```
p.status    = medium_out_of_domain      (expected 'ok')
hstar       = NaN                       (expected finite and > 0)
p.crit_star = NaN                       (expected finite)
p.slope0    = -14.1262769841743

p.medium_status  (n=33):  ok x29 , ref_denom_nonpositive x4
p.node_term_reason:       converged x29 , medium_out_of_domain x4
trc.nodes (n=34):         ok x30 , ref_denom_nonpositive x4
```

**The affected region is a bounded interior window — four consecutive nodes with solved `ok` nodes
on both sides:**

```
 idx   h                 m          medium_status            term_reason
  20   0.0021219622      4.16509    ok                       converged
  21   0.002633223       4.54081    ref_denom_nonpositive    medium_out_of_domain
  22   0.0032676658      4.81953    ref_denom_nonpositive    medium_out_of_domain
  23   0.0040549699      5.00369    ref_denom_nonpositive    medium_out_of_domain
  24   0.0050319653      5.11630    ref_denom_nonpositive    medium_out_of_domain
  25   0.0062443557      5.18629    ok                       converged
```

Nodes 1–20 rise smoothly `m: 0 → 4.17`, all `ok`; nodes 25–34 run `m: 5.19 → 5.44`, all `ok`.

### What this does and does not establish

**Established:** the current warm-started **and** cold-retry Picard attempts encounter
`1 + Sigma0 <= 0` in a bounded intermediate-`h` window. That is sufficient to fail the present
implementation and Gate 0 clauses (a) and (e).

**Not established:** that `1 + Sigma0` *passes through zero* as a property of a converged solution,
or that no positive-denominator fixed point exists at those fields. The flagged nodes are **not
converged self-consistent states**. Verified in source:

- `invz_ordered_node_solve.m:226-228` — the outer map `break`s the moment `medium_status` leaves
  `{'not_applicable','ok'}`, **before** `invz_lambdas` at `:232`. The docstring states this at
  `:113-115`.
- `:101-105` and `:165-171` — a non-accepted warm attempt is retried exactly once from a cold start
  (`seed_kind = 'cold_after_warm_fail'`), so two independent seeds were tried and both halted.

Proving a branch crossing, or the nonexistence of a valid fixed point, would require a separate
continuation/root study. Under prereg §3 that is diagnosis toward a future candidate — **never a
rescue of this one**.

### Coverage gap

Prereg §8 requires ordered fields `[0.05 0.25 0.5 1 2 2.5 2.9 3.0] T` at `T = 0.10 K`, plus separate
PM controls `[3.1 3.5] T`. G11 exercised **one of the eight** required ordered fields and **none** of
the PM controls. Prereg §3's "STOPS AT DIAGNOSIS" licenses halting at the first failure, but any
Gate-0 report must state the evaluated set explicitly rather than implying full coverage.

---

## 3. Gate 0 comparator — all three schemes on the G11 fixture

Run by the controller (not part of the committed tests), same fixture as above:

| `static_medium` | `status` | `hstar` | `crit_star` | `slope0` | `medium_status` | **max `omit_max`** |
|---|---|---|---|---|---|---|
| `resummed` | **node_failed** | NaN | NaN | NaN | `not_applicable` ×33 | — |
| `strict_1z_dyson_ref` | **medium_out_of_domain** | NaN | NaN | −14.126276984174282 | `ok` ×29, `ref_denom_nonpositive` ×4 | **1.20077** |
| `strict_1z_bare_ref` | **ok** | 6.6684674214151885e-06 | 11.341529084851707 | −7.5305398555933216 | `ok` ×33 | **16.4358** |

`m_star` under the bare reference: `0.0185025557511318`.
`omit_max` ranges: dyson `[9.86627e-06, 1.20077]`, bare `[9.39745e-06, 16.4358]`.

1. **The incumbent fails outright on real couplings.** `resummed` → `node_failed`. This is the
   defect the branch exists to repair, now measured at the anchor point.

2. **The bare reference is a comparator only, not positive evidence.** `strict_1z_bare_ref` shows
   that the code can reach a finite root when the Dyson reference denominator is omitted. It is not
   evidence of a controlled physical solution: its `omit_max = 16.4358` is 164× the frozen limit,
   and its root (`hstar ≈ 6.7e-6`, `m_star ≈ 0.0185`) is very small. Its inability to report a
   reference-denominator failure is an **expected limitation of the convention** — `denom ≡ 1` by
   construction, so it has nothing to screen (Task 16 finding A) — not a clean bill of health.

3. **The omitted-order diagnostic exceeds the frozen limit under both strict conventions**, by 12×
   and 164×. On the synthetic 24-point fixture the same quantity was `0.0866` prospective (Task 16)
   and `0.0350` solved — comfortably under. **"Well under `omit_promote`" was an artifact of the
   synthetic multiset.** On the real 16384-point BZ density the omitted next-order terms dominate
   the retained one.

   Mechanism: `omit_cubic = |2*mu2^2 − mu4| * Gref^2 / mu2`; with `mu4 = 2.389*mu2^2` this is
   `0.389 * mu2 * Gref^2`, which grows exactly where `|Gref|` is large — next to the near-singular
   nodes.

**Path-integral ratios**, reported by like-for-like comparison rather than a single figure:

| profile | `int_r_minus_1 / int_Sigma0` |
|---|---|
| bare-ref, real couplings | `25.4842` |
| strict, synthetic fixture | `0.248441` → real is **102.6×** |
| resummed, synthetic fixture | `0.103301` → real is **246.7×** |

The spec's binding caution against substituting `∫Sigma0 dh` for `∫(r−1) dh` is roughly two orders
of magnitude more severe on the real density than on the fixture that motivated it.

---

## 4. G1a–d — the identity holds; the frozen gate is incompatible with the frozen grid

| test | measured | frozen gate | over by |
|---|---|---|---|
| G1a `max abs(dm/dh + G0bare_avg)/max(1,abs(gb))` | `1.252839309096e-03` | `1e-6` | 1253× |
| G1b `max abs(dF/dh − crit_avg)/max(1,abs(cm))` | `3.682882251199726e-04` | `1e-6` | 368× |
| G1c `max abs(dF/dm − crit/chi_path)/…` | `1.738339996023436e-05` | `1e-6` | 17× |
| G1d `e(3)` at nH=129 | `3.682882251199726e-04` | `1e-6` | 368× |

**Measured nH scaling** (controller and Codex independently, numbers agree exactly):

| nH | G1a | G1b | G1c |
|---:|---:|---:|---:|
| 33 | 1.8216829803e-02 | 5.843122218586044e-03 | 2.2609487299e-04 |
| 65 | 4.8713648324e-03 | 1.472635188345395e-03 | 6.4840455639e-05 |
| 129 | 1.2528393091e-03 | 3.682882251199726e-04 | 1.7383399960e-05 |

G1b reductions: `3.967800`, `3.998594` → empirical order `p12 = 1.9883`, `p23 = 1.9995`.

**The continuous identity `F' = crit` is confirmed, not refuted.** The residual is O(h²) panel error
from comparing secant/trapezoidal quantities to differential identities. Supporting evidence:

- G1d's own first two assertions — `e(1) >= e(2) >= e(3)` — **passed**. Only its absolute
  `e(3) <= 1e-6` check failed.
- `test_G1d_nested_trapezoid_is_second_order` **passed** on a pure `trapz(exp)` fixture, ratio inside
  `[3.5, 4.5]`, independently confirming the quadrature machinery is second order.

### G1a and G1b are the same test, provably

`invz_hmf_ordered.m:397-401` builds `h0 = cumtrapz(hgrid, r)` and `F = h0 - J0eff*m`; `:539` sets
`crit = r + J0eff*G0bare`. Because `h0` is a cumulative trapezoid of `r`, the secant `diff(h0)./diff(h)`
is *exactly* the midpoint average `r_avg`. Hence

```
dF/dh     = r_avg - J0eff * dm/dh
crit_avg  = r_avg + J0eff * G0bare_avg
--------------------------------------------------
dF/dh - crit_avg = -J0eff * (dm/dh + G0bare_avg)        [exact, to reassociation]
```

Codex measured this relation holding to `5.7e-15`; the derivation above shows it is an algebraic
identity, not an empirical coincidence. **G1b is therefore `J0eff` times G1a, differing only in the
normalising denominator. It is not an independent gate.**

Consequence for remediation: prereg §5 lists all three G1 forms as though independent. Two of them
are the same check. That redundancy belongs in the amendment.

### Both dials are frozen — an amendment is required, not a tolerance edit

- Prereg **§5** freezes all three `1e-6` identity thresholds.
- Prereg **§8** freezes `nH = [33 65 129]`, and states: *"Changing these grids **after seeing strict
  output** requires a dated preregistration amendment."* Strict output now exists, so the amendment
  requirement is triggered by construction.

The `1e-6` target is **not mathematically impossible** — it is incompatible with the frozen
`nH = 129` grid. At second order it would need roughly `sqrt(368) × 129 ≈ 2500` nodes. Any
remedy — finer grid, restructured gate, or revised threshold — touches a frozen value and requires a
dated re-registration section, **never an edit in place** (prereg §9 amendment rule).

Note for whoever drafts it: a "local derivative test retaining the `1e-6` continuum target" inherits
the same discretisation problem unless it uses Richardson extrapolation or a much finer local
stencil. Relabelling the test as local does not by itself buy three orders of magnitude.

---

## 5. G2 — the closure identity is exact; the failing assertion is out of scope

Three of four comparisons are **exactly zero**:

```
p.slope0 = -0.163499403500027 ; ptp.crit = -0.163499403500027
|p.slope0 - ptp.crit|            = 0.000000e+00      (gate AbsTol 1e-6)

GrefO = -171.948499251132 ; GrefP = -171.948499251132
|GrefO-GrefP|/|GrefP|            = 0.000000e+00      (gate RelTol 1e-6)

KO = KP = 0.00299680289420946
|KO - KP|                        = 0.000000e+00      <-- BITWISE
```

Prereg §7.3 explicitly does **not** guarantee this: *"The two callers reach Gref through different
expressions (`G0_PM(0)/(1+Sigma0)` vs `G0bare0/(1+Sigma0)`, equal only because `G0el0 -> 0`), so
bitwise identity is not a property the design guarantees."* It held anyway at this fixture — a
stronger result than the design promised.

The single failing comparison:

```
p.K0_pm0  = 0.00299680289420946
ptp.K(1)  = 0.00299680289523058
|p.K0_pm0 - ptp.K(1)| = 1.021112e-12
Kgate                 = 1.600000e-14
                      = 1.000e-14 + 1.000e-12*max([2.996803e-03, 2.996803e-03, 6.000000e-03])
headroom              = 63.819528   (must be <= 1)
```

`K_atol/K_rtol` is the **strict block-B recomputation gate** (prereg §2, whose own text reads: *"A
one-shot formula recomputed from the exported state should agree to floating-point reassociation
only; this gate catches mis-wiring, not physics."*). It is designed to catch closure mis-wiring
*within one exported state*. `p.K0_pm0` and `ptp.K(1)` are outputs of **two independently terminated
outer iterations** — the ordered predictor's node solve and the PM leg's converged
`invz_emt_scalar`. `1.02e-12` absolute is `3.4e-10` relative: cross-solver iteration residue, which
this gate was not designed to bound.

**Remedy is not to drop full-caller coverage.** Amend G2 so it checks caller provenance/convergence
and the exact cross-leg closure inputs and outputs under `K_atol/K_rtol`; if a direct cross-solver
exported-`K` comparison remains useful, record it under a bound derived from the declared
outer-solver tolerance instead.

---

## 6. G13 — no PM-slot leak, but the asserted invariant is not one the code promises

```
state.K(1) = 0.0027214553174579
med.K(1)   = 0.00271674377166701
difference = 4.711546e-06   (relative 1.734262e-03)     <-- verifyNotEqual PASSED

state.lam vs lam_from_exported:
  abs diff = [-4.634917e-13  -2.865552e-12  -1.771636e-11]
  rel diff = [-1.959624e-10  -3.615619e-10  -4.705009e-10]
  max relative error = 4.705009e-10        (asserted RelTol = 1e-10)

leak-detectability margin max(abs(lam_leaked - state.lam)) = 2.974382e-05
```

The sentinel's primary assertion passed: the ordered `K(1)` differs materially from the PM slot.

**Cause of the residue, verified in source:** `state.lam` is computed **inside** the outer loop at
`invz_ordered_node_solve.m:232`, whereas `state.K(1)` is overwritten **after** the loop by the
post-loop static refresh at `:257` (`K(1) = K0s`). Exact agreement between `state.lam` and
`invz_lambdas(state.K, …)` is therefore **not a guaranteed invariant** — it compares a pre-refresh
quantity against a post-refresh one. The `4.7e-10` is structural, not iteration noise.

> **Correction to the previous revision of this record.** It reported the leak margin as "about
> 63,000× larger than the observed residue". That divided an **absolute** lambda displacement
> (`2.974382e-05`) by a **relative** residual (`4.705009e-10`) — dimensionally invalid. In
> consistent units the margin is `2.974382e-05 / 1.771636e-11` ≈ **1.7×10⁶**. The corrected figure
> is larger and strengthens the conclusion, but the original arithmetic was wrong.

**Remedy:** replace the `RelTol = 1e-10` equality with a behavioural discriminator —

- retain the assertion that ordered and PM static slots are distinguishable;
- require the observed lambdas to be closer to the ordered-slot reconstruction than to the
  PM-slot-leaked reconstruction;
- capture `info` from `invz_ordered_node_solve` and require `info.accepted`;
- exercise the residual/final-refresh route explicitly rather than discarding its solver result.

This preserves G13's purpose — detecting a PM-slot leak — without asserting a stale-versus-refreshed
floating-point equality.

---

## 7. What passed

- **G3** — `r = 1 + Sigma0` holds for every `K0 ∈ {0, 1e-3, 5e-3, 0.05}` at `RelTol 1e-12`.
  Characterisation showed the algebraic condition producing this is `G0el0 = 0`: with `tl.m = 0.3`
  and `G0el0` still `0`, `out.r` remains `1.25` at `RelTol 1e-12` for all four `K0`; with
  `G0el0 = 0.4` it deviates (`1.24875669110374`, `1.2485702033296`, `1.24834791474165`,
  `1.26823452951256`).

  **Interpretation.** `G0el0 = 0` is the relevant *consequence* of the physical m→0 limit — prereg
  §7.3 independently identifies it as the condition through which the two `Gref` expressions
  coincide. The `tl.m ≠ 0` construction demonstrates algebraic sufficiency of `G0el0 = 0`; it does
  **not** show that physical m→0 behaviour is irrelevant. (An earlier revision of this record said
  the test's name "overstates what it discriminates" — too strong, and withdrawn.)

- **G15** — fatal ids escape every layer. Required a fix: the brief's
  `invz_try_solver_call(@() error(...))` is the bare anonymous void form that the primitive's own
  docstring forbids. Measured: `nargout` returns `-1`, the guard cannot fire, and the call throws
  `MATLAB:maxlhs` instead of returning the identifier. Replaced with a named local thrower declaring
  an output, per that docstring — then `completed=0`, `rid='invz:degenerateDoublet'`.
- **G15b** — a domain outcome keeps its category. `ref_margin = 1e9` drove all 33 nodes through
  `ref_denom_small`, `p.status = medium_out_of_domain`, confirming the intended path fired.
- **G1d nested trapezoid** — second order confirmed independently on a smooth analytic fixture.

---

## 8. Review verification log

Every load-bearing claim in that Codex review was checked against source or the
frozen preregistration before adoption:

| claim | verified against | result |
|---|---|---|
| Gate 0 is a lettered predicate; (a)/(c)/(e) fail | prereg §3 | confirmed verbatim |
| Nodes halt on domain exit, before `invz_lambdas` | `invz_ordered_node_solve.m:226-228`, docstring `:113-115` | confirmed |
| A cold retry runs after a failed warm attempt | same file `:101-105`, `:165-171` | confirmed |
| `state.lam` pre-refresh vs `state.K(1)` post-refresh | same file `:232` vs `:257` | confirmed |
| `K_atol/K_rtol` is the block-B mis-wiring gate | prereg §2, §7.3 | confirmed verbatim |
| The three `1e-6` G1 thresholds are frozen | prereg §5 | confirmed |
| Ratios 103× / 247× | arithmetic on measured values | confirmed (102.6, 246.7) |
| The `63,000×` margin is dimensionally inconsistent | arithmetic | confirmed — corrected to 1.7×10⁶ |
| G1a/G1b not independent at panel level | derived from `invz_hmf_ordered.m:397-401`, `:539` | confirmed, and shown to be exact algebra rather than measurement |

Additions made beyond the review: the algebraic proof of the G1a↔G1b relation and its consequence
for prereg §5; the fact that **both** the threshold (§5) and the grid (§8) are frozen; and the
Gate-0 field-coverage gap (1 of 8 required ordered fields, 0 of 2 PM controls).

---

## 9. Process failure worth recording

The Task 17 implementer **backgrounded the full-suite MATLAB run and polled its log, then died
mid-poll**. Consequences:

- `.superpowers/sdd/task-17-report.md` was never written.
- Its final message was a stray progress thought, not a status.
- `.superpowers/sdd/suite-task17.txt` is **truncated mid-run** at `test_invz_task2_matrix` and is
  not an authoritative number.

Everything above was recovered from the leftover scratch outputs
(`task17_slow_run_output.txt`, `task17_characterize_output.txt`, `task17_g11_detail_output.txt`,
`task17_g3_mcheck_output.txt`) and re-measured by the controller.

**Rule going forward: no agent may background the suite. Foreground only, with the result inline.**

Commit hygiene was otherwise clean: `%B` body carries no authorship trailer, exactly two test files
changed, worktree back to 30 unrelated dirty paths and 11 pre-existing stashes.

---

## 10. Open items

1. **Task 18** writes the formal Gate-0 **FAIL** report against clauses (a), (c), (e). It must state
   the evaluated field set (1 of 8 required ordered fields) rather than implying full coverage.
2. **A dated preregistration amendment** is required before G1, G2 or G13 can be corrected, because
   every relevant value is frozen. Draft: `Task-17_prereg_amendment_DRAFT.md`. It must be appended
   as a new dated section, **never edited in place** (prereg §9).
3. **Optional diagnosis, user's call:** map the failure across the remaining seven required ordered
   fields. This is diagnosis, not rescue — but it should be commissioned explicitly so it cannot be
   mistaken for a search for a passing field.

Nothing has been changed pending these. No tolerance was widened, no assertion weakened, no scheme
switched, and no production path altered.
