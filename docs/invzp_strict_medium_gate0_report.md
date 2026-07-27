# Gate 0 verdict: strict-order static medium (`strict_1z_dyson_ref`)

Task 18. Branch `invzp-stage2c-diagnostic`. Frozen preregistration:
`docs/invzp_strict_medium_prereg.md` (frozen 2026-07-25). This document records the measured
Gate-0/G5/G17/G7 output on the real production coupling multiset and states the verdict. It
does not amend the preregistration.

**Evaluated field set (full coverage, unlike Task 17's G11 which covered 1 of 8 ordered fields
and 0 of 2 PM controls):** all 8 required ordered fields `[0.05 0.25 0.5 1 2 2.5 2.9 3.0] T` at
all 3 required `nH in [33 65 129]` (24/24 combinations), both PM controls `[3.1 3.5] T`, and the
exact `B=0` hard-domain control. Nothing was stopped early; every preregistered point was run.

## 1. The frozen predicate (prereg SS3, quoted verbatim)

> ## 3. Gate-0 negative-outcome rule (spec SS6.0(3), verbatim intent)
> Promotion FAILS if any of:
>   (a) any required solved-path node has a non-'ok' REFERENCE denominator status;
>   (b) any skipped or invalid node is unaccounted for in the coverage counters;
>   (c) max(omit_max) over the solved path exceeds omit_promote (below);
>   (d) any LOCAL Gstat denominator crossing at which r or crit is non-finite or discontinuous
>       (G17 covers the algebra; this covers the actual path);
>   (e) any required ordered field does not return status='ok', a finite nonzero root, and a
>       stable endpoint under the frozen crit/D_uni/Dq margins, or either PM control fails to
>       return a converged finite positive-mass PM state.
> A local Gstat crossing that satisfies G17 does NOT fail promotion (spec SS1: the singularity is
> removable in the integrand).
> On failure the run STOPS AT DIAGNOSIS. Carrying another moment, changing Gref, or truncating
> other Matsubara sectors is a NEW theory candidate requiring a new spec and fresh
> preregistration -- never an in-run fallback. Regularisation, broadening and tolerance widening
> remain forbidden.

Frozen constants used below (prereg SS1/SS2/SS4/SS5, user-approved 2026-07-25): `crit_tol =
1e-6`, `D_tol = Dq_tol = 1e-6*max(1,|Gstat|*Jscale)` (state-dependent, evaluated at the root),
`ref_margin = 1e-6`, `omit_promote = 0.10`, `pole_cont_tol = 1e-3`, `I_atol = 1e-10 meV` (G5).
None were widened, loosened, or re-derived from the output below.

## 2. Exact input fingerprint / provenance

Requested coupling tuple: `grid=[16 16 16], dpRng=30, dipole='bruteforce', cache=false` (no
grid-policy fields). Coupling build wall-clock: **55.890 s**.

```
info.grid ABSENT       isfield(info,'grid') == false   (confirmed -- the legacy/no-grid-policy route)
info.dipole.backend          = 'bruteforce'
info.dipole.ewald            = struct('alpha',[],'r_cut',[],'g_cut',[],'boundary','')  (sentinel-empty)
info.dipole.q_reduction      = 'bruteforce: q used directly as MF_dipole/exchange Miller
                                indices (q*geom.b); no canonical q-domain reduction applied'
info.dipole.primitive_schema = 'MF_dipole+exchange (legacy, unversioned)'
digest (invz_exact_numeric_digest(Jnu(:))) =
  ddb9532d11326458554b93b2ce09c80a3299cca9baa69202dc311f138b4fae17
  MATCH vs frozen preregistered digest = 1   (verified inside invz_gate0_report on EVERY
                                              invocation, before any solve, per prereg SS8)
info.Jcc0 (J0eff)      =  6.424435656e-3 meV
info.Jaa0              =  3.510446205e-3 meV
n                      =  16384  (4096 q x 4 branches)
Jbar  = mean_q J       =  1.207664433e-4 meV
mu2   = mean((J-Jbar)^2) = 5.482637653e-6
mu3   = mean((J-Jbar)^3) = -3.42227577e-11
mu4   = mean((J-Jbar)^4) = 7.182350058e-11
J_max, J_min           =  5.985138929e-3, -6.763100317e-3 meV
```

Every one of these matches the frozen prereg SS10 record to the last digit quoted there
(reproduced independently in this task, not copied from it).

## 3. Step 3 -- Gate-0/G5/G17 measured tables

### 3.1 Ordered fields, all (B, nH) combinations (24/24)

```
     B   nH status                      hstar      crit*     D_uni*    Dq_min* D_tol*  n_nod  n_ok n_mod  n_dd  n_sf n_unrec   min_refmar   maxomitL     int_Sig0      int_rm1 g17_r        g17_crit
  0.05   33 medium_out_of_domain          NaN        NaN        NaN        NaN      1e-06     34     11    23     0     0     0  -1.5593e+05     0.3278     -10.9003          NaN no_crossing  no_crossing
  0.05   65 medium_out_of_domain          NaN        NaN        NaN        NaN      1e-06     66     21    45     0     0     0  -1.5593e+05     0.3278     -10.7806          NaN no_crossing  no_crossing
  0.05  129 medium_out_of_domain          NaN        NaN        NaN        NaN      1e-06    130     40    90     0     0     0  -1.5593e+05     0.3278     -10.7527          NaN no_crossing  no_crossing
  0.25   33 medium_out_of_domain          NaN        NaN        NaN        NaN      1e-06     34     11    23     0     0     0      -951.89    0.33239    -0.790492          NaN no_crossing  no_crossing
  0.25   65 medium_out_of_domain          NaN        NaN        NaN        NaN      1e-06     66     21    45     0     0     0      -951.89    0.33239    -0.728301          NaN no_crossing  no_crossing
  0.25  129 medium_out_of_domain          NaN        NaN        NaN        NaN      1e-06    130     40    90     0     0     0      -954.19    0.33239    -0.703678          NaN no_crossing  no_crossing
   0.5   33 medium_out_of_domain          NaN        NaN        NaN        NaN      1e-06     34     22    12     0     0     0      -39.037     1.4857   -0.0628266          NaN no_crossing  no_crossing
   0.5   65 medium_out_of_domain          NaN        NaN        NaN        NaN      1e-06     66     42    24     0     0     0      -39.037     3.6055   -0.0629776          NaN no_crossing  no_crossing
   0.5  129 medium_out_of_domain          NaN        NaN        NaN        NaN      1e-06    130     80    48     0     2     0      -329.06     3.6055   -0.0714393          NaN no_crossing  no_crossing
     1   33 medium_out_of_domain          NaN        NaN        NaN        NaN      1e-06     34     30     4     0     0     0      -1.4155     1.2011   -0.0060487          NaN no_crossing  no_crossing
     1   65 medium_out_of_domain          NaN        NaN        NaN        NaN      1e-06     66     58     8     0     0     0      -1.4155     1.2822  -0.00615107          NaN no_crossing  no_crossing
     1  129 medium_out_of_domain          NaN        NaN        NaN        NaN      1e-06    130    113    17     0     0     0      -15.468     1.2822   -0.0123479          NaN no_crossing  no_crossing
     2   33 ok                      0.0297551     0.8047     0.8768     0.8856      1e-06     43     43     0     0     0     0      0.88735    0.66027 -0.000638912   0.00157468 jump_exceeded jump_exceeded
     2   65 ok                      0.0297862      0.805     0.8771     0.8858      1e-06     74     74     0     0     0     0      0.88654    0.66027 -0.000654474   0.00154521 ok           ok
     2  129 ok                      0.0297964     0.8051     0.8772     0.8859      1e-06    137    137     0     0     0     0       0.8864    0.66027 -0.000658348   0.00153787 ok           ok
   2.5   33 ok                      0.0278292     0.7145     0.7779     0.7939      1e-06     43     43     0     0     0     0      0.92783     0.3165   3.6715e-05   0.00118161 no_crossing  no_crossing
   2.5   65 ok                      0.0278414     0.7147     0.7781     0.7941      1e-06     74     74     0     0     0     0      0.92692     0.3165  1.92035e-05   0.00115672 no_crossing  no_crossing
   2.5  129 ok                      0.0278476     0.7147     0.7782     0.7942      1e-06    137    137     0     0     0     0      0.92692     0.3165   1.4844e-05   0.00115052 no_crossing  no_crossing
   2.9   33 ok                      0.0253359     0.6247     0.6754     0.6991      1e-06     43     43     0     0     0     0      0.94467    0.19766  0.000424218   0.00120891 no_crossing  no_crossing
   2.9   65 ok                      0.0253713     0.6254     0.6761     0.6997      1e-06     74     74     0     0     0     0      0.94467    0.19766  0.000408509   0.00118845 no_crossing  no_crossing
   2.9  129 ok                       0.025365     0.6253     0.6759     0.6996      1e-06    137    137     0     0     0     0      0.94467    0.19766  0.000404595   0.00118335 no_crossing  no_crossing
     3   33 ok                      0.0246288     0.6006     0.6474     0.6733      1e-06     43     43     0     0     0     0      0.94819    0.17782  0.000508085   0.00123209 no_crossing  no_crossing
     3   65 ok                      0.0246567     0.6012      0.648     0.6738      1e-06     74     74     0     0     0     0      0.94815    0.17782  0.000493232   0.00121303 no_crossing  no_crossing
     3  129 ok                       0.024654     0.6012      0.648     0.6738      1e-06    137    137     0     0     0     0      0.94794    0.17782   0.00048953   0.00120827 no_crossing  no_crossing
```
Columns: `n_nod`=`numel(trc.nodes)`; `n_ok/n_mod/n_dd/n_sf/n_unrec` = the four-category coverage
buckets (`ok` / `medium_out_of_domain` / `degenerate_doublet` / `solver_failed`) plus
`unrecognized`, always summing to `n_nod` (`n_accounted == numel(trc.nodes)` held on all 24
rows: `n_unrec` is 0 everywhere). `min_refmar` = minimum finite `ref_margin` (denom - floor)
over the ledger; `maxomitL` = `max(omit_max)` over ledger nodes labelled `medium_status='ok'`.

**The `1 T` row is the mandatory direct comparison with Task 17's G11 measurement**: `status =
medium_out_of_domain`, `hstar = NaN`, `n_nod = 34` at `nH=33`. `max(omit_max)` over the COMPLETE
ledger (34 nodes, predictor included, per this task's "read the complete trace ledger" mandate)
is `1.201054187`; restricted to the 33-node PROFILE alone (excluding the predictor, matching what
`prof.omit_max` -- the population G11 characterized -- would give) it is `1.200766324`, which
rounds to G11's recorded `1.20077` exactly. The small difference between the two is real and
understood, not display rounding: the predictor node (`h=0`) itself carries `omit_max=
1.201054187`, marginally the single largest value in the whole ledger at this field, so including
it (as this task's brief requires) versus excluding it (as the profile-only characterization
does) genuinely changes the reported maximum in the 4th significant figure. Both values obviously
still exceed `omit_promote=0.10` by more than 12x either way. The **`3 T`** row (the closest
required ordered field to the legacy 3.025 T PM boundary) is `status='ok'`
with a finite root, but still exceeds `omit_promote` (`max(omit_max) = 0.17782`), and its
`int_r_minus_1` (`0.0012`) is over **150,000 times** the low-field rows' `NaN` -- see SS5 below on
why that boundary-adjacent shift cannot bound the low-field integrals.

Note the `0.5 T, nH=129` row: `n_sf = 2` -- the only ledger entries in the entire sweep classified
`solver_failed` (a genuine convergence non-acceptance, distinct from every other failure mode
measured here, which is a reference/closure domain event or an endpoint tolerance miss).

### 3.2 PM controls (`invz_solve_point` only -- the Jensen HMF solver was never called on these)

```
     B converged         crit   crit_tol medium_status
   3.1         1    -0.884482      1e-06 ok
   3.5         1    -0.528219      1e-06 ok
```
Both PM points **converge** and report `medium_status='ok'`, yet both report a **NEGATIVE**
`crit` -- clause (e) requires "a converged finite **positive**-mass PM state", so both fire (e)
independently of every ordered-field failure below. This was cross-checked directly against
`invz_solve_point` outside the driver (bypassing all Gate-0 machinery): at `B=3.1 T` the
`resummed` scheme does not even converge (`pt.converged=false`, matching Task 17's "resummed
gives node_failed" note) and reports `crit=-0.860014`; at `B=3.5 T` `resummed` also fails to
converge but reports `crit=+0.374153` (the physically-expected sign for a field-polarized PM
state). The strict scheme's own outer loop reports `converged=true` at both fields while its
`crit` sign is wrong at both -- a new, independent PM-route symptom beyond the ordered-1T anchor
Task 17 measured, not an artifact of this driver (confirmed by the bypass call).

### 3.3 Exact `B=0` -- hard-domain control (excluded from clauses (a)/(e) by construction)

```
status=degenerate_doublet hstar=NaN n_nodes=34 n_accounted=34 expected_degenerate=1
```
Expected and confirmed: with no field at all (neither transverse `Bx` nor a nonzero longitudinal
`hz`), the zero-field doublet is exactly degenerate at the `h=0` predictor node, which forces
`prof.status='degenerate_doublet'` by `invz_hmf_status`'s binding precedence (spec SS5.3) even
though the rest of the 33-node profile is otherwise well-formed. This control was never passed
into `invz_gate0_aggregate` and therefore cannot contribute to `fail_a` or `fail_e` by
construction, not merely by the numbers happening to cooperate.

### 3.4 G5 -- path integral convergence (33->65 diagnostic, 65->129 frozen criterion)

```
     B |  d33->65(Sig0)   d33->65(rm1)     (diag) | d65->129(Sig0)  d65->129(rm1)    g5_pass   (tolS/tolR = I_atol + 1e-3*max(|I_fine|,|I_prev|))
  0.05 |       0.119727            NaN         -- |      0.0278231            NaN      false   (tolS=0.0108 tolR=NaN)
  0.25 |      0.0621912            NaN         -- |      0.0246231            NaN      false   (tolS=0.000728 tolR=NaN)
   0.5 |       0.000151            NaN         -- |     0.00846166            NaN      false   (tolS=7.14e-05 tolR=NaN)
     1 |    0.000102368            NaN         -- |     0.00619687            NaN      false   (tolS=1.23e-05 tolR=NaN)
     2 |    1.55618e-05    2.94724e-05         -- |    3.87346e-06    7.34317e-06      false   (tolS=6.58e-07 tolR=1.55e-06)
   2.5 |    1.75115e-05    2.48911e-05         -- |    4.35947e-06    6.19985e-06      false   (tolS=1.93e-08 tolR=1.16e-06)
   2.9 |    1.57089e-05    2.04554e-05         -- |    3.91356e-06    5.09724e-06      false   (tolS=4.09e-07 tolR=1.19e-06)
     3 |    1.48537e-05    1.90605e-05         -- |    3.70193e-06    4.75107e-06      false   (tolS=4.93e-07 tolR=1.21e-06)
```
`int_r_minus_1` is `NaN` at all three `nH` for all four low-field rows (`0.05, 0.25, 0.5, 1 T`) --
**the finest integral is missing** at those fields, exactly the condition prereg SS5 says Gate 0
cannot pass under (those fields already fail independently via (a)/(c)/(e), so this does not
require a sixth gating clause; see SS6 of this document). G5 is a **separate numerical-quality
gate** (prereg SS5, alongside G1), not one of the (a)-(e) Gate-0 clauses, and is reported here in
full rather than folded into `rep.pass`: note it fails its own 65->129 criterion on **every**
field, including the four with `status='ok'` -- a further, independent numerical symptom beyond
the (a)/(c)/(e) clauses, measured and stated plainly rather than absorbed into the pass/fail
Booleans it was never defined to be part of.

### 3.5 G17 -- actual-path pole continuity, nH=65 vs nH=129 (clause (d))

```
     B |       r@65      jump@65        |      r@129     jump@129        |    crit@65      jump@65        |   crit@129     jump@129
  0.05 | no_crossing          NaN        | no_crossing          NaN        | no_crossing          NaN        | no_crossing          NaN
  0.25 | no_crossing          NaN        | no_crossing          NaN        | no_crossing          NaN        | no_crossing          NaN
   0.5 | no_crossing          NaN        | no_crossing          NaN        | no_crossing          NaN        | no_crossing          NaN
     1 | no_crossing          NaN        | no_crossing          NaN        | no_crossing          NaN        | no_crossing          NaN
     2 |         ok    0.0008346        |         ok    0.0001594        |         ok    0.0004729        |         ok     5.04e-05
   2.5 | no_crossing          NaN        | no_crossing          NaN        | no_crossing          NaN        | no_crossing          NaN
   2.9 | no_crossing          NaN        | no_crossing          NaN        | no_crossing          NaN        | no_crossing          NaN
     3 | no_crossing          NaN        | no_crossing          NaN        | no_crossing          NaN        | no_crossing          NaN
```
Only `B=2 T` has an actual local-`Gstat`-denominator crossing on this multiset (`nH=33` alone
showed `jump_exceeded` on the coarse grid; both refinements resolve it `ok` with the jump
**decreasing** `0.0008346->0.0001594` (r) and `0.0004729->0.0000504` (crit), well inside
`pole_cont_tol=1e-3` and satisfying "must not increase" -- a genuine **removable** singularity in
the integrand, exactly the SS1 exemption, not a clause-(d) failure). No field showed an
unresolved crossing or a non-finite integrand on an `'ok'`-labelled node. **Clause (d) does not
fire anywhere in this sweep.**

## 4. Step 4 -- G7 scheme-jump measurement (non-gating, prereg SS6)

```
T=0.05  jump|K1s-K1r|=0.0001492  dispersion|K2-K1|=0.0002524  ratio=0.5913
T=0.10  jump|K1s-K1r|=0.0001307  dispersion|K2-K1|=0.000246  ratio=0.5311
T=0.31  jump|K1s-K1r|=8.617e-05  dispersion|K2-K1|=0.0001905  ratio=0.4523
T=1.00  jump|K1s-K1r|=5.955e-05  dispersion|K2-K1|=0.0003434  ratio=0.1734
```
`B=6 T`, `Ecut=40 meV`, same coupling fixture. The `omega=0` scheme jump `|K(1)_strict -
K(1)_resummed|` is **17%-59%** of the physical dispersion `|K(2)-K(1)|` across the four
temperatures -- not negligible relative to the coupling's own q-dependence. No pass threshold is
defined for G7 (prereg SS6); it is reported as a measurement.

## 5. The ~0.3% PM boundary shift bounds NEITHER path integral

Prereg SS4's binding caution states explicitly that the legacy ~0.3% PM boundary shift "bounds
neither integral Sigma0 dh nor integral (r-1) dh deep in the ordered phase." The measurements
above confirm this directly: at `B=0.05 T`, `int_Sigma0 = -10.75` (nH=129) -- **~3500x** larger in
magnitude than a 0.3% shift could plausibly bound -- and at every low field `int_r_minus_1` is
not even finite. A ~0.3%-scale boundary observation cannot be extrapolated into a bound on either
path integral away from the boundary, and this task did not attempt to use it as one.

## 6. Verdict

Per-clause Booleans, from SS3 above, evaluated over the complete 24-row ordered ledger, both PM
rows, with the `B=0` control excluded by construction:

```
fail_a = TRUE   (12 of 24 ordered rows: all three nH at B in {0.05, 0.25, 0.5, 1} T have >=1
                 ledger node with medium_status not in {'ok','not_applicable'})
fail_b = FALSE  (coverage identity holds on every row: n_accounted == numel(trc.nodes) always;
                 every predictor phase present; every finite hstar backed by a matching, value-
                 agreeing 'root' ledger entry; no dropped chunk)
fail_c = TRUE   (max(omit_max) over 'ok'-labelled ledger nodes exceeds omit_promote=0.10 on
                 ALL 24 rows -- the smallest measured value, 0.17782 at B=3 T nH=129, is still
                 78% over the frozen limit)
fail_d = FALSE  (the one genuine local Gstat crossing, at B=2 T, is resolved 'ok' at both nH=65
                 and nH=129 with a decreasing jump, well under pole_cont_tol; no unresolved
                 crossing, no increasing jump, no non-finite integrand on an 'ok' node, anywhere)
fail_e = TRUE   (12 of 24 ordered rows fail status='ok'/finite-root/stable-endpoint, exactly the
                 same B in {0.05,0.25,0.5,1} T rows that fail (a); AND, independently, both PM
                 controls report a converged but NEGATIVE crit)

rep.pass = ~(fail_a || fail_b || fail_c || fail_d || fail_e) = FALSE
```

## VERDICT: **FAIL** -- clauses (a), (c) and (e) fired; (b) and (d) held.

Per prereg SS3: **the run stops at diagnosis here.** No scheme switch, tolerance widening,
regularisation, or filtered node was applied anywhere above, and none is proposed. A revised
theory (a different retained moment, a different reference denominator, a truncated Matsubara
sector, or anything else) is a new theory candidate requiring its own spec and a fresh
preregistration -- never an in-run fallback from this one. Stage 4 (G6/G6d, G8, G10, G12, G14,
G16, the default-flip decision) is not worth planning against this scheme; the blind Stage-4
freeze in prereg SS9 remains exactly what it was: a blind freeze for a candidate that did not
survive Gate 0.
