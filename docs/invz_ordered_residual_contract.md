# Ordered jensen node — complete residual/scaling contract (Stage 2c, Task 1a)

**Status.** Pre-declared contract, committed BEFORE the shared solver (Task 1b) per the
Stage-2c plan of record (`invzp_stage3_plan.md` rev.3) and its independent review. This
document, plus the checker it specifies (`invz_projected/invz_ordered_residual.m`), is the
"measuring instrument": Task 1b will gate node acceptance on it. **Nothing in this document
changes any production acceptance gate.** The existing acceptance criteria in
`invz_solve_point_ordered.m` (`dS < tol_outer && sout.converged`, plus the post-loop
`final_resid < tol_outer` check) and `invz_hmf_ordered.m`'s `eval_node` (the duplicate of the
same loop) are untouched by this task.

## 0. Why this contract exists

Both ordered "jensen" node loops accept on a **pre-mix Sigma step** (`dS`) plus the static
closure's own convergence flag. Neither loop, nor the post-loop refresh in
`invz_solve_point_ordered.m:225-235`, ever recomputes the **dynamic** medium `K(2:end)` from
the **final, mixed** `Sigma`. So "accepted" today certifies only: (a) the Sigma step was
small, and (b) the static (`w=0`) sector closed. It does **not** certify that the exported
`(Sigma, K, lambda, K0s)` tuple is a genuine joint fixed point of the full coupled map. A small
step is not a small residual — this is precisely the gap `invz_ordered_residual.m` measures.

## 1. State vector

The nested formulation (decided; not re-litigated here) treats the node state as:

| Quantity | Role | Size | Units | Computed by |
|---|---|---|---|---|
| `Sigma` | **Outer unknown** (the only independent residual target) | `[nw x 1]`, `nw = numel(wn)` | dimensionless | `invz_sigma_ordered` |
| `K0s` | **Continuation seed** (NOT independent) | scalar | meV | `invz_emt_static_ordered`'s own Picard loop over `K0` |
| `K` | **Derived** each map evaluation | `[nw x 1]` | meV | `K(2:end)` from `invz_emt_scalar`; `K(1)` OVERWRITTEN by the closed `K0s` |
| `lam` | **Derived** each map evaluation | `[3 x 1]`, `p = [1 2 3]` | `lam(1)`: meV; `lam(2)`: dimensionless; `lam(3)`: meV⁻¹ | `invz_lambdas(K, g, wts, beta, [1 2 3])` |

`K0s` and `lam` are **not independently gated** as residual targets (per the plan's rule: a
derived quantity is never described as an independent residual component). They are inputs
to, or outputs of, the blocks below, never quantities whose own "correctness" is asserted in
isolation.

**Units derivation** (traced to the codebase's own conventions, not invented here):
- `G`/`chi` quantities are meV⁻¹ by the repo-wide convention `G = -chi`
  (`invz_gstat_ordered.m` docstring). `beta = 1/(kB T)` is meV⁻¹ (`invz_const.m`: `kB` in
  meV/K). `g(z) = 2 n01 Delta/(Delta^2 - z^2)` (`invz_g.m`) is meV⁻¹ (`Delta` in meV).
- `Sigma` is dimensionless: it enters as `D = 1 + Sigma` in `invz_emt_scalar.m:24`, and as
  `1 + Sigma0 - J0eff*chi0cc0` in `pt.crit` (`invz_solve_point_ordered.m:252`, `J0eff` meV
  times `chi0cc0` meV⁻¹ — dimensionless), so `Sigma` must be dimensionless for both sums to
  type-check.
- `K` is meV: in `invz_emt_scalar.m`, `A = mean(Jf./(D + Jf.*G0))` has units meV (Jf meV over a
  dimensionless denominator), and `K = A.*D./(1 - A.*G0)` keeps those units (`A*G0` meV·meV⁻¹ =
  dimensionless, so the denominator is dimensionless). Cross-checked by
  `D_uni = 1 + (J0eff - K0)*Gstat` (`invz_emt_static_ordered.m:60`): `J0eff` meV, so `K0` must
  be meV for `(J0eff-K0)` to type-check against `Gstat` meV⁻¹ producing a dimensionless sum.
- `lam(p) = (1/beta) sum wts.*K.*g.^p` (`invz_lambdas.m`): `(1/beta)` is meV; `K` is meV; `g^p`
  is meV⁻ᵖ. So `lam(p)` is meV^(2-p): `lam(1)` meV, `lam(2)` dimensionless, `lam(3)` meV⁻¹.
  Cross-checked against `invz_sigma_ordered.m`'s `alpha_m` sum, whose four additive terms
  (`lam(2)`; `g0*lam(1)`, `g0` meV⁻¹; `(4/g0)*lam(3)`, `4/g0` meV; and
  `(1-n01^2)(1+...)*K0*g0`, `K0` meV times `g0` meV⁻¹) are ALL dimensionless together, matching
  `lam(2)` alone — confirms the per-component units above.

## 2. Node input bundle

Constant across one node's outer loop (computed once in each caller's preamble —
`invz_solve_point_ordered.m:177-205`; `invz_hmf_ordered.m` `eval_node`, mirrored statement for
statement). All are **read-only inputs** to the map; none is mutated by the checker.

| Field | Meaning | Units |
|---|---|---|
| `tl` | Two-level struct (`invz_twolevel_ordered`): `m, M2, n01, g0, h0, hz, Delta, transverse_mf` | mixed (see `invz_twolevel_ordered.m`) |
| `G0` | Full electronuclear cc propagator, `-real(squeeze(c0(3,3,:)))`, `size(wn)` | meV⁻¹ |
| `g` | `real(invz_g(tl, 1i*wn))`, `size(wn)` | meV⁻¹ |
| `wts` | Matsubara quadrature weights (`[1; 2; 2; ...]`) | dimensionless |
| `wn` | Bosonic Matsubara grid, `n>=0` | meV |
| `beta` | `1/(kB T)` | meV⁻¹ |
| `J0eff` | Uniform cc coupling `J(q=0)` (ODD-shifted if applicable) | meV |
| `G0inel0` | Static inelastic-only propagator at `wn=0` | meV⁻¹ |
| `G0el0` | `G0bare0 - G0inel0`, the elastic+feedback static weight | meV⁻¹ |
| `G0bare0` | Path-consistent bare static response (`-(X(3,3)+fb)`) | meV⁻¹ (node-input only; **not** read by any of blocks A-D — it is the precursor `G0el0` was derived from in the caller's preamble, retained here for provenance/documentation, per the brief's explicit node-input list) |
| `eso` | `emt_static` opts struct threaded to `invz_emt_static_ordered` (carries `resid_tol`, default 1e-10 meV⁻¹; `warn`, forced `false` by both callers) | n/a |
| `eopts` | `emt` opts struct threaded to `invz_emt_scalar` (`K0`/`debug`/`freq_block`) | n/a |
| `Jnu_flat` | Static mode spectrum, `[nJ x 1]` (or `[nJ x nw]`, T2.1 retardation — not exercised by the synthetic control fixture) | meV |

`eopts.K0` (backward-compatible input to `invz_emt_scalar`) is **provably inert**:
`invz_emt_scalar.m`'s own docstring states "`tol`/`max_iter`/`mix`/`K0` are accepted for
backward compatibility but unused: the solve is direct". The checker therefore calls
`invz_emt_scalar(G0, Sigma, Jnu_flat, eopts)` without setting `.K0`, matching the function's
own documented behaviour exactly (setting it would change zero numerics; omitting it keeps the
checker's code simpler without altering fidelity to the production call).

## 3. The pure node map `F`

One call `Sigma_next = F(Sigma; lam_seed, K0s_seed | node)` reproduces **one full outer
iteration's body**, exactly as implemented (`invz_solve_point_ordered.m:206-221` /
`invz_hmf_ordered.m` `eval_node:313-328`, verbatim duplicates):

```
med  = invz_emt_scalar(G0, Sigma, Jnu_flat, eopts);            % (1) dynamic sector -- DIRECT
K    = med.K;                                                  %     solve, pure fn of Sigma
[K0s_new, Gstat, sout] = invz_emt_static_ordered(tl, lam_seed(1:2), Sigma(1), Jnu_flat, ...
                                    K0s_seed, beta, J0eff, G0inel0, G0el0, eso);  % (2) static
K(1) = K0s_new;                                                %     closure to ITS OWN tol
lam_next = invz_lambdas(K, g, wts, beta, [1 2 3]);             % (3) derived lambdas (FRESH,
sg   = invz_sigma_ordered(tl, lam_next, K, g, beta);           % (4) ordered Sigma map
Sigma_next = sg.Sigma;
```

**On the `lam_seed` input to step (2) — a necessary precision beyond the brief's own
one-line summary, recorded here explicitly rather than glossed over.** Step (2)'s static
closure needs *a* `lam(1:2)` value to evaluate `invz_gstat_ordered` internally, but `lam` is
only *produced* by step (3), which needs the *already-closed* `K` (i.e., `K(1)` from step (2)).
This is a genuine circularity in the map's own definition, resolved in production by
successive substitution **across outer iterations**: step (2) of iteration `i` uses the `lam`
value **left over from iteration `i-1`**'s step (3) (or the `[0;0;0]` cold start), not a
value freshly derived within the same pass. `K0s` plays the identical structural role for
step (2)'s `K0` seed. **Both `lam(1:2)` and `K0s` are therefore continuation seeds in the
identical operational sense** — this is why `state` (Sec. 1) bundles `lam` alongside `K0s`
rather than only `Sigma`: `F`, evaluated at the exported state, uses `state.lam(1:2)` and
`state.K0s` as the two seeds step (2) needs, mirroring exactly what the live loop would do if
fed this state as one more warm start. At a genuine fixed point (the loop has actually
converged) this lag is immaterial: consecutive iterations' `lam` values have themselves
converged, so `state.lam` and the step-(3) value it would produce agree to the same order as
the overall residual — confirmed empirically on the synthetic control (Sec. 9): `state.lam`
and a step-(3) recomputation from the converged `K` agree to machine precision.

`F` performs **both inner solves internally** (the dynamic EMT direct solve; the static `K0`
closure to its own `resid_tol`) — a caller of `F` never sees an intermediate `K`/`lam`.
`K`/`lam` are not read from `state` by `F` at all (`state.K` plays no role in evaluating `F`);
they are entirely regenerated fresh by steps (1)-(3). This is why a corruption of `state.K`
alone (leaving `Sigma`/`lam`/`K0s` untouched) is invisible to Block A (Sec. 4) — `F` never
looks at it — but is caught by Blocks C/D, which explicitly validate `state.K` against
independent recomputations (Sec. 9, Perturbation 2).

## 4. Residual blocks

Every block **recomputes independently from the exported `state`** — never reusing another
block's intermediate result, never reusing the original solve's cached values. All four blocks
are evaluated even if one fails (no short-circuit), so a caller always sees the complete
picture.

**NaN-propagation in the vector-difference residuals (`rA`, `rC`, `rD`).** MATLAB's `max`
**ignores individual `NaN` entries** by default (`max(abs([1 NaN 3])) == 3`, unlike `sum`/
`mean`, which propagate `NaN`). A bare `max(abs(a-b))` over a mostly-finite vector with one
corrupted entry would therefore silently report the max over the *remaining* finite entries,
masking that one component at the **block level** — even though the separate top-level
`.finite` flag (Sec. 7) would still catch it via `state`'s own fields directly. Reporting a
block as `pass = true` on a partially-`NaN` comparison would not be honest measurement, so the
checker computes these three residuals via an explicit NaN-propagating helper (`any
non-finite entry in either operand -> residual = NaN`) rather than relying on `max`'s default.
This changes nothing for any fully-finite comparison (identical arithmetic result whenever no
entry is non-finite — confirmed on every fixture in Sec. 9) and only strengthens the
non-finite case.

### Block A — outer Sigma map

```
rA = max| F(Sigma_exp; lam_exp, K0s_exp) - Sigma_exp |
```
**Recomputation point:** one full independent pass of steps (1)-(4), Sec. 3, seeded from
`state.lam`/`state.K0s`, starting from `state.Sigma`.
**Units:** dimensionless (Sigma is dimensionless, Sec. 1).
**Scale:** `scale_abs = scale_rel = tol_outer` (default `1e-8`, the caller's own outer
tolerance — `opts.tol_outer`, matching `invz_solve_point_ordered.m`'s `tolo` default).
**Justification (units-based):** `Sigma` is constructed to enter as `D = 1 + Sigma`
(`invz_emt_scalar.m:24`) — i.e. it is measured in units of "the bare inverse propagator",
which is exactly `1`. Its natural comparison scale is therefore `O(1)` **by construction**,
so an absolute tolerance in these units already **is** the natural relative tolerance; no
separate relative scale is needed (`scale_rel` is recorded identically to `scale_abs` purely
for schema uniformity with Blocks B/D, not because a distinct relative quantity was derived).
**Pass:** `isfinite(rA) && rA < scale_abs`.

### Block B — static EMT closure

```
[K0_b, Gstat_b, so] = invz_emt_static_ordered(tl, lam_exp(1:2), Sigma_exp(1), Jnu_flat, ...
                                              K0s_exp, beta, J0eff, G0inel0, G0el0, eso)
rB = so.resid   (= |mean(Gq) - Gstat_b|, the closure's own internal residual)
```
**Recomputation point:** a fresh call to `invz_emt_static_ordered` at the exported
`Sigma(1)`/`lam(1:2)`, seeded from the exported `K0s` (independent of any cached solve output
— this call is never reused from Block A's own internal step (2), even though the two calls
are analytically identical when both are fed the same, correct `state`; keeping every block's
recomputation self-contained means a bug in one block's wiring cannot silently hide behind
another's).
**Units:** meV⁻¹ (inverse-energy: `Gstat`/`Gq`, per the repo's `G=-chi` convention).
**Scale:** `rtol_B = getf(node.eso, 'resid_tol', 1e-10)` (meV⁻¹) — reused as **both**
`scale_abs = rtol_B` and `scale_rel = rtol_B`. **Combined gate:**
`rB < scale_abs + scale_rel * abs(Gstat_b)`.
**Justification (units-based, not tuned):** `rtol_B` is not a new number invented for this
task — it is `invz_emt_static_ordered`'s **own, already-declared** primary convergence
criterion (`invz_emt_static_ordered.m:38`, documented there as "meV⁻¹"). Reusing it, rather
than picking a fresh constant, is the units-based justification. A **bare absolute** gate
`rB < 1e-10` is wrong for this residual because `Gstat` is measured to be **O(300)** near the
synthetic control's own accepted node (Sec. 9) and is documented to reach **O(500)** at
genuinely problematic physical nodes (`invzp_stage3_plan.md` SS3.4: `Gstat` from -475 to -44
across three sampled real-coupling nodes) — demanding an *absolute* residual of `1e-10` on a
quantity of magnitude `100`-`500` asks for 12-14 significant digits of closure, which is not a
meaningful or achievable standard for a Picard iteration whose own documented conditioning is
`resid/|dK0| ~ 1.8e5` (`invz_emt_static_ordered.m` amendment-round-1 note). The **combined**
form `abs + rel*|Gstat|` is the standard atol+rtol convention (as in, e.g., ODE solvers):
it demands 10 correct digits when `Gstat = O(1)` (far from the pole) and relaxes
**proportionally** to `Gstat`'s own magnitude nearer the pole, without weakening the
*relative* precision demanded (still `1e-10` relative either way). **Both the abs and rel
component are the identical, single, pre-existing constant** — nothing here is tuned to any
run.
**Additional diagnostics exposed (never folded away):** `so.converged` (`res.blockB.converged`
— the closure's own PURE-ABSOLUTE flag, `so.resid < rtol_B`; this can legitimately be `false`
while `res.blockB.pass` is `true` for a large-`|Gstat|` node — that is the combined gate doing
its job, not a contradiction); `D_uni = so.D_uni` and `Dq = 1+(Jnu_flat-K0_b).*Gstat_b`
(`Dq_min`, `Dq_max`) are surfaced at the **top level** of `res` (Sec. 6), per the physics note
in `stage2c-context.md`: a sign-changing `Dq`/`D_uni` may be **physical instability**, not
noise, and must never be smoothed over or folded into a single pass/fail bit.
**Pass:** `isfinite(rB) && rB < scale_abs + scale_rel*abs(Gstat_b)`.

### Block C — Sigma self-consistency of the derived lambda/Sigma chain

```
lam_check = invz_lambdas(K_exp, g, wts, beta, [1 2 3])     % FRESH from K_exp, NOT state.lam
sg        = invz_sigma_ordered(tl, lam_check, K_exp, g, beta)
rC = max| sg.Sigma - Sigma_exp |
```
**Recomputation point:** `lam` and `sg` are recomputed **fresh from `state.K`**, exactly
mirroring production's own `final_resid` computation
(`invz_solve_point_ordered.m:232-234`) — `state.lam` is **not** read here (lam is derived, Sec.
1; this block is the one place production already treats it that way, and this contract keeps
that convention rather than inventing a new one). This block is the "keep it, name it" block
the brief specifies: it is the existing `final_resid`, given a name and independent-scale
statement, not new behaviour.
**Units:** dimensionless (same argument as Block A).
**Scale:** `scale_abs = scale_rel = tol_outer` (identical reasoning to Block A: this measures
the same physical Sigma-space quantity via a different recomputation path).
**Pass:** `isfinite(rC) && rC < scale_abs`.

### Block D — dynamic EMT identity (the block `final_resid` omits)

```
med_D = invz_emt_scalar(G0, Sigma_exp, Jnu_flat, eopts)     % FRESH from Sigma_exp
rD = max| K_exp(2:end) - med_D.K(2:end) |
```
**Recomputation point:** a fresh `invz_emt_scalar` call at the exported `Sigma`, independent of
Block A's internal step (1) (same non-sharing discipline as Block B). **`K(1)` is EXCLUDED by
design**, not by oversight: `med_D.K(1)` is the **ordinary-Dyson** static value, whereas the
exported `K(1)` is the **elastic-hybrid** `K0s` from `invz_gstat_ordered` — these are two
different physical quantities by construction
(`invz_emt_static_ordered.m` docstring, lines 3-4: "the elastic `G(0)` ... breaks the ordinary
Dyson structure, so the closed-form direct solve of INVZ_EMT_SCALAR does not apply at `w=0`
for `m != 0`"). This is the
same principle `test_invz_ordered_jensen.m` already pins in `G`-space
(`verifyGreaterThan(..., abs(pt.G(1) - G_dyson)/abs(G_dyson), 1e-6)`, `G_dyson` built from the
ordinary Dyson formula fed the SAME `K(1)`) — that test shows the two *formulas* disagree;
`test_invz_ordered_residual.m` independently confirms the same disagreement directly in
`K`-space on this fixture (`medD.K(1) = 4.03e-3` vs the exported `K(1) = 3.62e-3`, Sec. 9),
a related but distinct comparison. Comparing them here (in Block D) would be comparing two
different physical quantities, not measuring a residual.
**Units:** meV (`K` is meV, Sec. 1).
**Scale:** `scale_abs = tol_outer * Jscale`, `scale_rel = tol_outer`, where
`Jscale = max(abs(Jnu_flat))` (meV) — a **problem-native** physical scale, always available
from the node's own coupling spectrum, never invented per-run. **Combined gate:**
`rD < scale_abs + scale_rel * max(abs(med_D.K(2:end)))`.
**Justification (units-based, not tuned):** `K` is a `Jnu`-weighted effective coupling
(`invz_emt_scalar.m`: `K = A.*D./(1-A.*G0)`, `A` a `Jnu`-weighted average), so its own natural
magnitude scale is set by the coupling spectrum it was averaged from. Non-dimensionalizing by
`Jscale` turns `tol_outer` (already the codebase's own dimensionless closure standard, reused
verbatim from Blocks A/C, not a new number) into a K-space absolute floor with the correct
units (meV); the relative term additionally guards against a degenerate `Jscale -> 0` and
scales correctly if the coupling spectrum's overall magnitude changes between problems (e.g. a
future rescaled ion or a differently-normalized `Jnu`) — a pure `abs`-only gate would not.
**Finite requirement (block-local, per the brief's own Block-D text, which names exactly
these three quantities):** `all(isfinite(K_exp)) && all(isfinite(lam_exp)) &&
all(isfinite(Sigma_exp))` is folded into `res.blockD.pass` directly (this is narrower than,
and independent of, the global `res.finite` in Sec. 7, which additionally requires `K0s`
finite — the two checks are deliberately not identical; see Sec. 7).
**Degenerate-size guard:** if `numel(K_exp) <= 1` (no dynamic frequencies exist to check —
never the case for any Ecut-derived Matsubara grid used in this repo, but guarded rather than
left to crash), `rD` is defined as `0` and the block passes trivially (subject to the finite
clause and `med_D.converged`).
**Pass:** `isfinite(rD) && rD < scale_abs + scale_rel*max(abs(med_D.K(2:end)))
&& all(isfinite(K_exp)) && all(isfinite(lam_exp)) && all(isfinite(Sigma_exp))
&& med_D.converged`.

## 5. Why lambda is not a fifth, independently-gated block

`state.lam` is read as an input by Block B (Sec. 4) and by Block A's internal step (2), but no
block asserts `state.lam == invz_lambdas(state.K, ...)` directly. This is intentional, not an
oversight: per Sec. 1, `lam` is **derived**, and "do NOT describe a derived quantity as an
independent residual component" is the plan's explicit rule. Concretely: **Block B cannot
detect a `lam` corrupted independently of `K`/`Sigma`**, because `invz_emt_static_ordered`'s
own internal Picard loop treats `lam`/`Sigma0` as **fixed parameters** and converges `K0` to
self-consistency **for whatever `lam` it is given** — `so.resid` measures "did this closure
converge", not "was `lam` the physically correct value". This mirrors a pre-existing property
of the production code itself: `pt.lambda` (the loop's carried `lam`) is **never** compared
against the post-loop `lam_check` production computes for `final_resid`
(`invz_solve_point_ordered.m:232`) — the two are formally different variables there too. This
contract does not widen that scope. A future task may choose to add a diagnostic
`lam_check`-vs-`state.lam` comparison; it is explicitly **out of scope here** (Task 1a
implements exactly the four blocks the brief specifies).

## 6. Aggregate verdict

```
res.accepted = res.finite && blockA.pass && blockB.pass && blockC.pass && blockD.pass
```
Acceptance is a **logical AND of four independently-scaled pass/fail flags** — never a single
blended/weighted norm. A blended norm risks exactly the units-mismatch failure mode Block B's
own `O(100)` excursions warn against (a tiny absolute residual in one block's units could mask
a large relative failure in another's if combined into one number). Each block's raw residual,
scale, and pass flag are exposed separately (`res.blockA` ... `res.blockD`); `res.accepted` is
the only aggregate, and it is boolean, not a norm.

## 7. `finite` / `max_iter` / `stall` semantics

**`finite`** (global, `res.finite`): `all(isfinite(Sigma)) && all(isfinite(K)) &&
all(isfinite(lam)) && isfinite(K0s) && isfinite(rA) && isfinite(rB) && isfinite(rC) &&
isfinite(rD)` — **all four state components AND all four raw block residuals**. This is
deliberately broader than Block D's own block-local finite clause (Sec. 4), which names only
`K`/`lam`/`Sigma` (matching the brief's own, separately-scoped Block-D sentence) and omits
`K0s`. The asymmetry is faithful to two distinct statements in the brief, not an
inconsistency: Block D's own pass/fail is about the **dynamic** identity specifically (which
never uses `K0s`); the global `finite` is about the **whole exported state**, which does
include the continuation seed.

**`max_iter`** (a property of the CALLING LOOP's history, not of one exported-state snapshot —
**not computed by this checker**, which has no visibility into an outer-loop counter): defined
here for Task 1b's use. `max_iter := ` the calling loop's own outer index reached
`opts.max_outer` **without** `res.accepted` (as this contract defines it) ever having been
true. Task 1b's solver owns tracking this; it must not substitute `dS < tol_outer` or
`sout.converged` alone as a stand-in for `accepted` when deciding whether `max_iter` is the
honest verdict.

**`stall`** (diagnostic only, `res.stall`; **NEVER** used to set `res.accepted`): contrasts the
OLD, incomplete acceptance signal (`dS`, the live loop's own pre-mix step) against the NEW,
complete one (`res.accepted`). Requires the caller to supply `opts.dS` (the loop's own last
`dS`) — this is not reconstructable from the exported state alone (`dS` is an iteration-history
quantity, not part of `state`). If `opts.dS` is omitted (default `NaN`), `res.stall = NaN`
("undecidable without `dS`"). Otherwise:
```
res.stall = isfinite(opts.dS) && (opts.dS < tol_outer) && ~res.accepted
```
i.e. **stall = true** means "the OLD gate would have accepted this node (small `dS`), but the
COMPLETE residual says it should not have" — precisely the failure mode this task exists to
surface. `stall` never gates `res.accepted`; it is reported for Task 1b's diagnostics only.

## 8. Non-mutation rule (contract invariant)

`invz_ordered_residual(node, state, opts)` **never writes into any field of `node` or
`state`**. Every quantity it needs is copied into a local variable before use; the four blocks
recompute everything from those local copies. This is a hard invariant, not an optimization:
`isequaln(node_before, node_after)` and `isequaln(state_before, state_after)` must hold across
every call, including calls where one or more blocks fail or a helper throws. The test suite
(`invz_projected/tests/test_invz_ordered_residual.m`) asserts this directly.

## 9. Empirical confirmation on the synthetic control (not a tuning exercise — a check that
the declared, units-derived scales above are actually satisfiable by a genuinely converged
state, and that the perturbation tests are non-vacuous)

Fixture: `T=0.31`, `B=[2.85 0 0]` T, `Jnu=linspace(-2e-3,6e-3,24).'`, `J0eff=6.42e-3`,
`Jxx0=ion.Jxx0`, `ordered_mode='jensen'` (the Task-0 pinned synthetic control,
`test_invz_ordered_trace.m:89-93`). On this fixture (`invz_solve_point_ordered` converges,
`outer_iters=14`):

| Block | raw residual | combined gate | margin |
|---|---|---|---|
| A | `2.26e-9` | `1e-8` | ~4x |
| B | `6.82e-11` (`Gstat=-310`) | `1e-10 + 1e-10*310 = 3.11e-8` | ~450x |
| C | `2.25e-9` (`= pt.final_resid` exactly) | `1e-8` | ~4x |
| D | `1.26e-12` | `6.0e-11 + 2.7e-11 = 8.7e-11` | ~70x |

All four blocks pass; `res.accepted = true`. `med_D.K(1) = 4.03e-3` vs the exported
(elastic-hybrid) `K(1) = 3.62e-3` — confirmed to differ, as Sec. 4's Block-D design requires.

Two independent perturbations of this same accepted state (magnitudes chosen only to exceed
every gate above by several orders of magnitude — not tuned to graze any threshold):

- **`Sigma <- Sigma + 0.01`** (uniform shift, all `nw` entries): `blockA.pass=false`
  (`rA=1.13e-2`), `blockC.pass=false` (`rC=1.00e-2`), `blockD.pass=false` (`rD=8.1e-6`,
  collateral: `K` is Sigma-dependent through `invz_emt_scalar`'s closed-form solve, so a
  Sigma-space corruption legitimately cascades into a K-space disagreement too) — **but**
  `blockB.pass=true` still (`rB=9.0e-11`, well inside gate `3.0e-8`): Block B's own internal
  closure re-converges regardless of the (now-wrong) `Sigma(1)` it was handed, exactly as
  Sec. 4 predicts ("Block B measures internal self-consistency, not global correctness").
- **`K(2:end) <- K(2:end) + 1e-4`** (uniform shift, dynamic entries only, `Sigma`/`lam`/`K0s`
  untouched): `blockD.pass=false` (`rD=1.00e-4`, the direct target) and `blockC.pass=false`
  (`rC=8.44e-3`, collateral: `lam_check`/`sg` are recomputed from the now-perturbed `K`) —
  **but** `blockA.pass=true` and `blockB.pass=true`, **bit-identical** to the unperturbed
  values (`rA`/`rB` unchanged to the last printed digit): direct, empirical confirmation that
  neither block reads `state.K` at all (Sec. 3/4).

These two signatures are what `test_invz_ordered_residual.m` asserts.

## 10. Scope fence

This task implements the contract above and the checker that measures it. It does **not**:
change any acceptance gate in `invz_solve_point_ordered.m` or `invz_hmf_ordered.m`; merge the
two node loops; add seed rollback/cold-start retry; or alter any existing tolerance, sign, or
gate in production code. Task 1b consumes this contract to build the shared solver.
