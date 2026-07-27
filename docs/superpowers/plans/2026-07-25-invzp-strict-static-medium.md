# Strict-Order Static Medium Implementation Plan (stages 1–3)

**Status (Codex review, 2026-07-25): conditionally approved as an implementation layout, but
execution-blocked at Task 0.** The design choice and the three judgement-bearing constants
(`crit_tol`, `omit_promote`, `pole_cont_tol`) are approved. The revised, blind Stage-4 G8/G14/G16
protocol in §9 still requires explicit user approval, and Task 0 must record the baseline and exact
coupling fingerprint. No strict-mode production run may precede those remaining Task-0 steps.

This file is self-contained; it does not depend on a locally installed “superpowers” skill. Steps use
checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the beyond-order resummed `omega_n = 0` effective medium with a one-shot strict-`O(1/z)` moment closure in both the PM and ordered legs, behind a non-default opt-in flag, and measure whether the solved ordered `h`-path stays inside its domain.

**Architecture:** Shared primitives in `invz_common/` compute coupling moments, construct the
reference propagator, evaluate `K0 = Jbar - mu2*Gref`, resolve the scheme, classify recoverable
errors, and centralize the outer “catch only whitelisted signals” boundary. Both existing medium
leaves (`invz_emt_scalar` at `omega = 0`, `invz_emt_static_ordered`) gain a scheme-gated strict
branch that calls those primitives. The residual checker's block B is revised in place. Nothing
numerical changes on the default `'resummed'` path.

**Tech Stack:** MATLAB R2025a, `matlab.unittest` function-based tests (`functiontests(localfunctions)`), no toolboxes beyond base.

**Spec:** `docs/superpowers/specs/2026-07-25-invzp-ordered-solver-static-medium-design.md`. Section references below (`§4.1`, `§6.0`, `G17`, …) are to that spec.

**Out of scope (stage 4, a follow-on plan):** G6/G6d, G8, G10, the end-to-end part of G12,
G14, G16; the `one_field` un-nesting (§7 item 4, recommended deferred); the `[nJ,nw]`
**ordered-retarded** flattening fix (§7 item 5, split out — this plan rejects that ordered
combination). The PM leaf continues to support `[nJ,nw]`: its strict static slot uses column 1 and
leaves columns 2:end bit-identical.

## Global Constraints

- Conventions: `G = -chi` (meV^-1), ferromagnetic positive `J`. Copy verbatim into new docstrings.
- **Production default stays `static_medium = 'resummed'`.** No task in this plan flips a default.
- **The durable gotcha:** never broaden or regularize the static response, add a pole regularizer, flip a sign, or widen a tolerance as a convergence patch. Exact algebraic reassociation is permitted (Task 6); nothing else is.
- Error policy: `invz:*` identifiers only. New ids introduced here: `invz:staticMedium`, `invz:staticMediumConflict` is **NOT** introduced (§4.5 — no ingestion path exists).
- All moments use **population** normalization `mu_n = mean((J-Jbar).^n)`. MATLAB's `var(J)` (N-1) is forbidden: the difference is 4% at the `N = 24` synthetic fixtures.
- Tests live in `invz_projected/tests/`, named `test_<function>.m`, using the `functiontests(localfunctions)` + `setupOnce` addpath idiom of `invz_projected/tests/test_invz_ordered_residual.m:1-10`.
- Fast suite: `runtests('invz_projected/tests')`. At execution start, record a baseline from the
  **current dirty worktree**, not from the historical HEAD count. Every task must leave FAILED=0 and
  must not turn a previously passing test into failed/incomplete. New slow tests may add incompletes
  when `INVZ_SLOW` is absent; compare named tests as well as aggregate counts.
- MATLAB invocation: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "<cmd>"` from the repo root.
- The worktree contains unrelated user edits (`invz_tensor/*`, `invz_projected/invz_run_spectra.m`,
  READMEs, session docs, figures). **Never stage them.** Every commit stages exact paths only — never
  `git add -A`, never `git add .`.
- Do not add an authorship trailer for a model/person that did not author the implementation. The
  recent branch convention is no trailer.
- A task may not commit a knowingly failing test pending a later task. Where production units are
  mutually dependent, land them in one atomic task or defer the integration test until both sides are
  wired. “Expected to fail until Task N” is not compatible with the per-task FAILED=0 rule.

---

### Task 0: Freeze the preregistration and settle the open scope decisions

Documentation plus one shared primitive — the exact-byte digest that freezes the coupling fingerprint, created here because Task 0 is its first consumer. **This is a sign-off gate: no later task may run in strict mode until the user has approved this file.** The spec is deliberately unexecutable until it exists (§6.0).

**Files:**
- Create: `docs/invzp_strict_medium_prereg.md`
- Create: `invz_common/invz_exact_numeric_digest.m`
- Test: `invz_projected/tests/test_invz_exact_numeric_digest.m`
- Reference (read, do not modify): the spec §6.0, §7

**Interfaces:**
- Produces: the frozen constants every later task cites by name — `crit_tol`, `D_tol`, `Dq_tol`,
  `K_atol`, `K_rtol`, `ref_floor`, `omit_report`, `omit_promote`, `Delta_floor`,
  `pole_cont_tol`; the deferred-stage G8/G14/G16 bounds; and the five §7 decisions. Also
  produces `d = invz_exact_numeric_digest(x)` (a lowercase-hex char row vector), the shared
  exact-byte SHA-256 primitive used by Tasks 17 and 18.

- [ ] **Step 1: Write the preregistration document**

Draft `docs/invzp_strict_medium_prereg.md` from the proposal below, present the judgement-bearing
entries to the user, and replace `PROPOSED` by `FROZEN` only after explicit approval. Some scales are
inherited; the stability and asymptotic-control bounds are methodological choices and must be labelled
as such. Claiming that every number is “derived” would be false.

**Approval status as of 2026-07-25.** `crit_tol = 1e-6`, `omit_promote = 0.10` and
`pole_cont_tol = 1e-3` are FROZEN by explicit user approval and carry that label below. The §9
Stage-4 blind table is likewise populated below and requires its own explicit approval before Task 1
begins — it is the last outstanding item. Nothing else in this file may be changed by editing in
place; a value that needs revisiting gets a new dated re-registration section.

```markdown
# Strict-order static medium — FROZEN preregistration

Frozen: <DATE OF USER APPROVAL>. Branch invzp-stage2c-diagnostic. Spec:
docs/superpowers/specs/2026-07-25-invzp-ordered-solver-static-medium-design.md §6.0.
Changing any value below requires a new dated re-registration section, never an edit in place.

## Inherited scales (not new numbers)
tol_outer      = 1e-8      invz_ordered_node_solve.m:129 default
resid_tol      = 1e-10     invz_emt_static_ordered.m:38 default
Jscale         = max(abs(Jnu_flat))   invz_ordered_residual.m:88; = 6.7631e-3 meV on the
                                      production multiset (grid [16 16 16], dpRng 30, bruteforce)
dB_boundary    = 0.02 T    invz_critical.m:21 crossing tolerance -- the field resolution at
                           which this project already declares a boundary located
Delta_floor    = 1e-4 meV  invz_twolevel_ordered.m:18, reused verbatim

## 1. Endpoint stability tolerances (spec §1 two-tier verdict, G4)
crit_tol = 1e-6 FROZEN 2026-07-25 (user approval)
                           dimensionless. 100x tol_outer, so the SIGN of crit is resolved well
                           above outer-loop noise while staying far below crit's O(1) physical
                           scale away from the boundary.
D_tol    = 1e-6 * max(1, abs(Gstat)*Jscale)      state-dependent, evaluated per node.
Dq_tol   = D_tol           same construction and units; Dq >= D_uni holds on this multiset
                           (spec §B), so a common floor is consistent.
Rationale for the D_uni/Dq scaling: their absolute noise scales with |Gstat| (which reaches
O(1e2-1e3)), unlike crit which is O(1). A flat absolute floor would be meaningless at both ends.

## 2. Strict block-B gate and domain margins (spec §4.4, §5.2, G0)
K_atol     = 1e-14 meV
K_rtol     = 1e-12
           gate = K_atol + K_rtol*max(|K0|, |Kstrict|, Jscale)  ~= 1.0e-14 meV in production,
           i.e. ~1e4 x eps(Jscale). A one-shot formula recomputed from the exported state should
           agree to floating-point reassociation only; this gate catches mis-wiring, not physics.
           NOT eps(Jbar) (2.7e-20, fails on any reassociation) and NOT mu2*|Gref| (vanishes as
           Gref -> 0).
ref_floor  = 1e-6          1 + Sigma0 at or below this is 'ref_denom_small'. 1+Sigma0 is O(1)
                           because Sigma0 is O(1/z); 1e-6 means the Dyson reference is
                           degenerate, i.e. genuinely out of domain.
Boundary-indeterminate rule (spec §4.4 case 3): a PM probe is boundary_band when
           |crit_pm| <= crit_tol.
           Field resolution is represented separately by the returned Bc interval. Do NOT use
           |d crit/dB|*dB_boundary inside one_field: a parfor field point does not own a sweep
           derivative, and doing so would make phase classification grid-dependent.

## 3. Gate-0 negative-outcome rule (spec §6.0(3), verbatim intent)
Promotion FAILS if any of:
  (a) any required solved-path node has a non-'ok' REFERENCE denominator status;
  (b) any skipped or invalid node is unaccounted for in the coverage counters;
  (c) max(omit_max) over the solved path exceeds omit_promote (below);
  (d) any LOCAL Gstat denominator crossing at which r or crit is non-finite or discontinuous
      (G17 covers the algebra; this covers the actual path);
  (e) any required ordered field does not return status='ok', a finite nonzero root, and a
      stable endpoint under the frozen crit/D_uni/Dq margins, or either PM control fails to
      return a converged finite positive-mass PM state.
A local Gstat crossing that satisfies G17 does NOT fail promotion (spec §1: the singularity is
removable in the integrand).
On failure the run STOPS AT DIAGNOSIS. Carrying another moment, changing Gref, or truncating
other Matsubara sectors is a NEW theory candidate requiring a new spec and fresh
preregistration -- never an in-run fallback. Regularisation, broadening and tolerance widening
remain forbidden.

## 4. Omitted-order thresholds (spec §4.1, §B)
omit_report  = 0.05        flag/report above 5%. This is a reporting convention, not a node
                           rejection threshold.
omit_promote = 0.10 FROZEN 2026-07-25 (user approval)
                           candidate-promotion bound on max(omit_max), over the ACTUAL SOLVED
                           path (the bare scan is only prospective).
Rationale: 0.10 is an explicit asymptotic-control choice, not something derivable from the
0.3% PM boundary shift. The binding caution in the spec says that boundary shift bounds neither
integral Sigma0 dh nor integral (r-1) dh deep in the ordered phase. The earlier proposed 0.25
argument incorrectly multiplied that boundary shift by an omitted-term ratio and is rejected.
If the user prefers another bound, freeze it here before the first strict run and keep it fixed.
Both ratios (omit_mu3 AND omit_cubic) enter omit_max: mu3*Gref^2 is the FIRST omitted term,
and its near-zero value is a measured property of ONE multiset.

## 5. Numerical convergence tolerances (G1, G5)
G1 identity errors:
   |dm/dh + G0bare|                        <= 1e-6 * max(1, |G0bare|)
   |Delta F/Delta h - trapz_avg(crit)|     <= 1e-6 * max(1, |crit|)
   |dF/dm - crit/chi_path|                 <= 1e-6 * max(1, |crit/chi_path|)
   A synthetic smooth, nested-grid quadrature fixture must show second-order convergence.
   The adaptive/geometric production profile is required to converge monotonically and meet
   the finest-grid tolerance; do not require a [3,5] error ratio from non-nested adaptive grids.
G5 path integrals: int_Sigma0 and int_r_minus_1 stable under nH 33/65/129 using
   |I_fine-I_prev| <= I_atol + 1e-3*max(|I_fine|,|I_prev|), I_atol = 1e-10 meV.
   A relative-only rule is invalid when an integral crosses or approaches zero.
G17 actual-path crossing continuity:
   pole_cont_tol = 1e-3 FROZEN 2026-07-25 (user approval), dimensionless relative jump.
   Export d_local = 1+Sigma0+K0*G0inel0. For each sign change of d_local, estimate its
   zero by linear interpolation, extrapolate r and crit to that same h from the nearest
   two finite nodes on each side, and require
       |y_left-y_right|/max(1,|y_left|,|y_right|) <= pole_cont_tol
   EXACT-ZERO NODE: when d_local is exactly 0 at a grid node, that node is a direct evaluation
   AT the crossing and is therefore the strongest available evidence, not a special case to be
   waived. Build y_left and y_right from the nearest two finite nodes on each side EXCLUDING
   it, then require ALL THREE relative differences -- |y_left-y_right|, |y_node-y_left| and
   |y_node-y_right|, each normalised by max(1,|.|,|.|) of its own pair -- to be <= pole_cont_tol.
   Requiring only that y_node be finite would admit a finite but arbitrarily wrong value sitting
   between two mutually consistent extrapolants.
   for y in {r,crit}. The finest-grid jump must not increase from nH=65 to 129.
   If fewer than two finite nodes exist on either side, the crossing is unresolved and
   Gate 0 fails; it is not silently omitted.

## 6. Scheme comparison (G7, non-gating)
Report K(1)_strict - K(1)_resummed (the exact scheme jump, since omega_n != 0 is unchanged by
construction) against the physical dispersion K(2) - K(1), versus T. No pass threshold: this is
a measurement, and neither quantity is assumed negligible.

## 7. Settled scope decisions (spec §7)
7.1 (item 1) These values are the preregistration. Filled and frozen here.
7.2 (item 2) invz_static_domain_scan accepts an EXPLICIT hgrid. The initial geometric grid is
    factored into invz_hmf_grid (shared by HMF and the scanner). Solved-path margins are read
    off prof.hgrid AFTER adaptation. "Grid identity" = shared initial helper + complete
    accounting of every actual solved node, not two algorithms that agree in one test.
7.3 (item 3) G2 uses K_atol/K_rtol, NOT bitwise. The two callers reach Gref through different
    expressions (G0_PM(0)/(1+Sigma0) vs G0bare0/(1+Sigma0), equal only because G0el0 -> 0), so
    bitwise identity is not a property the design guarantees.
7.4 (item 4) DEFER the one_field un-nesting to its own change. Keep only the three-way
    dispatcher, which IS required (it is what stops solver availability from labelling a
    phase). The strict dispatcher handles both successful auto states (`phase` 1 and 2);
    what is deferred is running the 1/z leg when auto itself returns `phase=0`. Rationale:
    below Bc_1z the bare set orders so phase == 1 holds, and a separate transverse-PM branch
    already exists at invz_spectra_map.m:372+, so the remaining auto-failure un-nesting is
    robustness, not part of the masking fix.
7.5 (item 5) The resolver REJECTS strict mode combined with a [nJ,nw] retardation Jnu_flat on
    the ordered leg (invz:staticMedium). The PM leg REMAINS supported and uses Jnu_flat(:,1)
    plus moment column 1 for omega_n=0; K/G slots 2:end retain the existing matrix branch.
    The pre-existing invz_emt_static_ordered.m:43 all-column flattening is split out.

## 8. Frozen stage-3 grids (must be approved with the constants)
Coupling fixture: grid [16 16 16], dpRng 30, dipole='bruteforce', cache=false, with no
explicit grid-policy fields. On this legacy/absent-policy route `info.grid` is intentionally
ABSENT (`invz_bz_couplings.m:26`); record the requested tuple, assert that absence, record
`info.dipole`, and record invz_exact_numeric_digest(Jnu(:)) -- the shared exact-byte digest of
class + shape + data. Tasks 17 and 18 must call that same primitive, never a reimplementation.
Gate-0/G5: T=0.10 K; required ordered fields [0.05 0.25 0.5 1 2 2.5 2.9 3.0] T;
separate PM controls [3.1 3.5] T; exact B=0 is a labelled hard-domain test, not a required
ordered-path node. nH=[33 65 129], Ecut=40 meV.
G7: T=[0.05 0.10 0.31 1.0] K, B=6 T, Ecut=40 meV, same coupling fixture.
Changing these grids after seeing strict output requires a dated preregistration amendment.

## 9. Blind freeze for deferred Stage 4 (BLIND: fixed before any strict output exists)
The design spec §6.0 requires these values before the FIRST strict run, not when Stage 4 begins.
Stage-4 IMPLEMENTATION stays out of scope; only its acceptance bounds are frozen here, so they
cannot be chosen once Gate-0 output is visible. Every entry below is either INHERITED (a value
already in the code, cited file:line), DERIVED (follows from an inherited value), or CHOICE (a
methodological decision with its reasoning stated). No entry is derived from strict output,
because none exists yet.

### 9.0 Resolution scales these bounds are built from
dB_boundary = 0.02 T          INHERITED, invz_critical.m:21 -- the field resolution at which
                              this project already declares a boundary located.
eta_prod    = 5e-5 meV        INHERITED, invz_run_spectra.m:42. NOTE the library default is
                              5e-3 meV (invz_spectra_map.m:129, invz_chi_realaxis.m:22) but the
                              production driver forwards 5e-5 and overrides it. G16 freezes the
                              value PRODUCTION actually uses; a later eta change must not be
                              readable as physics.
w_res       = 2*eta_prod
            = 1.0e-4 meV      DERIVED. A Lorentzian of half-width eta has FWHM 2*eta. There is
                              NO linewidth/FWHM code anywhere in the repo (invz_peak_energy.m is
                              peak-position only), so the broadening INPUT is the only defensible
                              resolution scale; it is not measured from output.
dw          = 0.01 GHz
            = 4.136e-5 meV    DERIVED from invz_run_spectra.m:41 (w = (0:0.01:5.5).', GHz) and
                              the GHz->meV conversion. w_res spans only ~2.4 grid points, so the
                              peak is marginally resolved; G16 must report that, not hide it.
I_atol      = 1e-10 meV       CHOICE, shared with G5 -- a relative-only rule is invalid where an
                              integral crosses zero, and int_r_minus_1 can.
chi_atol    = 1e-8 meV^-1     CHOICE, a numerical-zero floor for spectral-vector comparisons.
                              It is used only as sqrt(nw)*chi_atol in an L2 gate and never as a
                              broadening, censor, or replacement for a physical response.
h_atol      = 1e-10 meV       CHOICE, the absolute floor for comparing roots near zero.
conv_rel    = 1e-3            CHOICE, the common relative refinement tolerance in G14.

Unless a G14 ladder explicitly overrides one factor, Stage 4 uses the §8 production coupling
fixture (`N=16`, dpRng=30, bruteforce, no grid-policy fields), `Ecut=40 meV`, `nH=33`,
`tol_root=1e-3`, `hyp=true`, `field_dir=[1 0 0]`, `transverse_mf='legacy_x'`, and `demag=0`.
These are the user-facing driver's intrinsic transverse benchmark settings, not values selected
from strict output.

### 9.1 G8 -- strict vs resummed, COMMON STABLE STATES ONLY
This gate is deliberately **not** evaluated on the unstable ordered low-h interval where the
resummed theory has no value. It uses a frozen common-stable PM response path, so "common stable"
cannot become a post-output choice:

  common scalar paths   production coupling fixture (§8), T=0.10 K, Ecut=40 meV,
                        B in {3.5, 6.0} T, and
                        hmax = 8*J0eff,
                        hgrid = [0, hmax*10.^linspace(-6,0,128)] meV.             CHOICE
                        The factor 8 is the Ho electronic J ceiling. Evaluate both schemes at
                        every identical (B,h) node. Every node must return an accepted,
                        finite state with `crit>crit_tol`, `D_uni>D_tol`, and
                        `Dq_min>Dq_tol` under both schemes; otherwise G8 is `uncontrolled`,
                        not silently shortened to the surviving tail.
  equilibrium spectra  B in {3.1, 3.5} T on the §9.3 frequency grid and eta_prod. CHOICE
                        Both PM points must be accepted under both schemes. The denominator-pole
                        tracker is the §9.3 tracker; invz_peak_energy is report-only.
  boundary              compute each scheme's PM mass root in window [2.5,3.5] T with
                        fieldstep=0.05 T and root interval <=0.005 T.             CHOICE

For every scalar formula below, suffixes `s` and `r` mean strict and resummed at the **same**
state. No unsuffixed quantity may be used as an implicit denominator:

  K0          |K0_s-K0_r| <= 1e-12
                            + omit_promote*mu2*max(|Gref_s|,|Gref_r|)
              CHOICE(floor)/DERIVED(relative). The 1e-12 meV floor is far above reassociation
              noise and far below the retained medium correction.
  Sigma0      |Sigma0_s-Sigma0_r| <= 1e-10
                                      + omit_promote*max(|Sigma0_s|,|Sigma0_r|)
  mass        |crit_s-crit_r| and |Duni_s-Duni_r| <=
                 crit_tol + omit_promote*
                 max(|r_s-1|,|r_r-1|,|Sigma0_s|,|Sigma0_r|)
              crit's scheme dependence enters through the medium-carrying part; the absolute
              floor is required where a mass approaches zero.
  integrals   for I in {int_Sigma0,int_r_minus_1},
                 |I_s-I_r| <= I_atol + omit_promote*max(|I_s|,|I_r|)
              Both integrals use the complete identical hgrid and endpoints above.
  Bc_1z       |Bc_s-Bc_r| <= dB_boundary = 0.02 T.
  spectra     (i)  |E_pole_s-E_pole_r| <= w_res;
              (ii) ||Imchi_s-Imchi_r||_2 <= sqrt(nw)*chi_atol
                       + 0.05*max(||Imchi_s||_2,||Imchi_r||_2).       CHOICE
              Both must hold. The absolute term makes the norm executable at numerical zero;
              the symmetric relative term does not privilege either scheme.

Ordered equilibrium spectra remain the subject of G16, not a G8 comparison: the resummed
ordered leg's failure is the defect under repair and cannot be converted into a "common state" by
dropping failed nodes. G8 must report its exact attempted/common counts and all rejection reasons.

### 9.2 G14 -- convergence ladders and acceptance
G14 is a one-factor-at-a-time study. Except for the factor being refined, hold all controls at
the finest values `N=24`, `nH=129`, `tol_root=1e-4`, `Ecut=80 meV`, boundary spacing
`dB=0.01 T`, the production bruteforce/dpRng=30 coupling convention, and eta_prod. Compare every
adjacent ladder pair; every pair containing a production value (`N=16`, `nH=33`,
`tol_root=1e-3`, or `Ecut=40`) is load-bearing, not merely printed. This protocol is a one-off
default-flip gate, so cost is not a reason to mix refinements.

  forward/reverse B   serially traverse the §8 ordered list ascending and descending. CHOICE
                      The first field in each direction is cold; every later field receives the
                      preceding accepted ordered root/profile as its continuation seed. Merely
                      reversing an independently solved/parfor list does NOT satisfy this gate.
                      If the public solver has no seed interface, Stage 4 adds a diagnostic-only
                      seed input without changing default dispatch. Run cold endpoint controls.
  field spacing       let F0=linspace(0,9,101), the INHERITED production grid, and define the
                      nested CHOICE ladder
                        F1=F0 union (2.5:0.04:3.5),
                        F2=F1 union (2.5:0.02:3.5),
                        F3=F2 union (2.5:0.01:3.5).
                      Sort each set and coalesce values within 1e-12 T, retaining an existing
                      coarser/F0 value on a tie. This avoids binary-colon near-duplicates and
                      makes the production-to-refinement comparison load-bearing while retaining
                      identical common fields. The 0.01-T finest spacing is finer than
                      dB_boundary; the former 0.045-T endpoint was rejected because it was still
                      coarser.
  nH                  {33,65,129}.                                      INHERITED/CHOICE
  tol_root            {1e-3,3e-4,1e-4}. The first value is inherited from
                      invz_hmf_ordered.m:51; the two refinements are CHOICE.
  Ecut / Matsubara    {40,60,80} meV. The first value is inherited from
                      invz_hmf_ordered.m:114; the two refinements are CHOICE.
  q-grid              N in {12,16,24}, dpRng=30, same legacy grid convention. CHOICE

At the §8 ordered and PM fields, and at common nodes of the nested boundary grids, apply the
following symmetric formulas to every adjacent pair `(a,b)`; treat the forward/reverse results as
one additional `(a,b)` pair at each common field:

  hstar       |hstar_a-hstar_b| <= h_atol
                                      + conv_rel*max(|hstar_a|,|hstar_b|)
  mstar       |mstar_a-mstar_b| <= 1e-8
                                      + conv_rel*max(|mstar_a|,|mstar_b|)
  Sigma0      |Sigma0_a-Sigma0_b| <= 1e-10
                                      + conv_rel*max(|Sigma0_a|,|Sigma0_b|)
  crit        |crit_a-crit_b| <= crit_tol
                                      + conv_rel*max(|crit_a|,|crit_b|)
  D_uni       |Duni_a-Duni_b| <= max(Dtol_a,Dtol_b)
                                      + conv_rel*max(|Duni_a|,|Duni_b|)
  Dq_min      |Dqmin_a-Dqmin_b| <= max(Dqtol_a,Dqtol_b)
                                      + conv_rel*max(|Dqmin_a|,|Dqmin_b|)
  integrals   each G5 integral obeys I_atol
                                      + conv_rel*max(|I_a|,|I_b|)
  Bc_1z       |Bc_a-Bc_b| <= dB_boundary
  pole        |Epole_a-Epole_b| <= w_res at every common G16 field
  spectrum    ||Imchi_a-Imchi_b||_2 <= sqrt(nw)*chi_atol
                                      + conv_rel*max(||Imchi_a||_2,||Imchi_b||_2)

The q-grid ladder is subject to **all** formulas, not only Bc_1z. `(Jbar,mu2,mu3,mu4)` are also
reported for each N but are not forced to be invariant: they are properties of the changing input
quadrature. At every solve in every ladder, `n_accounted == numel(trc.nodes)` and every required
output above must be finite. One unaccounted/missing value fails G14 rather than shrinking the
comparison set.

### 9.3 G16 -- end-to-end symptom gate
G16 has two field grids at T=0.10 K:

  context grid        linspace(0,9,101) T. INHERITED from invz_run_spectra.m:40. It verifies
                      the production panel and complete phase/reason coverage.
  critical grid       (2.50:0.01:3.50) T. CHOICE, blind-centered on the already known legacy
                      3.025-T PM boundary and fine enough to resolve dB_boundary.
  frequency grid      w=(0:0.01:5.5).' GHz. INHERITED from invz_run_spectra.m:41 and converted
                      by the driver's wMeV=w/eScale. Record eScale. The stale inline
                      "0-108 GHz" comment is not authoritative; the literal array is.
  frequency ladder    dw in {0.01, 0.005, 0.0025} GHz over the same 0-5.5 GHz span, with
                      **eta held at eta_prod** on every rung.                          CHOICE
                      In meV: dw = 4.136e-5, 2.068e-5, 1.034e-5, so w_res = 1.0e-4 meV spans
                      about 2.4, 4.8 and 9.7 bins. The coarsest rung is the production grid, so
                      the ladder certifies it rather than replacing it. eta is FIXED across the
                      ladder deliberately: refining dw at fixed eta separates grid resolution
                      from broadening, whereas shrinking both together would confound them and
                      could not distinguish a sampling artifact from a physical linewidth.
  broadening          eta=eta_prod=5e-5 meV, held fixed and recorded. INHERITED from the
                      production driver, not the 5e-3 library default.

**Boundary/mass refinement.** The fixed grids are for coverage and spectral shape, not for
pretending that a 0.01- or 0.09-T sample lands in a 1e-6 mass band. Starting from the nearest
accepted ordered/PM bracket on the critical grid, bisect/continue each physical branch until:

  Bc interval width <= 0.005 T = dB_boundary/4,                         CHOICE
  ordered anchor: 0 < D_ord <= D_tol(anchor),
  PM anchor:      0 < crit_pm <= crit_tol.

Use at most 60 refinements per side. Every refinement node is ledger-accounted. Failure to obtain
either stable anchor is a G16 failure; do not widen a mass tolerance or substitute the adjacent
regular-grid column. Report the two anchors, the final Bc interval, and their distances from it.

**Denominator-pole tracker.** `invz_peak_energy` remains an unchanged presentation diagnostic and
is never the G16 verdict. Stage 4 must expose the intrinsic uniform response denominator before
the demagnetization transform,

  Dresp(w,B) = 1 - J0eff*chitilde(w,B),

on the frozen complex-frequency evaluation `w+i*eta_prod`.

`Dresp` is DEFINED here, blindly, before any strict output exists — it is not calibrated from
unseen results, which is precisely what a preregistration is for. What it does do is expand Stage
4: that stage must add a **diagnostic-only export** of `Dresp` together with algebraic and unit
tests (limiting forms, units, and agreement with the existing response at states where both are
defined), **without changing production response arithmetic**. The export is new observable
plumbing, not a new physics path; if exposing it would require altering how any production
response is computed, that is a design change requiring a new spec, not a Stage-4 implementation
detail.

> `D_ord` is the static Matsubara endpoint-stability mass, whereas `Dresp(w)` is the retarded
> dynamical pole denominator. They need not be numerically equal away from the transition and must
> never be substituted for one another; G16 requires their independently defined static and
> dynamic softening limits to coincide at the boundary.

At each field enumerate every interior
local minimum of `|Dresp|^2`, refine its frequency by a local quadratic vertex, and retain it as a
pole candidate only when the corresponding `Im chi` is finite and positive. Treat `w=0`
separately: retain it without vertex refinement when it is a one-sided local minimum with finite
`Dresp` and finite nonnegative `Im chi`. This enumerated list, not the global maximum, is the
tracking input.

Anchor independently at the ordered and PM refined mass points by the lowest-frequency candidate
in `[0,10*w_res]`; break a frequency tie by smaller `|Dresp|`, then lower grid index. Absence of
such a candidate fails. Track **outward from the boundary** on each side of the critical grid.
For the first step, choose the nearest candidate within `5*w_res`.
Thereafter linearly predict from the preceding two pole energies and choose the closest candidate
within

  track_window = max(5*w_res, 2*|E_prev-E_prevprev|).                  CHOICE

Break a prediction-distance tie by smaller `|Dresp|`, then lower frequency. If none exists, mark
the branch `lost` and fail; never reacquire the global maximum. Store the candidate index,
`|Dresp|`, energy, and response weight at every field so branch switching is auditable.
Scope note on indices: within a SINGLE grid an index is a legitimate deterministic tie-break (as
above) and an audit record. It is never a cross-grid identity — the frequency-ladder matching in
the resolution classification below pairs poles BY ENERGY, because an index scales with 1/dw and
refinement can reorder the local minima.

**Coverage and exclusions.** Exact B=0 is the sole hard-domain exclusion and remains a labelled
`degenerate_doublet` control. A point inside the final 0.005-T Bc interval may be labelled
`boundary_indeterminate`, but it is listed and replaced only by the two refined anchors above.
Outside that interval every context/critical-grid column must have `phase_1z` in `{1,2}`, a finite
spectrum, and an accounted reason. `pm_probe_unknown`, `solver_failed`, an auto-overlay failure
that suppresses an otherwise valid 1/z solve, or any silent exclusion **fails G16**.

**Two-sided continuous-softening predicate.** Let K=5 be the five valid critical-grid pole points
nearest to, but outside, the final Bc interval on each side (CHOICE). Sort the ordered points by
increasing B toward Bc and the PM points by increasing B away from Bc. Require all of:

1. approaching Bc from the ordered side, `E(i+1) <= E(i)+w_res`; moving away from Bc on the PM
   side, `E(i+1) >= E(i)-w_res`;
2. fit `E^2=a*B+c` independently on each side, require `a_ordered<0` and `a_PM>0`, and require
   every fitted residual
   `|E_i^2-(a*B_i+c)| <= 2*E_i*w_res + w_res^2` (propagated frozen energy resolution);
3. for each fitted zero `B0=-c/a`, require
   `dist(B0,[Bc_lo,Bc_hi]) <= dB_boundary`, and require the ordered and PM fitted zeros to differ
   by at most dB_boundary;
4. the two refined-anchor energies are each `<=w_res` and differ by `<=w_res`, evaluated under
   the frequency-ladder classification below rather than on the production grid alone;
5. the ordered anchor satisfies its positive D_tol margin and the PM anchor its positive
   crit_tol margin as specified above.

This is the load-bearing symptom test: both masses and the **same denominator branch** must close
from their own sides. Ordered-only softening or two masses that close while the spectral branches
jump cannot pass.

**Resolution classification.** `w_res` is only ~2.4 production bins, so a predicate failure on the
production grid alone is ambiguous: it may mean the branch does not soften, or merely that the grid
cannot see that it does. Run the frozen frequency ladder (`dw` in {0.01, 0.005, 0.0025} GHz,
`eta = eta_prod` on every rung) and apply the following procedure **in this order** — convergence is
established BEFORE the physics predicate is interpreted, because a predicate evaluated on
unconverged grids means nothing either way.

*Step A — grid convergence, across the two finest rungs.* Both must hold:

  branch matching   Match tracked poles ACROSS RUNGS BY ENERGY, never by index. For each side,
                    pair each tracked pole of the coarser rung with the tracked pole of the finer
                    rung nearest in energy, requiring the pairing to be ONE-TO-ONE (every tracked
                    pole matched exactly once, in both directions) and every matched pair to lie
                    within `max(dw_fine, dw_prev)` in meV. An unmatched pole, a many-to-one
                    assignment, or a differing tracked-pole count fails.
                    A raw frequency index CANNOT identify a branch across grids: it scales with
                    1/dw, and refinement can resolve, merge or reorder local minima, so even the
                    local-minimum ordinal changes. Indices are recorded as PROVENANCE ONLY and
                    never used as an identity or matching criterion.
  quantity          Every load-bearing pole/fit quantity agrees across the two finest rungs:
  convergence         - each matched tracked pole energy `E_i`, both sides, to
                        `max(dw_fine, dw_prev)`;
                      - both refined-anchor energies, to `max(dw_fine, dw_prev)`;
                      - `sign(a_ordered)` and `sign(a_PM)` identical on both rungs;
                      - each fitted zero `B0_ordered`, `B0_PM`, to `dB_boundary`.
                    Any of these failing to converge fails Step A.

*Step B — the complete predicate.* Only if Step A holds, evaluate G16 in full: predicate items
1-5, the two refined boundary/mass anchors, the coverage-and-exclusions requirements, and the
per-column `phase_1z`/reason status requirements. This is the COMPLETE gate; no subset of it
may stand in for the whole.

The class is then exactly one of:

  resolution_unresolved   Step A fails -- one-to-one continuation/energy branch matching fails, or
                          any load-bearing pole/fit quantity does not converge across the two
                          finest grids. The grid, not the physics, is the limiting factor, and the
                          predicate is NOT interpreted.
  physics_fail            Step A holds but Step B fails: the grids converge and any complete G16
                          item 1-5, anchor, coverage or status requirement fails. This is evidence
                          against the strict-medium physics.
  pass                    Step A and Step B both hold in full.

`resolution_unresolved` **blocks promotion** exactly as `physics_fail` does — an unresolved gate
certifies nothing. But it carries its own label and must **never** be presented as evidence against
the strict-medium construction: it is a statement about the frequency grid. The report states the
class explicitly, together with the per-rung matched-pole energies, anchor energies, fit
coefficients and fitted zeros that produced it — and the grid indices as provenance only — so the
two failure modes can never be conflated after the fact.

Every §9 entry is frozen on user approval together with §§1-8. Changing one afterwards requires a
new dated re-registration section, never an edit in place.

Every entry above is a number or an executable formula, each labelled INHERITED / DERIVED /
CHOICE. None was selected from Gate-0 or any other strict output, because none exists yet -- that
is what makes this freeze blind, and it is the whole reason it happens now rather than when
Stage 4 begins.
```

- [ ] **Step 2: Create the shared exact-byte digest primitive**

This is the ONE definition of the exact-byte SHA-256 digest that three consumers rely on: this
task's own coupling fingerprint (Step 3 below), Task 17's G11 real-coupling anchor, and Task 18's
Gate-0 driver. Reproducing the algorithm from prose in three places would let a silent mismatch
destroy the fingerprint's whole purpose and surface as a spurious Gate-0 abort.

Create `invz_projected/tests/test_invz_exact_numeric_digest.m`:

```matlab
function tests = test_invz_exact_numeric_digest
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_digest_is_64char_lowercase_hex_row(testCase)
d = invz_exact_numeric_digest([1 2 3]);
verifyEqual(testCase, size(d), [1 64]);
verifyTrue(testCase, ischar(d));
verifyTrue(testCase, all(ismember(d, '0123456789abcdef')));
end

function test_digest_is_stable_across_repeated_calls(testCase)
x = [1 2 3; 4 5 6];
verifyEqual(testCase, invz_exact_numeric_digest(x), invz_exact_numeric_digest(x));
end

% NOT shape-invariant by design: a reshaped multiset is a different input.
function test_digest_is_shape_sensitive(testCase)
verifyNotEqual(testCase, invz_exact_numeric_digest([1 2 3]), ...
    invz_exact_numeric_digest([1 2 3].'));
end

% NOT order-invariant by design: a reordered multiset is a different input.
function test_digest_is_order_sensitive(testCase)
verifyNotEqual(testCase, invz_exact_numeric_digest([1 2 3]), ...
    invz_exact_numeric_digest([3 2 1]));
end

function test_digest_is_class_sensitive(testCase)
verifyNotEqual(testCase, invz_exact_numeric_digest(int32([1 2 3])), ...
    invz_exact_numeric_digest(double([1 2 3])));
end

% REGRESSION on the "exact-byte" claim itself. These two int64 values are distinct but collapse
% to the SAME double (2^53 is the last integer with an exact double neighbour), so any
% implementation that hashes double(x) instead of the original class bytes digests them
% identically -- same class, same shape, same converted data. Hashing x(:) directly separates
% them. Without this test "exact-byte" would be true only for doubles.
function test_digest_hashes_original_class_bytes_not_doubles(testCase)
a = int64(9007199254740992);        % 2^53
b = int64(9007199254740993);        % 2^53 + 1; double(a) == double(b)
verifyEqual(testCase, double(a), double(b), 'AbsTol', 0);   % the collision is real
verifyNotEqual(testCase, invz_exact_numeric_digest(a), invz_exact_numeric_digest(b));
end

function test_nonnumeric_and_sparse_inputs_raise_exactDigest(testCase)
verifyError(testCase, @() invz_exact_numeric_digest('abc'), 'invz:exactDigest');
verifyError(testCase, @() invz_exact_numeric_digest(sparse([1 0 2])), 'invz:exactDigest');
verifyError(testCase, @() invz_exact_numeric_digest([1 2i]), 'invz:exactDigest');
end
```

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_exact_numeric_digest.m'); disp(table(r))"`
Expected before the implementation exists: all 7 error with "Unrecognized function or variable
'invz_exact_numeric_digest'".

Create `invz_common/invz_exact_numeric_digest.m`:

```matlab
function d = invz_exact_numeric_digest(x)
%INVZ_EXACT_NUMERIC_DIGEST Exact-byte SHA-256 of a numeric array's class, shape and data.
% ONE definition, THREE consumers: Task 0 freezes the coupling fingerprint, Task 17's G11
% anchor asserts it, and Task 18's Gate-0 driver re-checks it before any solve. Reproducing
% this algorithm from prose in three places would let a silent mismatch destroy the
% fingerprint's whole purpose and surface as a spurious Gate-0 abort.
%
% Digest input, in this exact order: the class name, the size vector as doubles, then the
% ORIGINAL class-specific bytes of the data in column-major order. Deliberately NOT
% order-invariant and NOT shape-invariant: a reshaped or reordered multiset is a different
% input and must digest differently.
%
% The data bytes are taken with typecast(x(:),'uint8'), NEVER via double(x). Converting first
% would make the digest exact only for doubles: int64(2^53) and int64(2^53+1) are distinct
% inputs that collapse to the same double, so a double-converting digest would call them equal
% (pinned by test_digest_hashes_original_class_bytes_not_doubles). The size vector is genuinely
% double because size() returns doubles; that is a canonical encoding, not a data conversion.
% Contract: FULL, REAL, numeric x. Sparse and complex inputs are rejected rather than silently
% densified or split into re/im, since either choice would be an unstated encoding decision.
%
% JVM REQUIREMENT: uses java.security.MessageDigest, so it does not work under -nojvm. This
% matches the existing precedent and constraint documented for invz_bz_couplings.m's local
% qhash helper; this repository is plotting-oriented and -nojvm is not a supported mode.
if ~isnumeric(x) || ~isreal(x) || issparse(x)
    error('invz:exactDigest', ['x must be a full real numeric array; got %s%s.'], ...
          class(x), repmat(' (sparse)', 1, issparse(x)));
end
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(uint8(class(x)));
md.update(typecast(double(size(x)), 'uint8'));
md.update(typecast(x(:), 'uint8'));      % ORIGINAL class bytes -- never via double(x)
h = typecast(md.digest(), 'uint8');
d = lower(reshape(dec2hex(h, 2).', 1, []));
end
```

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_exact_numeric_digest.m'); disp(table(r))"`
Expected: 7 passed, 0 failed.

- [ ] **Step 3: Record the execution baseline and coupling fingerprint**

Run the full fast suite and record its named/aggregate result in the preregistration. Compute the
production coupling fixture with the explicit tuple in §8; record the requested grid tuple,
assert that `info.grid` is absent on this legacy route, record complete `info.dipole` provenance,
and compute `invz_exact_numeric_digest(Jnu(:))`. Do not invent a field the API intentionally omits.
This is a pre-strict input fingerprint, not a fitted output.

- [ ] **Step 4: Present to the user and get explicit approval**

Do not proceed to Task 1 until the user approves. `crit_tol = 1e-6`, `omit_promote = 0.10` and
`pole_cont_tol = 1e-3` were approved on 2026-07-25 and are labelled FROZEN in §§1, 4, 5 — do not
re-open them. What remains outstanding is the **§9 blind Stage-4 table**: present it as its own
approval, with each entry's INHERITED / DERIVED / CHOICE label visible, and separately from the
inherited scales, the five §7 decisions and the stage-3 grids. Explicitly call out the new choices:
the G8 common paths and 5% spectrum-shape bound; the G14 refinement ladders and `conv_rel`; and
G16's 0.01-T critical grid, 0.005-T adaptive boundary interval, denominator-candidate anchor/window,
and K=5 two-sided fit. Also call out `eta_prod=5e-5 meV`: production overrides the 5e-3 library
default and the FWHM spans only about 2.4 frequency bins — which is why §9.3 adds the frozen
frequency ladder `dw in {0.01, 0.005, 0.0025}` GHz at fixed `eta` and the three-way
`physics_fail` / `resolution_unresolved` / `pass` classification. State plainly that
`resolution_unresolved` blocks promotion but is never evidence against the strict-medium physics.
Finally, flag the Stage-4 obligation `Dresp` creates: a diagnostic-only export with algebraic and
unit tests, and no change to production response arithmetic.

- [ ] **Step 5: Commit**

```bash
git add docs/invzp_strict_medium_prereg.md invz_common/invz_exact_numeric_digest.m \
  invz_projected/tests/test_invz_exact_numeric_digest.m
git commit -m "docs(invzp): freeze strict-static-medium preregistration and scope decisions"
```

---

### Task 1: `invz_coupling_moments` — population central moments of the coupling multiset

**Files:**
- Create: `invz_common/invz_coupling_moments.m`
- Test: `invz_projected/tests/test_invz_coupling_moments.m`

**Interfaces:**
- Consumes: nothing.
- Produces: `mom = invz_coupling_moments(Jnu_flat)` returning `struct('Jbar','mu2','mu3','mu4','n')`. Scalars for a vector input; `1 x nw` row vectors for an `[nJ,nw]` matrix input (one column of moments per frequency). Used by Tasks 3, 7, 8, 10, 12, 14.

- [ ] **Step 1: Write the failing test**

Create `invz_projected/tests/test_invz_coupling_moments.m`:

```matlab
function tests = test_invz_coupling_moments
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% Population (N) normalization, NOT MATLAB's sample (N-1) var(): the difference is 4% at
% N = 24, exactly the size of the legacy synthetic fixtures.
function test_population_normalization_not_sample(testCase)
J = linspace(-2e-3, 6e-3, 24).';
mom = invz_coupling_moments(J);
Jb = mean(J);
verifyEqual(testCase, mom.Jbar, Jb, 'AbsTol', 0);
verifyEqual(testCase, mom.mu2, mean((J - Jb).^2), 'AbsTol', 0);
% and it is demonstrably NOT var()
verifyGreaterThan(testCase, abs(var(J) - mom.mu2)/mom.mu2, 0.03);
end

function test_third_and_fourth_moments(testCase)
J = [1e-3; -2e-3; 5e-3; 0; 3e-3];
Jb = mean(J);
mom = invz_coupling_moments(J);
verifyEqual(testCase, mom.mu3, mean((J - Jb).^3), 'AbsTol', 0);
verifyEqual(testCase, mom.mu4, mean((J - Jb).^4), 'AbsTol', 0);
verifyEqual(testCase, mom.n, numel(J));
end

% [nJ,nw] retardation interface: one moment set PER COLUMN, never one flattened multiset.
function test_matrix_input_is_per_column(testCase)
A = [1e-3; 2e-3; 3e-3];
B = [-4e-3; 0; 4e-3];
mom = invz_coupling_moments([A B]);
momA = invz_coupling_moments(A);
momB = invz_coupling_moments(B);
verifyEqual(testCase, size(mom.Jbar), [1 2]);
verifyEqual(testCase, mom.Jbar, [momA.Jbar momB.Jbar], 'AbsTol', 0);
verifyEqual(testCase, mom.mu2,  [momA.mu2  momB.mu2],  'AbsTol', 0);
verifyEqual(testCase, mom.n,    [momA.n    momB.n]);
% a flattened interpretation would give a single scalar equal to neither
flat = invz_coupling_moments([A; B]);
verifyNotEqual(testCase, flat.mu2, momA.mu2);
end

function test_rejects_nonfinite_and_complex(testCase)
verifyError(testCase, @() invz_coupling_moments([1e-3; NaN]), 'invz:couplingMoments');
verifyError(testCase, @() invz_coupling_moments([1e-3; Inf]), 'invz:couplingMoments');
verifyError(testCase, @() invz_coupling_moments([1e-3; 1i*1e-3]), 'invz:couplingMoments');
verifyError(testCase, @() invz_coupling_moments([]), 'invz:couplingMoments');
end

% Production-multiset regression, with its provenance asserted (spec §B: the numbers are valid
% ONLY for this tuple). INVZ_SLOW-gated like the other real-coupling anchors.
function test_production_multiset_moments(testCase)
if isempty(getenv('INVZ_SLOW')), assumeFail(testCase, 'set INVZ_SLOW=1 to run'); end
ion = invz_ion();
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [16 16 16], 'dpRng', 30, ...
    'dipole', 'bruteforce', 'cache', false));
mom = invz_coupling_moments(Jnu(:));
verifyEqual(testCase, mom.n, 16384);
verifyEqual(testCase, info.dipole.backend, 'bruteforce');
verifyFalse(testCase, isfield(info, 'grid'), ...
    'the frozen fixture uses the legacy absent-grid-policy route');
verifyEqual(testCase, info.Jcc0, 6.42444e-3, 'RelTol', 1e-4);
verifyEqual(testCase, mom.Jbar,  1.20766e-4, 'RelTol', 1e-4);
verifyEqual(testCase, mom.mu2,   5.48264e-6, 'RelTol', 1e-4);
verifyEqual(testCase, sqrt(mom.mu2), 2.3415e-3, 'RelTol', 1e-4);
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_coupling_moments.m'); disp(table(r))"`
Expected: every test errors with "Unrecognized function or variable 'invz_coupling_moments'".

- [ ] **Step 3: Write the implementation**

Create `invz_common/invz_coupling_moments.m`:

```matlab
function mom = invz_coupling_moments(Jnu_flat)
%INVZ_COUPLING_MOMENTS Population central moments of the coupling multiset (strict-order
% static medium, spec SS4.1). G = -chi, ferromagnetic positive J.
%   mom.Jbar = mean_q J,  mom.mu_n = mean((J - Jbar).^n) for n = 2,3,4,  mom.n = count
% NORMALIZATION IS POPULATION (divide by N), matching the BZ average mean_q that
% invz_emt_scalar / invz_emt_static_ordered actually take. MATLAB's default var() is
% sample-normalized (N-1) and is NOT interchangeable: the difference is 6e-5 relative at
% N = 16384 but 4% at the N = 24 synthetic test fixtures, i.e. largest exactly where it would
% go unnoticed.
%
% Jnu_flat: [nJ,1] vector -> scalar fields.
%           [nJ,nw] matrix (T2.1 retardation interface) -> 1 x nw row-vector fields, ONE
%           moment set PER FREQUENCY COLUMN. Static callers use index 1 only; flattening all
%           columns into one static multiset is never a valid interpretation (spec SS4.3).
% The exact q/branch weighting of the input is preserved: every entry counts once.
if ~isnumeric(Jnu_flat) || isempty(Jnu_flat) || ~isreal(Jnu_flat) || ~all(isfinite(Jnu_flat(:)))
    error('invz:couplingMoments', ...
        'Jnu_flat must be a nonempty real finite numeric array; got %s (%d elements).', ...
        class(Jnu_flat), numel(Jnu_flat));
end
if isvector(Jnu_flat), J = Jnu_flat(:); else, J = Jnu_flat; end
Jbar = mean(J, 1);
d    = J - Jbar;                       % implicit expansion over columns
mom  = struct('Jbar', Jbar, ...
              'mu2',  mean(d.^2, 1), ...
              'mu3',  mean(d.^3, 1), ...
              'mu4',  mean(d.^4, 1), ...
              'n',    repmat(size(J, 1), 1, size(J, 2)));
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_coupling_moments.m'); disp(table(r))"`
Expected: 4 passed, 1 incomplete (the `INVZ_SLOW` anchor). Then run it once with the anchor enabled:
`INVZ_SLOW=1 "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_coupling_moments.m'); disp(table(r))"`
Expected: 5 passed, 0 failed.

- [ ] **Step 5: Run the full suite for no regression**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests'); fprintf('P=%d F=%d I=%d T=%d\n', nnz([r.Passed]), nnz([r.Failed]), nnz([r.Incomplete]), numel(r))"`
Expected: `F=0`; compare the named and aggregate results with the Task-0 baseline rather than a
historical hard-coded count.

- [ ] **Step 6: Commit**

```bash
git add invz_common/invz_coupling_moments.m invz_projected/tests/test_invz_coupling_moments.m
git commit -m "feat(invzp): add invz_coupling_moments population-normalized coupling moments"
```

---

### Task 2: `invz_static_medium_reference` — construct `Gref` and own its denominator metadata

The closure leaf cannot diagnose a denominator it never sees (spec §4.1). This primitive owns `Gref` construction *and* the reference-denominator status.

**Files:**
- Create: `invz_common/invz_static_medium_reference.m`
- Test: `invz_projected/tests/test_invz_static_medium_reference.m`

**Interfaces:**
- Consumes: `ref_margin` from Task 0's prereg (default `1e-6`).
- Produces: `[Gref, ref] = invz_static_medium_reference(G0bare, Sigma0, scheme, opts)`.
  `Gref` is scalar (NaN when out of domain).
  `ref = struct('denom','floor','margin','status','scheme')`, where `margin = denom-floor`,
  `ref.status` in `{'ok','ref_denom_nonpositive','ref_denom_small','nonfinite'}`.
  `opts.ref_margin` (default `1e-6`). Used by Tasks 3, 7, 8.

- [ ] **Step 1: Write the failing test**

Create `invz_projected/tests/test_invz_static_medium_reference.m`:

```matlab
function tests = test_invz_static_medium_reference
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_dyson_ref_divides_by_one_plus_sigma(testCase)
[Gref, ref] = invz_static_medium_reference(-0.5, 0.25, 'strict_1z_dyson_ref');
verifyEqual(testCase, ref.denom, 1.25, 'AbsTol', 0);
verifyEqual(testCase, Gref, -0.5/1.25, 'AbsTol', 0);
verifyEqual(testCase, ref.status, 'ok');
verifyEqual(testCase, ref.scheme, 'strict_1z_dyson_ref');
end

function test_bare_ref_denominator_is_one(testCase)
[Gref, ref] = invz_static_medium_reference(-0.5, 0.25, 'strict_1z_bare_ref');
verifyEqual(testCase, ref.denom, 1, 'AbsTol', 0);
verifyEqual(testCase, Gref, -0.5, 'AbsTol', 0);
verifyEqual(testCase, ref.status, 'ok');
end

% The two conventions coincide only as Sigma0 -> 0; that is the O(1/z^2) difference the spec
% labels a Dyson-improved scheme choice, and it must be visible, not hidden.
function test_conventions_differ_at_finite_sigma_and_agree_at_zero(testCase)
gD = invz_static_medium_reference(-0.5, 0.25, 'strict_1z_dyson_ref');
gB = invz_static_medium_reference(-0.5, 0.25, 'strict_1z_bare_ref');
verifyGreaterThan(testCase, abs(gD - gB)/abs(gB), 0.1);
g0D = invz_static_medium_reference(-0.5, 0, 'strict_1z_dyson_ref');
g0B = invz_static_medium_reference(-0.5, 0, 'strict_1z_bare_ref');
verifyEqual(testCase, g0D, g0B, 'AbsTol', 0);
end

% Domain events RETURN a status; they never throw (spec §5.2).
function test_nonpositive_denominator_returns_status(testCase)
[Gref, ref] = invz_static_medium_reference(-0.5, -1.5, 'strict_1z_dyson_ref');
verifyEqual(testCase, ref.status, 'ref_denom_nonpositive');
verifyTrue(testCase, isnan(Gref));
verifyEqual(testCase, ref.denom, -0.5, 'AbsTol', 1e-15);
end

function test_small_denominator_returns_status(testCase)
[Gref, ref] = invz_static_medium_reference(-0.5, -1 + 1e-9, 'strict_1z_dyson_ref');
verifyEqual(testCase, ref.status, 'ref_denom_small');
verifyTrue(testCase, isnan(Gref));
verifyEqual(testCase, ref.floor, 1e-6, 'AbsTol', 0);
verifyLessThanOrEqual(testCase, ref.margin, 0);
end

function test_margin_is_configurable_and_boundary_is_inclusive(testCase)
[~, ref] = invz_static_medium_reference(-0.5, -1 + 1e-3, 'strict_1z_dyson_ref', ...
                                       struct('ref_margin', 1e-2));
verifyEqual(testCase, ref.status, 'ref_denom_small');
[~, ok] = invz_static_medium_reference(-0.5, -1 + 1e-1, 'strict_1z_dyson_ref', ...
                                       struct('ref_margin', 1e-2));
verifyEqual(testCase, ok.status, 'ok');
verifyGreaterThan(testCase, ok.margin, 0);
end

function test_nonfinite_input_returns_status(testCase)
[Gref, ref] = invz_static_medium_reference(NaN, 0.1, 'strict_1z_dyson_ref');
verifyEqual(testCase, ref.status, 'nonfinite');
verifyTrue(testCase, isnan(Gref));
end

% 'resummed' must never reach this primitive: that is a wiring error, not a domain event.
function test_resummed_scheme_is_a_wiring_error(testCase)
verifyError(testCase, @() invz_static_medium_reference(-0.5, 0.1, 'resummed'), ...
            'invz:staticMedium');
verifyError(testCase, @() invz_static_medium_reference(-0.5, 0.1, 'nonsense'), ...
            'invz:staticMedium');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_static_medium_reference.m'); disp(table(r))"`
Expected: all 8 error with "Unrecognized function or variable 'invz_static_medium_reference'".

- [ ] **Step 3: Write the implementation**

Create `invz_common/invz_static_medium_reference.m`:

```matlab
function [Gref, ref] = invz_static_medium_reference(G0bare, Sigma0, scheme, opts)
%INVZ_STATIC_MEDIUM_REFERENCE Reference propagator for the strict-O(1/z) static medium
% (spec SS4.1). G = -chi, ferromagnetic positive J.
%   'strict_1z_dyson_ref' : denom = 1 + Sigma0;  Gref = G0bare/denom   (SELECTED convention)
%   'strict_1z_bare_ref'  : denom = 1;           Gref = G0bare         (systematic comparator)
%
% WHY THIS IS A SEPARATE PRIMITIVE: invz_medium_moment_closure receives only the QUOTIENT
% Gref, so it cannot reconstruct the denominator or report its margin. Reference construction
% therefore cannot live inside the closure.
%
% The Dyson convention is O(1/z)-equivalent to the bare one (Sigma0 is itself O(1/z)), so it
% is a Dyson-IMPROVED scheme choice, not uniquely the literal first-order expansion. It is
% selected because it makes the PM-leg expression exactly G0/D with invz_emt_scalar's own
% D = 1 + Sigma, so both legs' expressions are textually the same object. Gref carries NO K0,
% lambda or xi dependence -- that is what makes the closure one-shot.
%
% Domain events RETURN a status and Gref = NaN; they never throw (spec SS5.2). A caller must
% not evaluate the closure on a non-'ok' reference. An unknown or 'resummed' scheme IS a
% wiring error and throws invz:staticMedium -- the resummed path must bypass this primitive
% entirely.
% opts.ref_margin (default 1e-6; named ref_floor in the preregistration): denom at or below
% this is 'ref_denom_small'. 1 + Sigma0 is O(1) because Sigma0 is O(1/z), so a
% denominator this small means the reference is degenerate.
if nargin < 4 || isempty(opts), opts = struct(); end
margin = getf(opts, 'ref_margin', 1e-6);
if ~(isscalar(margin) && isfinite(margin) && margin > 0)
    error('invz:staticMedium', 'opts.ref_margin must be a positive finite scalar.');
end
switch scheme
    case 'strict_1z_dyson_ref'
        denom = 1 + Sigma0;
    case 'strict_1z_bare_ref'
        denom = 1;
    otherwise
        error('invz:staticMedium', ['invz_static_medium_reference: scheme must be ' ...
            '''strict_1z_dyson_ref'' or ''strict_1z_bare_ref''; got ''%s''. The ' ...
            '''resummed'' path must not call this primitive.'], scheme);
end
ref = struct('denom', denom, 'floor', margin, 'margin', denom - margin, ...
             'status', 'ok', 'scheme', scheme);
if ~isfinite(G0bare) || ~isfinite(denom)
    ref.status = 'nonfinite';
elseif denom <= 0
    ref.status = 'ref_denom_nonpositive';
elseif denom <= margin
    ref.status = 'ref_denom_small';
end
if strcmp(ref.status, 'ok')
    Gref = G0bare/denom;
else
    Gref = NaN;
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_static_medium_reference.m'); disp(table(r))"`
Expected: 8 passed, 0 failed.

- [ ] **Step 5: Run the full suite**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests'); fprintf('P=%d F=%d I=%d T=%d\n', nnz([r.Passed]), nnz([r.Failed]), nnz([r.Incomplete]), numel(r))"`
Expected: `F=0`; no baseline-passing test changes state.

- [ ] **Step 6: Commit**

```bash
git add invz_common/invz_static_medium_reference.m invz_projected/tests/test_invz_static_medium_reference.m
git commit -m "feat(invzp): add invz_static_medium_reference with reference-denominator domain status"
```

---

### Task 3: `invz_medium_moment_closure` — the one-shot `K0 = Jbar - mu2*Gref`

**Files:**
- Create: `invz_common/invz_medium_moment_closure.m`
- Test: `invz_projected/tests/test_invz_medium_moment_closure.m`

**Interfaces:**
- Consumes: `mom` from Task 1 (`struct('Jbar','mu2','mu3','mu4','n')`, scalar fields); `Gref` from Task 2.
- Produces: `[K0, info] = invz_medium_moment_closure(Gref, mom, scheme)`.
  `info = struct('scheme','retained','Kstrict','omit_mu3','omit_cubic','omit_max','status')`,
  `info.retained = 'mu2'`, `info.status` in `{'ok','nonfinite'}`. Used by Tasks 7, 8, 9.

- [ ] **Step 1: Write the failing test**

Create `invz_projected/tests/test_invz_medium_moment_closure.m`:

```matlab
function tests = test_invz_medium_moment_closure
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function mom = fixture_mom()
% production-multiset values (spec §B), used as plain numbers so this stays a unit test
mom = struct('Jbar', 1.20766e-4, 'mu2', 5.48264e-6, 'mu3', -3.42228e-11, ...
             'mu4', 2.3894*5.48264e-6^2, 'n', 16384);
end

function test_one_shot_formula(testCase)
mom = fixture_mom();  Gref = -200;
[K0, info] = invz_medium_moment_closure(Gref, mom, 'strict_1z_dyson_ref');
verifyEqual(testCase, K0, mom.Jbar - mom.mu2*Gref, 'AbsTol', 0);
verifyEqual(testCase, info.Kstrict, K0, 'AbsTol', 0);   % identical for one-shot, by definition
verifyEqual(testCase, info.retained, 'mu2');
verifyEqual(testCase, info.status, 'ok');
verifyEqual(testCase, info.scheme, 'strict_1z_dyson_ref');
end

% Gref < 0 (G = -chi, chi > 0) so the mu2 term RAISES K0 above Jbar. A sign flip here would
% invert the whole medium correction, so pin the direction explicitly.
function test_correction_sign(testCase)
mom = fixture_mom();
K0 = invz_medium_moment_closure(-200, mom, 'strict_1z_dyson_ref');
verifyGreaterThan(testCase, K0, mom.Jbar);
end

% mu3*Gref^2 is the FIRST omitted term, before the cubic. Both ratios are always reported --
% mu3's near-zero value is a measured property of one multiset, never a general licence.
function test_both_omitted_ratios_reported(testCase)
mom = fixture_mom();  Gref = -200;
[~, info] = invz_medium_moment_closure(Gref, mom, 'strict_1z_dyson_ref');
verifyEqual(testCase, info.omit_mu3, abs(mom.mu3*Gref^2)/abs(mom.mu2*Gref), 'RelTol', 1e-12);
verifyEqual(testCase, info.omit_cubic, ...
    abs((2*mom.mu2^2 - mom.mu4)*Gref^3)/abs(mom.mu2*Gref), 'RelTol', 1e-12);
verifyEqual(testCase, info.omit_max, max(info.omit_mu3, info.omit_cubic), 'AbsTol', 0);
end

% A skewed multiset must make omit_mu3 DOMINATE, proving the ratio is not decoration.
function test_skewed_multiset_makes_mu3_dominate(testCase)
mom = fixture_mom();
mom.mu3 = 100*mom.mu2^1.5;                 % strongly skewed, unlike the production multiset
[~, info] = invz_medium_moment_closure(-200, mom, 'strict_1z_dyson_ref');
verifyGreaterThan(testCase, info.omit_mu3, info.omit_cubic);
verifyEqual(testCase, info.omit_max, info.omit_mu3, 'AbsTol', 0);
end

% Frozen thresholds (prereg §4) are the CALLER's gate, not this leaf's: the leaf never rejects.
function test_leaf_never_rejects_on_large_ratio(testCase)
mom = fixture_mom();
[K0, info] = invz_medium_moment_closure(-5000, mom, 'strict_1z_dyson_ref');
verifyGreaterThan(testCase, info.omit_max, 0.25);       % beyond omit_promote
verifyEqual(testCase, info.status, 'ok');               % still 'ok': the polynomial is defined
verifyTrue(testCase, isfinite(K0));
end

% Explicit zero convention (spec §4.1): 0/0 -> 0 only when the numerator vanishes too.
function test_zero_denominator_convention(testCase)
mom = struct('Jbar', 1e-4, 'mu2', 5e-6, 'mu3', 0, 'mu4', 0, 'n', 10);
[~, a] = invz_medium_moment_closure(0, mom, 'strict_1z_dyson_ref');
verifyEqual(testCase, a.omit_mu3, 0, 'AbsTol', 0);      % numerator also 0
verifyEqual(testCase, a.omit_cubic, 0, 'AbsTol', 0);
mom.mu3 = 1e-11;
[~, b] = invz_medium_moment_closure(0, mom, 'strict_1z_dyson_ref');
verifyEqual(testCase, b.omit_mu3, 0, 'AbsTol', 0);      % Gref = 0 => numerator 0 as well
mom2 = struct('Jbar', 1e-4, 'mu2', 0, 'mu3', 1e-11, 'mu4', 1e-14, 'n', 10);
[~, c] = invz_medium_moment_closure(-200, mom2, 'strict_1z_dyson_ref');
verifyEqual(testCase, c.omit_mu3, Inf);                 % mu2 = 0, numerator nonzero
end

function test_nonfinite_gref_returns_status(testCase)
mom = fixture_mom();
[K0, info] = invz_medium_moment_closure(NaN, mom, 'strict_1z_dyson_ref');
verifyEqual(testCase, info.status, 'nonfinite');
verifyTrue(testCase, isnan(K0));
end

function test_resummed_and_unknown_schemes_are_wiring_errors(testCase)
mom = fixture_mom();
verifyError(testCase, @() invz_medium_moment_closure(-200, mom, 'resummed'), ...
            'invz:staticMedium');
verifyError(testCase, @() invz_medium_moment_closure(-200, mom, 'nope'), 'invz:staticMedium');
end

% Non-scalar moment fields are a wiring error: the static slot must pass column 1, not the
% whole [nJ,nw] moment set (spec §4.3).
function test_nonscalar_moments_rejected(testCase)
mom = struct('Jbar', [1e-4 2e-4], 'mu2', [5e-6 6e-6], 'mu3', [0 0], 'mu4', [0 0], 'n', [10 10]);
verifyError(testCase, @() invz_medium_moment_closure(-200, mom, 'strict_1z_dyson_ref'), ...
            'invz:staticMedium');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_medium_moment_closure.m'); disp(table(r))"`
Expected: all 9 error with "Unrecognized function or variable 'invz_medium_moment_closure'".

- [ ] **Step 3: Write the implementation**

Create `invz_common/invz_medium_moment_closure.m`:

```matlab
function [K0, info] = invz_medium_moment_closure(Gref, mom, scheme)
%INVZ_MEDIUM_MOMENT_CLOSURE One-shot strict-O(1/z) static effective medium (spec SS4.1, SSB).
% G = -chi, ferromagnetic positive J.
%
%   K0 = Jbar - mu2*Gref                                          (ONE-SHOT, no feedback)
%
% Both legs call THIS function with an already-constructed Gref (invz_static_medium_reference),
% which is what makes same-order truncation structural rather than a coincidence of two edits:
% the exact closures of invz_emt_scalar and invz_emt_static_ordered are the SAME function of
% (local G; moments) --
%   K = Jbar - mu2*G + mu3*G^2 + (2*mu2^2 - mu4)*G^3 + (mu5 - 5*mu2*mu3)*G^4 + O(G^5)
% -- verified term-by-term to O(G^4). Truncating both at mu2 with the same Gref makes the
% m -> 0 cross-phase identity exact within the scheme.
%
% ORDER ACCOUNTING: under the high-density counting mu2 ~ 1/z, mu2*Gref is the O(1/z) medium
% correction. Solving K0 = Jbar - mu2*G(K0) self-consistently would re-admit denominator
% feedback of the same class that exceeds retained order, so it is NOT done here (the
% self-consistent quadratic is a separately named diagnostic comparator, not a scheme).
%
% info.omit_mu3   = |mu3*Gref^2|            / |mu2*Gref|   <-- the FIRST omitted term
% info.omit_cubic = |(2*mu2^2-mu4)*Gref^3|  / |mu2*Gref|
% info.omit_max   = max of the two.
% Both are ALWAYS reported. mu3 is near zero on the production multiset, but that is a measured
% property of one grid/cutoff/backend -- generalising it is the same inference error that
% produced the synthetic-Jnu defect this work repairs. Zero convention: if mu2*Gref == 0, a
% ratio is 0 when its own numerator is also 0, else Inf.
%
% This leaf NEVER rejects a node on a large ratio: the truncated polynomial stays defined, and
% the frozen omit_report/omit_promote thresholds (docs/invzp_strict_medium_prereg.md SS4) are the
% CALLER's promotion gate. A large ratio must never trigger a scheme switch.
% An unknown or 'resummed' scheme is a wiring error (invz:staticMedium); the resummed path
% bypasses this primitive entirely.
if ~any(strcmp(scheme, {'strict_1z_dyson_ref', 'strict_1z_bare_ref'}))
    error('invz:staticMedium', ['invz_medium_moment_closure: scheme must be ' ...
        '''strict_1z_dyson_ref'' or ''strict_1z_bare_ref''; got ''%s''.'], scheme);
end
req = {'Jbar', 'mu2', 'mu3', 'mu4'};
for k = 1:numel(req)
    if ~isfield(mom, req{k}) || ~isscalar(mom.(req{k}))
        error('invz:staticMedium', ['invz_medium_moment_closure: mom.%s must be a SCALAR. ' ...
            'For a [nJ,nw] retardation moment set, pass the static column (index 1) -- ' ...
            'never the whole set (spec SS4.3).'], req{k});
    end
end

corr = mom.mu2*Gref;                                   % the retained O(1/z) correction
num3 = abs(mom.mu3*Gref^2);
num4 = abs((2*mom.mu2^2 - mom.mu4)*Gref^3);
den  = abs(corr);
info = struct('scheme', scheme, 'retained', 'mu2', 'Kstrict', NaN, ...
              'omit_mu3', ratio(num3, den), 'omit_cubic', ratio(num4, den), ...
              'omit_max', NaN, 'status', 'ok');
info.omit_max = max(info.omit_mu3, info.omit_cubic);

if ~isfinite(Gref) || ~isfinite(corr)
    info.status = 'nonfinite';
    K0 = NaN;
else
    K0 = mom.Jbar - corr;
    if ~isfinite(K0), info.status = 'nonfinite'; end
end
info.Kstrict = K0;      % one-shot: the checked value IS the returned value, by construction
end

% ---------------------------------------------------------------------------------------------
function r = ratio(num, den)
%RATIO Omitted-term ratio with the explicit zero convention (spec SS4.1): 0/0 is 0, x/0 is Inf.
if den == 0
    if num == 0, r = 0; else, r = Inf; end
else
    r = num/den;
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_medium_moment_closure.m'); disp(table(r))"`
Expected: 9 passed, 0 failed.

- [ ] **Step 5: Run the full suite**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests'); fprintf('P=%d F=%d I=%d T=%d\n', nnz([r.Passed]), nnz([r.Failed]), nnz([r.Incomplete]), numel(r))"`
Expected: `F=0`; the new file passes and no baseline-passing test changes state.

- [ ] **Step 6: Commit**

```bash
git add invz_common/invz_medium_moment_closure.m invz_projected/tests/test_invz_medium_moment_closure.m
git commit -m "feat(invzp): add one-shot strict-1z moment closure with both omitted-order ratios"
```

---

### Task 4: `invz_check_static_medium` — sole scheme authority, stamps both legs

**Files:**
- Create: `invz_common/invz_check_static_medium.m`
- Test: `invz_projected/tests/test_invz_check_static_medium.m`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `[sm, eopts, eso] = invz_check_static_medium(opts, eopts, eso)`.
  `sm = struct('scheme','is_strict','ref_margin')`. The returned `eopts`/`eso` carry
  `.static_medium` and `.ref_margin` stamped identically. Used by Tasks 7, 8, 12, 13.

- [ ] **Step 1: Write the failing test**

Create `invz_projected/tests/test_invz_check_static_medium.m`:

```matlab
function tests = test_invz_check_static_medium
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_default_is_resummed_and_not_strict(testCase)
sm = invz_check_static_medium(struct());
verifyEqual(testCase, sm.scheme, 'resummed');
verifyFalse(testCase, sm.is_strict);
end

function test_accepts_the_three_canonical_ids(testCase)
for id = {'resummed', 'strict_1z_dyson_ref', 'strict_1z_bare_ref'}
    sm = invz_check_static_medium(struct('static_medium', id{1}));
    verifyEqual(testCase, sm.scheme, id{1});
end
verifyTrue(testCase, invz_check_static_medium(struct('static_medium', ...
    'strict_1z_dyson_ref')).is_strict);
verifyFalse(testCase, invz_check_static_medium(struct('static_medium', 'resummed')).is_strict);
end

function test_unknown_id_throws(testCase)
verifyError(testCase, @() invz_check_static_medium(struct('static_medium', 'strict')), ...
            'invz:staticMedium');
verifyError(testCase, @() invz_check_static_medium(struct('static_medium', 42)), ...
            'invz:staticMedium');
% the rejected self-consistent comparator is NOT a selectable production scheme (spec §B)
verifyError(testCase, @() invz_check_static_medium(struct('static_medium', ...
            'strict_1z_selfconsistent')), 'invz:staticMedium');
end

% The stamp is what makes the two legs unable to disagree (spec §4.2).
function test_stamps_both_leg_structs_identically(testCase)
[sm, eopts, eso] = invz_check_static_medium( ...
    struct('static_medium', 'strict_1z_dyson_ref'), struct('K0', 7), struct('warn', false));
verifyEqual(testCase, eopts.static_medium, sm.scheme);
verifyEqual(testCase, eso.static_medium,   sm.scheme);
verifyEqual(testCase, eopts.ref_margin,    sm.ref_margin);
verifyEqual(testCase, eso.ref_margin,      sm.ref_margin);
verifyEqual(testCase, eopts.K0, 7);              % pre-existing fields preserved
verifyFalse(testCase, eso.warn);
end

% A per-leg value that DISAGREES with the resolved scheme is a conflict: it would split the two
% sectors across different truncation orders.
function test_disagreeing_per_leg_values_are_conflicts(testCase)
verifyError(testCase, @() invz_check_static_medium( ...
    struct('emt', struct('static_medium', 'strict_1z_dyson_ref'))), 'invz:staticMedium');
verifyError(testCase, @() invz_check_static_medium(struct( ...
    'static_medium', 'strict_1z_dyson_ref', ...
    'emt_static', struct('static_medium', 'resummed'))), 'invz:staticMedium');
end

% IDEMPOTENCY: re-validating an ALREADY-STAMPED opts struct must not be a conflict.
% invz_solve_point_ordered forwards its full numerical context (including opts.emt /
% opts.emt_static) into invz_hmf_ordered, which validates again -- so a stamp that agrees with
% the resolved scheme has to be accepted, or the second pass would throw on its own output.
function test_restamping_an_already_stamped_opts_is_idempotent(testCase)
[sm1, eopts, eso] = invz_check_static_medium(struct('static_medium', 'strict_1z_dyson_ref'));
forwarded = struct('static_medium', sm1.scheme, 'emt', eopts, 'emt_static', eso);
[sm2, eopts2, eso2] = invz_check_static_medium(forwarded, eopts, eso);
verifyEqual(testCase, sm2.scheme, sm1.scheme);
verifyEqual(testCase, eopts2.static_medium, sm1.scheme);
verifyEqual(testCase, eso2.static_medium, sm1.scheme);
verifyEqual(testCase, sm2.ref_margin, sm1.ref_margin, 'AbsTol', 0);
end

function test_ref_margin_default_and_override(testCase)
verifyEqual(testCase, invz_check_static_medium(struct()).ref_margin, 1e-6, 'AbsTol', 0);
sm = invz_check_static_medium(struct('ref_margin', 1e-8));
verifyEqual(testCase, sm.ref_margin, 1e-8, 'AbsTol', 0);
end

% §7.5: strict mode + an [nJ,nw] ordered-retarded coupling matrix is rejected, not silently
% flattened (invz_emt_static_ordered.m:43 would average all frequency columns together).
function test_validator_does_not_guess_coupling_shape(testCase)
% Retarded compatibility is decided only after Jnu_flat has actually been resolved. The PM
% leg supports [nJ,nw]; the ordered solver rejects that matrix at its public entry.
sm = invz_check_static_medium(struct('static_medium', 'strict_1z_dyson_ref'));
verifyTrue(testCase, sm.is_strict);
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_check_static_medium.m'); disp(table(r))"`
Expected: all 8 error with "Unrecognized function or variable 'invz_check_static_medium'".

- [ ] **Step 3: Write the implementation**

Create `invz_common/invz_check_static_medium.m`:

```matlab
function [sm, eopts, eso] = invz_check_static_medium(opts, eopts, eso)
%INVZ_CHECK_STATIC_MEDIUM Sole public authority for the omega_n = 0 static-medium scheme
% (spec SS4.2), following the shared-validator idiom of invz_check_coupling_opts.m /
% invz_check_solve_opts.m / invz_common/invz_check_transverse_mf.m.
%
%   opts.static_medium : 'resummed'             (DEFAULT -- legacy, bit-identical)
%                        'strict_1z_dyson_ref'  (selected strict candidate, spec SS0.3)
%                        'strict_1z_bare_ref'   (systematic comparator)
%   opts.ref_margin    : reference-denominator floor (default 1e-6, FROZEN in
%                        docs/invzp_strict_medium_prereg.md SS2)
%
% Returns sm = struct('scheme','is_strict','ref_margin') and STAMPS the resolved scheme plus
% ref_margin into BOTH internal leg option structs -- eopts (invz_emt_scalar, the PM static
% slot) and eso (invz_emt_static_ordered, the ordered static sector). Stamping from one
% authority is what makes it impossible for the two legs to run different truncation orders,
% which spec SS0.2 requires: a formally O(1/z^2) PM/FM mismatch is dangerous at a continuous
% boundary, where the target mass is exactly zero.
%
% Setting the scheme directly on one leg (opts.emt.static_medium or
% opts.emt_static.static_medium) is a CONFLICT, not an override: it would split the sectors.
% 'strict_1z_selfconsistent' is deliberately NOT selectable -- the rejected self-consistent
% quadratic needs different inputs (G0, D, a branch choice) and lives in a separately named
% diagnostic comparator (spec SSB).
% This validator does not infer coupling shape from option names. The PM leaf supports an
% [nJ,nw] matrix by using column 1 for the strict static slot. The ordered public solver,
% after resolving ODD/retardation, rejects a non-vector Jnu_flat under strict mode rather
% than allowing invz_emt_static_ordered.m:43 to flatten all columns.
if nargin < 1 || isempty(opts), opts = struct(); end
if nargin < 2 || isempty(eopts), eopts = struct(); end
if nargin < 3 || isempty(eso),   eso   = struct(); end

valid = {'resummed', 'strict_1z_dyson_ref', 'strict_1z_bare_ref'};
scheme = getf(opts, 'static_medium', 'resummed');
if ~(ischar(scheme) || isstring(scheme)) || ~any(strcmp(char(scheme), valid))
    error('invz:staticMedium', ['opts.static_medium must be one of %s; got %s. ' ...
        '(''strict_1z_selfconsistent'' is a diagnostic comparator, not a selectable ' ...
        'production scheme -- see the spec SSB.)'], strjoin(valid, ' | '), ...
        local_describe(scheme));
end
scheme = char(scheme);

% A per-leg value that AGREES with the resolved scheme is this function's own stamp being
% re-read, so validation is IDEMPOTENT: invz_solve_point_ordered forwards its full numerical
% context (opts.emt / opts.emt_static included) into invz_hmf_ordered, which validates again.
% A value that DISAGREES is a genuine conflict -- it would split the two sectors across
% different truncation orders, which spec SS0.2 forbids.
for f = {'emt', 'emt_static'}
    if isfield(opts, f{1}) && isstruct(opts.(f{1})) && isfield(opts.(f{1}), 'static_medium')
        inner = opts.(f{1}).static_medium;
        if ~(ischar(inner) || isstring(inner)) || ~strcmp(char(inner), scheme)
            error('invz:staticMedium', ['opts.%s.static_medium = %s conflicts with the ' ...
                'resolved scheme ''%s''. Set the scheme ONCE via opts.static_medium and let ' ...
                'this function stamp both legs (spec SS4.2); a per-leg override would split ' ...
                'the sectors across different truncation orders.'], f{1}, ...
                local_describe(inner), scheme);
        end
    end
end

% Also validate the explicit internal structs passed as arguments; otherwise a disagreeing
% eopts/eso value would simply be overwritten and the conflict hidden.
for pair = {{'eopts', eopts}, {'eso', eso}}
    label = pair{1}{1}; innerStruct = pair{1}{2};
    if isfield(innerStruct, 'static_medium')
        inner = innerStruct.static_medium;
        if ~(ischar(inner) || isstring(inner)) || ~strcmp(char(inner), scheme)
            error('invz:staticMedium', '%s.static_medium conflicts with resolved scheme ''%s''.', ...
                  label, scheme);
        end
    end
end

is_strict = ~strcmp(scheme, 'resummed');
ref_margin = getf(opts, 'ref_margin', 1e-6);
if ~(isscalar(ref_margin) && isfinite(ref_margin) && ref_margin > 0)
    error('invz:staticMedium', 'opts.ref_margin must be a positive finite scalar.');
end
sm = struct('scheme', scheme, 'is_strict', is_strict, 'ref_margin', ref_margin);
eopts.static_medium = sm.scheme;   eopts.ref_margin = sm.ref_margin;
eso.static_medium   = sm.scheme;   eso.ref_margin   = sm.ref_margin;
end

% ---------------------------------------------------------------------------------------------
function s = local_describe(v)
%LOCAL_DESCRIBE Readable rendering of a rejected value (mirrors invz_check_coupling_opts.m).
if ischar(v) || isstring(v)
    s = sprintf('''%s''', char(v));
elseif isnumeric(v) && isscalar(v)
    s = sprintf('%g (%s)', v, class(v));
else
    s = sprintf('a %s', class(v));
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_check_static_medium.m'); disp(table(r))"`
Expected: 8 passed, 0 failed.

- [ ] **Step 5: Run the full suite**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests'); fprintf('P=%d F=%d I=%d T=%d\n', nnz([r.Passed]), nnz([r.Failed]), nnz([r.Incomplete]), numel(r))"`
Expected: `F=0`; no baseline-passing test changes state.

- [ ] **Step 6: Commit**

```bash
git add invz_common/invz_check_static_medium.m invz_projected/tests/test_invz_check_static_medium.m
git commit -m "feat(invzp): add invz_check_static_medium scheme authority stamping both legs"
```

---

### Task 5: shared recoverable-error classifier and outer-call boundary

Today at least four sites absorb *every* `invz:*` id: `invz_ordered_node_solve.m:213-216`,
`invz_ordered_residual.m:178-192` (`safe_eval`), `invz_spectra_map.m:297-300` and
`:325-328`, plus `invz_solve_auto`. Narrowing one only relocates the swallow. This task adds the
shared predicate for **outer dispatcher catches**. Tasks 9–10 remove the inner catches entirely:
after domain events become statuses, the node map and residual checker have no expected recoverable
exception.

**Files:**
- Create: `invz_common/invz_is_recoverable_solver_error.m`
- Create: `invz_common/invz_try_solver_call.m`
- Test: `invz_projected/tests/test_invz_is_recoverable_solver_error.m`

**Interfaces:**
- Produces: `tf = invz_is_recoverable_solver_error(id)`, logical scalar, and
  `[value,completed,error_id] = invz_try_solver_call(fn)`. The latter is the one outer catch
  boundary used by Task 15 and `invz_solve_auto`: it returns an empty value plus the exact id only
  for whitelisted signals and rethrows everything else.

- [ ] **Step 1: Write the failing test**

Create `invz_projected/tests/test_invz_is_recoverable_solver_error.m`:

```matlab
function tests = test_invz_is_recoverable_solver_error
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% The whitelist is exactly the two existing branch/domain signals (spec §5.1).
function test_whitelisted_ids_are_recoverable(testCase)
verifyTrue(testCase, invz_is_recoverable_solver_error('invz:orderedPhase'));
verifyTrue(testCase, invz_is_recoverable_solver_error('invz:degenerateDoublet'));
end

% Wiring errors must NOT be recoverable -- this is the whole point. A masked column caused by a
% swallowed wiring error is the exact failure mode that hid the original defect for a stage.
function test_wiring_ids_are_not_recoverable(testCase)
for id = {'invz:staticMedium', 'invz:emtJnu', 'invz:hzFixed', 'invz:nodeSolveNode', ...
          'invz:residualNode', 'invz:residualState', 'invz:couplingMoments', ...
          'invz:transverseMF', 'invz:hmfOpts', 'invz:oddArgs'}
    verifyFalse(testCase, invz_is_recoverable_solver_error(id{1}), ...
        sprintf('%s must NOT be recoverable', id{1}));
end
end

% Unclassified ids rethrow by DEFAULT -- including unseen invz:* ones. Adding to the whitelist
% is a reviewed contract change, not a response to a failing run.
function test_unknown_invz_id_is_not_recoverable(testCase)
verifyFalse(testCase, invz_is_recoverable_solver_error('invz:somethingBrandNew'));
verifyFalse(testCase, invz_is_recoverable_solver_error('MATLAB:nomem'));
verifyFalse(testCase, invz_is_recoverable_solver_error(''));
end

function test_non_char_input_is_not_recoverable(testCase)
verifyFalse(testCase, invz_is_recoverable_solver_error(42));
verifyFalse(testCase, invz_is_recoverable_solver_error([]));
end

function test_string_class_accepted(testCase)
verifyTrue(testCase, invz_is_recoverable_solver_error("invz:orderedPhase"));
end

% The classifier must be usable directly on an MException identifier.
function test_works_on_a_caught_exception(testCase)
try
    error('invz:degenerateDoublet', 'synthetic');
catch err
    verifyTrue(testCase, invz_is_recoverable_solver_error(err.identifier));
end
try
    error('invz:staticMedium', 'synthetic');
catch err2
    verifyFalse(testCase, invz_is_recoverable_solver_error(err2.identifier));
end
end

function test_try_call_returns_exact_recoverable_category(testCase)
[v, completed, id] = invz_try_solver_call(@() error('invz:degenerateDoublet','synthetic'));
verifyEmpty(testCase, v);
verifyFalse(testCase, completed);
verifyEqual(testCase, id, 'invz:degenerateDoublet');
end

function test_try_call_rethrows_fatal_and_returns_success(testCase)
verifyError(testCase, @() invz_try_solver_call(@() error('invz:staticMedium','synthetic')), ...
    'invz:staticMedium');
[v, completed, id] = invz_try_solver_call(@() 42);
verifyEqual(testCase, v, 42);
verifyTrue(testCase, completed);
verifyEqual(testCase, id, '');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_is_recoverable_solver_error.m'); disp(table(r))"`
Expected: classifier tests fail on the missing classifier and wrapper tests fail on the missing
outer-call helper.

- [ ] **Step 3: Write the implementation**

Create `invz_common/invz_is_recoverable_solver_error.m`:

```matlab
function tf = invz_is_recoverable_solver_error(id)
%INVZ_IS_RECOVERABLE_SOLVER_ERROR Single classifier for every catch on the strict-medium path
% (spec SS5.1). Returns true ONLY for identifiers that are genuine physics/branch/domain
% signals; everything else -- including unseen invz:* identifiers -- is treated as a wiring or
% programming error and must be rethrown.
%
% WHY A WHITELIST, NOT A PREFIX MATCH. invz_ordered_node_solve.m:213-216 currently absorbs
% every 'invz:*' id. Its docstring justifies that by the premise that the only error() site in
% its whole chain (invz_emt_scalar / invz_emt_static_ordered / invz_gstat_ordered /
% invz_lambdas / invz_sigma_ordered / invz_sigma) is invz:emtJnu. The strict scheme adds throw
% sites in that same chain, so a prefix match would silently downgrade a wiring error to "node
% not accepted" -- a masked column. That is exactly the failure mode that let the original
% masking defect hide for a whole stage. And there are at least four such absorbers on the path
% (invz_ordered_residual.safe_eval, invz_spectra_map x2, invz_solve_auto), so narrowing one
% only relocates the swallow: every catch must use THIS predicate.
%
% ADDING AN IDENTIFIER HERE IS A REVIEWED CONTRACT CHANGE, never a convenience response to a
% failing run. Strict-medium domain outcomes deliberately return STATUSES rather than throwing
% (spec SS5.2), so the inner node map and residual checker should have no expected recoverable
% throw at all.
%
% Whitelist rationale:
%   invz:orderedPhase      the strict-PM solver's m = 0 branch gate under a longitudinal tilt
%                          (invz_spectra_map.m:291) -- a branch signal, not a defect.
%   invz:degenerateDoublet the two-level domain floor Delta < 1e-4 meV
%                          (invz_twolevel_ordered.m:19) -- a domain signal, not a defect.
recoverable = {'invz:orderedPhase', 'invz:degenerateDoublet'};
tf = (ischar(id) || isstring(id)) && isscalar(string(id)) && ...
     any(strcmp(char(id), recoverable));
end
```

Create `invz_common/invz_try_solver_call.m`:

```matlab
function [value, completed, error_id] = invz_try_solver_call(fn)
%INVZ_TRY_SOLVER_CALL Outer dispatcher boundary: absorb only reviewed recoverable signals.
% Strict-medium domain outcomes are statuses and do not arrive here.
if ~isa(fn, 'function_handle')
    error('invz:solverCall', 'fn must be a function handle.');
end
try
    value = fn();
    completed = true;
    error_id = '';
catch err
    if ~invz_is_recoverable_solver_error(err.identifier)
        rethrow(err);
    end
    value = [];
    completed = false;
    error_id = err.identifier;       % exact category; never collapse to "nonconverged"
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_is_recoverable_solver_error.m'); disp(table(r))"`
Expected: every classifier/wrapper test passes.

- [ ] **Step 5: Run the full suite**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests'); fprintf('P=%d F=%d I=%d T=%d\n', nnz([r.Passed]), nnz([r.Failed]), nnz([r.Incomplete]), numel(r))"`
Expected: `F=0`; no baseline-passing test changes state. Stage-1 primitives complete; no
production file has been modified yet, so G9 (legacy bitwise) holds trivially.

- [ ] **Step 6: Commit**

```bash
git add invz_common/invz_is_recoverable_solver_error.m invz_common/invz_try_solver_call.m \
  invz_projected/tests/test_invz_is_recoverable_solver_error.m
git commit -m "feat(invzp): centralize whitelisted outer solver-call recovery"
```

---

## Stage 2 — opt-in wiring, default `'resummed'`

Every task below touches a production file. The default path must stay numerically unchanged; each task proves it.

---

### Task 6: `invz_gstat_ordered` — strict-gated exact reassociation (G17 without violating G9)

`1/Gtil0 = 1/Gstat - K0`, so `r -> -G0bare*K0` and `crit -> G0bare*(J0eff - K0)` as `Gstat -> ±Inf`: the divergence **cancels in the quantity Jensen's condition integrates**. But `:30-31` computes `Gtil0 = Gstat/(1 - K0*Gstat)` then `r = G0bare/Gtil0`, which at `Gstat = Inf` evaluates `Inf/(-Inf) = NaN`. This is a reassociation, **not** a regulariser — identical value, no broadening, no tolerance — and it is what makes the spec's removability claim true in floating point.

**Files:**
- Modify: `invz_common/invz_gstat_ordered.m:28-32`
- Test: `invz_projected/tests/test_invz_gstat_removable_pole.m` (new)
- Regression: `invz_projected/tests/test_invz_gstat_ordered.m` (existing, must stay green)

**Compatibility correction.** Applying the reassociation unconditionally changes last bits on the
default resummed path, contradicting G9's bitwise requirement. Extend the signature with an optional
eighth `opts` argument and use `opts.stable_form` (default `false`). Existing/resummed calls retain the
old arithmetic exactly; `invz_emt_static_ordered` passes `stable_form=true` only in strict mode.

**Interfaces:**
- Consumes: optional `opts.stable_form`.
- Produces: `[Gstat, out] = invz_gstat_ordered(..., opts)`; under the strict stable form,
  `out.r` and `out.Gtil0` are finite and continuous through a `Gstat` pole and
  `out.gstat_local_denom = 1+Sigma0+K0*G0inel0` exposes the signed local denominator used by
  G17/Gate 0. Existing seven-argument calls retain their arithmetic bit-identically; appending an
  output-struct field does not alter their numeric values. Relied on by Tasks 8, 10, 12.

- [ ] **Step 1: Measure the reassociation delta before changing anything**

The new arrangement is algebraically identical but not bitwise. Quantify it first, so any later test movement is explained rather than discovered.

Run:
```
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); \
rng_vals=[-1e-3 -1 -50 -475 -1e6]; K0=5e-3; gb=-4e-3; \
for g=rng_vals, a=gb/(g/(1-K0*g)); b=gb*(1/g-K0); \
fprintf('Gstat=%-12.4g as-written=%.17g rearranged=%.17g ulp=%.1f\n', g, a, b, abs(a-b)/eps(abs(a))); end"
```
Expected: `ulp` at most a few (single-digit) for every finite value. **If any value shows more than ~4 ulp, stop and report — that would mean the two forms are not the same expression and the change is not a reassociation.**

- [ ] **Step 2: Write the failing G17 test**

Create `invz_projected/tests/test_invz_gstat_removable_pole.m`:

```matlab
function tests = test_invz_gstat_removable_pole
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% A two-level struct whose hybrid denominator 1 + Sigma0 + K0*G0inel0 can be driven to zero,
% so Gstat sweeps through its own pole. m ~= 0 so the elastic xi term is live (the m = 0 case
% cannot pole this way).
function [tl, beta] = fixture_tl()
beta = 1/(0.0862*0.31);                     % kB*T at T = 0.31 K, meV
tl = struct('Delta', 0.02, 'M2', 0.8, 'm', 0.6, 'n01', tanh(0.02*beta/2), 'g0', 0.5);
tl.g0 = 2*tl.n01/tl.Delta;
end

% G17: r, Gtil0 and crit are finite and continuous through the Gstat crossing, and match the
% analytic limits. Sigma0 is chosen so the hybrid denominator crosses zero.
function test_r_and_crit_finite_through_the_crossing(testCase)
[tl, beta] = fixture_tl();
K0 = 5e-3;  G0inel0 = -300;  G0el0 = -20;  J0eff = 6.42444e-3;
lam = [0.01; 0.02; 0.005];
% denominator 1 + Sigma0 + K0*G0inel0 = 0 at Sigma0 = -1 - K0*G0inel0
S0 = -1 - K0*G0inel0;
eps_list = [-1e-6 -1e-9 0 1e-9 1e-6];
rv = nan(size(eps_list));  cv = nan(size(eps_list));
for k = 1:numel(eps_list)
[Gs, out] = invz_gstat_ordered(tl, lam, K0, S0 + eps_list(k), beta, G0inel0, G0el0, ...
                               struct('stable_form', true));
    verifyTrue(testCase, isfinite(out.r), ...
        sprintf('r must be finite at eps=%g (Gstat=%g)', eps_list(k), Gs));
    verifyTrue(testCase, isfinite(out.Gtil0), sprintf('Gtil0 finite at eps=%g', eps_list(k)));
    rv(k) = out.r;
    cv(k) = out.r + J0eff*out.G0bare;                                  % crit
end
% continuity: the spread across the crossing is small relative to |r| itself
verifyLessThan(testCase, (max(rv) - min(rv))/max(abs(rv)), 1e-3);
verifyLessThan(testCase, (max(cv) - min(cv))/max(abs(cv)), 1e-3);
end

% The analytic limits, hit exactly by driving Gstat to +-Inf. Under the OLD arrangement these
% are NaN, so this test is the direct guard on the reassociation.
function test_analytic_limits_at_infinite_gstat(testCase)
[tl, beta] = fixture_tl();
K0 = 5e-3;  G0inel0 = -300;  G0el0 = 0;   % G0el0 = 0 => Gstat is exactly the hybrid term
lam = [0.01; 0.02; 0.005];
S0inf = -1 - K0*G0inel0;                  % denominator exactly 0 => Gstat = +-Inf
[Gs, out] = invz_gstat_ordered(tl, lam, K0, S0inf, beta, G0inel0, G0el0, ...
                               struct('stable_form', true));
verifyTrue(testCase, isinf(Gs), 'fixture must actually drive Gstat to Inf');
verifyEqual(testCase, out.Gtil0, -1/K0, 'RelTol', 1e-12);
verifyEqual(testCase, out.r, -out.G0bare*K0, 'RelTol', 1e-12);
verifyEqual(testCase, out.r + J0eff*out.G0bare, ...
            out.G0bare*(J0eff-K0), 'RelTol', 1e-12);
verifyTrue(testCase, isfinite(out.r));
end

% Away from any pole the new arrangement must agree with the old expression to float noise.
function test_agrees_with_old_expression_away_from_pole(testCase)
[tl, beta] = fixture_tl();
K0 = 5e-3;  G0inel0 = -300;  G0el0 = -20;  lam = [0.01; 0.02; 0.005];
for S0 = [0, 0.1, -0.3]
    [Gs, out] = invz_gstat_ordered(tl, lam, K0, S0, beta, G0inel0, G0el0);
    old_Gtil0 = Gs/(1 - K0*Gs);
    verifyEqual(testCase, out.Gtil0, old_Gtil0, 'RelTol', 1e-12);
    verifyEqual(testCase, out.r, out.G0bare/old_Gtil0, 'RelTol', 1e-12);
end
end

% The m = 0 pinned identity is untouched by the arrangement (it holds for ANY K0).
function test_m_zero_identity_still_r_equals_one_plus_sigma(testCase)
[tl, beta] = fixture_tl();
tl.m = 0;
lam = [0.01; 0.02; 0.005];
for K0 = [0, 1e-3, 5e-3]
    [~, out] = invz_gstat_ordered(tl, lam, K0, 0.25, beta, -300, 0);
    verifyEqual(testCase, out.r, 1 + 0.25, 'RelTol', 1e-12);
end
end
```

- [ ] **Step 3: Run the new test to verify the limit cases fail**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_gstat_removable_pole.m'); disp(table(r))"`
Expected: `test_analytic_limits_at_infinite_gstat` FAILS (`out.r` is NaN, not `-G0bare*K0`), and `test_r_and_crit_finite_through_the_crossing` fails at `eps = 0`. The other two pass already.

- [ ] **Step 4: Apply the reassociation**

In `invz_common/invz_gstat_ordered.m`, add optional `opts` handling and replace the final arithmetic
with a compatibility branch:

```matlab
gstat_local_denom = 1 + Sigma0 + K0*G0inel0;
Gstat  = G0inel0/gstat_local_denom + xi*G0el0;
G0bare = G0inel0 + G0el0;
Gtil0  = Gstat/(1 - K0*Gstat);
out = struct('xi', xi, 'h0', h0, 'G0bare', G0bare, 'Gtil0', Gtil0, 'r', G0bare/Gtil0);
end
```

with:

```matlab
gstat_local_denom = 1 + Sigma0 + K0*G0inel0;
Gstat  = G0inel0/gstat_local_denom + xi*G0el0;
G0bare = G0inel0 + G0el0;
if nargin >= 8 && getf(opts, 'stable_form', false)
% EXACT REASSOCIATION (not a regulariser -- identical value, no broadening, no tolerance):
%   1/Gtil0 = 1/Gstat - K0   <=>   Gtil0 = Gstat/(1 - K0*Gstat)
% Written this way because Gstat has its own pole where 1 + Sigma0 + K0*G0inel0 -> 0, and the
% divergence CANCELS in the quantities Jensen's condition actually integrates:
%   Gtil0 -> -1/K0,   r -> -G0bare*K0,   crit = r + J0eff*G0bare -> G0bare*(J0eff - K0)
% all finite, with the same limit from both sides. The former arrangement
% (Gtil0 = Gstat/(1-K0*Gstat); r = G0bare/Gtil0) evaluates Inf/(-Inf) = NaN at the crossing and
% loses precision approaching it, which would turn a removable singularity into a node failure.
% Pinned by test_invz_gstat_removable_pole.m (spec G17).
    invGtil0 = 1/Gstat - K0;
    Gtil0    = 1/invGtil0;
    r        = G0bare*invGtil0;
else
    % G9 compatibility: preserve the historical arithmetic on every existing/resummed call.
    Gtil0 = Gstat/(1 - K0*Gstat);
    r     = G0bare/Gtil0;
end
out = struct('xi', xi, 'h0', h0, 'G0bare', G0bare, 'Gtil0', Gtil0, 'r', r, ...
             'gstat_local_denom', gstat_local_denom);
end
```

- [ ] **Step 5: Run the new test to verify it passes**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_gstat_removable_pole.m'); disp(table(r))"`
Expected: 4 passed, 0 failed.

- [ ] **Step 6: Run the existing `gstat` and downstream ordered tests for movement**

Run:
```
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); \
r=runtests({'invz_projected/tests/test_invz_gstat_ordered.m', \
            'invz_projected/tests/test_invz_emt_static_ordered.m', \
            'invz_projected/tests/test_invz_hmf_ordered.m', \
            'invz_projected/tests/test_invz_ordered_residual.m'}); disp(table(r))"
```
Expected: all existing tests pass **bitwise**, because seven-argument calls retain the historical
form. Add a direct G9 test comparing the seven-argument result before/after and a strict-form test
covering the finite limit. If a legacy assertion moves, stop; do not loosen it.

- [ ] **Step 7: Run the full suite**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests'); fprintf('P=%d F=%d I=%d T=%d\n', nnz([r.Passed]), nnz([r.Failed]), nnz([r.Incomplete]), numel(r))"`
Expected: `F=0`; no baseline-passing test changes state.

- [ ] **Step 8: Commit**

```bash
git add invz_common/invz_gstat_ordered.m invz_projected/tests/test_invz_gstat_removable_pole.m
git commit -m "fix(invzp): reassociate gstat r/Gtil0 so the removable Gstat pole stays finite"
```

---

### Task 7: `invz_emt_scalar` — strict `omega_n = 0` slot (PM leg)

Only the `omega = 0` slot changes. `K(2:end)`, `G(2:end)`, both the vector and `[nJ,nw]`
branches, and the `opts.debug` closure diagnostic are untouched. For a matrix input the strict slot
uses `Jnu_flat(:,1)` and moment column 1; rejecting the entire PM-retarded interface would contradict
spec §4.3.

**Files:**
- Modify: `invz_projected/invz_emt_scalar.m` (insert after `:52`, before `med.G = G;`)
- Test: `invz_projected/tests/test_invz_emt_scalar_strict.m` (new)
- Regression: `invz_projected/tests/test_invz_emt_scalar.m` (existing)

**Interfaces:**
- Consumes: `invz_coupling_moments` (Task 1), `invz_static_medium_reference` (Task 2), `invz_medium_moment_closure` (Task 3).
- Produces: `opts.static_medium` (default `'resummed'`), `opts.Jmom` (optional, scalar-field moment
  struct for the static column), `opts.ref_margin` (default `1e-6`). `med` gains `med.medium`
  (`struct('scheme','status','ref','closure')`), `med.medium_status` (char), and
  `med.dynamic_converged = all(isfinite(G(2:end))) && all(isfinite(K(2:end)))`. The existing
  `med.converged` keeps its PM meaning and includes slot 1. The separate dynamic flag is required
  because the ordered node discards/replaces PM slot 1. Used by Tasks 9, 12, 13, 14.

- [ ] **Step 1: Write the failing test**

Create `invz_projected/tests/test_invz_emt_scalar_strict.m`:

```matlab
function tests = test_invz_emt_scalar_strict
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function [G0, Sigma, Jnu] = fixture()
nw = 8;
G0    = -(0.5 ./ (1:nw)').^0 * 0.5;      % [nw x 1], O(1) negative (G = -chi)
G0(1) = -0.9;
Sigma = 0.05*ones(nw, 1);
Jnu   = linspace(-2e-3, 6e-3, 24).';
end

% Absent field => legacy path, BIT-IDENTICAL. This is the G9 guard at leaf level.
function test_absent_scheme_is_bit_identical(testCase)
[G0, Sigma, Jnu] = fixture();
a = invz_emt_scalar(G0, Sigma, Jnu, struct());
b = invz_emt_scalar(G0, Sigma, Jnu, struct('static_medium', 'resummed'));
verifyEqual(testCase, a.K, b.K, 'AbsTol', 0);
verifyEqual(testCase, a.G, b.G, 'AbsTol', 0);
verifyEqual(testCase, a.medium_status, 'not_applicable');
end

% Strict mode changes slot 1 ONLY.
function test_strict_touches_only_slot_one(testCase)
[G0, Sigma, Jnu] = fixture();
leg = invz_emt_scalar(G0, Sigma, Jnu, struct());
st  = invz_emt_scalar(G0, Sigma, Jnu, struct('static_medium', 'strict_1z_dyson_ref'));
verifyEqual(testCase, st.K(2:end), leg.K(2:end), 'AbsTol', 0);
verifyEqual(testCase, st.G(2:end), leg.G(2:end), 'AbsTol', 0);
verifyTrue(testCase, st.dynamic_converged);
verifyNotEqual(testCase, st.K(1), leg.K(1));
verifyEqual(testCase, st.medium_status, 'ok');
verifyEqual(testCase, st.medium.scheme, 'strict_1z_dyson_ref');
end

% The strict slot is exactly the primitive composition -- no re-derivation inside the leaf.
function test_strict_slot_equals_the_primitives(testCase)
[G0, Sigma, Jnu] = fixture();
st = invz_emt_scalar(G0, Sigma, Jnu, struct('static_medium', 'strict_1z_dyson_ref'));
mom  = invz_coupling_moments(Jnu);
Gref = invz_static_medium_reference(G0(1), Sigma(1), 'strict_1z_dyson_ref');
K0   = invz_medium_moment_closure(Gref, mom, 'strict_1z_dyson_ref');
verifyEqual(testCase, st.K(1), K0, 'AbsTol', 0);
verifyEqual(testCase, st.G(1), G0(1)/(1 + Sigma(1) + K0*G0(1)), 'AbsTol', 0);
end

% Supplying Jmom must give the identical answer to deriving it (Task 3 hot-path optimisation
% must not be a semantic change).
function test_supplied_jmom_matches_derived(testCase)
[G0, Sigma, Jnu] = fixture();
o = struct('static_medium', 'strict_1z_dyson_ref');
a = invz_emt_scalar(G0, Sigma, Jnu, o);
o.Jmom = invz_coupling_moments(Jnu);
b = invz_emt_scalar(G0, Sigma, Jnu, o);
verifyEqual(testCase, a.K, b.K, 'AbsTol', 0);
verifyEqual(testCase, a.G, b.G, 'AbsTol', 0);
end

% An out-of-domain reference is a STATUS, not a throw, and it is not reported as
% non-convergence-without-explanation (spec §5.2).
function test_out_of_domain_reference_is_a_status(testCase)
[G0, Sigma, Jnu] = fixture();
Sigma(1) = -1;                                   % 1 + Sigma0 = 0 exactly
med = invz_emt_scalar(G0, Sigma, Jnu, struct('static_medium', 'strict_1z_dyson_ref'));
verifyEqual(testCase, med.medium_status, 'ref_denom_nonpositive');
verifyTrue(testCase, isnan(med.K(1)));
verifyFalse(testCase, med.converged);
verifyEqual(testCase, med.medium.ref.status, 'ref_denom_nonpositive');
end

% The bare-reference comparator is selectable and differs at finite Sigma0.
function test_bare_ref_comparator_differs(testCase)
[G0, Sigma, Jnu] = fixture();
d = invz_emt_scalar(G0, Sigma, Jnu, struct('static_medium', 'strict_1z_dyson_ref'));
b = invz_emt_scalar(G0, Sigma, Jnu, struct('static_medium', 'strict_1z_bare_ref'));
verifyGreaterThan(testCase, abs(d.K(1) - b.K(1))/abs(b.K(1)), 1e-3);
end

% Strict PM + [nJ,nw] remains supported: slot 1 uses column 1 and slots 2:end stay legacy.
function test_strict_matrix_uses_static_column_only(testCase)
[G0, Sigma, Jnu] = fixture();
Jm = repmat(Jnu, 1, numel(G0));
Jm(:,2:end) = Jm(:,2:end) .* (1 + (1:numel(G0)-1));
leg = invz_emt_scalar(G0, Sigma, Jm, struct());
st  = invz_emt_scalar(G0, Sigma, Jm, struct('static_medium', 'strict_1z_dyson_ref'));
mom = invz_coupling_moments(Jm);
mom0 = structfun(@(x) x(1), mom, 'UniformOutput', false);
Gref = invz_static_medium_reference(G0(1), Sigma(1), 'strict_1z_dyson_ref');
K0 = invz_medium_moment_closure(Gref, mom0, 'strict_1z_dyson_ref');
verifyEqual(testCase, st.K(1), K0, 'AbsTol', 0);
verifyEqual(testCase, st.K(2:end), leg.K(2:end), 'AbsTol', 0);
verifyEqual(testCase, st.G(2:end), leg.G(2:end), 'AbsTol', 0);
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_emt_scalar_strict.m'); disp(table(r))"`
Expected: `test_absent_scheme_is_bit_identical` fails on the missing `medium_status` field; the strict tests fail because slot 1 is unchanged.

- [ ] **Step 3: Insert the strict slot**

In `invz_projected/invz_emt_scalar.m`, after line 52 (`G = G0 ./ (D + K .* G0);`) and **before** `med.G = G;  med.K = K;`, insert:

```matlab
% --- strict-O(1/z) static slot: omega_n = 0 ONLY (spec SS4.2, SSB) -------------------------
% The closed-form solve above resums mean_q 1/(D + J*G0) to all orders. That resummation's
% feedback into K exceeds retained order and carries the q-denominator that dies below Bc, so
% under a strict scheme slot 1 is replaced by the one-shot moment closure. K(2:end)/G(2:end)
% and both Jnu branches above are untouched: this is a single-slot substitution, and the
% resulting O(1/z^2) mismatch between K(1) and K(2) is a DOCUMENTED artifact measured by the
% G7 gate, not an assumption.
smid = getf(opts, 'static_medium', 'resummed');
medium = struct('scheme', smid, 'status', 'not_applicable', 'ref', [], 'closure', []);
if ~strcmp(smid, 'resummed')
    if isfield(opts, 'Jmom') && ~isempty(opts.Jmom)
        mom_all = opts.Jmom;                   % hot-path: threaded once per resolved point
    else
        mom_all = invz_coupling_moments(Jf);   % per-column for [nJ,nw]
    end
    mom = local_static_mom(mom_all);            % omega_n=0 is column/index 1, never flatten
    [Gref, refi] = invz_static_medium_reference(G0(1), Sigma(1), smid, ...
        struct('ref_margin', getf(opts, 'ref_margin', 1e-6)));
    [K0s, clo] = invz_medium_moment_closure(Gref, mom, smid);
    medium.ref = refi;  medium.closure = clo;
    if strcmp(refi.status, 'ok') && strcmp(clo.status, 'ok')
        medium.status = 'ok';
        K(1) = K0s;
        G(1) = G0(1)/(D(1) + K0s*G0(1));       % Dyson at the strict medium, same form as :52
    else
        if strcmp(refi.status, 'ok'), medium.status = clo.status;
        else,                         medium.status = refi.status;  end
        K(1) = NaN;  G(1) = NaN;               % domain event: reported, never regularised
    end
end
```

Then, after the existing `med.G = G;  med.K = K;` line, add:

```matlab
med.medium = medium;  med.medium_status = medium.status;
% Ordered callers replace slot 1 with the elastic-hybrid static sector. They must be able to
% validate the dynamic slots without letting the discarded PM slot vote on node acceptance.
med.dynamic_converged = all(isfinite(G(2:end))) && all(isfinite(K(2:end)));
```

Add a private `local_static_mom` helper at the end of the file that copies the first element of
`Jbar`, `mu2`, `mu3`, `mu4`, and `n` into a scalar-field struct, rejecting missing/empty fields with
`invz:staticMedium`. Do not use `structfun` in production: its field-order/output-shape semantics make
the contract less explicit.

- [ ] **Step 4: Run the test to verify it passes**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_emt_scalar_strict.m'); disp(table(r))"`
Expected: 7 passed, 0 failed.

- [ ] **Step 5: Verify the legacy leaf is untouched**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_emt_scalar.m'); disp(table(r))"`
Expected: all pass, unchanged count. The only new behaviour on the legacy path is two additional output fields.

- [ ] **Step 6: Run the full suite**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests'); fprintf('P=%d F=%d I=%d T=%d\n', nnz([r.Passed]), nnz([r.Failed]), nnz([r.Incomplete]), numel(r))"`
Expected: `F=0`; no baseline-passing test changes state.

- [ ] **Step 7: Update the docstring**

Add to `invz_emt_scalar.m`'s header block, after the `Jnu_flat` paragraph:

```matlab
% opts.static_medium ('resummed' default): 'strict_1z_dyson_ref' / 'strict_1z_bare_ref' replace
% the omega_n = 0 slot with the one-shot moment closure K0 = Jbar - mu2*Gref (spec SS4.2). Set it
% ONCE via the public opts.static_medium and let invz_check_static_medium stamp it here -- never
% by hand, or the two legs can diverge in truncation order. opts.Jmom (optional) is the
% invz_coupling_moments struct for the static column; derived here when absent. opts.ref_margin
% (1e-6) is the reference-denominator floor. med.medium/med.medium_status report the reference
% and closure outcome; a domain event leaves K(1)/G(1) = NaN and med.converged false, WITH a
% status that says why. med.dynamic_converged checks slots 2:end only for ordered callers that
% replace slot 1; PM callers continue to use med.converged.
```

- [ ] **Step 8: Commit**

```bash
git add invz_projected/invz_emt_scalar.m invz_projected/tests/test_invz_emt_scalar_strict.m
git commit -m "feat(invzp): add strict-1z omega=0 static slot to invz_emt_scalar (opt-in)"
```

---

### Task 8: `invz_emt_static_ordered` — strict mode removes the inner iteration entirely

Under a strict scheme there is **no fixed point**: `K0` is evaluated once. The damped Picard loop, its `mix`/`maxit`/`tol` knobs and the `invz:emtStatic` non-convergence warning all become inapplicable, which is why node non-convergence ceases to be a possible failure mode in this sector (spec §1).

**Files:**
- Modify: `invz_projected/invz_emt_static_ordered.m:43-66`
- Test: `invz_projected/tests/test_invz_emt_static_ordered_strict.m` (new)
- Regression: `invz_projected/tests/test_invz_emt_static_ordered.m` (existing)

**Interfaces:**
- Consumes: Tasks 1–3 primitives; Task 6's reassociated `out.r`.
- Produces: same signature `[K0, Gstat, out] = invz_emt_static_ordered(tl, lam, Sigma0, Jnu_flat, K0_seed, beta, J0eff, G0inel0, G0el0, opts)`. `opts` gains `.static_medium`, `.Jmom`, `.ref_margin`. `out` gains `.medium` (`struct('scheme','status','ref','closure')`), `.medium_status`, `.omit_mu3`, `.omit_cubic`, `.omit_max`, `.Dq_min`, `.Dq_max`, `.Dq_neg_count`, and `.gstat_local_denom`; under strict, `out.resid = |K0 - Kstrict|` and `out.iters = 0`. Used by Tasks 9, 10, 12.

- [ ] **Step 1: Write the failing test**

Create `invz_projected/tests/test_invz_emt_static_ordered_strict.m`:

```matlab
function tests = test_invz_emt_static_ordered_strict
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function [tl, args] = fixture()
beta = 1/(0.0862*0.31);
tl = struct('Delta', 0.02, 'M2', 0.8, 'm', 0.6, 'n01', tanh(0.02*beta/2), 'g0', 1);
tl.g0  = 2*tl.n01/tl.Delta;
args = struct('lam', [0.01; 0.02], 'Sigma0', 0.05, ...
              'Jnu', linspace(-2e-3, 6e-3, 24).', 'K0_seed', 0, 'beta', beta, ...
              'J0eff', 6.42444e-3, 'G0inel0', -300, 'G0el0', -20);
end

function [K0, Gstat, out] = call(tl, a, o)
[K0, Gstat, out] = invz_emt_static_ordered(tl, a.lam, a.Sigma0, a.Jnu, a.K0_seed, a.beta, ...
                                           a.J0eff, a.G0inel0, a.G0el0, o);
end

% Absent field => legacy iteration, bit-identical.
function test_absent_scheme_is_bit_identical(testCase)
[tl, a] = fixture();
[K1, G1] = call(tl, a, struct('warn', false));
[K2, G2] = call(tl, a, struct('warn', false, 'static_medium', 'resummed'));
verifyEqual(testCase, K1, K2, 'AbsTol', 0);
verifyEqual(testCase, G1, G2, 'AbsTol', 0);
end

% Strict mode does not iterate, and K0 is exactly the primitive composition.
function test_strict_is_one_shot(testCase)
[tl, a] = fixture();
[K0, ~, out] = call(tl, a, struct('warn', false, 'static_medium', 'strict_1z_dyson_ref'));
verifyEqual(testCase, out.iters, 0);
verifyEqual(testCase, out.medium_status, 'ok');
mom  = invz_coupling_moments(a.Jnu);
Gref = invz_static_medium_reference(a.G0inel0 + a.G0el0, a.Sigma0, 'strict_1z_dyson_ref');
Kexp = invz_medium_moment_closure(Gref, mom, 'strict_1z_dyson_ref');
verifyEqual(testCase, K0, Kexp, 'AbsTol', 0);
end

% The strict residual is the algebraic K0 check, identically zero for a correct one-shot call,
% and out.converged reflects DOMAIN validity, not iteration success.
function test_strict_residual_is_the_algebraic_check(testCase)
[tl, a] = fixture();
[~, ~, out] = call(tl, a, struct('warn', false, 'static_medium', 'strict_1z_dyson_ref'));
verifyEqual(testCase, out.resid, 0, 'AbsTol', 0);
verifyTrue(testCase, out.converged);
end

% The K0_seed is IGNORED under strict: a one-shot medium has no warm start to inherit, so a
% contaminated seed can no longer propagate between nodes.
function test_seed_is_ignored_under_strict(testCase)
[tl, a] = fixture();
o = struct('warn', false, 'static_medium', 'strict_1z_dyson_ref');
[Ka] = call(tl, a, o);
a.K0_seed = 0.05;
[Kb] = call(tl, a, o);
verifyEqual(testCase, Ka, Kb, 'AbsTol', 0);
end

% Domain event: a status, no throw, no warning flood.
function test_out_of_domain_reference_is_a_status(testCase)
[tl, a] = fixture();
a.Sigma0 = -1;                                   % 1 + Sigma0 = 0
[K0, ~, out] = call(tl, a, struct('warn', false, 'static_medium', 'strict_1z_dyson_ref'));
verifyEqual(testCase, out.medium_status, 'ref_denom_nonpositive');
verifyTrue(testCase, isnan(K0));
verifyFalse(testCase, out.converged);
end

% Strict mode emits no invz:emtStatic warning even with warn = true: nothing iterates.
function test_no_closure_warning_under_strict(testCase)
[tl, a] = fixture();
w = warning('off', 'all');  restore = onCleanup(@() warning(w));
lastwarn('');
call(tl, a, struct('warn', true, 'static_medium', 'strict_1z_dyson_ref'));
[~, id] = lastwarn();
verifyNotEqual(testCase, id, 'invz:emtStatic');
end

% Both omitted-order ratios are surfaced for the caller's promotion gate.
function test_omitted_ratios_exposed(testCase)
[tl, a] = fixture();
[~, ~, out] = call(tl, a, struct('warn', false, 'static_medium', 'strict_1z_dyson_ref'));
verifyTrue(testCase, isfinite(out.omit_mu3) && isfinite(out.omit_cubic));
verifyEqual(testCase, out.omit_max, max(out.omit_mu3, out.omit_cubic), 'AbsTol', 0);
end

% Dq / D_uni are still built from the physical Gstat and still reported in full (spec §0).
function test_collective_observables_still_reported(testCase)
[tl, a] = fixture();
[K0, Gstat, out] = call(tl, a, struct('warn', false, 'static_medium', 'strict_1z_dyson_ref'));
verifyEqual(testCase, out.D_uni, 1 + (a.J0eff - K0)*Gstat, 'RelTol', 1e-12);
Dq = 1 + (a.Jnu - K0).*Gstat;
verifyEqual(testCase, out.Dq_min, min(Dq), 'AbsTol', 0);
verifyEqual(testCase, out.Dq_max, max(Dq), 'AbsTol', 0);
verifyEqual(testCase, out.Dq_neg_count, nnz(Dq <= 0));
verifyTrue(testCase, isfinite(out.r));                       % Task 6 reassociation in effect
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_emt_static_ordered_strict.m'); disp(table(r))"`
Expected: the strict tests fail (the loop runs and `out.iters >= 1`; `medium_status` missing).

- [ ] **Step 3: Add the strict branch**

In `invz_projected/invz_emt_static_ordered.m`, replace lines 43–56:

```matlab
Jf = Jnu_flat(:);
K0 = K0_seed;
for it = 1:maxit
    Gs = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0);
    Gq = Gs ./ (1 + (Jf - K0).*Gs);
    Gbar = mean(Gq);
    if abs(Gbar - Gs) < rtol, break; end % closed at the CURRENT K0 -- exported as-is
    K0_new = mean(Jf .* Gq) / Gbar;
    dK = abs(K0_new - K0);
    if dK < max(tol, 4*eps(abs(K0)))     % TRUE stall: no representable progress possible
        break;
    end
    K0 = K0 + mix*(K0_new - K0);
end
```

with:

```matlab
Jf = Jnu_flat(:);
smid = getf(opts, 'static_medium', 'resummed');
medium = struct('scheme', smid, 'status', 'not_applicable', 'ref', [], 'closure', []);
if strcmp(smid, 'resummed')
    K0 = K0_seed;
    for it = 1:maxit
        Gs = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0);
        Gq = Gs ./ (1 + (Jf - K0).*Gs);
        Gbar = mean(Gq);
        if abs(Gbar - Gs) < rtol, break; end % closed at the CURRENT K0 -- exported as-is
        K0_new = mean(Jf .* Gq) / Gbar;
        dK = abs(K0_new - K0);
        if dK < max(tol, 4*eps(abs(K0)))     % TRUE stall: no representable progress possible
            break;
        end
        K0 = K0 + mix*(K0_new - K0);
    end
else
    % --- STRICT O(1/z): ONE-SHOT, NO ITERATION (spec SS1, SSB) -----------------------------
    % K0 = Jbar - mu2*Gref with Gref = G0bare/(1+Sigma0), a K0/lambda/xi-INDEPENDENT
    % reference. There is no fixed point here, so mix/maxit/tol and the invz:emtStatic
    % non-convergence warning are all inapplicable, and K0_seed is deliberately IGNORED: a
    % one-shot medium has no warm start, so a contaminated seed cannot propagate between nodes.
    % The resummed q-average this replaces carries the denominator that dies below Bc; its
    % feedback into K0 is what exceeds retained order.
    if ~isvector(Jnu_flat)
        error('invz:staticMedium', ['invz_emt_static_ordered: static_medium ''%s'' does not ' ...
            'support a [nJ,nw] coupling matrix in this phase (spec SS7.5): the pre-existing ' ...
            'Jnu_flat(:) flattening would average every frequency column into the static ' ...
            'q-average.'], smid);
    end
    if isfield(opts, 'Jmom') && ~isempty(opts.Jmom)
        mom = opts.Jmom;                          % threaded once per resolved point / node
    else
        mom = invz_coupling_moments(Jf);          % compatibility fallback
    end
    [Gref, refi] = invz_static_medium_reference(G0inel0 + G0el0, Sigma0, smid, ...
        struct('ref_margin', getf(opts, 'ref_margin', 1e-6)));
    [K0, clo] = invz_medium_moment_closure(Gref, mom, smid);
    medium.ref = refi;  medium.closure = clo;
    if strcmp(refi.status, 'ok') && strcmp(clo.status, 'ok')
        medium.status = 'ok';
    elseif ~strcmp(refi.status, 'ok')
        medium.status = refi.status;
    else
        medium.status = clo.status;
    end
    it = 0;                                       % nothing iterated
end
```

Then replace lines 57–66 (the export tail):

```matlab
[Gstat, go] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0);
Gq = Gstat ./ (1 + (Jf - K0).*Gstat);
out = go;
out.D_uni = 1 + (J0eff - K0)*Gstat;
out.resid = abs(mean(Gq) - Gstat);
out.iters = it;
out.converged = out.resid < rtol;        % measured on the EXPORTED tuple
if warn && ~out.converged
    warning('invz:emtStatic', 'static closure not converged after %d iterations: resid = %.3g', it, out.resid);
end
end
```

with:

```matlab
strict_mode = ~strcmp(medium.status, 'not_applicable');
if strict_mode && ~strcmp(medium.status, 'ok')
    % Do not feed an invalid reference/K0 into Jensen's local denominators. The caller will
    % stop before lambdas/Sigma consume it.
    Gstat = NaN;
    go = struct('xi', NaN, 'h0', NaN, 'G0bare', G0inel0 + G0el0, ...
                'Gtil0', NaN, 'r', NaN, 'gstat_local_denom', NaN);
else
    if strict_mode
        [Gstat, go] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0, ...
                                         struct('stable_form', true));
    else
        [Gstat, go] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0);
    end
end
Gq = Gstat ./ (1 + (Jf - K0).*Gstat);
out = go;
out.D_uni = 1 + (J0eff - K0)*Gstat;      % collective observable, built from the PHYSICAL Gstat
out.Dq_min = min(1 + (Jf - K0).*Gstat);
out.Dq_max = max(1 + (Jf - K0).*Gstat);
out.Dq_neg_count = nnz(1 + (Jf - K0).*Gstat <= 0);
out.iters = it;
out.medium = medium;  out.medium_status = medium.status;
out.omit_mu3 = NaN;  out.omit_cubic = NaN;  out.omit_max = NaN;
if strcmp(medium.status, 'not_applicable')
    out.resid = abs(mean(Gq) - Gstat);   % resummed: the q-average closure residual
    out.converged = out.resid < rtol;    % measured on the EXPORTED tuple
    if warn && ~out.converged
        warning('invz:emtStatic', 'static closure not converged after %d iterations: resid = %.3g', it, out.resid);
    end
else
    % STRICT: the load-bearing residual is the ALGEBRAIC K0 check, not the closure of the
    % discarded resummation (spec SS4.4). Running that discarded inner solve here would restore
    % the very iteration and pole exposure this scheme removes, so it is not computed.
    % out.converged reports DOMAIN validity: there is no iteration to converge.
    out.resid = abs(K0 - medium.closure.Kstrict);
    out.converged = strcmp(medium.status, 'ok');
    out.omit_mu3 = medium.closure.omit_mu3;
    out.omit_cubic = medium.closure.omit_cubic;
    out.omit_max = medium.closure.omit_max;
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_emt_static_ordered_strict.m'); disp(table(r))"`
Expected: 8 passed, 0 failed.

- [ ] **Step 5: Verify the legacy sector is untouched, then the full suite**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_emt_static_ordered.m'); disp(table(r))"`
Expected: all pass.
Run the full suite: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests'); fprintf('P=%d F=%d I=%d T=%d\n', nnz([r.Passed]), nnz([r.Failed]), nnz([r.Incomplete]), numel(r))"`
Expected: `F=0`; no baseline-passing test changes state.

- [ ] **Step 6: Update the docstring**

Add to the header block of `invz_emt_static_ordered.m`, after the `opts` paragraph:

```matlab
% opts.static_medium ('resummed' default): under 'strict_1z_dyson_ref'/'strict_1z_bare_ref' this
% function does NOT iterate. K0 = Jbar - mu2*Gref once; the selected Dyson reference is
% Gref = (G0inel0+G0el0)/(1+Sigma0), while the explicit bare comparator omits that division.
% mix/maxit/tol and the invz:emtStatic warning are inapplicable and K0_seed is IGNORED.
% out.resid becomes the algebraic |K0 - Kstrict| check (zero for a correct call), out.iters = 0,
% and out.converged reports DOMAIN validity via out.medium_status
% ('ok' | 'ref_denom_nonpositive' | 'ref_denom_small' | 'nonfinite'). D_uni and Dq are still
% built from the physical Gstat and reported through Dq_min/Dq_max/Dq_neg_count.
% Strict mode requests invz_gstat_ordered's stable-form reassociation; resummed mode preserves
% the historical seven-argument arithmetic bitwise. opts.Jmom / opts.ref_margin as in
% invz_emt_scalar.
```

- [ ] **Step 7: Commit**

```bash
git add invz_projected/invz_emt_static_ordered.m invz_projected/tests/test_invz_emt_static_ordered_strict.m
git commit -m "feat(invzp): one-shot strict-1z ordered static medium (no inner iteration)"
```

---

### Task 9: `invz_ordered_residual` — block B in place, `res.stability` separate, inner catch removed

Three changes, all in one file because a reviewer cannot accept one without the others: block B's
load-bearing residual becomes the algebraic `K0` check; the stability masses move into their own
field that **never** feeds `res.accepted`; and `safe_eval` is removed. Physics non-convergence and
strict-medium domain outcomes already return values/statuses. Any exception inside the independent
checker is therefore a wiring/programming error and must escape.

**Files:**
- Modify: `invz_projected/invz_ordered_residual.m` (`:58-63` req_node, `:106-126` block B,
  `:162-175` aggregate; delete `:178-192` safe_eval and call all four blocks directly)
- Modify: `docs/invz_ordered_residual_contract.md` (block-B section)
- Test: `invz_projected/tests/test_invz_ordered_residual_strict.m` (new)
- Regression: `invz_projected/tests/test_invz_ordered_residual.m` (existing)

**Interfaces:**
- Consumes: Task 8's `out.medium`.
- Produces: `res.blockB` gains `.status`, `.scheme`, `.ref_denom`, `.omit_mu3`, `.omit_cubic`; `res.blockB.resid_resummed` present only when `opts.debug_resummed`; new `res.stability = struct('crit','D_uni','Dq_min','class','pass')` with `.class` in `{'stable','unstable','boundary_band','undefined'}`. `opts` gains `.K_atol` (1e-14), `.K_rtol` (1e-12), `.crit_tol` (1e-6), `.debug_resummed` (false). Used by Tasks 10, 12.
  Block D uses `med.dynamic_converged`, never `med.converged`, because it compares only
  `K(2:end)` and the ordered map replaces slot 1.

- [ ] **Step 1: Write the failing test**

Create `invz_projected/tests/test_invz_ordered_residual_strict.m`:

```matlab
function tests = test_invz_ordered_residual_strict
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% Build a fixed-h node directly. Do NOT call invz_solve_point_ordered here: that public solver
% does not become scheme-aware until Task 14, and committing a test that knowingly fails until
% then violates the per-task green-suite rule.
function [node, state] = build_strict_fixture()
ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];  hz = 0.02;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
J0eff = 6.42e-3;  Jxx0 = ion.Jxx0;  tmf = 'legacy_x';  Ecut = 40;
[wn, wts, beta] = invz_matsubara(T, Ecut);
si = invz_single_ion(ion, T, Bx, struct('hyp', true, 'hz_fixed', hz, 'Jxx0', Jxx0, ...
                                        'transverse_mf', tmf));
tl = invz_twolevel_ordered(ion, T, Bx, hz, struct('Jxx0', Jxx0, 'transverse_mf', tmf));
c0  = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
G0  = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si, T, 1i*wn(1), struct('elastic', false));
G0inel0 = -real(c0i(3,3,1));
X  = real(c0(:, :, 1));
fb = X(3,1) * (Jxx0 / (1 - Jxx0*X(1,1))) * X(1,3);
G0bare0 = -(X(3,3) + fb);
G0el0   = G0bare0 - G0inel0;
g = real(invz_g(tl, 1i*wn));
mom = invz_coupling_moments(Jnu);
eso = struct('warn', false, 'static_medium', 'strict_1z_dyson_ref', 'Jmom', mom);
eopts = struct('static_medium', 'strict_1z_dyson_ref', 'Jmom', mom);
node = struct('tl', tl, 'G0', G0, 'g', g, 'wts', wts, 'wn', wn, 'beta', beta, ...
              'J0eff', J0eff, 'G0inel0', G0inel0, 'G0el0', G0el0, 'G0bare0', G0bare0, ...
              'eso', eso, 'eopts', eopts, 'Jnu_flat', Jnu, 'Jmom', mom);
[state, info] = invz_ordered_node_solve(node, [], struct('trace', false));
assert(info.accepted, 'invz:testFixture', ...
       'strict residual fixture did not reach an accepted state; choose a deterministic algebra fixture');
end

% Block B is revised IN PLACE: same field name, new load-bearing residual. No fifth block.
function test_blockB_is_the_algebraic_k0_check(testCase)
[node, state] = build_strict_fixture();
res = invz_ordered_residual(node, state);
verifyTrue(testCase, isfield(res, 'blockB'));
verifyFalse(testCase, isfield(res, 'strict'));        % NOT a parallel block
verifyEqual(testCase, res.blockB.scheme, 'strict_1z_dyson_ref');
verifyEqual(testCase, res.blockB.status, 'ok');
verifyLessThan(testCase, res.blockB.resid, res.blockB.scale_abs + ...
    res.blockB.scale_rel*max([abs(state.K0s), 6.0e-3]));
verifyTrue(testCase, res.blockB.pass);
end

% The gate must not be scaled by a vanishing correction or by eps(Jbar) (prereg §2).
function test_gate_has_a_problem_native_floor(testCase)
[node, state] = build_strict_fixture();
res = invz_ordered_residual(node, state);
verifyGreaterThan(testCase, res.blockB.scale_abs + res.blockB.scale_rel*max(abs(node.Jnu_flat)), 1e-15);
end

% Mis-wiring is exactly what this residual exists to catch.
function test_mis_wired_k0_fails_blockB(testCase)
[node, state] = build_strict_fixture();
state.K0s = state.K0s * 1.01;                          % 1% off the strict formula
state.K(1) = state.K0s;
res = invz_ordered_residual(node, state);
verifyFalse(testCase, res.blockB.pass);
verifyFalse(testCase, res.accepted);
end

% The discarded resummed closure is opt-in ONLY, and never feeds .finite/.accepted.
function test_resummed_diagnostic_is_opt_in(testCase)
[node, state] = build_strict_fixture();
res = invz_ordered_residual(node, state);
verifyFalse(testCase, isfield(res.blockB, 'resid_resummed'));
dbg = invz_ordered_residual(node, state, struct('debug_resummed', true));
verifyTrue(testCase, isfield(dbg.blockB, 'resid_resummed'));
verifyEqual(testCase, dbg.accepted, res.accepted);     % diagnostic changes no verdict
verifyEqual(testCase, dbg.finite, res.finite);
end

% THE TWO-TIER SEPARATION (spec §1, G4). Stability is computed but never gates acceptance --
% intermediate path nodes ARE the unstable Landau interval by construction.
function test_stability_is_separate_and_never_gates_acceptance(testCase)
[node, state] = build_strict_fixture();
res = invz_ordered_residual(node, state);
verifyTrue(testCase, isfield(res, 'stability'));
for f = {'crit', 'D_uni', 'Dq_min', 'class', 'pass'}
    verifyTrue(testCase, isfield(res.stability, f{1}), f{1});
end
verifyTrue(testCase, any(strcmp(res.stability.class, ...
    {'stable', 'unstable', 'boundary_band', 'undefined'})));
% forcing the stability verdict negative must NOT change res.accepted
node2 = node;  node2.J0eff = 10*node.J0eff;            % drives crit/D_uni negative
res2 = invz_ordered_residual(node2, state);
verifyEqual(testCase, res2.stability.class, 'unstable');
verifyEqual(testCase, res2.accepted, res.accepted, ...
    'stability must not feed acceptance (spec §1 two-tier verdict)');
end

% crit is the dimensionless mass r + J0eff*G0bare, not an inverse susceptibility.
function test_crit_definition(testCase)
[node, state] = build_strict_fixture();
res = invz_ordered_residual(node, state);
[~, ~, so] = invz_emt_static_ordered(node.tl, state.lam(1:2), state.Sigma(1), node.Jnu_flat, ...
    state.K0s, node.beta, node.J0eff, node.G0inel0, node.G0el0, node.eso);
verifyEqual(testCase, res.stability.crit, so.r + node.J0eff*node.G0bare0, 'RelTol', 1e-12);
end

% Wiring errors must escape the residual checker, not become a failed block (spec §5.1).
function test_wiring_error_is_not_absorbed(testCase)
[node, state] = build_strict_fixture();
node.eso.static_medium = 'not_a_scheme';
verifyError(testCase, @() invz_ordered_residual(node, state), 'invz:staticMedium');
end

% Missing Jmom is a wiring error under strict, and harmless under resummed.
function test_jmom_required_only_under_strict(testCase)
[node, state] = build_strict_fixture();
bad = rmfield(node, 'Jmom');
verifyError(testCase, @() invz_ordered_residual(bad, state), 'invz:residualNode');
leg = bad;  leg.eso.static_medium = 'resummed';  leg.eopts.static_medium = 'resummed';
r = invz_ordered_residual(leg, state);
verifyTrue(testCase, isstruct(r));                     % no throw on the legacy path
end

% A strict reference-domain event returns a complete nonaccepted residual WITHOUT evaluating
% A/C/D on an invalid medium.
function test_domain_preflight_short_circuits_checker(testCase)
[node, state] = build_strict_fixture();
state.Sigma(1) = -1;                                  % 1+Sigma0 = 0
res = invz_ordered_residual(node, state);
verifyEqual(testCase, res.blockB.status, 'ref_denom_nonpositive');
verifyFalse(testCase, res.accepted);
verifyFalse(testCase, res.finite);
verifyEqual(testCase, res.stability.class, 'undefined');
verifyTrue(testCase, all(isnan([res.blockA.resid,res.blockC.resid,res.blockD.resid])));
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_ordered_residual_strict.m'); disp(table(r))"`
Expected before the edit: failures because `res.stability` is absent and block B still demands the
resummed closure. The direct fixed-h fixture deliberately avoids any dependency on the later public
point-solver wiring.

- [ ] **Step 3: Add the conditional `Jmom` requirement**

In `invz_ordered_residual.m`, after the existing `req_node` loop (`:58-63`), insert:

```matlab
% Jmom is required only once a STRICT scheme is selected (spec SS4.3): legacy/resummed
% direct-node fixtures may omit it without changing their numerical path.
smid_node = getf(node.eso, 'static_medium', 'resummed');
if ~strcmp(smid_node, 'resummed') && (~isfield(node, 'Jmom') || isempty(node.Jmom))
    error('invz:residualNode', ['node.Jmom is required under static_medium ''%s'' ' ...
        '(invz_coupling_moments of the static coupling column).'], smid_node);
end
```

- [ ] **Step 4: Replace block B with the scheme-aware version**

Before evaluating any of Blocks A/C/D, add a strict-only independent reference preflight using
`invz_static_medium_reference(node.G0bare0,Sigma(1),smid_node,eso)` and
`invz_medium_moment_closure(...,node.Jmom,smid_node)`. If either status is non-`ok`, return a complete
residual schema immediately: Block B carries the exact reference/closure status and denominator;
A/C/D have `resid=NaN`, `pass=false`, `err=''`; stability is `undefined`; `finite=false`,
`accepted=false`. Do not call `local_F`, lambdas, Sigma, or either EMT leaf after that preflight.
This is an independent recomputation from the exported state, not trust in the live loop's `sout`.

Replace `invz_ordered_residual.m:106-126` (the whole Block B section) with:

```matlab
% =========================================================================================
% Block B -- static medium (contract Sec. 4, REVISED IN PLACE for the strict scheme).
%   'resummed'          : unchanged -- the q-average closure residual from a fresh
%                         invz_emt_static_ordered at the exported Sigma(1)/lam(1:2).
%   strict_1z_*         : the load-bearing residual becomes the ALGEBRAIC check
%                         |K0s - Kstrict(Gref)|, independently recomputed from the exported
%                         state. It is identically zero for a correctly wired one-shot call,
%                         so it exists to catch MIS-WIRING, which is exactly the class of
%                         defect a prefix-matching error catch used to hide.
% The discarded resummed closure is NOT run under a strict scheme unless
% opts.debug_resummed: doing so would restore the inner iteration and pole exposure this
% design removes, and the analytic-continuation path deliberately crosses that pole.
% =========================================================================================
outB = local_blockB(tl, lam, Sigma, Jnu_flat, K0s, beta, J0eff, ...
                    G0inel0, G0el0, eso);
strictB = ~strcmp(smid_node, 'resummed');
if strictB
    scaleB_abs = getf(opts, 'K_atol', 1e-14);
    scaleB_rel = getf(opts, 'K_rtol', 1e-12);
else
    rtolB = getf(eso, 'resid_tol', 1e-10);
    scaleB_abs = rtolB;  scaleB_rel = rtolB;
end
statusB = 'nonfinite';  omit3 = NaN;  omit4 = NaN;  refdenB = NaN;
Gstat_b = outB.Gstat;
D_uni   = outB.so.D_uni;
Dq      = 1 + (Jnu_flat(:) - outB.K0) .* Gstat_b;
Dq_min  = min(Dq);  Dq_max = max(Dq);
convB   = outB.so.converged;
if strictB
    statusB = outB.so.medium_status;
    omit3 = outB.so.omit_mu3;  omit4 = outB.so.omit_cubic;
    refdenB = outB.so.medium.ref.denom;
    if strcmp(statusB, 'ok')
        Kstrict = outB.so.medium.closure.Kstrict;
        rB = abs(K0s - Kstrict);
        gate = scaleB_abs + scaleB_rel*max([abs(K0s), abs(Kstrict), Jscale]);
        passB = isfinite(rB) && (rB < gate);
    else
        rB = NaN;  passB = false;  convB = false;    % domain status, not an exception
    end
else
    statusB = 'not_applicable';
    rB = outB.so.resid;
    passB = isfinite(rB) && (rB < scaleB_abs + scaleB_rel*abs(Gstat_b));
end
res.blockB = struct('resid', rB, 'scale_abs', scaleB_abs, 'scale_rel', scaleB_rel, ...
                     'pass', passB, 'converged', convB, 'err', '', 'status', statusB, ...
                     'scheme', smid_node, 'ref_denom', refdenB, ...
                     'omit_mu3', omit3, 'omit_cubic', omit4);
if strictB && getf(opts, 'debug_resummed', false)
    esoR = eso;  esoR.static_medium = 'resummed';  esoR.warn = false;
    outR = local_blockB(tl, lam, Sigma, Jnu_flat, K0s, beta, J0eff, ...
                        G0inel0, G0el0, esoR);
    res.blockB.resid_resummed = outR.so.resid;
end
```

- [ ] **Step 5: Add `res.stability` and adopt the classifier**

Replace the aggregate block (`:162-175`) so `res.stability` is added and acceptance is untouched by it:

```matlab
% ---- top-level diagnostics / finite / stall / aggregate (contract Sec. 6-7) -------------
res.D_uni  = D_uni;
res.Dq_min = Dq_min;
res.Dq_max = Dq_max;
% ---- STABILITY TIER (spec SS1): computed for every node, folded into .accepted for NONE.
% Intermediate path nodes may legitimately be unstable -- they are the analytic continuation
% through the unstable Landau interval, and requiring per-node positivity would re-mask the
% ordered phase. Only an ENDPOINT root is held to this tier, by the caller.
crit_tol = getf(opts, 'crit_tol', 1e-6);
if all(isfinite([D_uni, Dq_min]))
    critv = outB.so.r + J0eff*node.G0bare0;
    D_tol = 1e-6*max(1, abs(Gstat_b)*Jscale);        % prereg SS1: noise scales with |Gstat|
    Dq_tol = D_tol;
    if ~isfinite(critv)
        cls = 'undefined';
    elseif critv > crit_tol && D_uni > D_tol && Dq_min > D_tol
        cls = 'stable';
    elseif critv < -crit_tol || D_uni < -D_tol || Dq_min < -Dq_tol
        cls = 'unstable';
    else
        cls = 'boundary_band';
    end
    res.stability = struct('crit', critv, 'D_uni', D_uni, 'Dq_min', Dq_min, ...
                           'class', cls, 'pass', strcmp(cls, 'stable'));
else
    res.stability = struct('crit', NaN, 'D_uni', NaN, 'Dq_min', NaN, ...
                           'class', 'undefined', 'pass', false);
end
res.finite = all(isfinite(Sigma)) && all(isfinite(K)) && all(isfinite(lam)) && isfinite(K0s) ...
             && isfinite(rA) && isfinite(rB) && isfinite(rC) && isfinite(rD);
if isnan(dS_in)
    res.stall = NaN;
else
    res.stall = isfinite(dS_in) && (dS_in < tol_outer) && ...
                ~(passA && passB && passC && passD && res.finite);
end
res.accepted = res.finite && passA && passB && passC && passD;
end
```

Remove `safe_eval` and evaluate Blocks A–D directly. Preserve each block's `.err` field as `''` for
schema compatibility, but do not catch. Replace every `[out, err] = safe_eval(@() fn(...))` by
`out = fn(...); err = '';`. The opt-in resummed diagnostic follows the same rule. This is deliberate:
`invz:degenerateDoublet` must already have been converted to a return-mode domain record before a
node reaches this checker, and an `invz:orderedPhase` signal cannot arise inside it. Absorbing either
here would hide a wiring defect.

```matlab
% No safe_eval local function remains. Exceptions are fatal at this layer.
```

In Block D replace both uses of `medD.converged` with `medD.dynamic_converged`. This is not a
tolerance relaxation: Block D already excludes `K(1)` by contract, and the ordered loop overwrites
that slot before lambdas. Letting the discarded PM static slot make Block D fail would reintroduce a
static-sector veto through a nominally dynamic residual.

- [ ] **Step 6: Run the existing residual test for regression, then the full suite**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_ordered_residual.m'); disp(table(r))"`
Expected: all pass. `res.blockB.status` is `'not_applicable'` on the legacy path and every existing assertion is unchanged.
Then the full suite: `F=0`; the new strict residual file must pass now. Do not commit a
failing-by-design test.

- [ ] **Step 7: Update the binding contract document**

In `docs/invz_ordered_residual_contract.md`, add a dated subsection to the Block-B section:

```markdown
### Block B under a strict static medium (added 2026-07-25)

Block B is revised IN PLACE; there is no fifth block. Under
`static_medium` in {`strict_1z_dyson_ref`, `strict_1z_bare_ref`}:

- load-bearing residual: `|K0s - Kstrict(Gref)|`, with `Kstrict` independently recomputed from
  the exported `Sigma(1)`/`lam(1:2)` via `invz_static_medium_reference` +
  `invz_medium_moment_closure`. Gate `K_atol + K_rtol*max(|K0s|,|Kstrict|,Jscale)`
  (`1e-14`/`1e-12`, frozen in docs/invzp_strict_medium_prereg.md §2). Identically zero for a
  correct one-shot call: it exists to catch mis-wiring.
- `blockB.pass` additionally requires `blockB.status == 'ok'`, the reference-denominator domain
  verdict. An out-of-domain node is not accepted, but the CALLER must distinguish that from
  non-convergence (`prof.status` `'medium_out_of_domain'` vs `'node_failed'`).
- a strict reference-domain failure is preflighted independently from the exported state and
  returns a complete nonaccepted schema without evaluating A/C/D on an invalid medium.
- the resummed q-average closure residual is NOT computed unless `opts.debug_resummed`, and is
  reported as the nullable diagnostic `blockB.resid_resummed`. It never enters `.finite` or
  `.accepted`: the analytic-continuation path deliberately crosses its pole.
- `res.stability` (`crit`, `D_uni`, `Dq_min`, `class`, `pass`) is computed for every node and
  folded into `res.accepted` for NONE. Requiring per-node stability would re-mask the ordered
  phase, since intermediate nodes are the unstable Landau interval by construction.
- the checker contains no exception absorber. All exceptions escape; domain and numerical
  non-convergence are represented by statuses/residuals before this layer.
- Block D checks `med.dynamic_converged` (slots `2:end`), not whole-PM `med.converged`, because
  ordered callers replace the discarded PM static slot before lambdas.
```

- [ ] **Step 8: Commit**

```bash
git add invz_projected/invz_ordered_residual.m docs/invz_ordered_residual_contract.md invz_projected/tests/test_invz_ordered_residual_strict.m
git commit -m "feat(invzp): revise residual block B in place for strict medium, split stability tier"
```

---

### Task 10: `invz_ordered_node_solve` — thread `Jmom`, stop on a domain event, remove the inner catch

**Files:**
- Modify: `invz_projected/invz_ordered_node_solve.m` (`:119-125` req_node, `:163-180` locals,
  `:181-221` loop and catch removal, `:228-239` term_reason/info)
- Test: `invz_projected/tests/test_invz_ordered_node_solve_strict.m` (new)
- Regression: `invz_projected/tests/test_invz_ordered_node_solve.m` (existing)

**Interfaces:**
- Consumes: Task 8's `out.medium_status`.
- Produces: `node` accepts a 14th field `Jmom`; `info` gains `.medium` (the winning attempt's `sout.medium`) and `.medium_status`; `info.term_reason` gains `'medium_out_of_domain'`. Used by Task 12.

- [ ] **Step 1: Write the failing test**

Create `invz_projected/tests/test_invz_ordered_node_solve_strict.m`:

```matlab
function tests = test_invz_ordered_node_solve_strict
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% Minimal real node built from public calls at a field where the bare set orders.
function node = build_node(scheme)
ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];  hz = 0.02;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
J0eff = 6.42e-3;  Jxx0 = ion.Jxx0;  tmf = 'legacy_x';
[wn, wts, beta] = invz_matsubara(T, 40);
si  = invz_single_ion(ion, T, Bx, struct('hyp', true, 'hz_fixed', hz, 'Jxx0', Jxx0, ...
                                         'transverse_mf', tmf));
tl  = invz_twolevel_ordered(ion, T, Bx, hz, struct('Jxx0', Jxx0, 'transverse_mf', tmf));
c0  = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
G0  = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si, T, 1i*wn(1), struct('elastic', false));
G0inel0 = -real(c0i(3,3,1));
X  = real(c0(:, :, 1));
fb = X(3,1) * (Jxx0 / (1 - Jxx0*X(1,1))) * X(1,3);
G0bare0 = -(X(3,3) + fb);
node = struct('tl', tl, 'G0', G0, 'g', real(invz_g(tl, 1i*wn)), 'wts', wts, 'wn', wn, ...
    'beta', beta, 'J0eff', J0eff, 'G0inel0', G0inel0, 'G0el0', G0bare0 - G0inel0, ...
    'G0bare0', G0bare0, 'eso', struct('warn', false, 'static_medium', scheme), ...
    'eopts', struct('static_medium', scheme), 'Jnu_flat', Jnu, ...
    'Jmom', invz_coupling_moments(Jnu));
end

function test_strict_node_solves_and_reports_medium(testCase)
node = build_node('strict_1z_dyson_ref');
[state, info] = invz_ordered_node_solve(node, [], struct('trace', false));
verifyEqual(testCase, info.medium_status, 'ok');
verifyEqual(testCase, info.medium.scheme, 'strict_1z_dyson_ref');
verifyTrue(testCase, isfinite(state.K0s));
verifyTrue(testCase, isfield(info.res, 'stability'));
end

% Jmom must actually be threaded into BOTH leaves, not silently re-derived.
function test_jmom_is_threaded_to_both_leaves(testCase)
node = build_node('strict_1z_dyson_ref');
bad = node;
bad.Jmom.mu2 = node.Jmom.mu2 * 2;          % a deliberately wrong moment
[sa] = invz_ordered_node_solve(node, [], struct('trace', false));
[sb] = invz_ordered_node_solve(bad,  [], struct('trace', false));
verifyNotEqual(testCase, sa.K0s, sb.K0s, ...
    'node.Jmom must reach the static leaf; identical K0s means it was re-derived');
end

% Missing Jmom under strict is a wiring error; under resummed it is harmless.
function test_missing_jmom_under_strict_throws(testCase)
node = rmfield(build_node('strict_1z_dyson_ref'), 'Jmom');
verifyError(testCase, @() invz_ordered_node_solve(node, [], struct()), 'invz:nodeSolveNode');
leg = rmfield(build_node('resummed'), 'Jmom');
[~, info] = invz_ordered_node_solve(leg, [], struct());
verifyEqual(testCase, info.medium_status, 'not_applicable');
end

% A domain event stops the attempt BEFORE lambdas/Sigma consume an invalid reference, and is
% reported as its own term_reason -- never as generic max_iter.
function test_domain_event_stops_before_lambdas(testCase)
node = build_node('strict_1z_dyson_ref');
node.eso.ref_margin = 1e9;                 % force every reference denominator out of domain
[state, info] = invz_ordered_node_solve(node, [], struct('trace', false));
verifyEqual(testCase, info.term_reason, 'medium_out_of_domain');
verifyEqual(testCase, info.medium_status, 'ref_denom_small');
verifyFalse(testCase, info.accepted);
verifyEqual(testCase, state.lam, [0; 0; 0], 'AbsTol', 0);   % lambdas never ran
end

% Wiring errors escape; they are not scored as a failed node.
function test_wiring_error_escapes(testCase)
node = build_node('strict_1z_dyson_ref');
node.eso.static_medium = 'not_a_scheme';
verifyError(testCase, @() invz_ordered_node_solve(node, [], struct()), 'invz:staticMedium');
end

% Legacy path unchanged.
function test_resummed_path_untouched(testCase)
node = build_node('resummed');
[~, info] = invz_ordered_node_solve(node, [], struct('trace', false));
verifyTrue(testCase, any(strcmp(info.term_reason, ...
    {'accepted', 'loop_converged_not_accepted', 'max_iter'})));
verifyEqual(testCase, info.medium_status, 'not_applicable');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_ordered_node_solve_strict.m'); disp(table(r))"`
Expected: fails on the missing `info.medium_status`, `info.medium`, and the `'medium_out_of_domain'` term reason.

- [ ] **Step 3: Add the conditional `Jmom` requirement**

In `invz_ordered_node_solve.m`, after the `req_node` loop (`:121-125`), insert:

```matlab
% Jmom is required only once a STRICT scheme is selected (spec SS4.3).
smid_node = getf(node.eso, 'static_medium', 'resummed');
if ~strcmp(smid_node, 'resummed') && (~isfield(node, 'Jmom') || isempty(node.Jmom))
    error('invz:nodeSolveNode', ['node.Jmom is required under static_medium ''%s'' ' ...
        '(invz_coupling_moments of the static coupling column).'], smid_node);
end
```

- [ ] **Step 4: Thread `Jmom` and stop on a domain event inside `run_attempt`**

In `run_attempt`, after the existing local copies (`:173`, `eopts = node.eopts;`), insert:

```matlab
if isfield(node, 'Jmom') && ~isempty(node.Jmom)
    eopts.Jmom = node.Jmom;                 % strict slot reads it; legacy path ignores it
    eso_local  = node.eso;  eso_local.Jmom = node.Jmom;
else
    eso_local  = node.eso;
end
medium_status = 'not_applicable';  medium = struct('scheme', getf(node.eso, ...
    'static_medium', 'resummed'), 'status', 'not_applicable', 'ref', [], 'closure', []);
```

Extend the pre-loop `med` placeholder with `dynamic_converged=false`. The ordered map may retain a
non-finite PM slot 1 until it is replaced by `K0s/Gstat`; neither the live loop nor the residual
checker may use whole-PM `med.converged` as an ordered acceptance condition. Slots `2:end` remain
guarded by `med.dynamic_converged` in Block D.

Change the in-loop static call (`:188-189`) to pass `eso_local` instead of `node.eso`, and insert the domain-event break immediately after it, **before** `K(1) = K0s;`:

```matlab
        [K0s, Gstat_it, sout] = invz_emt_static_ordered(node.tl, lam(1:2), Sigma(1), ...
            Jnu_flat, K0s, node.beta, node.J0eff, node.G0inel0, node.G0el0, eso_local);
        % Domain event: stop BEFORE invz_lambdas / invz_sigma_ordered consume an invalid
        % reference (spec SS4.4). Propagating a non-finite K0 through the outer map would turn a
        % reportable domain outcome into an unexplained non-convergence.
        medium_status = getf(sout, 'medium_status', 'not_applicable');
        medium        = getf(sout, 'medium', medium);
        if ~any(strcmp(medium_status, {'not_applicable', 'ok'}))
            outer_used = outer;
            break;
        end
```

Do **not** run the post-loop refresh after a domain event. That would feed the invalid reference into
the static function a second time and contradict “stop before consumption.” Use:

```matlab
    if any(strcmp(medium_status, {'not_applicable', 'ok'}))
        [K0s, Gstat, so] = invz_emt_static_ordered(node.tl, lam(1:2), Sigma(1), ...
            Jnu_flat, K0s, node.beta, node.J0eff, node.G0inel0, node.G0el0, eso_local);
        medium_status = getf(so, 'medium_status', medium_status);
        medium        = getf(so, 'medium', medium);
        if any(strcmp(medium_status, {'not_applicable', 'ok'})), K(1) = K0s; end
    else
        so = sout;  Gstat = Gstat_it;       % preserve the classified domain record
    end
```

- [ ] **Step 5: Remove the catch and add the new term reason**

Delete the `try/catch` surrounding the node map. No exception is recoverable here: the two-level
domain check occurs before entry, and medium-domain events return statuses. Any thrown
`invz:emtJnu`, `invz:staticMedium`, or future identifier is a wiring/programming failure and must
escape.

```matlab
% No catch remains in run_attempt.
```

Replace the `term_reason` selection (`:228-234`) with:

```matlab
if res.accepted
    term_reason = 'accepted';
elseif ~any(strcmp(medium_status, {'not_applicable', 'ok'}))
    term_reason = 'medium_out_of_domain';   % distinct from a convergence failure
elseif loop_converged
    term_reason = 'loop_converged_not_accepted';
else
    term_reason = 'max_iter';
end
```

And extend the `info` struct (`:238-239`):

```matlab
info = struct('res', res, 'loop_converged', loop_converged, 'so', so_out, 'med', med, ...
              'outer_iters', outer_used, 'term_reason', term_reason, 'iters', iters, ...
              'medium', medium, 'medium_status', medium_status);
```

- [ ] **Step 6: Run the new and existing node-solve tests**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests({'invz_projected/tests/test_invz_ordered_node_solve_strict.m','invz_projected/tests/test_invz_ordered_node_solve.m'}); disp(table(r))"`
Expected: the new file 6 passed; the existing file all pass unchanged.

- [ ] **Step 7: Run the full suite**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests'); fprintf('P=%d F=%d I=%d T=%d\n', nnz([r.Passed]), nnz([r.Failed]), nnz([r.Incomplete]), numel(r))"`
Expected: `F=0`; no baseline-passing test changes state.

- [ ] **Step 8: Commit**

```bash
git add invz_projected/invz_ordered_node_solve.m invz_projected/tests/test_invz_ordered_node_solve_strict.m
git commit -m "feat(invzp): thread Jmom, halt on medium domain events, remove node-solve catch"
```

---

### Task 11: `invz_twolevel_ordered` — single-evaluation `domain_policy` for the `Delta` floor

`Delta` is evaluated at the node's own molecular field, and the geometric `h`-grid clusters at 0, so the predictor and lowest nodes are at risk whenever `Bx` is small. Today the throw escapes `invz_hmf_ordered` entirely (`:281` is outside any try/catch) and masks the column indistinguishably from a solver failure. The fix is a return mode with **one** diagonalization — not a pre-screen that duplicates it and then calls the constructor again.

**Files:**
- Modify: `invz_common/invz_twolevel_ordered.m:16-21`
- Test: `invz_projected/tests/test_invz_twolevel_ordered_domain.m` (new)

**Interfaces:**
- Consumes: nothing new.
- Produces: `opts.domain_policy` in `{'throw','return'}` (default `'throw'`, the current behaviour). `tl.valid` (logical) is now always present; `tl.Delta` is always set. In return mode an invalid `tl` carries **only** `valid` and `Delta` — the two-level fields are deliberately absent, so any consumer that ignores `valid` fails loudly. Used by Task 12.

- [ ] **Step 1: Read the file to locate its closing structure**

Run: `sed -n '1,40p' invz_common/invz_twolevel_ordered.m`
Confirm the guard is at `:18-21` and note which line assigns the last `tl.*` field. The edit below only touches `:16-21`, so nothing downstream needs relocating.

- [ ] **Step 2: Write the failing test**

Create `invz_projected/tests/test_invz_twolevel_ordered_domain.m`:

```matlab
function tests = test_invz_twolevel_ordered_domain
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% Default is unchanged: exactly B = 0, hz = 0 still throws, deliberately (README.html:208).
function test_default_still_throws_at_zero_field(testCase)
ion = invz_ion();
verifyError(testCase, @() invz_twolevel_ordered(ion, 0.31, [0 0 0], 0, ...
    struct('Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x')), 'invz:degenerateDoublet');
end

function test_valid_flag_present_on_the_normal_path(testCase)
ion = invz_ion();
tl = invz_twolevel_ordered(ion, 0.31, [2.85 0 0], 0.02, ...
    struct('Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x'));
verifyTrue(testCase, tl.valid);
verifyGreaterThan(testCase, tl.Delta, 1e-4);
verifyTrue(testCase, isfield(tl, 'M2') && isfield(tl, 'g0'));
end

% Return mode reports the domain outcome instead of throwing, and reports the Delta it measured.
function test_return_mode_reports_instead_of_throwing(testCase)
ion = invz_ion();
o = struct('Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x', 'domain_policy', 'return');
tl = invz_twolevel_ordered(ion, 0.31, [0 0 0], 0, o);
verifyFalse(testCase, tl.valid);
verifyLessThan(testCase, tl.Delta, 1e-4);
% the two-level fields are deliberately ABSENT: a consumer ignoring .valid must fail loudly
verifyFalse(testCase, isfield(tl, 'g0'));
end

% Return mode is behaviour-neutral where the doublet is resolved.
function test_return_mode_identical_when_valid(testCase)
ion = invz_ion();
base = struct('Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x');
a = invz_twolevel_ordered(ion, 0.31, [2.85 0 0], 0.02, base);
b = invz_twolevel_ordered(ion, 0.31, [2.85 0 0], 0.02, ...
    setfield(base, 'domain_policy', 'return'));  %#ok<SFLD>
verifyEqual(testCase, b, a);
end

function test_unknown_policy_is_a_wiring_error(testCase)
ion = invz_ion();
verifyError(testCase, @() invz_twolevel_ordered(ion, 0.31, [2.85 0 0], 0.02, ...
    struct('Jxx0', ion.Jxx0, 'domain_policy', 'maybe')), 'invz:twolevelDomainPolicy');
end
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_twolevel_ordered_domain.m'); disp(table(r))"`
Expected: the `valid`-flag and return-mode tests fail; the default-throw test passes already.

- [ ] **Step 4: Add the domain policy**

In `invz_common/invz_twolevel_ordered.m`, replace lines 16–21:

```matlab
si = invz_single_ion(ion, T, invz_field_vec(Bx), struct('hyp', false, 'hz_fixed', hz, 'Jxx0', Jxx0, 'transverse_mf', tmf));
tl.Delta = si.E(2) - si.E(1);
if tl.Delta < 1e-4
    error('invz:degenerateDoublet', ...
        'Doublet splitting %.2e meV too small for the two-level Sigma.', tl.Delta);
end
```

with:

```matlab
si = invz_single_ion(ion, T, invz_field_vec(Bx), struct('hyp', false, 'hz_fixed', hz, 'Jxx0', Jxx0, 'transverse_mf', tmf));
tl.Delta = si.E(2) - si.E(1);
tl.valid = true;
% Domain policy (spec SS5.3). Delta is measured at THIS node's molecular field hz, so the
% geometric h-grid clustered at 0 puts the predictor and lowest nodes at risk whenever Bx is
% small -- not only at exactly Bx = 0. 'throw' (DEFAULT) is the historical behaviour and stays
% the contract for every existing caller. 'return' reports the outcome instead, using the SAME
% single diagonalization already performed above: it returns before g0 or any division that
% assumes a resolved doublet, so an invalid tl carries only .valid and .Delta and a consumer
% that ignores .valid fails loudly rather than silently computing with a degenerate doublet.
% Do NOT pre-screen by repeating this diagonalization and then calling the constructor again.
% The 1e-4 meV floor is unchanged, not re-tuned.
policy = getf(opts, 'domain_policy', 'throw');
if ~any(strcmp(policy, {'throw', 'return'}))
    error('invz:twolevelDomainPolicy', ...
        'opts.domain_policy must be ''throw'' or ''return''; got ''%s''.', policy);
end
if tl.Delta < 1e-4
    if strcmp(policy, 'return')
        tl.valid = false;
        return;                              % ONLY .Delta and .valid are defined
    end
    error('invz:degenerateDoublet', ...
        'Doublet splitting %.2e meV too small for the two-level Sigma.', tl.Delta);
end
```

- [ ] **Step 5: Run the new test, then every `twolevel` consumer**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_twolevel_ordered_domain.m'); disp(table(r))"`
Expected: 5 passed.
Then the consumers (a new always-present `tl.valid` field could break an `isequal`-style pin):
Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests({'invz_projected/tests/test_invz_matsubara_twolevel.m','invz_projected/tests/test_invz_gstat_ordered.m','invz_projected/tests/test_invz_hmf_ordered.m','invz_projected/tests/test_invz_deltaF.m'}); disp(table(r))"`
Expected: all pass. If any test compares a whole `tl` struct with `isequal`/`verifyEqual`, it will now see the extra field — update that test to compare the fields it actually cares about, and say so in the commit message.

- [ ] **Step 6: Run the full suite**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests'); fprintf('P=%d F=%d I=%d T=%d\n', nnz([r.Passed]), nnz([r.Failed]), nnz([r.Incomplete]), numel(r))"`
Expected: `F=0`; no baseline-passing test changes state.

- [ ] **Step 7: Commit**

```bash
git add invz_common/invz_twolevel_ordered.m invz_projected/tests/test_invz_twolevel_ordered_domain.m
git commit -m "feat(invzp): add single-evaluation domain_policy to invz_twolevel_ordered"
```

---

### Task 12: `invz_hmf_ordered` — grid helper, moment threading, per-node diagnostics, path integrals

The measurement the spec's binding caution demands (`integral (r-1) dh`, not just `integral Sigma0 dh`) is free: it reads a profile that already exists. `crit_star` needs one discarded output recovered at `:243`.

**Files:**
- Create: `invz_projected/invz_hmf_grid.m`
- Modify: `invz_projected/invz_hmf_ordered.m` (`:66-71` prof init, `:114-116` opts, `:144-145` grid, `:152-156` integrals, `:190-198` re-densify + export, `:243-247` root, `:262-333` eval_node)
- Test: `invz_projected/tests/test_invz_hmf_grid.m`, `invz_projected/tests/test_invz_hmf_ordered_diagnostics.m` (new)
- Regression: `invz_projected/tests/test_invz_hmf_ordered.m` (existing)

**Interfaces:**
- Consumes: Tasks 1, 4, 8, 10.
- Produces: `[hgrid, ratio] = invz_hmf_grid(hmax, nH, hmin_frac)`. `opts.Jmom`
  (optional). `prof` gains predictor seeds `r_pm0`, `G0bare_pm0`; per-node `crit`,
  `r_minus_1`, `Delta`, `medium_status`, `node_term_reason`, `Dq_min`, `ref_denom`,
  `ref_margin`, `gstat_local_denom`, `omit_mu3`, `omit_cubic`, `omit_max`; and endpoint/integral fields
  `G0bare_star`, `crit_star`, `Dq_min_star`, `int_Sigma0`, `int_r_minus_1`. These solved-path
  fields are load-bearing for Gate 0; a prospective scan cannot substitute for them. Used by
  Tasks 13, 14, 16, 17.

- [ ] **Step 1: Write the failing grid-helper test**

Create `invz_projected/tests/test_invz_hmf_grid.m`:

```matlab
function tests = test_invz_hmf_grid
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% Byte-identical to the expression it replaces (invz_hmf_ordered.m:144-145).
function test_matches_the_inline_expression(testCase)
hmax = 0.37;  nH = 33;  hfrac = 1e-3;
ratio_ref = hfrac^(1/(nH-1));
hgrid_ref = hmax * ratio_ref.^((nH-1):-1:0);
[hgrid, ratio] = invz_hmf_grid(hmax, nH, hfrac);
verifyEqual(testCase, hgrid, hgrid_ref, 'AbsTol', 0);
verifyEqual(testCase, ratio, ratio_ref, 'AbsTol', 0);
end

function test_shape_and_endpoints(testCase)
[hgrid, ~] = invz_hmf_grid(0.5, 17, 1e-4);
verifyEqual(testCase, numel(hgrid), 17);
verifyEqual(testCase, hgrid(end), 0.5, 'RelTol', 1e-15);
verifyEqual(testCase, hgrid(1), 0.5*1e-4, 'RelTol', 1e-12);
verifyTrue(testCase, all(diff(hgrid) > 0));            % ascending, clustered at 0
end

function test_rejects_bad_arguments(testCase)
verifyError(testCase, @() invz_hmf_grid(0, 33, 1e-3), 'invz:hmfGrid');
verifyError(testCase, @() invz_hmf_grid(0.5, 1, 1e-3), 'invz:hmfGrid');
verifyError(testCase, @() invz_hmf_grid(0.5, 33, 0), 'invz:hmfGrid');
verifyError(testCase, @() invz_hmf_grid(0.5, 33, 1.5), 'invz:hmfGrid');
end
```

- [ ] **Step 2: Run it to verify it fails, then write the helper**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_hmf_grid.m'); disp(table(r))"`
Expected: all 3 error on the missing function.

Create `invz_projected/invz_hmf_grid.m`:

```matlab
function [hgrid, ratio] = invz_hmf_grid(hmax, nH, hmin_frac)
%INVZ_HMF_GRID The initial geometric H_MF profile grid, clustered at 0 (P1-4).
% Extracted VERBATIM from invz_hmf_ordered.m:144-145 so the prospective Gate-0 scanner
% (invz_static_domain_scan) and the solver build the SAME initial grid from one definition
% rather than two implementations that happen to agree (spec SS7.2).
%
% The scanner must NOT reproduce invz_hmf_ordered's adaptive extension or redensification: this
% helper covers the INITIAL grid only, and solved-path margins are read off prof.hgrid AFTER
% adaptation.
%   hgrid = hmax * ratio.^((nH-1):-1:0),  ratio = hmin_frac^(1/(nH-1))
% ascending, hgrid(end) = hmax, hgrid(1) = hmax*hmin_frac.
if ~(isscalar(hmax) && isfinite(hmax) && hmax > 0)
    error('invz:hmfGrid', 'hmax must be a positive finite scalar; got %s', mat2str(hmax));
end
if ~(isscalar(nH) && nH == round(nH) && nH >= 2)
    error('invz:hmfGrid', 'nH must be an integer >= 2; got %s', mat2str(nH));
end
if ~(isscalar(hmin_frac) && isfinite(hmin_frac) && hmin_frac > 0 && hmin_frac < 1)
    error('invz:hmfGrid', 'hmin_frac must be in (0,1); got %s', mat2str(hmin_frac));
end
ratio = hmin_frac^(1/(nH-1));
hgrid = hmax * ratio.^((nH-1):-1:0);
end
```

Run the test again. Expected: 3 passed.

- [ ] **Step 3: Write the failing diagnostics test**

Create `invz_projected/tests/test_invz_hmf_ordered_diagnostics.m`:

```matlab
function tests = test_invz_hmf_ordered_diagnostics
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function [ion, T, Bx, Jnu, o] = fixture()
ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
           'static_medium', 'strict_1z_dyson_ref');
end

% Every new per-node array is present and the same length as hgrid, on ANY exit path.
function test_per_node_arrays_are_present_and_aligned(testCase)
[ion, T, Bx, Jnu, o] = fixture();
[~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
n = numel(prof.hgrid);
for f = {'crit', 'r_minus_1', 'Delta', 'Dq_min', 'ref_denom', 'ref_margin', ...
         'gstat_local_denom', ...
         'omit_mu3', 'omit_cubic', 'omit_max'}
    verifyEqual(testCase, numel(prof.(f{1})), n, f{1});
end
verifyEqual(testCase, numel(prof.medium_status), n);
verifyEqual(testCase, numel(prof.node_term_reason), n);
verifyTrue(testCase, iscell(prof.medium_status));
end

% crit is the dimensionless mass r + J0eff*G0bare. The predictor identity is tested at h=0
% directly; do not compare it to the first nonzero geometric node with an arbitrary 5% band.
function test_crit_definition_and_slope0_consistency(testCase)
[ion, T, Bx, Jnu, o] = fixture();
[~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyEqual(testCase, prof.crit, prof.r + o.J0eff*prof.G0bare, 'RelTol', 1e-12);
verifyEqual(testCase, prof.slope0, prof.r_pm0 + o.J0eff*prof.G0bare_pm0, 'RelTol', 1e-12);
end

% r_minus_1 is r - 1 -- NOT Sigma0. They coincide only at m = 0 (spec §A Consequence 2), so a
% finite-moment node must show them differing, or the whole (r-1) measurement is vacuous.
function test_r_minus_1_is_not_sigma0_at_finite_moment(testCase)
[ion, T, Bx, Jnu, o] = fixture();
[~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyEqual(testCase, prof.r_minus_1, prof.r - 1, 'AbsTol', 0);
k = find(abs(prof.m) > 1e-3 & isfinite(prof.Sigma0) & isfinite(prof.r), 1, 'last');
verifyNotEmpty(testCase, k, 'fixture produced no finite-moment converged node');
verifyGreaterThan(testCase, abs(prof.r_minus_1(k) - prof.Sigma0(k)), 1e-9, ...
    'r - 1 must differ from Sigma0 at finite m (spec §A Consequence 2)');
end

% Both path integrals use the SAME first-panel seeding as h0, so the three are
% quadrature-consistent by construction.
function test_path_integrals_match_the_h0_quadrature(testCase)
[ion, T, Bx, Jnu, o] = fixture();
[~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyEqual(testCase, prof.status, 'ok', 'fixture must be a deterministic accepted profile');
verifyEqual(testCase, prof.int_r_minus_1, ...
    trapz([0 prof.hgrid], [prof.r_pm0 - 1, prof.r - 1]), ...
    'RelTol', 1e-9, 'int_r_minus_1 must use the actual h=0 r seed');
verifyEqual(testCase, prof.int_Sigma0, ...
    trapz([0 prof.hgrid], [prof.Sigma0_pm0, prof.Sigma0]), 'RelTol', 1e-9);
end

% crit_star requires the root call to stop discarding G0bare.
function test_crit_star_at_the_root(testCase)
[ion, T, Bx, Jnu, o] = fixture();
[hstar, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyTrue(testCase, isfinite(hstar), 'fixture must produce the preregistered root');
verifyTrue(testCase, isfinite(prof.G0bare_star));
verifyEqual(testCase, prof.crit_star, prof.r_star + o.J0eff*prof.G0bare_star, 'RelTol', 1e-12);
% the accepted root must be the INCREASING crossing of F, i.e. a Landau minimum
verifyGreaterThan(testCase, prof.crit_star, 0);
end

% Supplying Jmom must not change any result (hot-path optimisation, not a semantic change).
function test_supplied_jmom_is_behaviour_neutral(testCase)
[ion, T, Bx, Jnu, o] = fixture();
[h1, p1] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
o.Jmom = invz_coupling_moments(Jnu);
[h2, p2] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyEqual(testCase, h2, h1, 'AbsTol', 0);
verifyEqual(testCase, p2.r, p1.r, 'AbsTol', 0);
verifyEqual(testCase, p2.crit, p1.crit, 'AbsTol', 0);
end
```

- [ ] **Step 4: Run it to verify it fails**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_hmf_ordered_diagnostics.m'); disp(table(r))"`
Expected: all fail on the missing `prof` fields.

- [ ] **Step 5: Extend the `prof` initializer**

In `invz_hmf_ordered.m`, replace the `prof = struct(...)` block (`:66-71`) with:

```matlab
prof = struct('hgrid', [], 'r', [], 'h0', [], 'm', [], 'Sigma0', [], 'K0', [], ...
              'D_uni', [], 'G0bare', [], 'Gstat', [], 'node_conv', [], 'F', [], ...
              'crit', [], 'r_minus_1', [], 'Delta', [], 'Dq_min', [], ...
              'ref_denom', [], 'ref_margin', [], 'gstat_local_denom', [], ...
              'omit_mu3', [], 'omit_cubic', [], ...
              'omit_max', [], 'medium_status', {{}}, 'node_term_reason', {{}}, ...
              'slope0', NaN, 'r_pm0', NaN, 'G0bare_pm0', NaN, ...
              'Sigma0_pm0', NaN, 'K0_pm0', NaN, 'J0eff', J0eff, ...
              'n_extend', 0, 'hmin_initial', NaN, 'status', 'no_bare_order', ...
              'redensified', false, 'int_Sigma0', NaN, 'int_r_minus_1', NaN, ...
              'm_star', NaN, 'D_uni_star', NaN, 'r_star', NaN, 'Gstat_star', NaN, ...
              'G0bare_star', NaN, 'crit_star', NaN, 'Dq_min_star', NaN, ...
              'static_medium', sm.scheme);
```

- [ ] **Step 6: Resolve the moments once and use the grid helper**

Resolve `eopts` and `eso` near the top of the public function, before the bare-order early return:

```matlab
eopts = getf(opts, 'emt', struct());
eso   = getf(opts, 'emt_static', struct());
[sm, eopts, eso] = invz_check_static_medium(opts, eopts, eso);
eso.warn = false;
```

Do not resolve `eopts` a second time near the Matsubara block. Then insert:

```matlab
% Coupling moments, resolved ONCE per call (spec SS4.3): recomputing them inside the outer
% iteration of every node would repeat an O(nJ) pass up to max_outer x nNodes times per field.
% A caller that already resolved the point's coupling spectrum (invz_solve_point_ordered) passes
% opts.Jmom straight through; a direct call derives it here as a compatibility fallback.
Jmom = getf(opts, 'Jmom', []);
if isempty(Jmom), Jmom = invz_coupling_moments(Jnu_flat); end
if sm.is_strict && ~isvector(Jnu_flat)
    error('invz:staticMedium', 'strict ordered HMF does not support [nJ,nw] Jnu_flat.');
end
```

Replace `:144-145`:

```matlab
ratio = hfrac^(1/(nH-1));
hgrid = hmax * ratio.^((nH-1):-1:0);                 % geometric, clustered at 0 (P1-4)
```

with:

```matlab
[hgrid, ratio] = invz_hmf_grid(hmax, nH, hfrac);     % geometric, clustered at 0 (P1-4)
```

And in the re-densify block replace `:186-187`:

```matlab
    ratio2 = hfrac_eff^(1/(nH-1));
    hgrid = hmax * ratio2.^((nH-1):-1:0);
```

with:

```matlab
    hgrid = invz_hmf_grid(hmax, nH, hfrac_eff);
```

- [ ] **Step 7: Return a per-node record; do not grow the positional-output list**

The earlier three-extra-output sketch was insufficient for Gate 0 and would turn `eval_node` into a
fragile 18-output interface. Change the nested function to:

```matlab
function [rec, Sigma, K0s] = eval_node(hp, Sigma, K0s)
```

with a fixed record schema present on every exit:

```matlab
rec = struct('r',NaN,'m',NaN,'Sigma0',NaN,'K0',NaN,'D_uni',NaN,'Dq_min',NaN, ...
    'G0bare',NaN,'Gstat',NaN,'accepted',false,'crit',NaN,'Delta',NaN, ...
    'medium_status','not_applicable','term_reason','not_evaluated', ...
    'ref_denom',NaN,'ref_margin',NaN,'gstat_local_denom',NaN, ...
    'omit_mu3',NaN,'omit_cubic',NaN,'omit_max',NaN);
```

Populate it from the direct diagonalisation plus `state/info`; add `Jmom` to the node passed to
`invz_ordered_node_solve`. `ref_margin` is the actual distance-to-floor
`info.medium.ref.margin`, not the denominator or configured floor. `Dq_min` comes from
`info.res.stability.Dq_min`; omitted ratios come from `info.so`; `term_reason` comes from `info`.
`gstat_local_denom` comes from `info.so.gstat_local_denom`.
The `fbare` and Task-13 degenerate exits populate the same record. `run_sweep` returns a struct array
of these records plus the continuation carriers. Update initial, extension, redensification,
bisection, and root call sites together. A schema test must cover normal, domain, degenerate, and
bare-shortcut records. When `opts.trace=true`, route every exit through one nested
`append_trace_node(rec,info)` finalizer that copies the same record fields into each `trc.nodes`
entry (alongside its existing id/phase/seed metadata). No domain/bare early return may append a
shorter schema by hand. Gate 0 uses that trace ledger to account for the predictor,
extension/redensification nodes, bisection iterates, and final root—not just the final integration
grid.

- [ ] **Step 8: Export every solved-path array and compute the two path integrals**

```matlab
prof.hgrid = hgrid;  prof.r = rv;  prof.h0 = h0;  prof.m = mv;
prof.Sigma0 = S0v;   prof.K0 = K0v;  prof.D_uni = Dv;  prof.node_conv = cnv;  prof.F = F;
prof.G0bare = Gbv;   prof.Gstat = Gsv;
prof.crit = [nodes.crit];  prof.r_minus_1 = rv - 1;  prof.Delta = [nodes.Delta];
prof.Dq_min = [nodes.Dq_min];  prof.ref_denom = [nodes.ref_denom];
prof.ref_margin = [nodes.ref_margin];
prof.gstat_local_denom = [nodes.gstat_local_denom];
prof.omit_mu3 = [nodes.omit_mu3];
prof.omit_cubic = [nodes.omit_cubic];  prof.omit_max = [nodes.omit_max];
prof.medium_status = {nodes.medium_status};  prof.node_term_reason = {nodes.term_reason};
% Path corrections (spec SSA binding caution, gate G5). BOTH are needed: at finite ordered
% moment r - 1 is NOT Sigma0 -- the hybrid elastic factor xi makes r depend on K0 and
% lambda(1:2) -- so integral Sigma0 dh is a component diagnostic, not the whole correction.
% The ~0.3% PM boundary shift bounds NEITHER of them deep in the ordered phase.
% Same first-panel seeding as h0 above (:152), so all three integrals are quadrature-consistent.
if ok0
    prof.int_Sigma0    = trapz([0 hgrid], [prof.Sigma0_pm0, S0v]);
    prof.int_r_minus_1 = trapz([0 hgrid], [r0n - 1, rv - 1]);
end
```

Derive `rv`, `mv`, `S0v`, `K0v`, `Dv`, `Gbv`, `Gsv`, and `cnv` from the final record array
at one explicit point; do not maintain a second parallel set that can become misaligned. Store the
predictor values when `ok0`: `prof.r_pm0=r0n`, `prof.G0bare_pm0=Gb0`,
`prof.Sigma0_pm0=S0pm`, `prof.K0_pm0=K0pm`. Recompute local `h0/F` whenever the grid changes;
export `prof` and compute both integrals once from the **final** adapted grid. This avoids briefly
publishing integrals that describe a superseded grid.

- [ ] **Step 9: Export the complete root record**

```matlab
[root, ~, ~] = eval_node(hmf_star, Sigma, K0s);
if ~root.accepted
    prof.status = 'node_failed';  hmf_star = NaN;  return;
end
prof.m_star = root.m;  prof.D_uni_star = root.D_uni;  prof.Dq_min_star = root.Dq_min;
prof.r_star = root.r;  prof.Gstat_star = root.Gstat;  prof.G0bare_star = root.G0bare;
prof.crit_star = root.crit;
```

- [ ] **Step 10: Run the new tests, the existing HMF test, and the full suite**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests({'invz_projected/tests/test_invz_hmf_ordered_diagnostics.m','invz_projected/tests/test_invz_hmf_ordered.m','invz_projected/tests/test_invz_ordered_trace.m','invz_projected/tests/test_invz_deltaF.m'}); disp(table(r))"`
Expected: the new file 6 passed; all three existing files pass unchanged. `invz_ordered_trace` and `invz_deltaF` both consume `prof`, so they are the direct blast-radius check for the schema change.
Then the full suite: `F=0`; no baseline-passing test changes state.

- [ ] **Step 11: Commit**

```bash
git add invz_projected/invz_hmf_grid.m invz_projected/invz_hmf_ordered.m invz_projected/tests/test_invz_hmf_grid.m invz_projected/tests/test_invz_hmf_ordered_diagnostics.m
git commit -m "feat(invzp): hmf profile gains crit/r-1/Delta diagnostics and both path integrals"
```

---

### Task 13: `invz_hmf_ordered` — `Delta` domain screen and status precedence

Two new `prof.status` values are otherwise unreachable: an out-of-domain or degenerate node is *also* a non-accepted node, so `any(~cnv)` would claim it as `'node_failed'` first. The two new cases are **prepended** to the existing chain, leaving the legacy relative order of `'unresolved'` vs `'node_failed'` byte-identical.

**Files:**
- Create: `invz_projected/invz_hmf_status.m` (pure precedence reducer)
- Modify: `invz_projected/invz_hmf_ordered.m` (`:124` predictor call, `:200-215` status block, `:276-281` eval_node's two-level construction)
- Test: `invz_projected/tests/test_invz_hmf_ordered_status.m` (new)

**Interfaces:**
- Consumes: Task 11's `domain_policy`; Task 12's per-node arrays.
- Produces: `prof.status` may now be `'degenerate_doublet'` or `'medium_out_of_domain'`. Used by Tasks 14, 15.

- [ ] **Step 1: Write the failing test**

Create `invz_projected/tests/test_invz_hmf_ordered_status.m`:

```matlab
function tests = test_invz_hmf_ordered_status
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% B = 0 no longer escapes as a throw: it is a LABELLED domain outcome, distinguishable from a
% solver failure. Production sweeps start at B = 0, so this path is hit on every run.
function test_zero_field_is_degenerate_not_a_throw(testCase)
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
[hstar, prof] = invz_hmf_ordered(ion, 0.31, [0 0 0], Jnu, o);   % must not throw
verifyTrue(testCase, isnan(hstar));
verifyEqual(testCase, prof.status, 'degenerate_doublet');
verifyTrue(testCase, isfinite(prof.Delta(1)) || isnan(prof.Delta(1)));
end

% A medium domain event is its own status, NOT 'node_failed'.
function test_medium_domain_status_is_distinct(testCase)
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
           'static_medium', 'strict_1z_dyson_ref', ...
           'ref_margin', 1e9);                          % public authority: force domain event
[hstar, prof] = invz_hmf_ordered(ion, 0.31, [2.85 0 0], Jnu, o);
verifyTrue(testCase, isnan(hstar));
verifyEqual(testCase, prof.status, 'medium_out_of_domain');
end

% Legacy statuses keep their exact relative order: the two new cases are prepended only.
function test_legacy_status_order_unchanged(testCase)
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
[~, prof] = invz_hmf_ordered(ion, 0.31, [2.85 0 0], Jnu, o);
verifyTrue(testCase, any(strcmp(prof.status, ...
    {'ok', 'unresolved', 'node_failed', 'no_bare_order'})));
end

% Integration anchors for the two new statuses and the pre-existing no-bare-order exit.
function test_new_and_no_bare_statuses_are_reachable(testCase)
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
base = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
seen = {};
[~, p] = invz_hmf_ordered(ion, 0.31, [0 0 0], Jnu, base);          seen{end+1} = p.status;
o2 = base;  o2.static_medium = 'strict_1z_dyson_ref';
o2.ref_margin = 1e9;
[~, p] = invz_hmf_ordered(ion, 0.31, [2.85 0 0], Jnu, o2);         seen{end+1} = p.status;
[~, p] = invz_hmf_ordered(ion, 0.31, [9 0 0], Jnu, base);          seen{end+1} = p.status;
verifyTrue(testCase, ismember('degenerate_doublet', seen));
verifyTrue(testCase, ismember('medium_out_of_domain', seen));
verifyTrue(testCase, ismember('no_bare_order', seen));             % 9 T: bare does not order
end

% Pure precedence pin, including the mixed case node_failed > unresolved.
function test_status_precedence_all_five_cases(testCase)
base = struct('accepted', true, 'term_reason', 'accepted', 'medium_status', 'ok');
pred = base;  nodes = repmat(base, 1, 2);  F = [0.1 0.2];
verifyEqual(testCase, invz_hmf_status(pred, nodes, -1, F), 'unresolved');
nodes(1).accepted = false; nodes(1).term_reason = 'max_iter';
verifyEqual(testCase, invz_hmf_status(pred, nodes, -1, F), 'node_failed');
nodes(1).medium_status = 'ref_denom_small';
verifyEqual(testCase, invz_hmf_status(pred, nodes, -1, F), 'medium_out_of_domain');
nodes(2).term_reason = 'degenerate_doublet';
verifyEqual(testCase, invz_hmf_status(pred, nodes, -1, F), 'degenerate_doublet');
nodes = repmat(base, 1, 2);
verifyEqual(testCase, invz_hmf_status(pred, nodes, 1, [-0.1 0.2]), 'ok');
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_hmf_ordered_status.m'); disp(table(r))"`
Expected: the `B = 0` test errors with `invz:degenerateDoublet` escaping (proving the current behaviour), and the domain-status test reports `'node_failed'`.

- [ ] **Step 3: Use the return-mode two-level constructor in `eval_node`**

In `invz_hmf_ordered.m`, replace the two-level call inside the Task-12 record-returning
`eval_node`:

```matlab
    tl = invz_twolevel_ordered(ion, T, Bx, hp, struct('Jxx0', Jxx0, 'transverse_mf', tmf));
```

with:

```matlab
    % Delta domain screen (spec SS5.3), SINGLE evaluation: return mode reuses the one
    % diagonalization the constructor already performs, instead of pre-screening with a
    % duplicate one and then calling it again. Delta is measured at THIS node's molecular field
    % hp, and the geometric grid clusters at 0, so the predictor and lowest nodes are the ones
    % at risk whenever Bx is small -- not only at exactly Bx = 0. Previously the constructor's
    % throw escaped this function entirely and the column masked as a solver failure.
    tl = invz_twolevel_ordered(ion, T, Bx, hp, struct('Jxx0', Jxx0, 'transverse_mf', tmf, ...
                                                      'domain_policy', 'return'));
    if ~tl.valid
        rec.m = si.Jexp(3);  rec.Delta = tl.Delta;
        rec.accepted = false;  rec.term_reason = 'degenerate_doublet';
        if tracing, append_trace_node(rec, []); end  % Task-12 single schema finalizer
        return;
    end
```

- [ ] **Step 4: Capture the predictor node's new outputs**

Capture the predictor as a record:

```matlab
[pred, Sigma, K0s] = eval_node(0, Sigma, K0s);
r0n = pred.r;  S0pm = pred.Sigma0;  K0pm = pred.K0;  Gb0 = pred.G0bare;
ok0 = pred.accepted;
```

- [ ] **Step 5: Replace the status block with the precedence chain**

Replace `:200-215`:

```matlab
if ~ok0                          % predictor never converged: h0/F are undefined (NaN)
    prof.status = 'node_failed'; % above -- report the honest verdict now that the grid's
    return;                      % own (convergence-independent) diagnostics are exported;
end                              % NEVER fall through to the F-based search on NaN data
if slope_pred < 0 && all(F >= 0)                      % floor hit without a bracket:
    prof.status = 'unresolved';                       % NEVER silently PM (round-3 P0-3)
    warning('invz:hmfUnresolved', ...
        'ordering predicted (slope_pred = %.3g) but no negative F above hmin_abs = %.3g', ...
        slope_pred, hmin_abs);
    return;                                           % hmf_star stays NaN; the jensen solver
end                                                   % must return converged = false here
if any(~cnv)                                          % round-4 P1-C: status must be truthful
    prof.status = 'node_failed';                      % on node failure -- never 'ok'
    return;
end
prof.status = 'ok';
```

with:

```matlab
prof.status = invz_hmf_status(pred, nodes, slope_pred, F);
if strcmp(prof.status, 'unresolved')
    warning('invz:hmfUnresolved', ...
        'ordering predicted (slope_pred = %.3g) but no negative F above hmin_abs = %.3g', ...
        slope_pred, hmin_abs);
end
if ~strcmp(prof.status, 'ok'), return; end
```

Create the pure helper used above:

```matlab
function status = invz_hmf_status(pred, nodes, slope_pred, F)
%INVZ_HMF_STATUS Reduce Jensen node records with binding reason precedence.
% degenerate > reference-domain > node failure > unresolved > ok.
allnodes = [pred, nodes];
reasons = {allnodes.term_reason};
medium  = {allnodes.medium_status};
if any(strcmp(reasons, 'degenerate_doublet'))
    status = 'degenerate_doublet';
elseif any(~cellfun(@(s) any(strcmp(s, {'ok','not_applicable'})), medium))
    status = 'medium_out_of_domain';
elseif any(~[allnodes.accepted])
    status = 'node_failed';
elseif slope_pred < 0 && all(F >= 0)
    status = 'unresolved';
else
    status = 'ok';
end
end
```

Use the same helper before **every** early return caused by a predictor, extension, redensification,
bisection, or final-root record. Pass the records evaluated so far; domain/degenerate reasons must
not revert to `node_failed` merely because they occur after the initial profile. This order matches
the binding precedence `degenerate_doublet > medium_out_of_domain > node_failed > unresolved > ok`;
the pure mixed fixture pins `node_failed > unresolved`.

- [ ] **Step 6: Run the new test and the HMF regression, then the full suite**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests({'invz_projected/tests/test_invz_hmf_ordered_status.m','invz_projected/tests/test_invz_hmf_ordered.m','invz_projected/tests/test_invz_hmf_ordered_diagnostics.m'}); disp(table(r))"`
Expected: new file 5 passed; both existing files pass unchanged. **If an existing test asserted that
`B = 0` throws `invz:degenerateDoublet` through `invz_hmf_ordered`, that expectation is now wrong by
design** — update it to assert `prof.status == 'degenerate_doublet'` and note the contract change in
the commit message.
Full suite: `F=0`; no baseline-passing test changes state.

- [ ] **Step 7: Commit**

```bash
git add invz_projected/invz_hmf_status.m invz_projected/invz_hmf_ordered.m \
  invz_projected/tests/test_invz_hmf_ordered_status.m
git commit -m "feat(invzp): label degenerate-doublet and medium-domain profiles distinctly"
```

---

### Task 14: `invz_solve_point` / `invz_solve_point_ordered` — resolve once, thread, expose provenance

**Files:**
- Modify: `invz_projected/invz_solve_point.m` (opts block; `pt` export)
- Modify: `invz_projected/invz_solve_point_ordered.m` (`:63-76` opts, `:130-135` hopts, `:187-214` jensen block, `:246-265` jensen exports)
- Test: `invz_projected/tests/test_invz_solve_point_strict.m` (new)

**Interfaces:**
- Consumes: Tasks 1, 4, 12, 13.
- Produces: both point solvers accept `opts.static_medium` / `opts.ref_margin`.
  `invz_solve_point` gains `.static_medium`, `.Jmom`, `.medium_status`, `.medium_denom`,
  `.medium_margin`. `invz_solve_point_ordered` additionally gains `.stable_1z`, `.crit_1z`,
  `.Dq_min`, endpoint `.omit_mu3/.omit_cubic/.omit_max`, and `.path_omit_max`. Used by Task 15.

- [ ] **Step 1: Write the failing test**

Create `invz_projected/tests/test_invz_solve_point_strict.m`:

```matlab
function tests = test_invz_solve_point_strict
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function [ion, T, Bx, Jnu, o] = fx()
ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
end

% Provenance is mandatory on BOTH legs: an unlabelled result cannot be compared across schemes.
function test_pm_leg_carries_provenance(testCase)
[ion, T, Bx, Jnu, o] = fx();
pt = invz_solve_point(ion, T, Bx, Jnu, o);
verifyEqual(testCase, pt.static_medium, 'resummed');
verifyTrue(testCase, isfield(pt, 'Jmom') && isfield(pt, 'medium_margin'));
o.static_medium = 'strict_1z_dyson_ref';
pts = invz_solve_point(ion, T, Bx, Jnu, o);
verifyEqual(testCase, pts.static_medium, 'strict_1z_dyson_ref');
verifyEqual(testCase, pts.Jmom.mu2, invz_coupling_moments(Jnu).mu2, 'AbsTol', 0);
end

% Default absent => legacy numbers unchanged (G9 at point level).
function test_absent_scheme_reproduces_legacy_numbers(testCase)
[ion, T, Bx, Jnu, o] = fx();
a = invz_solve_point(ion, T, Bx, Jnu, o);
b = invz_solve_point(ion, T, Bx, Jnu, setfield(o, 'static_medium', 'resummed'));  %#ok<SFLD>
verifyEqual(testCase, b.Sigma0, a.Sigma0, 'AbsTol', 0);
verifyEqual(testCase, b.crit,   a.crit,   'AbsTol', 0);
end

% Strict changes the PM mass -- it must, since K(0) changed (spec §0.2).
function test_strict_shifts_the_pm_mass(testCase)
[ion, T, Bx, Jnu, o] = fx();
a = invz_solve_point(ion, T, Bx, Jnu, o);
b = invz_solve_point(ion, T, Bx, Jnu, setfield(o, 'static_medium', 'strict_1z_dyson_ref'));  %#ok<SFLD>
verifyNotEqual(testCase, b.crit, a.crit);
verifyTrue(testCase, isfinite(b.crit));
end

% The ordered leg exposes the two tiers SEPARATELY, and converged keeps its old meaning.
function test_ordered_leg_two_tier_fields(testCase)
[ion, T, Bx, Jnu, o] = fx();
o.ordered_mode = 'jensen';  o.static_medium = 'strict_1z_dyson_ref';
pt = invz_solve_point_ordered(ion, T, Bx, Jnu, o);
verifyTrue(testCase, isfield(pt, 'stable_1z') && isfield(pt, 'crit_1z'));
verifyTrue(testCase, isfield(pt, 'omit_max'));
verifyEqual(testCase, pt.static_medium, 'strict_1z_dyson_ref');
if pt.is_ordered && isfinite(pt.crit_1z)
    Dtol = 1e-6*max(1, abs(pt.G(1))*max(abs(Jnu)));
    expected = pt.converged && pt.crit_1z > 1e-6 && ...
               pt.D_uni > Dtol && pt.Dq_min > Dtol;
    verifyEqual(testCase, pt.stable_1z, expected);
end
end

% A per-leg override is rejected at the public entry, so the sectors cannot split.
function test_per_leg_override_rejected(testCase)
[ion, T, Bx, Jnu, o] = fx();
o.emt = struct('static_medium', 'strict_1z_dyson_ref');
verifyError(testCase, @() invz_solve_point(ion, T, Bx, Jnu, o), 'invz:staticMedium');
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_solve_point_strict.m'); disp(table(r))"`
Expected: fails on the missing provenance fields.

- [ ] **Step 3: Wire `invz_solve_point`**

Resolve the scheme after reading `eopts`, but compute moments only **after** the entire ODD branch has
resolved the actual `Jnu_flat`. The previous sketch computed moments immediately after `eopts`, where
ODD still has `Jnu_flat=[]`; that would fail and describe the wrong coupling.

```matlab
% Static-medium scheme: resolved ONCE here and stamped into the leg option struct, so
% invz_emt_scalar's omega = 0 slot and the ordered sector can never run different truncation
% orders (spec SS4.2). Absent => 'resummed' => numerically identical to the pre-strict path.
[sm, eopts] = invz_check_static_medium(opts, eopts);
```

Immediately before the Matsubara/single-ion block, after ODD/retardation:

```matlab
Jmom = invz_coupling_moments(Jnu_flat);
eopts.Jmom = Jmom;                            % matrix moments allowed; PM slot selects column 1
```

Inside the PM outer loop, immediately after `invz_emt_scalar`, stop before lambdas/Sigma when
`med.medium_status` is neither `'ok'` nor `'not_applicable'`. Return a non-converged point with that
exact domain status; do not iterate NaNs to `max_outer`.

On every full-solve return path, before returning `pt`, add:

```matlab
pt.static_medium = sm.scheme;
pt.Jmom = Jmom;
ref = getf(getf(med, 'medium', struct()), 'ref', struct());
pt.medium_status = getf(med, 'medium_status', 'not_applicable');
pt.medium_denom  = getf(ref, 'denom', NaN);
pt.medium_margin = getf(ref, 'margin', NaN);   % distance to floor, NOT the denominator
```

(For early returns that build `pt` from a fixed field set, add `pt.static_medium = sm.scheme; pt.Jmom = Jmom; pt.medium_margin = NaN;` so no caller ever probes a missing member — the same discipline `early_return` already follows in the ordered solver.)

- [ ] **Step 4: Wire `invz_solve_point_ordered`**

After `eopts = getf(opts, 'emt', struct());` (`:73`), insert:

```matlab
eso_pub = getf(opts, 'emt_static', struct());
[sm, eopts, eso_pub] = invz_check_static_medium(opts, eopts, eso_pub);
```

Replace the ODD/`Jnu_flat` resolution tail so the moments are taken **after** the coupling spectrum is final (the ODD branch rebuilds `Jnu_flat`), i.e. immediately before `[wn, wts, beta] = invz_matsubara(T, Ecut);` (`:110`) insert:

```matlab
% Moments AFTER the point's coupling spectrum is resolved: the ODD branch above rebuilds
% Jnu_flat from odd_blocks + deltaJ, so taking them earlier would describe the wrong multiset.
Jmom = invz_coupling_moments(Jnu_flat);
eopts.Jmom = Jmom;  eso_pub.Jmom = Jmom;
if sm.is_strict && ~isvector(Jnu_flat)
    error('invz:staticMedium', ['strict ordered/Jensen mode does not support the [nJ,nw] ' ...
        'retarded coupling matrix in this phase; PM strict mode remains supported.']);
end
```

Also reject `opts.odd_retarded` / `opts.odd_retarded_exact` explicitly on the strict ordered entry;
there is no invented `ordered_retarded` option. The check occurs after actual option resolution and
before HMF starts.

In the jensen dispatch (`:130-135`), forward the resolved context:

```matlab
    hopts = opts;                                    % FULL numerical context (P1-6) ...
    hopts.J0eff = J0eff;                             % ... with the ODD-shifted coupling
    hopts.Jmom = Jmom;                               % ... the resolved moments (no re-derive)
    hopts.static_medium = sm.scheme;                 % ... and the resolved scheme
    hopts.emt = eopts;  hopts.emt_static = eso_pub;  % ... stamped (validation is idempotent)
    for f = {'ordered_mode', 'forced_moment'}        % ... and mode fields stripped
        if isfield(hopts, f{1}), hopts = rmfield(hopts, f{1}); end
    end
```

In the jensen solve block, replace `eso = getf(opts, 'emt_static', struct());` (`:187`) with `eso = eso_pub;` and add `'Jmom', Jmom` to the `node` struct (`:209-211`).

Finally, in the jensen export block (`:246-265`), after `pt.D_uni = info.so.D_uni;` add:

```matlab
    % TWO-TIER EXPORT (spec SS1): pt.converged keeps its existing meaning (info.accepted, the
    % consistency tier) so no existing consumer shifts meaning under it. Endpoint stability is
    % a SEPARATE field -- collapsing them would re-mask the ordered phase, because intermediate
    % path nodes are the unstable Landau interval by construction.
    pt.crit_1z = hprof.crit_star;
    pt.Dq_min = info.res.stability.Dq_min;
    pt.stable_1z = pt.converged && info.res.stability.pass && ...
                   isfinite(hprof.crit_star) && ...
                   abs(hprof.crit_star-info.res.stability.crit) <= ...
                       getf(opts,'crit_tol',1e-6);
    pt.omit_mu3 = getf(info.so, 'omit_mu3', NaN);
    pt.omit_cubic = getf(info.so, 'omit_cubic', NaN);
    pt.omit_max = getf(info.so, 'omit_max', NaN);
    finite_omit = hprof.omit_max(isfinite(hprof.omit_max));
    if isempty(finite_omit), pt.path_omit_max = NaN;
    else,                    pt.path_omit_max = max(finite_omit); end
    ref = getf(getf(info, 'medium', struct()), 'ref', struct());
    pt.medium_status = getf(info, 'medium_status', 'not_applicable');
    pt.medium_denom = getf(ref, 'denom', NaN);
    pt.medium_margin = getf(ref, 'margin', NaN);
```

The actual endpoint stability classification must use the frozen `crit_tol`, `D_tol`, and `Dq_tol`
inside `res.stability`; do not recreate it with raw `>0` checks. Set all provenance/two-tier fields
on **every** return path, including `early_return`, so callers never probe a missing member.

- [ ] **Step 5: Run the new test and the ordered regressions, then the full suite**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests({'invz_projected/tests/test_invz_solve_point_strict.m','invz_projected/tests/test_invz_ordered_residual_strict.m','invz_projected/tests/test_invz_ordered_jensen.m','invz_projected/tests/test_invz_qcp_closure.m'}); disp(table(r))"`
Expected: the new file 5 passed; the strict-residual file remains green; the two existing Jensen
files pass unchanged. Full suite: `F=0`.

- [ ] **Step 6: Commit**

```bash
git add invz_projected/invz_solve_point.m invz_projected/invz_solve_point_ordered.m invz_projected/tests/test_invz_solve_point_strict.m
git commit -m "feat(invzp): thread static-medium scheme and moments through both point solvers"
```

---

### Task 15: `invz_spectra_map` — three-way dispatcher (strict-gated) and per-column reasons

**Gated on a strict scheme being active.** Under `'resummed'` the historical dispatch is preserved exactly: the PM probe below `Bc_1z` fails for the same pole reason, so an ungated three-way rule would classify it *unknown* and mask **more** than today — and would break G9.

Also adopts the shared classifier at this file's two catches and at `invz_solve_auto`'s (spec §5.1: narrowing one site only relocates the swallow).

**Files:**
- Modify: `invz_projected/invz_spectra_map.m` (`:280` init, `:288-333` dispatch, `:297-300` and `:325-328` catches, the per-column output assembly, `Bc_1z` reduction)
- Modify: `invz_projected/invz_solve_auto.m` (its `invz:*` catch)
- Modify: `invz_projected/invz_check_solve_opts.m` (reserve driver-owned
  `static_medium`/`ref_margin` inside `solve_opts`)
- Test: `invz_projected/tests/test_invz_spectra_map_phase_reasons.m` (new)

**Interfaces:**
- Consumes: Tasks 4, 5, 13, 14.
- Produces: `S.static_medium` (requested scalar scheme); per-column `S.static_medium_used`,
  `S.crit_pm`, `S.pm_probe_status`, `S.pm_probe_error_id`, `S.stability_1z`,
  `S.phase_1z_reason`; `S.Bc_1z_interval`, `S.Bc_1z_status`. `phase_1z` keeps its
  `{0,1,2}` enum. A bare escape-hatch column uses `'n/a_bare_escape'` in
  `S.static_medium_used` rather than corrupting the scalar requested provenance.

- [ ] **Step 1: Locate the per-column output assembly**

Run: `rg -n "phase_1z|Bc_1z|one_field|Sigma0" invz_projected/invz_spectra_map.m`
Record the line numbers where `one_field`'s outputs are collected into `S` and where `Bc_1z` is reduced from the `phase_1z` column. Every insertion below goes at those points. Do not restructure `one_field` — the un-nesting is deferred (prereg §7.4).

- [ ] **Step 2: Write the failing test**

Create `invz_projected/tests/test_invz_spectra_map_phase_reasons.m`:

```matlab
function tests = test_invz_spectra_map_phase_reasons
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function o = base_opts()
o = struct('grid', [4 4 4], 'dpRng', 12, 'cache', false, 'verbose', false);
end
function w = wgrid()
w = (0.02:0.04:0.42).';
end

% Provenance is scalar and mandatory.
function test_scheme_provenance_present(testCase)
if isempty(getenv('INVZ_SLOW')), assumeFail(testCase, 'set INVZ_SLOW=1'); end
ion = invz_ion();
S = invz_spectra_map(ion, 0.31, 8, wgrid(), base_opts());
verifyEqual(testCase, S.static_medium, 'resummed');
end

% Every masked column carries a reason. This is the honesty gate: no phase_1z = 0 without one.
function test_no_masked_column_without_a_reason(testCase)
if isempty(getenv('INVZ_SLOW')), assumeFail(testCase, 'set INVZ_SLOW=1'); end
ion = invz_ion();
S = invz_spectra_map(ion, 0.31, [0 3 8], wgrid(), base_opts());
verifyEqual(testCase, numel(S.phase_1z_reason), numel(S.phase_1z));
for k = 1:numel(S.phase_1z)
    verifyFalse(testCase, isempty(S.phase_1z_reason{k}), sprintf('column %d has no reason', k));
    if S.phase_1z(k) == 0
        verifyTrue(testCase, any(strcmp(S.phase_1z_reason{k}, ...
            {'unstable_endpoint', 'medium_out_of_domain', 'degenerate_doublet', ...
             'solver_failed', 'pm_probe_unknown', 'boundary_indeterminate', ...
             'not_attempted_longitudinal', 'bare_not_ordered'})), S.phase_1z_reason{k});
    end
end
end

% B = 0 is labelled, never a silent proxy and never "solver_failed".
function test_zero_field_column_is_labelled_degenerate(testCase)
if isempty(getenv('INVZ_SLOW')), assumeFail(testCase, 'set INVZ_SLOW=1'); end
ion = invz_ion();
o = base_opts();  o.ordered_1z = 'jensen';
S = invz_spectra_map(ion, 0.31, [0 8], wgrid(), o);
verifyTrue(testCase, any(strcmp(S.phase_1z_reason{1}, ...
    {'degenerate_doublet', 'bare_not_ordered', 'pm'})));
verifyNotEqual(testCase, S.phase_1z_reason{1}, 'solver_failed');
end

% The three-way rule is GATED: under 'resummed' the historical dispatch is byte-preserved.
function test_dispatcher_is_gated_on_strict(testCase)
if isempty(getenv('INVZ_SLOW')), assumeFail(testCase, 'set INVZ_SLOW=1'); end
ion = invz_ion();
a = invz_spectra_map(ion, 0.31, [0 3 8], wgrid(), base_opts());
b = invz_spectra_map(ion, 0.31, [0 3 8], wgrid(), ...
    setfield(base_opts(), 'static_medium', 'resummed'));  %#ok<SFLD>
verifyEqual(testCase, b.phase_1z, a.phase_1z);
verifyEqual(testCase, b.Sigma0, a.Sigma0, 'AbsTol', 0);
end

% An unknown PM probe must never vote 'ordered' under strict.
function test_unknown_pm_probe_cannot_label_ordered(testCase)
S = struct('crit_pm', NaN, 'pm_probe_status', {{'nonconverged'}});
verdict = invz_pm_verdict(NaN, false, 1e-6);            % helper under test
verifyEqual(testCase, verdict, 'unknown');
verifyEqual(testCase, invz_pm_verdict(1e-3, true, 1e-6), 'pm');
verifyEqual(testCase, invz_pm_verdict(-1e-3, true, 1e-6), 'ordered_eligible');
verifyEqual(testCase, invz_pm_verdict(1e-9, true, 1e-6), 'unknown');   % inside the band
end
```

- [ ] **Step 3: Add the verdict helper**

Create `invz_common/invz_pm_verdict.m` — a pure three-way classifier, so the rule is testable without running a map:

```matlab
function v = invz_pm_verdict(crit_pm, pm_ok, crit_band)
%INVZ_PM_VERDICT Three-way PM-probe verdict for the 1/z phase dispatcher (spec SS4.4).
%   'pm'                converged finite PM with crit_pm >  crit_band
%   'ordered_eligible'  converged finite PM with crit_pm < -crit_band
%   'unknown'           PM non-convergence / non-finite / recoverable error, OR
%                       |crit_pm| <= crit_band (boundary-indeterminate)
%
% A FAILED PM PROBE IS NOT EVIDENCE FOR ORDER. The historical dispatch ran the ordered solver
% whenever the PM probe was not valid and then labelled a converging jensen result 'ordered' --
% which makes solver availability the phase criterion. 'unknown' may run the ordered solver for
% DIAGNOSTICS but must never emit phase_1z = 1 without a separately validated free-energy /
% branch-selection rule.
% crit_band is the frozen crit_tol. Field-resolution uncertainty is represented separately
% by S.Bc_1z_interval; a parfor point never owns dcrit/dB.
if ~pm_ok || ~isfinite(crit_pm)
    v = 'unknown';
elseif crit_pm > crit_band
    v = 'pm';
elseif crit_pm < -crit_band
    v = 'ordered_eligible';
else
    v = 'unknown';
end
end
```

- [ ] **Step 4: Resolve the public scheme, then wire the gated dispatcher**

At the driver entry, after validating `solve_opts`:

```matlab
sm = invz_check_static_medium(opts);
% static_medium/ref_margin are driver-owned at this level.
% Extend invz_check_solve_opts so either field inside solve_opts throws invz:solveOpts.
sopts.static_medium = sm.scheme;
sopts.ref_margin = sm.ref_margin;
```

Pass `sm` (or an explicit `is_strict` scalar and canonical scheme) into `one_field`; the original
sketch referenced an out-of-scope `sm` variable. Keep the new per-column outputs sliced for `parfor`.
The implementation boundary is:

```matlab
phase_1z = 0;  crit_pm = NaN;  m_1z = NaN;  D_ord = NaN;
pm_probe_status = 'not_attempted';  pm_probe_error_id = '';
stability_1z = false;  phase_1z_reason = 'solver_failed';
static_medium_used = sm.scheme;
```

Use a top-level strict gate inside `one_field` so the complete existing resummed body remains in its
original order and arithmetic. Do not sprinkle verdict conditions through the legacy body:

```matlab
if sm.is_strict
    % new strict-only control flow described below
    ...
    return;
end
% EXISTING RESUMMED BODY, unchanged
```

The strict-only control flow first calls `invz_solve_auto` for the independent auto/RPA overlay.
It then handles the transverse 1/z leg once for **both successful auto states**:

1. If `phase==2` and `pt` is a converged strict PM point, reuse it as `ptp`; otherwise call
   `invz_solve_point`. Thus the PM probe is neither skipped on an auto-PM state nor solved twice.
2. Set `crit_pm`, `pm_probe_status`, and
   `verdict=invz_pm_verdict(crit_pm,pm_consistent,crit_tol)`.
3. `verdict='pm'`: set `phase_1z=2`, reason `'pm'`, and evaluate `chiz` from `ptp`.
4. `verdict='ordered_eligible'`: call `invz_solve_point_ordered` with
   `ordered_mode='jensen'`. Emit `phase_1z=1` only when its consistency tier and
   `stable_1z` both pass; otherwise keep the column masked with the exact
   `unstable_endpoint` / `medium_out_of_domain` / `degenerate_doublet` / `solver_failed`
   reason.
5. `verdict='unknown'`: the Jensen call is diagnostic only. It can refine the failure reason but
   can never emit `phase_1z=1`; preserve `'boundary_indeterminate'` for a finite PM point inside
   the band and `'pm_probe_unknown'` for non-convergence/non-finite/recoverable PM outcomes.

This fixes the previous draft's hole where `phase==2` plus a negative PM mass never reached Jensen.
The explicitly deferred robustness item is narrower: when `invz_solve_auto` itself returns
`phase==0`, this stage keeps the strict column masked as `'solver_failed'` rather than fully un-nesting
the 1/z leg from auto availability. Stage-4 G16 must expose that limitation; it must not be described
as a successful strict solve.

Compute `chirpa` from the auto state with the same branch-specific arithmetic as today, but never let
its success/failure vote on `phase_1z`. For longitudinal fields and `ordered_1z='bare'`, keep the
documented escape hatch, set reason `'bare_escape_hatch'`, and set
`static_medium_used='n/a_bare_escape'`.

Add the local helper at the end of the file:

```matlab
function s = local_pm_status(ptp, crit_pm, crit_tol)
%LOCAL_PM_STATUS Narrow PM-probe vocabulary (spec SS4.4), so a masked column's cause is
% recorded rather than inferred. Fatal/unclassified errors never reach here: they rethrow.
if isempty(ptp),                          s = 'not_attempted';
elseif ~isfield(ptp, 'converged') || ~ptp.converged, s = 'nonconverged';
elseif ~isfinite(crit_pm),                s = 'nonfinite';
elseif abs(crit_pm) <= crit_tol,           s = 'boundary_band';
elseif crit_pm > crit_tol,                 s = 'stable';
else,                                      s = 'unstable';
end
end
```

Every recoverable catch must store the exact identifier in `pm_probe_error_id` and set
`pm_probe_status='recoverable_error'`; a non-empty error id must never be collapsed into
`'nonconverged'`.

- [ ] **Step 5: Route every outer solver/response call through the shared boundary**

Search both files with `rg -n "try|catch|strncmp"`. Replace every outer one-output solver/response
catch reached by strict mode with `invz_try_solver_call(@() ...)`. In `invz_solve_auto`, use the
returned `completed/error_id` to populate `di.ordered_err` or `di.para_err`; in
`invz_spectra_map`, use it to populate the exact per-column error id/status. Do not leave a second
hand-written classifier branch beside the helper. The resummed numerical body is unchanged; only the
old error-prefix policy is narrowed.

After editing, `rg -n "strncmp\\(.*invz:|catch" invz_projected/invz_solve_auto.m
invz_projected/invz_spectra_map.m` must show no broad `invz:*` absorber. Any remaining `catch` must
be justified as non-solver cleanup and must rethrow by default.

- [ ] **Step 6: Collect the new columns and widen the `Bc_1z` reduction**

At the assembly point found in Step 1, add `S.static_medium = sm.scheme;` plus all per-column arrays.
Do not use the earlier `min/max(all unknown fields)` reduction: it referenced nonexistent `S.B` and
could widen a boundary with an unrelated unknown column on the far side of the sweep. Add a pure
`invz_boundary_interval(fields, orderedMask, pmMask, unknownMask)` helper that locates the last
ordered and first stable-PM anchors, includes only intervening unknown/boundary columns, and returns
`valid | widened | unbracketed | invalid`. Require finite strictly increasing `fields`; callers that
want another presentation order must reorder before this reduction rather than make “last/first”
depend on input order.

```matlab
unk = strcmp(S.phase_1z_reason, 'pm_probe_unknown') | ...
      strcmp(S.phase_1z_reason, 'boundary_indeterminate');
ord = strcmp(S.phase_1z_reason, 'ordered');       % excludes bare_escape_hatch
pm  = strcmp(S.phase_1z_reason, 'pm');
[S.Bc_1z_interval, S.Bc_1z_status, S.Bc_1z] = invz_boundary_interval( ...
    S.fields, ord, pm, unk);
```

Run that replacement only under `sm.is_strict`. Under `'resummed'`, preserve the historical scalar
`S.Bc_1z` reduction exactly for G9; populate the new interval/status fields diagnostically without
feeding them back into any pre-existing field.

Unit-test valid, widened, unknown-outside-bracket, unbracketed, and all-unknown cases without running
a spectrum. Emit at most one sweep-level summary warning with counts of domain/degenerate/unknown
columns; never warn per node.

- [ ] **Step 7: Assert the coupling cache is untouched (spec §4.5)**

Run, before and after a map call:
```
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); \
n0=numel(dir('invz_projected/cache/jq5_*.mat')); ion=invz_ion(); \
invz_spectra_map(ion, 0.31, 8, (0.02:0.04:0.42).', \
struct('grid',[4 4 4],'dpRng',12,'cache',false,'verbose',false, \
'static_medium','strict_1z_dyson_ref')); \
n1=numel(dir('invz_projected/cache/jq5_*.mat')); fprintf('cache %d -> %d\n', n0, n1)"
```
Expected: unchanged count. Moments are derived at call time, so no cache key, prefix or schema changes.

- [ ] **Step 8: Run the new test and the spectra regressions, then the full suite**

Run: `INVZ_SLOW=1 "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests({'invz_projected/tests/test_invz_spectra_map_phase_reasons.m','invz_projected/tests/test_invz_spectra_map.m','invz_projected/tests/test_invz_ewald_integration.m'}); disp(table(r))"`
Expected: new file 5 passed; both existing files pass unchanged (`test_invz_ewald_integration` is the Gate-C7 provenance matrix and is the direct guard that `S`'s schema growth did not disturb cache identity).
Full suite: `F=0`.

- [ ] **Step 9: Commit**

```bash
git add invz_common/invz_pm_verdict.m invz_common/invz_boundary_interval.m \
  invz_projected/invz_check_solve_opts.m invz_projected/invz_spectra_map.m \
  invz_projected/invz_solve_auto.m invz_projected/tests/test_invz_spectra_map_phase_reasons.m
git commit -m "feat(invzp): strict-gated three-way phase dispatch with per-column reasons"
```

---

### Task 16: `invz_static_domain_scan` — the prospective half of Gate 0

Bare diagonalizations only: no solve, no iteration. Accepts an **explicit** `hgrid` and never reproduces the adaptive extension/redensification (prereg §7.2).

**Files:**
- Create: `invz_projected/invz_static_domain_scan.m`
- Test: `invz_projected/tests/test_invz_static_domain_scan.m`

**Interfaces:**
- Consumes: Tasks 1, 2, 3, 11, 12 (`invz_hmf_grid`).
- Produces: `scan = invz_static_domain_scan(ion, T, Bx, Jnu_flat, opts)` with fields
  `predictor` (the separately required `h=0` record), `hgrid`, `Delta`, `valid`, `G0bare`,
  `Gref`, `ref_denom`, `ref_margin`, `ref_status` (cell), `omit_mu3`, `omit_cubic`,
  `omit_max`, `n_nodes`, `n_required`, `n_valid`, `n_skipped`, `n_out_of_domain`,
  `n_degenerate`, `scheme`, `grid_source`. Used by Task 17.

- [ ] **Step 1: Write the failing test**

Create `invz_projected/tests/test_invz_static_domain_scan.m`:

```matlab
function tests = test_invz_static_domain_scan
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function [ion, Jnu, o] = fx()
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'static_medium', 'strict_1z_dyson_ref');
end

function test_explicit_hgrid_is_honoured(testCase)
[ion, Jnu, o] = fx();
o.hgrid = [1e-4 1e-3 1e-2];
scan = invz_static_domain_scan(ion, 0.31, [2.85 0 0], Jnu, o);
verifyEqual(testCase, scan.hgrid, o.hgrid, 'AbsTol', 0);
verifyEqual(testCase, scan.grid_source, 'explicit');
verifyEqual(testCase, scan.n_nodes, 3);
verifyEqual(testCase, scan.n_required, 4);     % explicit profile nodes plus h=0 predictor
end

function test_default_grid_comes_from_the_shared_helper(testCase)
[ion, Jnu, o] = fx();
scan = invz_static_domain_scan(ion, 0.31, [2.85 0 0], Jnu, o);
verifyEqual(testCase, scan.grid_source, 'invz_hmf_grid');
verifyEqual(testCase, numel(scan.hgrid), 33);
verifyTrue(testCase, all(diff(scan.hgrid) > 0));
end

% Every array is aligned and every counter accounted for -- no silent caps (spec §5.5).
function test_arrays_aligned_and_counters_complete(testCase)
[ion, Jnu, o] = fx();
scan = invz_static_domain_scan(ion, 0.31, [2.85 0 0], Jnu, o);
n = scan.n_nodes;
for f = {'Delta', 'valid', 'G0bare', 'Gref', 'ref_denom', 'ref_margin', ...
         'omit_mu3', 'omit_cubic', 'omit_max'}
    verifyEqual(testCase, numel(scan.(f{1})), n, f{1});
end
verifyEqual(testCase, numel(scan.ref_status), n);
verifyEqual(testCase, scan.n_valid + scan.n_degenerate + ...
            scan.n_out_of_domain + scan.n_skipped, scan.n_required, ...
            'every predictor/profile node must be accounted for exactly once');
end

% Small transverse field: the lowest nodes are where Delta dips below the floor.
function test_degenerate_nodes_are_flagged_at_small_field(testCase)
[ion, Jnu, o] = fx();
scan = invz_static_domain_scan(ion, 0.31, [0 0 0], Jnu, o);
verifyGreaterThan(testCase, scan.n_degenerate, 0);
verifyEqual(testCase, scan.predictor.status, 'degenerate_doublet');
end

function test_resummed_scheme_is_rejected(testCase)
[ion, Jnu, o] = fx();
o.static_medium = 'resummed';
verifyError(testCase, @() invz_static_domain_scan(ion, 0.31, [2.85 0 0], Jnu, o), ...
            'invz:staticMedium');
end
```

- [ ] **Step 2: Run it to verify it fails, then write the scanner**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_static_domain_scan.m'); disp(table(r))"`
Expected: all 5 error on the missing function.

Create `invz_projected/invz_static_domain_scan.m`:

```matlab
function scan = invz_static_domain_scan(ion, T, Bx, Jnu_flat, opts)
%INVZ_STATIC_DOMAIN_SCAN Prospective (stage-a) half of Gate 0: does the strict-order static
% medium have a valid reference along the H_MF profile grid, and how large are the omitted
% moment terms there? G = -chi, ferromagnetic positive J.
%
% BARE DIAGONALIZATIONS ONLY -- no Sigma<->medium solve, no iteration, no convergence
% dependence. Sigma0 is taken as 0 (the kinematic proxy): Sigma0 is O(1/z), so this bounds the
% reference denominator 1 + Sigma0 to leading order. The SOLVED-path margins are the stage-b
% half of Gate 0 and are read off prof.medium_status / prof.Delta after the real solve -- this
% function does NOT predict them (spec SS7.2).
%
% GRID OWNERSHIP (prereg SS7.2): opts.hgrid is honoured verbatim when supplied. Otherwise the
% INITIAL geometric grid is built from the shared invz_hmf_grid helper with the same hmax rule
% invz_hmf_ordered uses (hmax_fac * |bare ordered hz|, or the exact opts.hmax_abs override).
% This function deliberately does NOT reproduce invz_hmf_ordered's adaptive extension or
% redensification: two implementations that agree in one test are not grid identity.
%
% opts: J0eff (required), Jxx0, transverse_mf, hyp, nH, hmax_fac, hmax_abs, hmin_frac, hgrid,
%       static_medium (must be a strict scheme), ref_margin, Jmom, Ecut.
% scan: hgrid, Delta, valid, G0bare, Gref, ref_denom, ref_status{}, omit_mu3, omit_cubic,
%       omit_max, predictor, n_nodes, n_required, n_valid, n_skipped, n_out_of_domain,
%       n_degenerate, scheme, grid_source.
% The h=0 predictor is required by Jensen but is not part of the nonzero geometric hgrid, so
% it is evaluated and accounted separately. Counters satisfy
% n_valid+n_degenerate+n_out_of_domain+n_skipped == n_required.
if nargin < 5, opts = struct(); end
J0eff  = opts.J0eff;
Jxx0   = getf(opts, 'Jxx0', ion.Jxx0);
tmf    = getf(opts, 'transverse_mf', 'legacy_x');
hyp    = getf(opts, 'hyp', true);
nH     = getf(opts, 'nH', 33);
hfrac  = getf(opts, 'hmin_frac', 1e-3);
Ecut   = getf(opts, 'Ecut', 40);
scheme = getf(opts, 'static_medium', '');
if ~any(strcmp(scheme, {'strict_1z_dyson_ref', 'strict_1z_bare_ref'}))
    error('invz:staticMedium', ['invz_static_domain_scan requires a strict scheme ' ...
        '(''strict_1z_dyson_ref'' or ''strict_1z_bare_ref''); got ''%s''. The resummed medium ' ...
        'has no reference denominator to scan.'], scheme);
end
refopt = struct('ref_margin', getf(opts, 'ref_margin', 1e-6));
Jmom = getf(opts, 'Jmom', []);
if isempty(Jmom), Jmom = invz_coupling_moments(Jnu_flat); end
if ~isvector(Jnu_flat)
    error('invz:staticMedium', 'ordered-domain scan requires a static vector Jnu_flat.');
end

sibase = struct('hyp', hyp, 'Jxx0', Jxx0, 'transverse_mf', tmf);
if isfield(opts, 'hgrid') && ~isempty(opts.hgrid)
    hgrid = opts.hgrid(:).';  grid_source = 'explicit';
else
    sibo = sibase;  sibo.order = true;  sibo.J0z = J0eff;
    sib = invz_single_ion(ion, T, Bx, sibo);
    if isfield(opts, 'hmax_abs'), hmax = opts.hmax_abs;
    else,                          hmax = getf(opts, 'hmax_fac', 1.25) * abs(sib.hz);  end
    if ~(isfinite(hmax) && hmax > 0)
        error('invz:hmfGrid', ['no positive bracket ceiling: the bare set does not order at ' ...
            'T = %g, |Bx| = %g. Supply opts.hgrid to scan anyway.'], T, norm(Bx));
    end
    hgrid = invz_hmf_grid(hmax, nH, hfrac);  grid_source = 'invz_hmf_grid';
end

n = numel(hgrid);
empty = struct('status','not_evaluated','Delta',NaN,'G0bare',NaN,'Gref',NaN, ...
    'ref_denom',NaN,'ref_margin',NaN,'omit_mu3',NaN,'omit_cubic',NaN,'omit_max',NaN);
nodes = repmat(empty, 1, n);
for k = 1:n
    nodes(k) = eval_proxy(hgrid(k));
end
predictor = eval_proxy(0);                 % the same routine, not a second implementation

node_status = {nodes.status};
statuses = [{predictor.status}, node_status];
known = strcmp(statuses, 'ok') | strcmp(statuses, 'degenerate_doublet') | ...
        strcmp(statuses, 'ref_denom_nonpositive') | strcmp(statuses, 'ref_denom_small') | ...
        strcmp(statuses, 'nonfinite');
n_valid         = nnz(strcmp(statuses, 'ok'));
n_degenerate    = nnz(strcmp(statuses, 'degenerate_doublet'));
n_skipped       = nnz(strcmp(statuses, 'not_evaluated'));
n_out_of_domain = nnz(known) - n_valid - n_degenerate;
n_required      = n + 1;
if nnz(known) + n_skipped ~= n_required
    error('invz:staticDomainScan', 'unclassified prospective node status.');
end

scan = struct('hgrid', hgrid, 'Delta', [nodes.Delta], ...
    'valid', strcmp({nodes.status}, 'ok'), 'G0bare', [nodes.G0bare], ...
    'Gref', [nodes.Gref], 'ref_denom', [nodes.ref_denom], ...
    'ref_margin', [nodes.ref_margin], 'ref_status', {node_status}, ...
    'predictor', predictor, 'omit_mu3', [nodes.omit_mu3], ...
    'omit_cubic', [nodes.omit_cubic], 'omit_max', [nodes.omit_max], ...
    'n_nodes', n, 'n_required', n_required, 'n_valid', n_valid, ...
    'n_skipped', n_skipped, 'n_out_of_domain', n_out_of_domain, ...
    'n_degenerate', n_degenerate, 'scheme', scheme, 'grid_source', grid_source);

    function rec = eval_proxy(hp)
    % One bare, non-iterative kinematic proxy. Every return preserves the fixed schema.
    rec = empty;
    sio = sibase;  sio.hz_fixed = hp;
    si = invz_single_ion(ion, T, Bx, sio);
    tl = invz_twolevel_ordered(ion, T, Bx, hp, ...
        struct('Jxx0', Jxx0, 'transverse_mf', tmf, 'domain_policy', 'return'));
    rec.Delta = tl.Delta;
    if ~tl.valid
        rec.status = 'degenerate_doublet';
        return;
    end
    wn0 = invz_matsubara(T, Ecut);
    c0 = invz_chi0z(si, T, 1i*wn0(1), struct('elastic', true));
    X = real(c0(:, :, 1));
    switch tmf
        case 'none',      fb = 0;
        case 'legacy_x',  fb = X(3,1) * (Jxx0/(1-Jxx0*X(1,1))) * X(1,3);
        case 'vector_ab', t = [1 2];
                          fb = X(3,t) * (Jxx0*((eye(2)-Jxx0*X(t,t))\X(t,3)));
        otherwise, error('invz:transverseMF', 'unknown transverse_mf ''%s''', tmf);
    end
    rec.G0bare = -(X(3,3) + fb);
    [rec.Gref, ref] = invz_static_medium_reference(rec.G0bare, 0, scheme, refopt);
    rec.ref_denom = ref.denom;  rec.ref_margin = ref.margin;  rec.status = ref.status;
    if ~strcmp(ref.status, 'ok'), return; end
    [~, clo] = invz_medium_moment_closure(rec.Gref, Jmom, scheme);
    rec.omit_mu3 = clo.omit_mu3;  rec.omit_cubic = clo.omit_cubic;
    rec.omit_max = clo.omit_max;  rec.status = clo.status;
    end
end
```

- [ ] **Step 3: Run the test, then the full suite**

Run: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests('invz_projected/tests/test_invz_static_domain_scan.m'); disp(table(r))"`
Expected: 5 passed.
Full suite: `F=0`.

- [ ] **Step 4: Commit**

```bash
git add invz_projected/invz_static_domain_scan.m invz_projected/tests/test_invz_static_domain_scan.m
git commit -m "feat(invzp): add prospective strict-medium domain scanner (Gate 0 stage a)"
```

---

### Task 17: build-blocking gate tests — G1, G2, G3, G11, G13, G15

Two test-only files. Nothing in production changes, so a reviewer can reject these without touching stage 2.

**Files:**
- Test: `invz_projected/tests/test_invz_strict_identities.m` (G1, G2, G3)
- Test: `invz_projected/tests/test_invz_strict_contracts.m` (G11, G13, G15)

**Interfaces:**
- Consumes: everything from Tasks 1–16, plus `invz_exact_numeric_digest` (Task 0). Produces: no
  production interface.

- [ ] **Step 1: Write the identity gates (G1, G2, G3)**

Create `invz_projected/tests/test_invz_strict_identities.m`:

```matlab
function tests = test_invz_strict_identities
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function [ion, T, Bx, Jnu, o] = fx()
ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
           'static_medium', 'strict_1z_dyson_ref');
end

% G1a: dm/dh = -G0bare (J 2.31), panelwise on the profile.
function test_G1a_dm_dh_equals_minus_G0bare(testCase)
[ion, T, Bx, Jnu, o] = fx();
o.nH = 129;
[~, p] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyEqual(testCase, p.status, 'ok');
dm = diff(p.m)./diff(p.hgrid);
gb = 0.5*(p.G0bare(1:end-1) + p.G0bare(2:end));
ok = isfinite(dm) & isfinite(gb);
verifyLessThanOrEqual(testCase, max(abs(dm(ok) + gb(ok))./max(1, abs(gb(ok)))), 1e-6, ...
    'dm/dh = -G0bare must hold to panel order');
end

% G1b: Delta F/Delta h against the trapezoidal average of crit -- i.e. F' = crit (spec §A).
function test_G1b_F_prime_equals_crit(testCase)
[ion, T, Bx, Jnu, o] = fx();
o.nH = 129;
[~, p] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyEqual(testCase, p.status, 'ok');
dF = diff(p.F)./diff(p.hgrid);
cm = 0.5*(p.crit(1:end-1) + p.crit(2:end));
ok = isfinite(dF) & isfinite(cm);
verifyLessThanOrEqual(testCase, max(abs(dF(ok) - cm(ok))./max(1, abs(cm(ok)))), 1e-6);
end

% G1c: dF/dm = crit/chi_path, chi_path = -G0bare.
function test_G1c_dF_dm_equals_crit_over_chi_path(testCase)
[ion, T, Bx, Jnu, o] = fx();
o.nH = 129;
[~, p] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyEqual(testCase, p.status, 'ok');
dFdm = diff(p.F)./diff(p.m);
cm   = 0.5*(p.crit(1:end-1) + p.crit(2:end));
chp  = -0.5*(p.G0bare(1:end-1) + p.G0bare(2:end));
ok = isfinite(dFdm) & isfinite(cm) & chp > 0;
verifyLessThanOrEqual(testCase, ...
    max(abs(dFdm(ok) - cm(ok)./chp(ok))./max(1, abs(cm(ok)./chp(ok)))), 1e-6);
end

% G1d: second-order convergence under nH refinement (prereg §5).
function test_G1d_second_order_convergence(testCase)
[ion, T, Bx, Jnu, o] = fx();
e = nan(1, 3);  nHs = [33 65 129];
for k = 1:3
    ok_ = o;  ok_.nH = nHs(k);
    [~, p] = invz_hmf_ordered(ion, T, Bx, Jnu, ok_);
    verifyEqual(testCase, p.status, 'ok');
    dF = diff(p.F)./diff(p.hgrid);
    cm = 0.5*(p.crit(1:end-1) + p.crit(2:end));
    m_ = isfinite(dF) & isfinite(cm);
    e(k) = max(abs(dF(m_) - cm(m_)));
end
verifyGreaterThanOrEqual(testCase, e(1), e(2));
verifyGreaterThanOrEqual(testCase, e(2), e(3));
verifyLessThanOrEqual(testCase, e(3), 1e-6);
% A separate pure nested-uniform-grid test pins second-order trapezoid convergence. The
% adaptive production grid is not required to show a [3,5] ratio.
end

% G1d companion: second-order is pinned on a genuinely nested smooth quadrature fixture,
% separately from the non-nested adaptive production grids.
function test_G1d_nested_trapezoid_is_second_order(testCase)
n = [17 33 65];  e = nan(size(n));
Iexact = exp(1) - 1;
for k = 1:numel(n)
    x = linspace(0, 1, n(k));
    e(k) = abs(trapz(x, exp(x)) - Iexact);
end
ratio = e(1:end-1)./e(2:end);
verifyGreaterThan(testCase, min(ratio), 3.5);
verifyLessThan(testCase, max(ratio), 4.5);
end

% G2: the two legs coincide at m -> 0, WITHIN the frozen K tolerances (prereg §7.3 -- not
% bitwise: the callers reach Gref through different expressions).
function test_G2_onset_coincidence_at_m_zero(testCase)
[ion, T, Bx, Jnu, o] = fx();
Jscale = max(abs(Jnu));  K_atol = 1e-14;  K_rtol = 1e-12;
ptp = invz_solve_point(ion, T, Bx, Jnu, o);                     % PM leg, m = 0
[~, p] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyTrue(testCase, isfinite(p.slope0) && isfinite(ptp.crit));
% the dimensionless masses agree at onset
verifyEqual(testCase, p.slope0, ptp.crit, 'AbsTol', 1e-6);
% and so do the references and the media
G0pm = -ptp.chi0cc0;
GrefP = invz_static_medium_reference(G0pm, ptp.Sigma0, 'strict_1z_dyson_ref');
GrefO = invz_static_medium_reference(p.G0bare_pm0, p.Sigma0_pm0, 'strict_1z_dyson_ref');
verifyEqual(testCase, GrefO, GrefP, 'RelTol', 1e-6);
mom = invz_coupling_moments(Jnu);
KP = invz_medium_moment_closure(GrefP, mom, 'strict_1z_dyson_ref');
KO = invz_medium_moment_closure(GrefO, mom, 'strict_1z_dyson_ref');
Kgate = K_atol + K_rtol*max([abs(KO), abs(KP), Jscale]);
verifyLessThanOrEqual(testCase, abs(KO - KP), Kgate);
verifyLessThanOrEqual(testCase, abs(p.K0_pm0-ptp.K(1)), Kgate);
end

% G3: the pinned m = 0 identity r = 1 + Sigma0 survives under strict, for ANY K0.
function test_G3_r_equals_one_plus_sigma_under_strict(testCase)
beta = 1/(0.0862*0.31);
tl = struct('Delta', 0.02, 'M2', 0.8, 'm', 0, 'n01', tanh(0.02*beta/2), 'g0', 1);
tl.g0 = 2*tl.n01/tl.Delta;
for K0 = [0, 1e-3, 5e-3, 0.05]
    [~, out] = invz_gstat_ordered(tl, [0.01; 0.02], K0, 0.25, beta, -300, 0);
    verifyEqual(testCase, out.r, 1.25, 'RelTol', 1e-12, sprintf('K0 = %g', K0));
end
end
```

- [ ] **Step 2: Write the contract gates (G11, G13, G15)**

Create `invz_projected/tests/test_invz_strict_contracts.m`:

Before saving the test, replace `<FROZEN_HASH_FROM_TASK_0>` with the exact digest recorded in the
approved preregistration. Treat a remaining angle-bracket token as a build error, not as an expected
first-run failure.

```matlab
function tests = test_invz_strict_contracts
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% G11: a REAL-coupling ordered anchor, with its full provenance asserted. The original masking
% defect survived a whole stage because no test fed the jensen leg real bz_couplings densities.
function test_G11_real_coupling_ordered_anchor(testCase)
if isempty(getenv('INVZ_SLOW')), assumeFail(testCase, 'set INVZ_SLOW=1'); end
ion = invz_ion();
prov = struct('grid', [16 16 16], 'dpRng', 30, 'dipole', 'bruteforce', 'cache', false);
[Jnu, info] = invz_bz_couplings(ion, prov);
mom = invz_coupling_moments(Jnu(:));
% provenance first: the hard-coded moments are valid ONLY for this tuple (spec §B)
verifyEqual(testCase, info.dipole.backend, 'bruteforce');
verifyFalse(testCase, isfield(info, 'grid'), ...
    'no grid-policy field means the bit-identical legacy route, whose info.grid is absent');
verifyEqual(testCase, invz_exact_numeric_digest(Jnu(:)), '<FROZEN_HASH_FROM_TASK_0>');
verifyEqual(testCase, mom.n, 16384);
verifyEqual(testCase, mom.Jbar, 1.20766e-4, 'RelTol', 1e-4);
verifyEqual(testCase, mom.mu2,  5.48264e-6, 'RelTol', 1e-4);
% Build-blocking real anchor: the diagnosed 1 T production point must now have a solved root.
o = struct('J0eff', info.Jcc0, 'Jxx0', ion.Jxx0, 'hyp', true, ...
           'static_medium', 'strict_1z_dyson_ref', 'Jmom', mom, 'trace', true);
[hstar, p, trc] = invz_hmf_ordered(ion, 0.1, [1 0 0], Jnu(:), o);
verifyEqual(testCase, p.status, 'ok');
verifyTrue(testCase, isfinite(hstar) && hstar > 0);
verifyTrue(testCase, all(strcmp(p.medium_status, 'ok')));
verifyTrue(testCase, trc.enabled && ~isempty(trc.nodes));
verifyTrue(testCase, all(strcmp({trc.nodes.medium_status}, 'ok')));
verifyTrue(testCase, isfinite(p.crit_star));
% Dq diagnostics still need the FULL multiset, not just the two moments
verifyEqual(testCase, numel(Jnu(:)), 16384);
end

% G13: behavioural sentinel -- the PM slot must not leak into ordered lambdas. Not a
% source-text-order assertion (those are brittle; a prior test-regex bug on this branch came
% from exactly that style).
function test_G13_pm_slot_does_not_leak_into_ordered_lambdas(testCase)
ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];  hz = 0.02;
Jnu = linspace(-2e-3, 6.0e-3, 24).';  J0eff = 6.42e-3;
[wn, wts, beta] = invz_matsubara(T, 40);
si = invz_single_ion(ion, T, Bx, struct('hyp', true, 'hz_fixed', hz, 'Jxx0', ion.Jxx0, ...
                                        'transverse_mf', 'legacy_x'));
tl = invz_twolevel_ordered(ion, T, Bx, hz, struct('Jxx0', ion.Jxx0, ...
                                                  'transverse_mf', 'legacy_x'));
c0  = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
G0  = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si, T, 1i*wn(1), struct('elastic', false));
G0inel0 = -real(c0i(3,3,1));
X  = real(c0(:, :, 1));
fb = X(3,1) * (ion.Jxx0/(1 - ion.Jxx0*X(1,1))) * X(1,3);
G0bare0 = -(X(3,3) + fb);
node = struct('tl', tl, 'G0', G0, 'g', real(invz_g(tl, 1i*wn)), 'wts', wts, 'wn', wn, ...
    'beta', beta, 'J0eff', J0eff, 'G0inel0', G0inel0, 'G0el0', G0bare0 - G0inel0, ...
    'G0bare0', G0bare0, 'eso', struct('warn', false, 'static_medium', 'strict_1z_dyson_ref'), ...
    'eopts', struct('static_medium', 'strict_1z_dyson_ref'), 'Jnu_flat', Jnu, ...
    'Jmom', invz_coupling_moments(Jnu));
[state, ~] = invz_ordered_node_solve(node, [], struct('trace', false));
% the exported K(1) must be the ORDERED static value, and the lambdas must be the ones derived
% from it -- recompute both and compare against invz_emt_scalar's own slot 1
med = invz_emt_scalar(node.G0, state.Sigma, node.Jnu_flat, node.eopts);
verifyNotEqual(testCase, state.K(1), med.K(1), ...
    'ordered K(1) must be the elastic-hybrid static value, not the PM slot');
lam_from_exported = invz_lambdas(state.K, node.g, node.wts, node.beta, [1 2 3]);
verifyEqual(testCase, state.lam, lam_from_exported, 'RelTol', 1e-10, ...
    'lambdas must derive from the exported K WITH the ordered K(1) substituted');
Kleak = state.K;  Kleak(1) = med.K(1);
lam_leaked = invz_lambdas(Kleak, node.g, node.wts, node.beta, [1 2 3]);
verifyGreaterThan(testCase, max(abs(lam_leaked - state.lam)), 1e-12, ...
    'the sentinel must be able to detect a leak at all');
end

% G15: fatal ids escape every layer; recoverable/domain outcomes keep their exact category.
function test_G15_fatal_ids_escape_every_layer(testCase)
verifyFalse(testCase, invz_is_recoverable_solver_error('invz:staticMedium'));
verifyTrue(testCase,  invz_is_recoverable_solver_error('invz:degenerateDoublet'));
[~, completed, rid] = invz_try_solver_call(@() error('invz:degenerateDoublet','synthetic'));
verifyFalse(testCase, completed);
verifyEqual(testCase, rid, 'invz:degenerateDoublet');
% node solver
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
           'static_medium', 'not_a_scheme');
verifyError(testCase, @() invz_hmf_ordered(ion, 0.31, [2.85 0 0], Jnu, o), ...
            'invz:staticMedium');
% point solver
verifyError(testCase, @() invz_solve_point(ion, 0.31, [2.85 0 0], Jnu, o), ...
            'invz:staticMedium');
% auto dispatcher: inject a fatal error INSIDE its ordered/PM try blocks, not at its own entry.
obad = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
              'static_medium', 'strict_1z_dyson_ref');
Jbad = complex(Jnu, 1e-12);
verifyError(testCase, @() invz_solve_auto(ion, 0.31, [2.85 0 0], Jbad, obad), ...
            'invz:couplingMoments');
% spectra outer boundary: the same injected fatal must cross one_field/parfor unchanged.
sp = struct('Jnu', Jbad, 'info', struct('Jcc0', 6.42e-3), ...
    'static_medium', 'strict_1z_dyson_ref', 'parallel', false, 'verbose', false);
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 2.85, (0.02:0.04:0.42).', sp), ...
            'invz:couplingMoments');
end

% G15b: a domain outcome must never be reported as generic node_failed.
function test_G15b_domain_outcome_keeps_its_category(testCase)
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
           'static_medium', 'strict_1z_dyson_ref', 'ref_margin', 1e9);
[~, p] = invz_hmf_ordered(ion, 0.31, [2.85 0 0], Jnu, o);
verifyEqual(testCase, p.status, 'medium_out_of_domain');
verifyNotEqual(testCase, p.status, 'node_failed');
end
```

- [ ] **Step 3: Run both files**

Run: `INVZ_SLOW=1 "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); r=runtests({'invz_projected/tests/test_invz_strict_identities.m','invz_projected/tests/test_invz_strict_contracts.m'}); disp(table(r))"`
Expected: every named test passes. **A G1 failure is a real signal, not a tolerance problem** — it
means `crit` as implemented is not `F'`, i.e. the construction in §A is not what the code computes.
Do not loosen the tolerance; report the measured error and its `nH` scaling.

- [ ] **Step 4: Run the full suite and commit**

Full suite: `F=0`; compare named and aggregate states with the Task-0 baseline.

```bash
git add invz_projected/tests/test_invz_strict_identities.m invz_projected/tests/test_invz_strict_contracts.m
git commit -m "test(invzp): add strict-medium identity and contract gates (G1-G3, G11, G13, G15)"
```

---

## Stage 3 — Gate 0 and the measurements

### Task 18: run Gate 0, G5 and G7, and write the verdict

**This is the decision point.** Gate 0 can fail, and its failure is a legitimate outcome: the run stops at diagnosis, and any revised theory starts a new spec and a fresh preregistration (prereg §3). Stage 4 is only worth planning if this task returns a pass.

**Files:**
- Create: `invz_projected/invz_gate0_report.m` (diagnostic driver, not a test)
- Create: `invz_projected/invz_pole_continuity.m` (pure actual-path G17 evaluator)
- Test: `invz_projected/tests/test_invz_pole_continuity.m`
- Create: `docs/invzp_strict_medium_gate0_report.md` (the verdict)

**Interfaces:**
- Consumes: Tasks 12, 14, 16 and the Task-0 frozen values.
- Produces:
  `rep = invz_gate0_report(ion,T,ordered_fields,pm_fields,Jnu_flat,opts)`, with
  `rep.ordered` (one row per field and `nH`), `rep.pm`, `rep.g5`, `rep.g17`, the exact
  coupling fingerprint/provenance, Boolean `fail_a` … `fail_e`, and `rep.pass`.

- [ ] **Step 1: Implement and test the pure actual-path pole check**

`invz_pole_continuity(h,d,y,tol)` takes aligned finite vectors, where
`d = prof.gstat_local_denom` and `y` is either `prof.r` or `prof.crit`. A crossing is an exact zero
or adjacent sign change of `d`. For every crossing:

1. locate `h_cross` by linear interpolation in `d`;
2. use the closest two finite nodes strictly on each side to form independent linear
   extrapolants of `y(h_cross)`;
3. return each relative jump
   `abs(yL-yR)/max(1,abs(yL),abs(yR))`;
4. return `status='ok'` only when every jump is finite and `<=tol`;
5. return `status='no_crossing'` when there is none, and `'unresolved'` when either side lacks
   two nodes. An exact-zero grid node is excluded from both fits — build `yL` and `yR` from the
   nearest two finite nodes strictly on each side — but it is then **checked against both of
   them**, not merely required to be finite: return `|y_node-yL|` and `|y_node-yR|` alongside
   `|yL-yR|`, each normalised by `max(1,|.|,|.|)` of its own pair, and require all three
   `<= tol`. That node is a direct evaluation *at* the crossing and so is the strongest evidence
   available; a finiteness-only rule would admit a finite but arbitrarily wrong value sitting
   between two mutually consistent extrapolants.

Unit tests cover a smooth removable crossing, an injected jump, no crossing, an exact-zero node
whose value agrees with both extrapolants, an exact-zero node whose value is finite but
displaced (must FAIL — this is the case a finiteness-only rule would have passed), and
insufficient side coverage. Do not use `interp1` across duplicate or exact-zero abscissae.

- [ ] **Step 2: Implement the driver against the complete solved-node ledger**

The driver must validate a strict scalar-vector configuration and must not contain a broad catch.
An unclassified/wiring error aborts the report. It must compute the coupling digest by calling the
shared `invz_exact_numeric_digest` primitive (Task 0), compare it against the preregistered digest
before any solve, and error on mismatch.
For every preregistered ordered field and each
`nH in [33 65 129]`:

1. Run `invz_static_domain_scan` on the shared initial grid. This is prospective metadata only;
   it never substitutes for a solved-path pass.
2. Run `[hstar,prof,trc] = invz_hmf_ordered(..., trace=true)`.
3. Treat **every** `trc.nodes` entry as required: predictor, initial sweep, extension,
   redensification, bisection, and root. Count exactly one of
   `ok`, `medium_out_of_domain`, `degenerate_doublet`, `solver_failed` for each ledger entry.
   `n_accounted == numel(trc.nodes)` is the coverage identity. An empty/unrecognised status,
   a missing trace phase, or a missing final root ledger entry fires (b).
4. Evaluate (a) from actual ledger `medium_status`, not from `scan`: every required strict node
   must be `'ok'`.
5. Evaluate (c) from `max(omit_max)` over the actual ledger. A missing/non-finite omitted-order
   ratio at a node otherwise labelled `ok` fails; it is not dropped by `isfinite` filtering.
6. Evaluate (d) using `invz_pole_continuity` on the final adapted `prof.hgrid` for both `r` and
   `crit`. Require finite integrands everywhere. If crossings exist, their maximum jump at
   `nH=129` must be within `pole_cont_tol` and must not increase relative to `nH=65`.
7. Evaluate (e): `prof.status='ok'`, finite `hstar>0`, and endpoint stability passes the frozen
   `crit_tol`, `D_tol`, and `Dq_tol`. Cross-check the root trace record against
   `prof.crit_star`, `prof.D_uni_star`, and `prof.Dq_min_star`.

For the separate PM controls `[3.1 3.5] T`, run only `invz_solve_point`; do **not** call the Jensen
HMF solver on a PM field. Each must converge with finite `crit > crit_tol`, strict provenance, and
an `'ok'` reference status, or (e) fires. Evaluate exact `B=0` separately with the return-mode
two-level/domain path and record the expected `degenerate_doublet`; it is a hard-domain control,
not a required ordered-path failure.

Compute G5 for each ordered field from the three `nH` runs. Both `int_Sigma0` and
`int_r_minus_1` must satisfy the frozen absolute-plus-relative 65→129 criterion. Record 33→65 as
the approach-to-convergence diagnostic. Gate 0 cannot pass if the finest integral is missing.

The report's pass predicate is exactly:

```matlab
rep.pass = ~(rep.fail_a || rep.fail_b || rep.fail_c || rep.fail_d || rep.fail_e);
```

Add a unit fixture for the aggregation logic so each of (a)–(e) can be made to fire
independently. This prevents a Boolean that is merely printed but never load-bearing.

- [ ] **Step 3: Run Gate 0/G5 on the frozen real multiset**

Use the exact Task-0 tuple and lists—no cache, no Ewald substitution, no omitted low-field rows:

```text
grid=[16 16 16], dpRng=30, dipole='bruteforce', cache=false
T=0.10 K
ordered_fields=[0.05 0.25 0.5 1 2 2.5 2.9 3.0] T
pm_fields=[3.1 3.5] T
nH=[33 65 129], Ecut=40 meV
```

The run command must print and the Markdown report must preserve:

- the requested coupling tuple, `info.dipole`, intentional absence of `info.grid`, exact
  class/shape/data SHA-256, `Jcc0`, and all four moments;
- every ordered `(B,nH)` row with root/status, endpoint margins, trace coverage counts,
  minimum reference margin, solved-ledger maximum omitted ratio, both integrals, and G17
  crossing/jump status;
- both PM-control rows and the B=0 domain-control result;
- the G5 33→65 and 65→129 differences and tolerances;
- separate Booleans (a)–(e), followed by the single overall verdict.

The `1 T` and `3 T` rows are mandatory direct comparisons with the diagnosed masked nodes.

- [ ] **Step 4: Run the G7 scheme-jump measurement**

Run:
```
INVZ_SLOW=1 "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath(genpath('invz_projected')); addpath('invz_common'); \
ion=invz_ion(); [Jnu,~]=invz_bz_couplings(ion,struct('grid',[16 16 16],'dpRng',30, \
'dipole','bruteforce','cache',false)); \
Jf=Jnu(:); mom=invz_coupling_moments(Jf); \
for T=[0.05 0.1 0.31 1.0], [wn,~,~]=invz_matsubara(T,40); \
si=invz_single_ion(ion,T,[6 0 0],struct('hyp',true,'Jxx0',ion.Jxx0,'transverse_mf','legacy_x')); \
c0=invz_chi0z(si,T,1i*wn,struct('elastic',true)); G0=-real(squeeze(c0(3,3,:))); S=zeros(size(wn)); \
a=invz_emt_scalar(G0,S,Jf,struct()); \
b=invz_emt_scalar(G0,S,Jf,struct('static_medium','strict_1z_dyson_ref','Jmom',mom)); \
fprintf('T=%-5.2f jump|K1s-K1r|=%.4g  dispersion|K2-K1|=%.4g  ratio=%.4g\n', T, \
abs(b.K(1)-a.K(1)), abs(a.K(2)-a.K(1)), abs(b.K(1)-a.K(1))/abs(a.K(2)-a.K(1))); end"
```
The scheme jump at `omega = 0` is *exactly* `K(1)_strict - K(1)_resummed`, because `omega_n != 0` is unchanged by construction. Reporting it against the physical dispersion `K(2) - K(1)` separates artifact from physics; a mixed `|K(iw_1) - K(0)|` would conflate them.

- [ ] **Step 5: Write the verdict document**

Create `docs/invzp_strict_medium_gate0_report.md` containing, in this order: the frozen predicate
quoted from prereg §3; the exact input fingerprint/provenance; the Step-3 ordered, PM, B=0, G5,
and G17 tables verbatim; the Step-4 G7 table verbatim; the explicit statement that the ~0.3% PM
boundary shift bounds neither path integral; and a one-line `PASS` / `FAIL` naming which of
(a)–(e) fired.

If **FAIL**: stop. Do not carry another moment, change `Gref`, or truncate other Matsubara sectors — those are new theory candidates needing a new spec and fresh preregistration (prereg §3), and regularisation, broadening and tolerance widening remain forbidden. Report the verdict and the diagnosis.

If **PASS**: stage 4 (G6/G6d, G8, G10, G12, G14, G16, then the default-flip decision) becomes worth planning, as its own plan.

- [ ] **Step 6: Commit**

```bash
git add invz_projected/invz_gate0_report.m invz_projected/invz_pole_continuity.m \
  invz_projected/tests/test_invz_pole_continuity.m docs/invzp_strict_medium_gate0_report.md
git commit -m "diag(invzp): Gate-0 domain/omitted-order report plus G5 and G7 measurements"
```

---

## Self-review

**Spec coverage.** Every in-scope §0–§5 requirement maps to a task: §0.3 one-shot → T6–T8; §1
construction/two-tier/removable pole → T6–T9, T12, T14, T18; §2 components → T1–T5, T9, T15,
T16; §4.1–4.3 → T1–T4, T7–T8, T12, T14; §4.4 → T9, T13–T15; §4.5 no-cache-change → T15
Step 7; §4.6 `K(1)` ordering → T17/G13; §5.1 classifier at all strict-path catches → T5,
T9–T10, T15, T17/G15; §5.2 → T8, T10, T13; §5.3 → T11, T13; §5.5 → T12 trace ledger,
T16 prospective counters, T18 solved-ledger counters; §6.0 → T0; G0/G5/actual G17 → T18;
G1/G2/G3 → T17; G4 → T9/T14; G7 → T18; G9 → every Stage-2 compatibility check; G11/G13/G15
→ T17. **Deliberately absent:** G6/G6d, G8, G10, full end-to-end G12, G14, G16—stage 4;
`one_field` un-nesting and ordered `[nJ,nw]` flattening—separate scoped changes. PM `[nJ,nw]`
support is retained in T7.

**Execution placeholders.** `<FROZEN_HASH_FROM_TASK_0>` is an instruction token, not literal test
data: when Task 17 writes the test it must copy the digest frozen by Task 0, and the token must not
remain in the created `.m` file. The three judgement-bearing values (`crit_tol`, `omit_promote`,
`pole_cont_tol`) are FROZEN by user approval on 2026-07-25, and the §9 blind Stage-4 table is now
populated with real numbers and formulas, each labelled INHERITED / DERIVED / CHOICE — it awaits its
own approval and is the last item blocking Task 1. `<DATE OF USER APPROVAL>` in the §0 header of the
generated prereg is the only remaining fill-in token, resolved when Task 0 Step 5 commits.
Task 15 Step 1 re-locates its assembly site with `rg` because
line numbers will move as preceding tasks land; it is a location step, not an unspecified behavior.

**Type and failure consistency.** Task 4 validation is idempotent for agreeing stamps and rejects
disagreement whether it arrives through `opts.emt`/`opts.emt_static` or explicit `eopts`/`eso`.
Task 7 retains the PM matrix interface and selects its static column; ordered strict mode rejects the
unsupported matrix only after actual coupling resolution. Task 9 uses a direct fixed-h fixture, so it
has no knowingly failing dependency on Task 14. Task 18 has no broad catch and reads the complete
trace ledger rather than filtering failed nodes out of the gate. The exact-byte coupling digest is
one shared primitive, `invz_exact_numeric_digest` (Task 0), with three consumers -- the Task-0
fingerprint, the Task-17 G11 anchor, and the Task-18 Gate-0 driver -- rather than three
hand-reproduced copies of the same algorithm.

**Known numerical risk, made testable.** Task 6's stable form is algebraically reassociated and can
move last bits, so it is strict-gated; existing seven-argument/resummed calls retain their historical
arithmetic. G17 tests the exact limit synthetically, while Task 18 tests finite continuity on the
actual solved path with a preregistered estimator.

---

## Review disposition and execution handoff

**Conditionally approved as an implementation plan.** No production implementation is authorised by
this review. The three judgement-bearing constants are approved and frozen. Task 0 remains a hard
gate only for the input fingerprint/baseline and explicit user approval of the revised blind
G8/G14/G16 protocol plus the settled scope/grids. Commit the approved preregistration before Task 1.
If any later task exposes a theory change rather than an implementation correction, stop and amend
the design/preregistration instead of widening this plan in flight.
