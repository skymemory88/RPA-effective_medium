# Strict-order static medium — FROZEN preregistration

Frozen: 2026-07-25. Branch invzp-stage2c-diagnostic. Spec:
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

## 10. Measured coupling fingerprint (recorded 2026-07-25, pre-strict baseline)
MEASURED, not a new tolerance or threshold decision: this section records the §8 production
fixture's actual computed output, per §8's "record the requested tuple ... record
invz_exact_numeric_digest(Jnu(:))" instruction. Appended rather than edited into §8 in place,
per this document's own amendment rule. Recomputed via
invz_bz_couplings(invz_ion(), struct('grid',[16 16 16],'dpRng',30,'dipole','bruteforce',
'cache',false)) on branch invzp-stage2c-diagnostic, HEAD 2ee310b.

info.grid    ABSENT -- isfield(info,'grid') == false, asserted. Confirms the legacy/
             no-grid-policy route (invz_bz_couplings.m:26); not invented.
info.dipole  backend          = 'bruteforce'
             ewald            = struct('alpha',[],'r_cut',[],'g_cut',[],'boundary','')  (sentinel-empty)
             q_reduction      = 'bruteforce: q used directly as MF_dipole/exchange Miller indices
                                 (q*geom.b); no canonical q-domain reduction applied'
             primitive_schema = 'MF_dipole+exchange (legacy, unversioned)'
n            = 16384 (4096 q x 4 branches)
digest       invz_exact_numeric_digest(Jnu(:)) =
             ddb9532d11326458554b93b2ce09c80a3299cca9baa69202dc311f138b4fae17
             (64 lowercase-hex chars, stable across two successive calls in the same run)

info.Jcc0 (J0eff)         =  6.424435656e-3 meV
Jaa0 (3rd output)         =  3.510446205e-3 meV
Jbar  = mean_q J          =  1.207664433e-4 meV
mu2   = mean((J-Jbar).^2) =  5.482637653e-6
mu3   = mean((J-Jbar).^3) = -3.42227577e-11        (skewness mu3/mu2^1.5 = -0.0026658)
mu4   = mean((J-Jbar).^4) =  7.182350058e-11        (= 2.389392*mu2^2)
J_max, J_min              =  5.985138929e-3, -6.763100317e-3 meV
(J_max-Jbar)/sqrt(mu2)    =  2.504533 sigma

mu2/mu3/mu4 use the population normalization mu_n = mean((J-Jbar).^n), NOT var() -- computed
INLINE here, not via invz_common/invz_coupling_moments.m, because that primitive is Task 1's
deliverable (design spec line 315) and does not exist yet; this record does not depend on it.
These values agree with the design spec's independently computed "Measured multiset" table
(§B) to the last digit quoted there (J0eff 6.42444e-3, Jbar 1.20766e-4, sqrt(mu2) 2.3415e-3,
mu2 5.48264e-6, mu3 -3.42228e-11, mu4/mu2^2 2.3894, J_max/J_min 5.98514e-3/-6.7631e-3,
(J_max-Jbar)/sqrt(mu2) 2.5045 sigma) -- an independent cross-check that this is the same
fixture computed the same way, not a coincidence of rounding.
