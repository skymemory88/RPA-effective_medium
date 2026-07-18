function A = invzt_anchors()
%INVZT_ANCHORS A0 preflight digit anchors for the invz_tensor CORE suite.
%   Plain fixture function (NOT a test; runtests on tests/ is non-recursive and
%   does not collect this file). CORE tests (e.g. Task 3's
%   test_onaxis_smallq_odd_decay) read these fields instead of hard-coding
%   unpinned digits. Do NOT edit a pinned value to make a test pass: a mismatch
%   is a finding to escalate (Global Constraints, anchors fixture rule).
%
%   Interop tests may instead read invz_odd_blocks / invz_odd_anchors()
%   directly (invz_projected/tests/invz_odd_anchors.m, superseded P0 preflight)
%   -- this file is the TENSOR-owned copy, measured independently on THIS tree
%   so invz_tensor/tests never depends on invz_projected/tests at run time.
%
%   PROVENANCE
%     Measured 2026-07-17 on branch invz-1z-lihof4 at git 8b0de0d (clean tree;
%     HEAD after Task 1's invz_common/ move -- pure `git mv`, behavior-neutral,
%     so these numbers are bit-identical to the superseded P0.3 preflight
%     anchors in invz_projected/tests/invz_odd_anchors.m, confirmed by direct
%     comparison during measurement). MATLAB R2025a. Produced by the read-only
%     exploratory script committed alongside:
%       invz_tensor/tests/exploratory/explore_tensor_blocks.m
%     which calls invz_odd_blocks (NOT raw MF_dipole) at dpRng = 30 (primary;
%     dpRng in {20, 40} additionally for the (v) sensitivity sweep). Full A0
%     ODD-LOG entry: docs/ODD-LOG.md, section A0 (controller-appended).
%
%   All literals are %.17g round-trips of the exact doubles emitted by the
%   script (invz_tensor/tests/exploratory/explore_tensor_blocks.m stdout,
%   captured in .superpowers/sdd/task-2-report.md).

% --- on-axis small-q ODD decay: max|J^{ca}(q)| along [q 0 0], dpRng 30 -------
% LINEAR decay (ratio ~10 per decade -> 0 at Gamma, C2-about-c symmetry). The
% residual at q = 1e-3 is 1.86e-5 meV = 2.894e-3 * Jcc0 -- an O(q) tail, NOT an
% absolute <= 1e-6*Jcc0 bound (that bound is unachievable at q = 1e-3 given
% linear-in-q decay; see the superseded P0.3 preflight note). [q: r.l.u.; maxca: meV]
A.odd_onaxis_smallq.q     = [1e-1 1e-2 1e-3];
A.odd_onaxis_smallq.maxca = [ 0.0017881297636183286 ...
                              0.00018564848786664122 ...
                              1.8590576478846918e-05 ];
% Bonus context (not gated by any test): the SAME q magnitudes along [0 0 q]
% (c*-axis). max|J^{ca}| is machine-zero (~1e-17) here by symmetry: with q along
% c the a<->b tetragonal rotation forces this ODD component to vanish exactly,
% confirming "ODD blocks strictly c<->(a,b), never within the (a,b) plane
% itself along high-symmetry c*" is not violated. Do NOT read the tiny
% dpRng-sensitivity of these numbers as physical (see dpRng_sensitivity below).
A.odd_onaxis_smallq.maxca_00q = [ 1.6919090887918477e-17 ...
                                   1.3678682610649618e-17 ...
                                   1.2513741335047746e-17 ];

% --- tilted-ray limit: q*[1 0 1]/sqrt(2) at q = 1e-3, dpRng 30 --------------
% NON-decaying direction-dependent macroscopic dipole limit (in contrast to the
% on-axis linear decay above): max|J^{ca}| stays a sizeable fraction of the
% shape scale 4*pi*C.gfac/ion.Vc as q -> 0 along this direction, so small-q
% decay assertions must stay ON-AXIS ONLY (never grid-extrapolate q -> 0 off
% high-symmetry rays; ODD pitfall (ii), Global Constraints). [meV]
A.odd_tilted_limit.q          = 1e-3;
A.odd_tilted_limit.dir         = [1 0 1]/sqrt(2);
A.odd_tilted_limit.maxca       = 1.3515062966663823e-05;
A.odd_tilted_limit.macroscale  = 0.0036613320164578279;   % 4*pi*C.gfac/ion.Vc (~3.66 ueV)
A.odd_tilted_limit.ratio       = 0.0036912967482634999;   % maxca / macroscale at q=1e-3

% --- generic q = [0.31 0.17 0.09], dpRng 30 ---------------------------------
% max|J^{ca}| element = 4.09e-3 meV = 0.637 * Jcc0: the off-diagonal ODD
% coupling is an O(Jcc0) effect at generic q (vanishing only near high-symmetry
% directions, per the on-axis/tilted contrast above). [meV]
A.odd_generic_q_max = 0.0040909010697032459;

% --- Gamma info from invz_odd_blocks (dpRng 30) -----------------------------
% Match the published anchors: Jcc0 6.421e-3 (here 6.4244e-3), Jcc0_dipole
% 6.821e-3 (here 6.8244e-3), Jaa0_dipole 3.912e-3 (here 3.9104e-3). jaa0 below
% is the FULL info.Jaa0 = Jaa0_dipole + 4*ion.J12 (NOT the dipole-only
% published anchor -- do not conflate the two). [meV]
A.jcc0 = 0.0064244356557014957;
A.jaa0 = 0.0035104462050649012;

% --- (v) dpRng sensitivity of the smallest 16^3-grid shell, a*/c* -----------
% q = [1/16 0 0] (a*) and [0 0 1/16] (c*), dpRng in {20, 30, 40}. The a* shell
% is well converged: relative spread (dpRng 20 vs 40) ~ 1.82e-4 (0.018%), i.e.
% dpRng = 30 (let alone 20) is already deep in the converged regime for this
% ODD residual. The c* shell values are ALL machine-noise-floor (~1e-17 to
% 1e-18 meV, i.e. exact zero by the same c*-axis symmetry null noted above) --
% their relative spread is dominated by floating-point noise and carries NO
% physical meaning; report the a* column as the sensitivity result. [meV]
A.dpRng_sensitivity.dpRng     = [20 30 40];
A.dpRng_sensitivity.maxca_a16 = [ 0.0011428204693998804 ...
                                   0.0011423739697444429 ...
                                   0.0011426126292293662 ];
A.dpRng_sensitivity.maxca_c16 = [ 2.3462557588702163e-18 ...
                                   1.0830376638918350e-17 ...
                                   1.2224671454486072e-17 ];   % noise floor, not physical
A.dpRng_sensitivity.Jcc0      = [ 0.0064230405010920833 ...
                                   0.0064244356557014957 ...
                                   0.0064218498276640218 ];
A.dpRng_sensitivity.Jaa0      = [ 0.0035111437823696070 ...
                                   0.0035104462050649012 ...
                                   0.0035117391190836412 ];

% =============================================================================
% TASK-14 HEADLINE ANCHORS (reproducibility pinning) -------------------------
% =============================================================================
% PROVENANCE: appended 2026-07-18 on branch invz-1z-lihof4 at git 856eeab
% (Tasks 1-13 complete; A0-A4 all landed; docs/ODD-LOG.md SSA1/SSA3/SSA4). These
% pin the CURRENT measured values -- the physics EVOLVED during execution, so
% these are NOT the plan's original assumptions (in particular, Task 12's
% original "lambda-slope >= 2.3" protocol was reframed mid-execution to the
% matched-truncation ratio+collapse+band gate below; see docs/ODD-LOG.md SSA3
% and the coordinator RESOLUTION in .superpowers/sdd/task-12-report.md -- the
% superseded slope>=2.3 number is deliberately NOT pinned here). Same read-only
% convention as the rest of this file: a mismatch on a future rerun is a
% finding to escalate, never silently updated to match new code.

% --- A1 proxy-Tc (SSA1, Task 7 commit 357e24c: invzt_critical.m / ---------
%     invzt_tc_pm_extrap.m; reproducibility-gated by ------------------------
%     test_headline_reproducibility_slow, tests/interop/ ----------------------
%     test_invzt_critical_parity.m) ---------------------------------------
% invzt_tc_pm_extrap on the SMALL-Bx proxy (0.05 T, ODD on), 16^3
% legacy_inclusive grid, dpRng 30, T-grid 1.30:1/30:1.75, mode 'a1'/nlevels
% 'std' (invzt_solve_point). T7 originally measured this to 1.5599 K
% (commit 357e24c); the full-precision literal below is a Task-14 controller
% re-run of the IDENTICAL protocol on THIS tree (ad hoc script, format long g,
% 2026-07-18, 562 s), reproducing 1.5598631081277523 K -- agreement to 5
% decimal places across the T8-T13 intervening commits (A1 untouched by A2-A4
% work, as expected). Comparator: invz_odd_zero_field('full', 16^3, dpRng 30)
% grid-matched projected closed form, NOT the production Richardson (12^3,24^3)
% anchor 1.509 K (see the T7 report's grid-match nuance). [K]
A.task14.a1_proxy_tc.Tc_tensor_16cubed          = 1.5598631081277523;
A.task14.a1_proxy_tc.Tc_projected_16cubed       = 1.5442359495209526;
A.task14.a1_proxy_tc.gap_tensor_minus_projected = 1.5598631081277523 - 1.5442359495209526;
A.task14.a1_proxy_tc.grid     = 'legacy_inclusive';
A.task14.a1_proxy_tc.n        = 16;
A.task14.a1_proxy_tc.dpRng    = 30;
A.task14.a1_proxy_tc.proxy_Bx_T = 0.05;

% --- A3 SS11.8 emergence, MATCHED-TRUNCATION form (SSA3, Task 12 commit -----
%     b22dd70: invzt_threestate.m / invzt_sigma_tensor.m; -------------------
%     tests/test_invzt_a3_threestate.m) ------------------------------------
% Two components, BOTH re-run at full precision by a Task-14 controller script
% replicating the committed tests verbatim (format long g, 2026-07-18):
%  (1) EXACT rho->0 scalar-compatibility identity (test_a3_scalar_compat_exact_
%      rho0, T=1.6 K/Bx=0.5 T, 6^3 halfopen/dpRng 10): A3's Sigma_cc_equiv(1)
%      vs the scalar invz_emt_scalar+invz_sigma chain on the SAME cc branch
%      spectrum. Originally reported |diff|=3.24e-11.
%  (2) MATCHED-TRUNCATION ODD-sector emergence (test_a3_emergence_matched_
%      truncation, STABLE PM anchor T=2.0 K/Bx=0.5 T, same grid): rd(1) =
%      dominant-dress-A3/A1 crit-shift ratio (matches E1 up to the O(1/z^2)
%      constraint-8 resummation band); rf(1) = full-A3/A1 ratio (REPORT: the
%      genuine beyond-E1 transverse-spectator dressing); collapse =
%      |rd(1)-1|/|rf(1)-1| (transverse-dressing removal fraction); slope is a
%      REPORT only (O(1/z^2)-capped, NOT the superseded >=2.3 gate).
%      Originally reported rd(1)=1.0159, rf(1)=1.1132, collapse=0.140.
A.task14.a3_emergence.rho0_A3_sigma_cc_equiv = 0.36403693825388755;
A.task14.a3_emergence.rho0_scalar_sigma      = 0.364036938221526;
A.task14.a3_emergence.rho0_absdiff           = 3.23615578778913e-11;
A.task14.a3_emergence.rho0_T_K  = 1.6;
A.task14.a3_emergence.rho0_Bx_T = 0.5;
A.task14.a3_emergence.rd1 = 1.0158938460775333;
A.task14.a3_emergence.rf1 = 1.1131532562138944;
A.task14.a3_emergence.collapse             = 0.14046300220904873;   % |rd1-1|/|rf1-1|
A.task14.a3_emergence.band_abs_rd1_minus1  = 0.015893846077533302;  % O(1/z^2) resummation band
A.task14.a3_emergence.slope_report         = 2.0563526186836825;    % REPORT only
A.task14.a3_emergence.resum_spread_crit    = -0.039327890102603319; % crit(dyson)-crit(additive)
A.task14.a3_emergence.anchor_T_K  = 2.0;
A.task14.a3_emergence.anchor_Bx_T = 0.5;
A.task14.a3_emergence.grid  = 'halfopen';
A.task14.a3_emergence.n     = 6;
A.task14.a3_emergence.dpRng = 10;

% --- A4 basis-defined ladder climb (SSA4, Task 13 commits 920f440 -----------
%     [feat] + 9414425 [N-adaptive fix]; invzt_run_ladder.m) ----------------
% Full-A3/A1 crit-shift ratio rf at the STABLE PM emergence anchor (2.0 K,
% 0.5 T; same anchor as a3_emergence above), across rungs three (N=3 toy) ->
% e3 (N=3 real CF basis) -> e6 (N=6 real CF basis, ~14 min/point, the largest
% rung proven end-to-end): the beyond-E1 transverse-spectator-dressing share
% (rf-1) CONVERGES 11.3% -> 2.6% as CF content grows, landing near the
% projected Tier-2 share (~2.8%, REPORT comparator, never tuned). Virtual-
% completeness deficit (chi0 rung vs full-136, DIAGNOSTIC not a bound) drops
% monotonically e3 0.041 -> e6 0.002. crit_shift_odd (+ = ODD lowers Tc): the
% e6 value is at THIS SAME (2.0 K, 0.5 T) anchor; the three/e3 values are at
% the SEPARATE LOCKED validation point (1.6 K, 0.1 T, sub-critical regime, NOT
% the clean-emergence point -- see the T13 report's regime note) since e6 was
% not re-run there (budget: a single e6 point is ~14-15 min; T13 measured both
% points for e6 already). These are PINNED AT THE PRECISION MEASURED AND
% RECORDED in commits 920f440/9414425 / docs/ODD-LOG.md SSA4 (NOT re-derived
% at higher precision here -- a fresh e6 rung re-run is production-ladder-scale
% compute, out of this task's cheap-spot-check budget; see the Task-14 report
% for the explicit time-budget rationale).
A.task14.a4_ladder.anchor_T_K  = 2.0;
A.task14.a4_ladder.anchor_Bx_T = 0.5;
A.task14.a4_ladder.rf_three = 1.1132;
A.task14.a4_ladder.rf_e3    = 1.1140;
A.task14.a4_ladder.rf_e6    = 1.0263;
A.task14.a4_ladder.vdef_e3  = 0.041;
A.task14.a4_ladder.vdef_e6  = 0.002;
A.task14.a4_ladder.validation_T_K  = 1.6;
A.task14.a4_ladder.validation_Bx_T = 0.1;
A.task14.a4_ladder.crit_shift_odd_e6_anchor           = 0.0547;
A.task14.a4_ladder.crit_shift_odd_three_validation     = 0.02812;
A.task14.a4_ladder.crit_shift_odd_e3_validation        = 0.03093;
end
