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
end
