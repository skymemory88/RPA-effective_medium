function A = invz_odd_anchors()
%INVZ_ODD_ANCHORS Controller-verified P0 preflight digit anchors (ODD main body).
%   Plain fixture function (NOT a test; runtests on tests/ is non-recursive and
%   does not collect this file). Later ODD tests reference these fields instead of
%   hard-coding unpinned digits. Do NOT edit a pinned value to make a test pass:
%   a mismatch is a finding to escalate (Global Constraints, anchors fixture).
%
%   PROVENANCE
%     Measured 2026-07-16 on branch invz-1z-lihof4 at git 360dfab (clean tree;
%     360dfab differs from the fast-baseline commit 5f4ff92 only by two plan .md
%     docs -- no module code). Numbers produced by the read-only exploratory
%     scripts committed alongside:
%       invz_projected/tests/exploratory/explore_chiperp.m    (P0.2, chi_perp)
%       invz_projected/tests/exploratory/explore_odd_blocks.m (P0.3, ODD blocks)
%     MATLAB R2025a. ODD lattice sums at dpRng = 30 (production default).
%     Full ODD-LOG entry: docs/ODD-LOG.md, section P0.
%
%   All literals are %.17g round-trips of the exact doubles emitted by the scripts.

% --- P0.2: transverse single-ion susceptibility chi_perp ---------------------
% Symmetrized (a,b)x(a,b) block of the full electronuclear invz_chi0z(si,T,0,
% struct('elastic',true)) at T = 1.53 K, B = [0 0 0], hyp = true, default
% transverse_mf ('legacy_x'). Van Vleck-dominated; chi_aa = chi_bb to ~4e-14 at
% Bx = 0 (C4); off-diagonal / antisymmetric parts are machine-zero (gyrotropic-
% free). Value 17.638 meV^-1 sits at/just above the 16-17 meV^-1 full-CF band
% (far from the 11 truncation and x2 convention-slip failure modes). [meV^-1]
A.chiperp_1p53K_0T = [ 17.63784561529863      -8.5583683405846167e-16 ; ...
                      -8.5583683405846167e-16  17.637845615298673     ];

% Max abs antisymmetric (gyrotropic) part of the same (a,b) block -- machine zero.
A.chiperp_asym_1p53K = 6.8157632679088228e-16;

% Share of chi_aa carried by the elastic (Curie/beta) term = (Xp - Xp_noelastic)/Xp
% at (1.53 K, 0 T): 0.064% -> chi_perp is Van Vleck (inelastic) dominated.
A.chiperp_elastic_share_1p53K = 0.00063874219711189612;

% Bx = 0:6 T sweep at T = 0.31 K (elastic on, hyp, default MF). chi_aa (response
% along the field a-axis) PEAKS at Bx = 1 T (21.03) then collapses -- the known
% (0.31 K, 1-2 T) island; every point is MF-converged (mf_residual < 1e-12), so
% the peak/collapse is physical, not a convergence artifact. NOTE (escalated in
% ODD-LOG): the max relative step is 0.51 (Bx 1->2 T), so a "< 0.25 no-jumps"
% gate on this curve is not achievable -- Task 3's smoothness threshold must be
% recalibrated to this measured curve, which IS the pinned truth. [meV^-1]
A.chiperp_0p31K_Bx.Bx     = 0:6;
A.chiperp_0p31K_Bx.chi_aa = [ 17.908870144743222  21.026941222931679 ...
                              10.366505439821841    5.1057826081515003 ...
                               3.0794604377402943   2.1173208498251412 ...
                               1.5855606520750753 ];
A.chiperp_0p31K_Bx.chi_bb = [ 17.908870144743165  20.226001745573051 ...
                              18.890579447920953   15.888029236787382 ...
                              13.683247121988876   12.147124794054848 ...
                              11.045629276479632 ];

% --- P0.3: ODD off-diagonal dipole blocks J^{c,alpha}(q) = -gfac*dip(3,alpha,:,:) --
% On-axis ray [q 0 0], dpRng 30: max |J^{ca}| element of the 4x4 sublattice block
% at q = {1e-1, 1e-2, 1e-3}. DECAYS LINEARLY in q (ratio ~10 per decade -> 0 at
% Gamma, C2-about-c). IMPORTANT (escalated in ODD-LOG): linear decay means the
% residual at q = 1e-3 is 1.86e-5 meV = 2.9e-3 * Jcc0, NOT <= 1e-6 * Jcc0. The
% plan/Task-2 absolute bound "m(3) <= 1e-6*Jcc0" is unachievable at q = 1e-3 and
% must be recalibrated to this measured linear-decay residual. [meV]
A.odd_onaxis_smallq.q     = [1e-1 1e-2 1e-3];
A.odd_onaxis_smallq.maxca = [ 0.0017881297636183286 ...
                              0.00018564848786664122 ...
                              1.8590576478846918e-05 ];

% Generic q = [0.31 0.17 0.09], dpRng 30: max |J^{ca}| element = 4.09e-3 meV =
% 0.637 * Jcc0. Off-diagonal ODD coupling is an O(Jcc0) effect at generic q
% (vanishing only near high-symmetry directions). [meV]
A.odd_generic_q_max = 0.0040909010697032459;

% --- Reference couplings from invz_jq_modes (dpRng 30) -----------------------
% Match the published anchors: Jcc0 6.421e-3 (here 6.4244e-3), Jcc0_dipole
% 6.821e-3 (here 6.8244e-3), Jaa0_dipole 3.912e-3 (here 3.9104e-3). [meV]
A.jcc0 = 0.0064244356557014957;
A.jaa0 = 0.0035104462050649012;

% --- T1.5: zero-field ODD headline (invz_odd_zero_field, mode = 'full') -------
% Richardson(12^3, 24^3, dpRng 30) on the Sigma_c-benchmark generator mesh
% (qVec_generator range [-0.5 0.5], Gamma-excluded), GOVERNING condition/Sigma-space
% split (controller adjudication 2026-07-17). Measured 2026-07-17 on branch
% invz-1z-lihof4. Tc(0) 1.743 K (off) -> 1.50937 K (ODD): dTc(0) = 0.2341 K (13.4%);
% Sigma_c(0) 0.29798 (off) -> 0.38880 (ODD); d(Tc) = 0.483 ueV (in the 0.3-0.5 ueV
% report band). Full-precision doubles emitted by invz_odd_zero_field; the headline
% slow test test_odd_headline_slow reproduces these at RelTol 1%. [K, meV, -]
A.odd_Tc_rich  = 1.509370677196421;       % Richardson Tc(0) with ODD (mode 'full')
A.odd_d_at_Tc  = 0.00048311966308299265;  % E5 uniform reduction d at Tc_rich (meV)
A.odd_Sc_rich  = 0.38879543801229982;     % Richardson Sigma_c(0) with ODD
end
