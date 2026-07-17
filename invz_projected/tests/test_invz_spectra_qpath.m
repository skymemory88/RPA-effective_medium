function tests = test_invz_spectra_qpath
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_gamma_approach_regression(testCase)
% Regression: the truncated MF_dipole sum used to collapse and jump when approaching
% the Gamma-equivalent point (2,0,0). The direction-aware snap must give a smooth,
% cutoff-stable approach; for this in-plane path (khat_z=0) the directional limit is
% the uniform-mode Lorentz value, so the endpoint equals info.Jcc0 by construction.
ion = invz_ion();
hs = [1.90 1.96 1.98 1.99 1.999 2.0].';
qpath = [hs zeros(numel(hs), 2)];
P20 = invz_jq_path(ion, qpath, struct('dpRng', 20, 'cache', false));
P40 = invz_jq_path(ion, qpath, struct('dpRng', 40, 'cache', false));
J20 = P20.Jnu(:, 4);  J40 = P40.Jnu(:, 4);          % max branch (ascending sort)
[~, iref] = invz_jq_modes(ion, [0 0 0], struct('dpRng', 40, 'cache', false));
verifyEqual(testCase, J40(end), iref.Jcc0, 'RelTol', 1e-9);   % endpoint = uniform mode
% no dip/jump: every point within 3% of the endpoint on this fine approach
verifyGreaterThan(testCase, min(J40)/J40(end), 0.97);
verifyGreaterThan(testCase, min(J20)/J20(end), 0.97);
% cutoff stability: dpRng 20 vs 40 agree pointwise within 2%
verifyEqual(testCase, J20, J40, 'RelTol', 0.02);
% the guard actually engaged near the endpoint, and reports it
verifyTrue(testCase, any(P40.snapped) && P40.snapped(end));
% dual path coordinates: index (r.l.u.) and Cartesian (Ang^-1)
verifyEqual(testCase, P40.s(end), 0.10, 'AbsTol', 1e-12);
verifyEqual(testCase, P40.s_cart(end), 2*pi*0.10/ion.a(1,1), 'RelTol', 1e-9);
end

function test_qpath_structure_and_gamma_limit(testCase)
% Structural contract + physics anchors of the q-path spectrum (default = uniform FM mode):
% shapes are [nw x nq]; (2,0,0) is Gamma-equivalent so S.Jq there equals info.Jcc0; and the
% uniform-mode coupling grows monotonically h=1->2, so Epeak softens monotonically toward
% Gamma (R 2007 Fig 3). Window must reach ~0.9 meV since the zone-edge (1,0,0) mode is highest.
ion = invz_ion();
T = 0.31;  B = 5.5;                        % paramagnetic side: fast, well-converged
w = (0.02:0.02:0.9).';
info = struct('Jcc0', 6.4e-3);             % synthetic medium (as in test_invz_spectra_map)
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
qpath = [1 0 0; 1.5 0 0; 2 0 0];
S = invz_spectra_qpath(ion, T, B, qpath, w, ...
                       struct('Jnu', Jnu, 'info', info, 'dpRng', 10));

verifySize(testCase, S.chiz,   [numel(w) 3]);
verifySize(testCase, S.chirpa, [numel(w) 3]);
verifyEqual(testCase, S.s, [0 0.5 1], 'AbsTol', 1e-12);
% plot coordinate: single-axis path -> the ACTUAL Miller component (R2007 Fig 3 axis,
% h = 1..2), NOT the distance-from-start s (which reads 0..1 for every unit-length path)
verifyEqual(testCase, S.x, [1 1.5 2], 'AbsTol', 1e-12);
verifyEqual(testCase, S.xlab, 'Q = (h, 0, 0) (r.l.u.)');
verifyEqual(testCase, S.phase, 2);
verifyTrue(testCase, all(isfinite(S.Epeak)));

[~, iref] = invz_jq_modes(ion, [0 0 0], struct('dpRng', 10, 'cache', false));
verifyEqual(testCase, S.Jq(3), iref.Jcc0, 'RelTol', 1e-9);   % (2,0,0) = uniform info.Jcc0
verifyGreaterThan(testCase, S.Jq(2), S.Jq(1));               % coupling grows toward Gamma...
verifyGreaterThan(testCase, S.Jq(3), S.Jq(2));               % ...monotonically
% ...so the mode softens monotonically toward Gamma (Epeak: 0.68 -> 0.45 -> 0.40 meV)
verifyGreaterThan(testCase, S.Epeak(1), S.Epeak(2));
verifyGreaterThan(testCase, S.Epeak(2), S.Epeak(3));
end

function test_qpath_peak_censoring(testCase)
% Review finding (blocker): a maximum in the first/last usable frequency bin means
% the true peak lies outside the window and must be censored (NaN), not reported.
ion = invz_ion();
T = 0.31;  B = 5.5;
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
qpath = [1 0 0; 2 0 0];
o = struct('Jnu', Jnu, 'info', info, 'dpRng', 10);
% (a) clipped window: the mode at 5.5 T sits above 0.15 meV at every path point
wclip = (0.02:0.02:0.15).';
Sc = invz_spectra_qpath(ion, T, B, qpath, wclip, o);
verifyTrue(testCase, all(isnan(Sc.Epeak)));
% (b) peak_wmin at/above the window top: nothing usable, all censored, no error
w = (0.02:0.02:0.6).';
o2 = o;  o2.peak_wmin = 1.0;
Sm = invz_spectra_qpath(ion, T, B, qpath, w, o2);
verifyTrue(testCase, all(isnan(Sm.Epeak)));
end

function test_qpath_no_ordering_channel_shape_leak(testCase)
% With the transverse channel held identical (synthetic info carries no Jaa0, so
% both runs fall back to ion.Jxx0), the path spectra must be bit-identical for
% demag on vs off: the guard and the couplings carry no ordering-channel shape
% dependence. (The REAL demag pathway into q-path spectra is info.Jaa0 only.)
ion0 = invz_ion();
ionS = invz_ion();  ionS.demag = 1;  ionS.alpha = 1;
T = 0.31;  B = 5.5;  w = (0.05:0.05:0.5).';
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
qpath = [1 0 0; 2 0 0];
o = struct('Jnu', Jnu, 'info', info, 'dpRng', 10);
S0 = invz_spectra_qpath(ion0, T, B, qpath, w, o);
SS = invz_spectra_qpath(ionS, T, B, qpath, w, o);
verifyEqual(testCase, SS.chiz, S0.chiz);
verifyEqual(testCase, SS.Jq,   S0.Jq);
end

function test_qpath_plot_coordinate_fallback(testCase)
% A path varying along more than one Miller axis has no single q coordinate to plot
% against: S.x must fall back to the distance-from-start S.s (with the matching label).
ion = invz_ion();
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
w = (0.05:0.05:0.3).';
qpath = [1 0 0; 1.5 0 0.5; 2 0 1];                 % h and l both vary
S = invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, struct('Jnu', Jnu, 'info', info, 'dpRng', 10));
verifyEqual(testCase, S.x, S.s);
verifyTrue(testCase, startsWith(S.xlab, 's along path'));
end

function test_qpath_fm_mode_monotonic(testCase)
% R 2007 Fig 3: along (1,0,0)->(2,0,0) the FM-mode coupling rises MONOTONICALLY as the
% mode softens toward the zone centre. The physical dispersion is the UNIFORM FM-mode
% projection v'*Jcc*v (P.Juni), not the max-eigenvalue branch P.Jnu(:,4), which instead
% mirrors about h=1.5 (wrong sublattice branch for h<1.5) -- the regression this guards against.
ion = invz_ion();
hs = linspace(1, 2, 21).';                     % 0.05 steps: hits 1.90,1.95,2.00 (avoids the
qpath = [hs zeros(numel(hs), 2)];              % lone h=1.96 truncation wiggle in the raw sum)
P = invz_jq_path(ion, qpath, struct('dpRng', 30, 'cache', false));
verifySize(testCase, P.Juni, [21 1]);
% physical FM-mode coupling: monotonically non-decreasing across the whole path
verifyGreaterThanOrEqual(testCase, min(diff(P.Juni)), -1e-9);
% and it genuinely rises: weak/negative at the (1,0,0) edge, strongest at Gamma (2,0,0)
verifyLessThan(testCase, P.Juni(1), 0);                          % ~ -5.2 ueV at (1,0,0)
[~, iref] = invz_jq_modes(ion, [0 0 0], struct('dpRng', 30, 'cache', false));
verifyEqual(testCase, P.Juni(end), iref.Jcc0, 'RelTol', 1e-9);   % (2,0,0) = uniform info.Jcc0
% the max-eig branch is demonstrably NON-monotonic (mirror artifact): it has a clear
% downhill step somewhere (the mirror falls from ~6.4 ueV at h=1.05 to ~2.9 at h=1.5),
% whereas a monotonic dispersion would have min(diff) >= 0. This is WHY the default
% follows Juni, not P.Jnu(:,4).
verifyLessThan(testCase, min(diff(P.Jnu(:,4))), -1e-4);
end

function test_qpath_out_of_plane_gamma_limit(testCase)
% The direction-aware guard is ANISOTROPIC: approaching Gamma along c (khat_z=1) the
% nonanalytic broadcast shifts the uniform branch by 4*gfac*(4pi/Vc) BELOW the in-plane
% limit, leaving the three non-uniform branches unchanged. Pins the kz2 coefficient and
% its sign, which the in-plane tests alone cannot see.
ion = invz_ion();  C = invz_const();
Px = invz_jq_path(ion, [1.99 0 0], struct('dpRng', 10, 'cache', false));  % in-plane, snapped
Pz = invz_jq_path(ion, [0 0 3.99], struct('dpRng', 10, 'cache', false));  % along c, snapped
verifyTrue(testCase, Px.snapped(1) && Pz.snapped(1));
d = 4*C.gfac*(4*pi/ion.Vc);                     % uniform-branch shift between the two limits
verifyEqual(testCase, Pz.Jnu(1,:), sort([Px.Jnu(1,1:3), Px.Jnu(1,4) - d]), ...
            'RelTol', 1e-9, 'AbsTol', 1e-15);
% a single-point path exactly AT G has no local direction: defaults to the
% in-plane (uniform-mode) convention, i.e. the max branch equals info.Jcc0
P1 = invz_jq_path(ion, [2 0 0], struct('dpRng', 10, 'cache', false));
[~, iref] = invz_jq_modes(ion, [0 0 0], struct('dpRng', 10, 'cache', false));
verifyEqual(testCase, P1.Jnu(1,4), iref.Jcc0, 'RelTol', 1e-9);
end
