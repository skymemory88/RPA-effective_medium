function tests = test_invz_spectra_qpath
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_gamma_approach_regression(testCase)
% Review finding (blocker): the sharply truncated MF_dipole sum collapses on the
% approach to the Gamma-equivalent (2,0,0) (max branch fell to ~0.0016 meV at
% h = 1.999, dpRng 30, vs the correct ~0.0064) and then jumped at the endpoint.
% The direction-aware snap must give a smooth, cutoff-stable approach: for this
% in-plane path (khat_z = 0) the directional limit is the uniform-mode Lorentz
% value, so the endpoint equals info.Jcc0 by construction.
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
% Structural contract + physics anchors of the q-path spectrum:
%  - shapes [nw x nq]; index path coordinate starts at 0;
%  - (2,0,0) IS Gamma-equivalent for the 4-site basis (structure factor 4): the max
%    coupling branch there equals the uniform-mode info.Jcc0 of the same dpRng;
%  - (1,0,0) is NOT Gamma-equivalent (structure factor 0): weaker coupling, so the
%    mode sits HIGHER in energy and the dispersion falls toward the zone centre.
ion = invz_ion();
T = 0.31;  B = 5.5;                        % paramagnetic side: fast, well-converged
w = (0.02:0.02:0.6).';
info = struct('Jcc0', 6.4e-3);             % synthetic medium (as in test_invz_spectra_map)
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
qpath = [1 0 0; 1.5 0 0; 2 0 0];
S = invz_spectra_qpath(ion, T, B, qpath, w, ...
                       struct('Jnu', Jnu, 'info', info, 'dpRng', 10));

verifySize(testCase, S.chiz,   [numel(w) 3]);
verifySize(testCase, S.chirpa, [numel(w) 3]);
verifyEqual(testCase, S.s, [0 0.5 1], 'AbsTol', 1e-12);
verifyEqual(testCase, S.phase, 2);
verifyTrue(testCase, all(isfinite(S.Epeak)));

[~, iref] = invz_jq_modes(ion, [0 0 0], struct('dpRng', 10, 'cache', false));
verifyEqual(testCase, S.Jq(3), iref.Jcc0, 'RelTol', 1e-9);
verifyGreaterThan(testCase, iref.Jcc0, S.Jq(1));         % (1,0,0) couples more weakly
verifyGreaterThan(testCase, S.Epeak(1), S.Epeak(3));     % mode softens toward Gamma
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

function test_qpath_out_of_plane_gamma_limit(testCase)
% The direction-aware guard is ANISOTROPIC: approaching a Gamma-equivalent point
% along c (khat_z = 1, e.g. (0,0,4-)) the nonanalytic broadcast is
% gfac*(4pi/Vc)*(1/3 - 1), so the uniform branch sits 4*gfac*(4pi/Vc) BELOW the
% in-plane limit while the three non-uniform branches are unchanged (a scalar
% sublattice broadcast only moves the uniform mode). Pins the kz2 coefficient
% and its sign -- the in-plane tests alone cannot see them.
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
