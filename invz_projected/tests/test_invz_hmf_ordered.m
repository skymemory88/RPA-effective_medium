function tests = test_invz_hmf_ordered
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function o = syn_opts(extra)
o = struct('J0eff', 6.4e-3, 'Jxx0', [], 'hyp', true, 'transverse_mf', 'legacy_x');
ion = invz_ion();  o.Jxx0 = ion.Jxx0;
if nargin, f = fieldnames(extra); for k = 1:numel(f), o.(f{k}) = extra.(f{k}); end, end
end

function test_bare_hook_reproduces_bare_mf_root(testCase)
% GATE 5 (wiring): with force_bare the relation collapses to h0 = hmf, so the root
% must satisfy hmf = J0eff*m(hmf) -- the bare MF fixed point invz_single_ion finds.
% Also pins the fixed-field contract: every node held its imposed hz (P0-1).
ion = invz_ion();  T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
[hstar, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, syn_opts(struct('force_bare', true)));
verifyTrue(testCase, isfinite(hstar) && hstar > 0);
si = invz_single_ion(ion, T, Bx, struct('hyp', true, 'order', true, ...
        'J0z', 6.4e-3, 'Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x'));
verifyEqual(testCase, hstar, si.hz, 'RelTol', 1e-3);        % same bare fixed point
verifyEqual(testCase, prof.r, ones(size(prof.r)), 'AbsTol', 1e-14);
end

function test_onset_coincides_with_pm_crit_zero(testCase)
% GATE 6 (THE closure theorem, SS5): the nonzero root exists at 2.85 T (PM crit
% -0.058) and NOT at 3.30 T (PM crit +0.087) -- the same Bc_1z the PM mass defines.
% Anchors from the stage-1 window probe (re-verify if the synthetic model changes).
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
[h_lo, p_lo] = invz_hmf_ordered(ion, T, [2.85 0 0], Jnu, syn_opts());
[h_hi, ~]    = invz_hmf_ordered(ion, T, [3.30 0 0], Jnu, syn_opts());
verifyTrue(testCase, isfinite(h_lo) && h_lo > 0, 'ordered root missing below Bc_1z');
verifyTrue(testCase, isnan(h_hi), 'spurious ordered root above Bc_1z');
verifyTrue(testCase, all(p_lo.node_conv), 'non-converged H_MF nodes at 2.85 T');
% small-field slope diagnostic ~ crit: negative below Bc_1z (SS5)
verifyLessThan(testCase, p_lo.slope0, 0);
% fluctuation suppression: the 1/z root sits BELOW the bare root, with a smaller moment
si = invz_single_ion(ion, T, [2.85 0 0], struct('hyp', true, 'order', true, ...
        'J0z', 6.4e-3, 'Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x'));
verifyLessThan(testCase, h_lo, si.hz);
verifyLessThan(testCase, p_lo.m_star, si.Jexp(3));
% pole observable at the root: ordered static inverse response is small and positive
% approaching closure from the ordered side (vanishes AT Bc_1z)
verifyTrue(testCase, isfinite(p_lo.D_uni_star));
end

function test_moment_vanishes_continuously_toward_onset(testCase)
% GATE 7: m(B) decreases toward the onset -- the continuity stage 1 lacked.
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
Bs = [2.85 2.95 3.05];  ms = nan(size(Bs));
for k = 1:numel(Bs)
    [h, p] = invz_hmf_ordered(ion, T, [Bs(k) 0 0], Jnu, syn_opts());
    if isfinite(h), ms(k) = p.m_star; end
end
fin = isfinite(ms);
verifyGreaterThanOrEqual(testCase, nnz(fin), 2);
verifyTrue(testCase, all(diff(ms(fin)) < 0));               % monotone decrease toward Bc_1z
end

function test_grid_convergence(testCase)
% GATE 8 (P1-4): the refined root is grid-converged -- doubling nH moves hmf_star by
% less than 1e-2 relative. This is the promised nH-refinement probe, as an actual test.
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
[h1, ~] = invz_hmf_ordered(ion, T, [2.85 0 0], Jnu, syn_opts(struct('nH', 21)));
[h2, ~] = invz_hmf_ordered(ion, T, [2.85 0 0], Jnu, syn_opts(struct('nH', 42)));
verifyEqual(testCase, h1, h2, 'RelTol', 1e-2);
end

function test_pm_limit_reproduces_solve_point(testCase)
% GATE 6b (round-2 P0-A, LOAD-BEARING): with hyp = true, the ordered machinery's
% DEDICATED h = 0 predictor node must land on the ACTUAL PM solver's fixed point --
% K(0), Sigma0, and the inverse mass -- not an isolated two-level reference. This is
% what makes the ordered onset coincide with the code's Bc_1z (SS5). PM-side field.
ion = invz_ion();  T = 0.31;  Bx = [3.30 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
ptpm = invz_solve_point(ion, T, Bx, Jnu, o);
verifyTrue(testCase, ptpm.converged);
[~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, ...
    struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true, 'nH', 5));
verifyEqual(testCase, prof.K0_pm0,     ptpm.K(1),   'RelTol', 1e-4);   % same medium
verifyEqual(testCase, prof.Sigma0_pm0, ptpm.Sigma0, 'RelTol', 1e-4);   % same self-energy
% same inverse mass: the h = 0 predictor slope equals the PM crit (SS5). NOTE: the
% cross response X_zt vanishes at hz = 0, so the path correction is inert here --
% this gate is sector-identity, not path-convention (that is the FD gate below).
verifyEqual(testCase, prof.slope0, ptpm.crit, 'RelTol', 1e-3);
end

function test_path_derivative_gate(testCase)
% Production FD gate (round-3 P0-2, SS4a; round-4 P1-A: ALL THREE transverse-MF modes):
% at a FINITE-moment node the chain-rule G0bare must equal the centered FD of the
% ACTUAL node-path moments under the SAME mode -- 'none' has NO feedback (fb = 0),
% 'legacy_x' the x channel, 'vector_ab' the {x,y} block. Small nH keeps this fast.
ion = invz_ion();  T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
for tmfc = {'none', 'legacy_x', 'vector_ab'}
    tmf = tmfc{1};
    [~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, ...
        syn_opts(struct('transverse_mf', tmf, 'nH', 7)));
    k = find(prof.hgrid > 0.5*max(prof.hgrid), 1);    % a well-ordered node
    hp = prof.hgrid(k);  d = 1e-6;
    sib = struct('hyp', true, 'Jxx0', ion.Jxx0, 'transverse_mf', tmf);
    sp = sib;  sp.hz_fixed = hp + d;  sm = sib;  sm.hz_fixed = hp - d;
    mp = invz_single_ion(ion, T, Bx, sp);  mm = invz_single_ion(ion, T, Bx, sm);
    fd = (mp.Jexp(3) - mm.Jexp(3)) / (2*d);           % path derivative by FD
    verifyEqual(testCase, -prof.G0bare(k), fd, 'RelTol', 1e-3, ...
        sprintf('path-derivative mismatch under transverse_mf = %s', tmf));
end
end

function test_failure_contract(testCase)
% Round-5 P1-A hook: the failure paths must be EXERCISED, not only documented. With an
% absurdly small outer budget the Sigma loop cannot converge, so every node fails and
% the profile must report 'node_failed' with hmf_star = NaN -- never an 'ok' status or
% a root built from unconverged nodes.
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
[h, p] = invz_hmf_ordered(ion, T, [2.85 0 0], Jnu, ...
    syn_opts(struct('max_outer', 1, 'nH', 5)));
verifyTrue(testCase, isnan(h));
verifyEqual(testCase, p.status, 'node_failed');
end

function test_near_boundary_root_not_missed(testCase)
% P1-C / round-3 P0-3 (blind interval, NON-CIRCULAR): refine Bc_1z from the PM
% crit = 0 side (INDEPENDENT of the root detector under test), pick a field just
% below it, and prove the ADAPTIVE EXTENSION actually ran: a deliberately coarse
% hmin_frac must miss the tiny root in its initial grid (n_extend >= 1) and still
% find it after extension. Several profiles -> INVZ_SLOW-gated.
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW=1 to run (multi-profile).');
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
% independent Bc_1z: PM crit interpolated to zero on a small field grid (crit is the
% PM solver's own mass -- NOT invz_hmf_ordered)
Bg = 2.90:0.05:3.20;  cg = nan(size(Bg));
for k = 1:numel(Bg)
    pt = invz_solve_point(ion, T, [Bg(k) 0 0], Jnu, o);
    if pt.converged, cg(k) = pt.crit; end
end
ok = isfinite(cg);
Bc_pm = interp1(cg(ok), Bg(ok), 0, 'linear', 'extrap');
Btest = Bc_pm - 0.01;                                   % just below the PM-refined boundary
% round-4 P2-D (DETERMINISTIC forcing): first obtain a fine-grid REFERENCE root, then
% CONSTRUCT the coarse lower limit to sit provably above it, so n_extend >= 1 is a
% guaranteed consequence of correct adaptivity, not luck.
[h_ref, p_ref] = invz_hmf_ordered(ion, T, [Btest 0 0], Jnu, syn_opts());
verifyTrue(testCase, isfinite(h_ref) && h_ref > 0, 'no reference root just below Bc_pm');
hmax_ref = max(p_ref.hgrid);
frac = min(0.5, 4*h_ref/hmax_ref);                      % initial lower node > 4x the root
[h, p] = invz_hmf_ordered(ion, T, [Btest 0 0], Jnu, ...
    syn_opts(struct('hmin_frac', frac, 'nH', 17)));
verifyLessThan(testCase, p.slope0, 0);                  % ordering predicted (independent)
verifyGreaterThan(testCase, p.hmin_initial, h_ref);     % coarse limit provably above the root
verifyGreaterThanOrEqual(testCase, p.n_extend, 1);      % the adaptive path actually ran
verifyTrue(testCase, isfinite(h) && h > 0 && h < p.hmin_initial, ...
    sprintf('root missed at B=%.3f (Bc_pm=%.3f): h=%.3g, hmin_initial=%.3g, n_extend=%d', ...
            Btest, Bc_pm, h, p.hmin_initial, p.n_extend));
verifyEqual(testCase, h, h_ref, 'RelTol', 5e-2);        % same root as the fine reference
end

function test_hmax_abs_common_cutoff(testCase)
% hmax_abs (stage-2 task 6b, step 1; third SS7 review P1): a positive finite value is
% an EXACT OVERRIDE of the geometric grid's ceiling -- NEVER a min cap (a cap would not
% give a COMMON cutoff across states whose natural bare-hz ceilings differ). Two states
% with DIFFERENT bare hz must report IDENTICAL max(prof.hgrid) under the same hmax_abs.
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
si1 = invz_single_ion(ion, T, [2.85 0 0], struct('hyp', true, 'order', true, ...
        'J0z', 6.4e-3, 'Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x'));
si2 = invz_single_ion(ion, T, [2.95 0 0], struct('hyp', true, 'order', true, ...
        'J0z', 6.4e-3, 'Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x'));
verifyGreaterThan(testCase, abs(si1.hz - si2.hz), 1e-6, ...
    'the two states'' bare hz must actually differ for this test to be meaningful');
hmax_abs = 5.0;
[~, p1] = invz_hmf_ordered(ion, T, [2.85 0 0], Jnu, syn_opts(struct('hmax_abs', hmax_abs, 'nH', 9)));
[~, p2] = invz_hmf_ordered(ion, T, [2.95 0 0], Jnu, syn_opts(struct('hmax_abs', hmax_abs, 'nH', 9)));
verifyEqual(testCase, max(p1.hgrid), hmax_abs, 'AbsTol', 1e-12);
verifyEqual(testCase, max(p2.hgrid), hmax_abs, 'AbsTol', 1e-12);
verifyEqual(testCase, max(p1.hgrid), max(p2.hgrid), 'AbsTol', 1e-12);
% nonpositive/nonfinite hmax_abs is a caller error (invz:hmfOpts), never silently used
verifyError(testCase, ...
    @() invz_hmf_ordered(ion, T, [2.85 0 0], Jnu, syn_opts(struct('hmax_abs', -1))), ...
    'invz:hmfOpts');
verifyError(testCase, ...
    @() invz_hmf_ordered(ion, T, [2.85 0 0], Jnu, syn_opts(struct('hmax_abs', 0))), ...
    'invz:hmfOpts');
verifyError(testCase, ...
    @() invz_hmf_ordered(ion, T, [2.85 0 0], Jnu, syn_opts(struct('hmax_abs', NaN))), ...
    'invz:hmfOpts');
verifyError(testCase, ...
    @() invz_hmf_ordered(ion, T, [2.85 0 0], Jnu, syn_opts(struct('hmax_abs', Inf))), ...
    'invz:hmfOpts');
end
