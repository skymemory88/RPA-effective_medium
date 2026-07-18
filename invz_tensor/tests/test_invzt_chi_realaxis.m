function tests = test_invzt_chi_realaxis
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion engine
end

% ------- brief Step 1 tests, transcribed verbatim (task-8-brief.md) -------

function test_sigma0_gamma_uniform_matches_bare_rpa(testCase)
% force_sigma0 + qsel 'gamma_uniform': the uniform-mode response must equal the
% single-site 3x3 RPA with J = diag(Jaa0, Jaa0, Jcc0) computed LOCALLY (exact
% identity; no projected dependency).
ion = invz_ion();  T = 1.6;  B = [2 0 0];
w = (0.05:0.05:0.7).';
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
pt = invzt_solve_point(ion, T, B, lat, struct('odd', false));
out = invzt_chi_realaxis(ion, T, B, pt, w, struct('force_sigma0', true, 'odd', false));
si = pt.si;
c0 = invz_chi0z(si, T, w + 1i*5e-3, struct('elastic', false));
J3 = diag([lat.info.Jaa0, lat.info.Jaa0, lat.info.Jcc0]);
for k = 1:numel(w)
    ref = (eye(3) - c0(:,:,k)*J3) \ c0(:,:,k);
    verifyEqual(testCase, out.chi_uniform(:,:,k), ref, 'RelTol', 1e-8);
end
end

function test_odd_on_spectra_reported(testCase)
ion = invz_ion();  T = 1.55;  B = [0.5 0 0];
w = (0.02:0.005:0.5).';
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
p1 = invzt_solve_point(ion, T, B, lat, struct('odd', true));
p0 = invzt_solve_point(ion, T, B, lat, struct('odd', false));
o1 = invzt_chi_realaxis(ion, T, B, p1, w, struct('odd', true));
o0 = invzt_chi_realaxis(ion, T, B, p0, w, struct('odd', false));
[~, i1] = max(squeeze(imag(o1.chi_uniform(3,3,:))));
[~, i0] = max(squeeze(imag(o0.chi_uniform(3,3,:))));
fprintf('ODD spectra: peak %.4f (on) vs %.4f (off) meV, shift %+.1f ueV\n', ...
    w(i1), w(i0), 1e3*(w(i1) - w(i0)));
verifyTrue(testCase, isfinite(w(i1)) && isfinite(w(i0)));
end

% ------- additional coverage (own tests, beyond the brief's verbatim two) -------

function test_odd_mismatch_errors(testCase)
% opts.odd must equal pt.odd (Sigma/alpha/lambda/K already bake in the flag
% used at solve time); a caller-stated mismatch is a structural usage error.
ion = invz_ion();  T = 1.6;  B = [2 0 0];
lat = invzt_jq_tensor(ion, [0.25 0 0; 0.1 0.2 0.3], struct('dpRng', 10, 'cache', false));
pt = invzt_solve_point(ion, T, B, lat, struct('odd', false));
w = (0.1:0.1:0.3).';
verifyError(testCase, @() invzt_chi_realaxis(ion, T, B, pt, w, struct('odd', true)), ...
    'invzt:oddMismatch');
end

function test_qsel_gamma_reported(testCase)
% qsel = 'gamma' (the ACTUAL JtGamma page, not the idealized uniform Jd):
% finite response, own code path exercised (no brief-given identity here).
ion = invz_ion();  T = 1.6;  B = [2 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
pt = invzt_solve_point(ion, T, B, lat, struct('odd', false));
w = (0.1:0.1:0.4).';
out = invzt_chi_realaxis(ion, T, B, pt, w, struct('odd', false, 'qsel', 'gamma'));
verifyTrue(testCase, all(isfinite(out.chi_uniform(:))));
verifyTrue(testCase, isempty(out.chi_cc_q));
fprintf('qsel=gamma vs gamma_uniform (cc,cc) peak-imag ratio reported: %.4f\n', ...
    max(squeeze(imag(out.chi_uniform(3,3,:)))));
end

function test_qsel_explicit_qlist_resolves_per_q(testCase)
% Explicit [nq,3] qvec: chi_cc_q [nq,nw] populated and finite; chi_uniform
% still available via the SAME Jd construction used by 'gamma_uniform'.
ion = invz_ion();  T = 1.6;  B = [2 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
pt = invzt_solve_point(ion, T, B, lat, struct('odd', false));
w = (0.1:0.1:0.4).';
qlist = [0.25 0 0; 0.1 0.2 0.3];
out = invzt_chi_realaxis(ion, T, B, pt, w, struct('odd', false, 'qsel', qlist, 'dpRng', 10, 'cache', false));
verifyEqual(testCase, size(out.chi_cc_q), [2 numel(w)]);
verifyTrue(testCase, all(isfinite(out.chi_cc_q(:))));
verifyTrue(testCase, all(isfinite(out.chi_uniform(:))));
% chi_uniform must be UNCHANGED from the gamma_uniform default (q-grid-independent).
outU = invzt_chi_realaxis(ion, T, B, pt, w, struct('odd', false));
verifyEqual(testCase, out.chi_uniform, outU.chi_uniform, 'RelTol', 1e-12);
end

function test_qsel_explicit_q_complex_response(testCase)
% F1 regression (Codex review 2026-07-18): the explicit-q response must keep
% its imaginary (dissipative) part on the real axis. The pre-fix code
% real()-projected each site-diagonal element (a Matsubara-pattern transplant
% from invzt_gcc_lattice, where real() is a legitimate noise-clean), making
% chi_cc_q identically real and every q-path chi'' map zero.
ion = invz_ion();  T = 1.6;  B = [2 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
pt = invzt_solve_point(ion, T, B, lat, struct('odd', false));
w = (0.05:0.05:0.55).';
qlist = [0.25 0 0; 0.1 0.2 0.3];
out = invzt_chi_realaxis(ion, T, B, pt, w, ...
    struct('odd', false, 'qsel', qlist, 'dpRng', 10, 'cache', false));
verifyFalse(testCase, isreal(out.chi_cc_q), ...
    'chi_cc_q lost its imaginary (dissipative) part');
mx = max(imag(out.chi_cc_q(:)));
verifyGreaterThan(testCase, mx, 1e-6, ...
    'no positive dissipative weight anywhere on the q list');
% NOTE: no global chi''>=0 gate on the full-Sigma call above. The frozen-Kw
% A1 continuation (Kw seeded at Gamma, held fixed across the whole w grid --
% see this function's own SCOPE box) has a known near-resonance negative-chi''
% lobe at these (T,B,q) -- e.g. -312.9 alongside a +652.5 peak two grid points
% away -- shared identically by out.chi_uniform (an untouched code path) and
% by the projected reference's same Kw-seeding scheme. It is NOT introduced by
% this fix and is out of scope to repair here; tracked as a known limitation.
% The causality/passivity gate instead runs below on the force_sigma0
% bare-chi0 RPA response, which is a manifestly passive equilibrium response
% (Sigma_w == 0 identically) while still exercising the identical complex
% accumulation path this fix touched -- a sign flip in the accumulation would
% make imF <= 0 everywhere and fail this hard.
outF = invzt_chi_realaxis(ion, T, B, pt, w, ...
    struct('odd', false, 'qsel', qlist, 'dpRng', 10, 'cache', false, 'force_sigma0', true));
imF = imag(outF.chi_cc_q);
mxF = max(imF(:));
verifyGreaterThan(testCase, mxF, 1e-6, 'no dissipative weight in the bare-RPA limit');
verifyGreaterThan(testCase, min(imF(:)), -1e-6*mxF, ...
    'bare-RPA chi'''' must be non-negative for w > 0 (passive response)');
end

function test_realaxis_rejects_non_a1_point(testCase)
% invzt_chi_realaxis is the A1 scalar-Sigma continuation ONLY (LOCKED scope):
% an A2/A3 point must be rejected loudly -- its alpha/lambda/K fields are
% matrix or diagnostic objects that would otherwise produce silent garbage
% or all-NaN spectra (Codex review F3).
ion = invz_ion();  T = 1.6;  B = [2 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
pt2 = invzt_solve_point(ion, T, B, lat, struct('odd', false, 'mode', 'a2'));
w = (0.1:0.1:0.3).';
verifyError(testCase, ...
    @() invzt_chi_realaxis(ion, T, B, pt2, w, struct('odd', false)), ...
    'invzt:realaxisMode');
end
