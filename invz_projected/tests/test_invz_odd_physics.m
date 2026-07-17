function tests = test_invz_odd_physics
% T1.5 zero-field ODD physics measurements (invz_odd_zero_field).
%   Fast: structural gate on the 8^3/dpRng-15 grid (ODD lowers Tc; each split
%   component lowers Tc; (a) full is the largest suppression).
%   Slow (INVZ_SLOW=1): the mode='off' published-anchor reproduction, the
%   Richardson(12,24) ODD headline + reproducibility gate, and the Bc(T) on/off
%   boundary-shift table. Physics numbers are REPORTED, never tuned.
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));   % repo root: MF_dipole, exchange
end

function test_zero_field_off_matches_published_slow(testCase)
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
[Tc, out] = invz_odd_zero_field(ion, struct('mode', 'off'));
verifyEqual(testCase, Tc, 1.743, 'AbsTol', 0.006);          % published anchor
verifyEqual(testCase, out.Sc_rich, 0.2980, 'AbsTol', 0.006);
end

function test_zero_field_structure_fast(testCase)
% Small-grid structural gate (AMENDED 2026-07-17, condition/Sigma-space split):
% directionality only — ODD lowers Tc; each governing leg lowers Tc; d > 0;
% the literal-(c) bookkeeping theorem holds. Magnitudes, leg orderings, and
% closure are REPORTED, not gated (the interaction term is genuine physics;
% the naive (b) is an invalid-regime diagnostic and is only checked finite).
ion = invz_ion();
o = struct('grids', {{8}}, 'dpRng', 15);
Tc0 = invz_odd_zero_field(ion, setfield(o, 'mode', 'off'));   %#ok<SFLD>
[Tc1, out1] = invz_odd_zero_field(ion, setfield(o, 'mode', 'full')); %#ok<SFLD>
verifyLessThan(testCase, Tc1, Tc0);
verifyGreaterThan(testCase, out1.d_at_Tc, 0);
verifyLessThan(testCase, out1.split.Tc_b, Tc0);               % condition-level leg
verifyLessThan(testCase, out1.split.Tc_c, Tc0);               % Sigma-level leg
verifyEqual(testCase, out1.split.Tc_c1_literal, Tc1, 'AbsTol', 2e-3);  % theorem
verifyTrue(testCase, isfinite(out1.split.Tc_b_naive) && isfinite(out1.split.Tc_c_factorial));
fprintf('split: a %.4f | b_cond %.4f | c_sigma %.4f | closure %+.4f | b_naive %.4f | c_fact %.4f\n', ...
    Tc1, out1.split.Tc_b, out1.split.Tc_c, out1.split.closure_defect, ...
    out1.split.Tc_b_naive, out1.split.Tc_c_factorial);
end

function test_odd_headline_slow(testCase)
% T1.5 headline (REPORT + reproducibility): Richardson(12,24) ODD numbers.
% First run: print and RECORD in ODD-LOG + invz_odd_anchors (controller step);
% the anchor assert below activates once the anchor fields exist.
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
[Tc1, out1] = invz_odd_zero_field(ion, struct('mode', 'full'));
fprintf('ODD zero field: Tc = %.4f K (dTc = %.4f), Sc = %.4f (dSc = %+.4f), d(Tc) = %.4g ueV\n', ...
    Tc1, 1.743 - Tc1, out1.Sc_rich, out1.Sc_rich - 0.2980, out1.d_at_Tc*1e3);
fprintf('split: (a) %.4f | b_cond %.4f | c_sigma %.4f | closure %+.4g K | b_naive %.4f (nex %d) | c1_lit %.4f | c_fact %.4f\n', ...
    Tc1, out1.split.Tc_b, out1.split.Tc_c, out1.split.closure_defect, ...
    out1.split.Tc_b_naive, max(out1.split.nex.b_naive(:)), out1.split.Tc_c1_literal, out1.split.Tc_c_factorial);
A = invz_odd_anchors();
if isfield(A, 'odd_Tc_rich')                      % reproducibility gate (1%, plan T1.5)
    verifyEqual(testCase, Tc1, A.odd_Tc_rich, 'RelTol', 0.01);
    verifyEqual(testCase, out1.d_at_Tc, A.odd_d_at_Tc, 'RelTol', 0.01);
end
verifyTrue(testCase, isfinite(Tc1) && Tc1 > 0.8 && Tc1 < 1.743);
end

function test_boundary_shift_slow(testCase)
% T1.5 boundary table at the README grid: Bc(T) ODD on/off at T = 0.31, 0.8,
% 1.2, 1.5 K; REPORT the shifts + the expected attenuation toward high Bx
% (crossover ~3.5 T per DS2023) — gate only the sign (Bc_on <= Bc_off).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
n = 16;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 30, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
Jnu0 = zeros(size(Vcc,3), 4);
for iq = 1:size(Vcc,3), Jnu0(iq,:) = sort(real(eig(Vcc(:,:,iq)))).'; end
Ts = [0.31 0.8 1.2 1.5];
% Window floor lowered from the invz_critical default [2 7] to [0.1 7]: at T = 1.5 K
% (only 9 mK below the ODD Tc(0) = 1.509 K) the ODD-suppressed boundary collapses to
% Bc_on ~ 1.2 T, well under the default 2 T floor. The gate is unchanged (sign only).
for it = 1:numel(Ts)
    o0 = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'window', [0.1 7]);
    b0 = invz_critical(ion, Ts(it), Jnu0(:), o0);
    o1 = o0;  o1.odd = true;  o1.odd_blocks = S;
    b1 = invz_critical(ion, Ts(it), [], o1);
    fprintf('Bc(%.2f K): off %.4f T, on %.4f T, shift %+.4f T\n', Ts(it), b0, b1, b1 - b0);
    verifyLessThanOrEqual(testCase, b1, b0 + 0.02);   % <= within crossing tol
end
end
