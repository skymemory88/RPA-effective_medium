function tests = test_invzt_solve_point
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion engine
end

function test_zero_field_guard(testCase)
ion = invz_ion();
lat = invzt_jq_tensor(ion, [0.25 0 0], struct('dpRng', 10, 'cache', false));
verifyError(testCase, @() invzt_solve_point(ion, 1.6, [0 0 0], lat, struct()), ...
    'invzt:a1ZeroField');
end

function test_odd_on_converges_crit_direction_reported(testCase)
% v3 (review Other 8): the EXACT PSD Schur relation (dJpre = xperp*(Vca*Vca' +
% Vcb*Vcb') >= 0) is the HARD gate, and it lives at FIXED chi0 in Task 4
% (test_schur_complement_equals_E1_direct / test_full_schur_enhancement_reported).
% Here the FULL odd-on solve also moves the self-consistent medium and scalar
% Sigma, so monotonicity of the final pt.crit does NOT follow from the Gaussian
% Schur identity alone. Hard-assert convergence + sum rule; MEASURE the crit
% direction (expected p1.crit < p0.crit; a positive shift is a finding to
% investigate, not an automatic pass/fail here).
ion = invz_ion();  T = 1.6;  B = [0.1 0 0];
g = invzt_qgrid(8, 'halfopen');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
p1 = invzt_solve_point(ion, T, B, lat, struct('odd', true));
p0 = invzt_solve_point(ion, T, B, lat, struct('odd', false));
verifyTrue(testCase, p1.converged && p0.converged);
verifyLessThan(testCase, p1.sumrule_rel, 0.10);
verifyTrue(testCase, isfinite(p1.crit) && isfinite(p0.crit));
fprintf('A1 odd on/off (crit direction REPORTED): crit %.5f / %.5f (d=%+.2e), Sigma0 %.4f / %.4f, diag4 spread %.2e\n', ...
    p1.crit, p0.crit, p1.crit - p0.crit, p1.Sigma0, p0.Sigma0, p1.diag4_spread);
end

function test_chi_rest_toggle_reported(testCase)
ion = invz_ion();  T = 1.6;  B = [0.5 0 0];
g = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 10, 'cache', true));
pa = invzt_solve_point(ion, T, B, lat, struct('chi_rest', true));
pb = invzt_solve_point(ion, T, B, lat, struct('chi_rest', false));
verifyTrue(testCase, pa.converged && pb.converged);
fprintf('chi_rest on/off: crit %.5f / %.5f, Sigma0 %.4f / %.4f\n', pa.crit, pb.crit, pa.Sigma0, pb.Sigma0);
end

function test_no_odd_vs_frozen_projected_baseline(testCase)
% PHYSICAL no-ODD comparison against the FROZEN baseline (core suite: no
% projected path). Category: physical full-tensor no-ODD mode; residual sources:
% cross-Cartesian chi0 elements, dominant-mask vs whole-cc division. Gate 2e-3.
fx = jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')), 'fixtures', 'projected_baseline.json')));
% NOTE (Task-6 impl): the brief's test code referenced fx.point, but the
% committed Task-2 fixture (invzt_projected_baseline_v1) names this entry
% solve_point_1p6K_0p5T. Same frozen numbers, same assertions/tolerances; only
% the JSON key is reconciled here (flagged in the Task-6 report).
ref = fx.solve_point_1p6K_0p5T;
ion = invz_ion();
g = invzt_qgrid(8, 'legacy_inclusive');                  % SAME convention as the baseline
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
pt = invzt_solve_point(ion, ref.T, [ref.Bx 0 0], lat, ...
    struct('odd', false, 'chi_rest', false));
verifyTrue(testCase, pt.converged);
verifyEqual(testCase, pt.Sigma0, ref.Sigma0, 'AbsTol', 2e-3);
verifyEqual(testCase, sign(pt.crit), sign(ref.crit));
fprintf('A1 no-ODD vs frozen projected: dSigma0 = %.3e\n', pt.Sigma0 - ref.Sigma0);
end
