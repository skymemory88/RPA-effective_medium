function tests = test_invz_twolevel_avg
%TEST_INVZ_TWOLEVEL_AVG T3.2 Gauss-Hermite-dressed doublet (invz_twolevel_avg).
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));   % repo root: MF_dipole, exchange
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_C0_bitwise(testCase)
ion = invz_ion();
tl  = invz_twolevel(ion, 1.6, 2, struct());
tla = invz_twolevel_avg(ion, 1.6, 2, zeros(2), struct());
fn = fieldnames(tl);
for i = 1:numel(fn), verifyEqual(testCase, tla.(fn{i}), tl.(fn{i})); end   % bitwise per shared field
end

function test_gh_weights_and_monotonicity(testCase)
% AMENDED (Task-9 adjudication, 2026-07-17). Level repulsion: Delta_eff grows
% with ||C||. The source plan's "M2_eff decreases monotonically" encodes the
% near-B=0 intuition and is REFUTED at the 2 T fixture: bare M2(Bx) is CONVEX
% along a (M2'' = +0.751/T^2), so the quenched average RAISES M2 there (the
% b-leg alone lowers it; verified == 0.5*M2''*sigma^2). The criticality-
% relevant variable-moments statement is the averaged static two-level
% susceptibility chi(0) = M2_eff * 2*n01_eff/Delta_eff, which IS monotonically
% suppressed — gate that; REPORT the M2 direction. fit_resid: g(0)/tail legs
% are exact by construction; the residual is the sum-rule truncation
% structure, second-order in ||C|| (measured band 7e-6..5e-4) — gate 1e-3.
ion = invz_ion();
[wn, ~, ~] = invz_matsubara(1.6, 40);
s = [1e-5 8e-5 6.4e-4];
[D, M2, X0] = deal(zeros(3,1));
for i = 1:3
    [~, avg] = invz_twolevel_avg(ion, 1.6, 2, s(i)*eye(2), struct('wn', wn));
    D(i) = avg.Delta_eff;  M2(i) = avg.M2_eff;
    X0(i) = avg.M2_eff * 2*avg.n01_eff/avg.Delta_eff;
    verifyLessThan(testCase, avg.fit_resid, 1e-3);
end
tl = invz_twolevel(ion, 1.6, 2, struct());
verifyTrue(testCase, all(diff(D) > 0));
verifyGreaterThanOrEqual(testCase, D(1), tl.Delta - 1e-12);   % Delta_eff >= Delta
verifyTrue(testCase, all(diff(X0) < 0));                      % variable moments suppress chi(0)
verifyLessThan(testCase, X0(1), tl.M2*tl.g0 + 1e-12);
fprintf('T3.2 ray: Delta_eff %s | chi0_2l %s | M2_eff %s (direction REPORTED; convexity note)\n', ...
    mat2str(D.', 6), mat2str(X0.', 6), mat2str(M2.', 6));
end

function test_sumrule_survives_averaging(testCase)
ion = invz_ion();
[wn, ~, ~] = invz_matsubara(1.6, 40);
[~, avg] = invz_twolevel_avg(ion, 1.6, 2, 2e-4*eye(2), struct('wn', wn));
verifyEqual(testCase, avg.sumrule_avg, 1, 'AbsTol', 2e-2);    % truncation-level, not drift
end

function test_zero_field_origin_node_limit(testCase)
% T3.4: B = 0 with C > 0 must not throw invz:degenerateDoublet; origin node
% handled by its h -> 0 limit, off-origin nodes lift the degeneracy.
ion = invz_ion();
[wn, ~, ~] = invz_matsubara(1.55, 40);
[tla, avg] = invz_twolevel_avg(ion, 1.55, 0, 2e-4*eye(2), struct('wn', wn));
verifyTrue(testCase, isfinite(tla.Delta) && tla.Delta > 0);
verifyGreaterThanOrEqual(testCase, avg.n_degenerate, 1);      % the origin node
fprintf('T3.4: B = 0 dressed doublet Delta_eff = %.4g meV (bare degenerate)\n', tla.Delta);
end

function test_avg_mode_comparison_reported(testCase)
ion = invz_ion();
[wn, ~, ~] = invz_matsubara(1.6, 40);
[~, ir] = invz_twolevel_avg(ion, 1.6, 2, 2e-4*eye(2), struct('wn', wn));
[~, ip] = invz_twolevel_avg(ion, 1.6, 2, 2e-4*eye(2), struct('wn', wn, 'avg', 'params'));
fprintf('T3.2 avg-mode Delta_eff: response %.6g vs params %.6g meV\n', ir.Delta_eff, ip.Delta_eff);
verifyTrue(testCase, isfinite(ir.Delta_eff) && isfinite(ip.Delta_eff));
end

function test_gh_node_convergence(testCase)
% 5 vs 7 vs 9 nodes at a realistic C: Delta_eff differences shrink (V4.2 sweep seed).
ion = invz_ion();
[wn, ~, ~] = invz_matsubara(1.6, 40);
D = zeros(3,1);  ns = [5 7 9];
for i = 1:3
    [~, avg] = invz_twolevel_avg(ion, 1.6, 2, 2e-4*eye(2), struct('wn', wn, 'ngh', ns(i)));
    D(i) = avg.Delta_eff;
end
verifyLessThan(testCase, abs(D(3) - D(2)), abs(D(2) - D(1)) + 1e-15);
fprintf('GH 5/7/9: Delta_eff = %.8g / %.8g / %.8g meV\n', D);
end

function test_G0_optin(testCase)
% AMENDED (Task-9 adjudication, 2026-07-17): avg.G0/chi0cc0 are OPT-IN via opts.G0
% (Task 10's consumer path) — opts.wn alone no longer triggers the electronuclear
% node solves. Cheap ngh 3 (9 nodes); G0 without wn is an argument error.
ion = invz_ion();
[wn, ~, ~] = invz_matsubara(1.6, 40);
[~, avg] = invz_twolevel_avg(ion, 1.6, 2, 2e-4*eye(2), struct('wn', wn, 'ngh', 3, 'G0', true));
verifyEqual(testCase, size(avg.G0), [numel(wn) 1]);
verifyTrue(testCase, all(isfinite(avg.G0)));
verifyGreaterThan(testCase, avg.chi0cc0, 0);                  % static cc susceptibility
verifyEqual(testCase, avg.chi0cc0, -avg.G0(1), 'AbsTol', 1e-12);   % chi = -G at wn = 0
verifyError(testCase, ...
    @() invz_twolevel_avg(ion, 1.6, 2, 2e-4*eye(2), struct('avg', 'params', 'G0', true)), ...
    'invz:tlavgArgs');
end
