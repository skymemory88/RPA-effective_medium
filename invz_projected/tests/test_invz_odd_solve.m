function tests = test_invz_odd_solve
%TEST_INVZ_ODD_SOLVE T1.4 ODD wiring: solve_point/_ordered flag + guards,
% T0field numeric-or-handle generalization, critical_T odd anchor guard.
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));   % repo root: MF_dipole, exchange
end

function test_solve_point_flag_off_bitwise(testCase)
% Default-opts call must take the identical code path (non-negotiable i).
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';               % synthetic fixture (house convention)
p1 = invz_solve_point(ion, 1.6, 0.5, Jnu, struct('J0eff', 6.4e-3));
p2 = invz_solve_point(ion, 1.6, 0.5, Jnu, struct('J0eff', 6.4e-3, 'odd', false));
verifyTrue(testCase, isequaln(p1, p2));
end

function test_solve_point_odd_args_guard(testCase)
ion = invz_ion();
verifyError(testCase, @() invz_solve_point(ion, 1.6, 0.5, [1e-3; 2e-3], ...
    struct('odd', true)), 'invz:oddArgs');          % blocks missing
S = struct('Vca', zeros(4,4,2), 'Vcb', zeros(4,4,2), 'Vcc', zeros(4,4,2), 'Jcc0', 6.4e-3);
verifyError(testCase, @() invz_solve_point(ion, 1.6, 0.5, [1e-3; 2e-3], ...
    struct('odd', true, 'odd_blocks', S)), 'invz:oddArgs');   % Jnu_flat not []
end

function test_solve_point_odd_pm_point(testCase)
% Hard gate at a GUARANTEED-PM point (1.80 K > no-ODD Tc(0) = 1.743 K): the ODD
% plan's own T1.4 acceptance point (1.60 K, 0.1 T) sits on the ORDERED side of
% the no-ODD baseline, so it cannot be a convergence gate — it is measured and
% REPORTED below instead (whether it now converges is itself the ODD Tc-shift
% signal). Gates: convergence, crit INCREASES with ODD on (sign contract),
% overhead <= 20%.
ion = invz_ion();
n = 8;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 15, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
Jnu0 = invz_odd_modes(Vcc, []);
o0 = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0);
t0 = tic;  p0 = invz_solve_point(ion, 1.80, 0.1, Jnu0(:), o0);            t_off = toc(t0);
o1 = o0;  o1.odd = true;  o1.odd_blocks = S;
t1 = tic;  p1 = invz_solve_point(ion, 1.80, 0.1, [], o1);                 t_on = toc(t1);
verifyTrue(testCase, p0.converged && p1.converged);
verifyGreaterThan(testCase, p1.crit, p0.crit);                            % ODD suppresses ordering
verifyGreaterThan(testCase, p1.odd.d, 0);
verifyLessThan(testCase, p1.sumrule_rel, 0.10);
fprintf('odd overhead: %.2fs vs %.2fs (%.0f%%), crit %.4f -> %.4f, d = %.3g ueV\n', ...
    t_on, t_off, 100*(t_on/t_off - 1), p0.crit, p1.crit, p1.odd.d*1e3);
verifyLessThan(testCase, t_on, 1.2*t_off + 0.5);                          % +0.5 s absolute floor for timer noise
% REPORT the plan's (1.60 K, 0.1 T) point, both flags — no convergence gate:
r0 = invz_solve_point(ion, 1.60, 0.1, Jnu0(:), o0);
r1 = invz_solve_point(ion, 1.60, 0.1, [], o1);
fprintf('plan point (1.60 K, 0.1 T): off converged=%d crit=%.4g; odd converged=%d crit=%.4g\n', ...
    r0.converged, r0.crit, r1.converged, r1.crit);
end

function test_t0field_handles_backcompat(testCase)
% Numeric args byte-identical; constant handles reproduce numerics exactly.
ion = invz_ion();
Tc1 = invz_critical_T0field(ion, 0.30, 6.4e-3);
Tc2 = invz_critical_T0field(ion, @(T) 0.30, @(T) 6.4e-3);
verifyEqual(testCase, Tc2, Tc1);                    % same bisection path
end

function test_critical_T_requires_anchor_with_odd(testCase)
ion = invz_ion();
S = struct('Vca', zeros(4,4,2), 'Vcb', zeros(4,4,2), 'Vcc', zeros(4,4,2), 'Jcc0', 6.4e-3);
verifyError(testCase, @() invz_critical_T(ion, 2.0, [], ...
    struct('odd', true, 'odd_blocks', S)), 'invz:oddTc0');
end
