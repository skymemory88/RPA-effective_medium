function tests = test_invzt_critical
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion engine
end

function test_tc_extrap_handle_math(testCase)
% Synthetic crit(T): linear through Tc = 1.5 with a non-converged hole.
critfun = @(T) deal(0.4*(T - 1.5), T > 1.52 || T > 1.6);
f = @(T) wrap2(critfun, T);   % local helper returning [c, ok]
Tg = 1.4:1/30:1.9;
tc = invzt_tc_pm_extrap(f, Tg);
verifyEqual(testCase, tc, 1.5, 'AbsTol', 1e-10);
end

function test_bc_finder_structure(testCase)
ion = invz_ion();  T = 1.6;
g = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 10, 'cache', true));
[Bc, out] = invzt_critical(ion, T, lat, [0.5 5], struct('odd', false, 'tol', 0.05));
verifyTrue(testCase, isfinite(Bc) && Bc > 0.5 && Bc < 5);
verifyTrue(testCase, out.iters(end).pt.converged);
end

% ------------------------------------------------------------------------------
function [c, ok] = wrap2(critfun, T)
% Convert the two-output anonymous critfun into the [c, ok] handle contract that
% invzt_tc_pm_extrap consumes (deal returns both outputs when both are requested).
[c, ok] = critfun(T);
end
