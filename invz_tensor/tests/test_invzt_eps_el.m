function tests = test_invzt_eps_el
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_eps_el_is_elastic_only(tc)
% Constraint 7: eps_el must use the ELASTIC-ONLY static cc weight c_d, not the
% full equal-time variance JzJz_fluct (review P1-4: the old formula was an
% upper bound). The three-state toy has zero diagonal Mz in the doublet ->
% c_d ~ 0 while JzJz_fluct = O(M^2) > 0: the two definitions differ by orders
% of magnitude, so the old formula CANNOT pass this test.
ion = invz_ion();  T = 2.0;  Bx = 0.5;
g   = invzt_qgrid(8, 'legacy_inclusive');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
pt  = invzt_solve_point(ion, T, [Bx 0 0], lat, struct('mode', 'a3', 'nlevels', 'three'));
verifyTrue(tc, pt.converged);
[wn, ~, beta] = invz_matsubara(T, 40);
c0_el   = invz_chi0z(pt.si, T, 1i*wn(1), struct('elastic', true));
c0_inel = invz_chi0z(pt.si, T, 1i*wn(1), struct('elastic', false));
c_d_ref = real(c0_el(3,3,1) - c0_inel(3,3,1)) / beta;
verifyEqual(tc, pt.c_d, c_d_ref, 'AbsTol', 1e-12, 'RelTol', 1e-10);
verifyEqual(tc, pt.eps_el, beta*abs(pt.K(1))*c_d_ref, 'AbsTol', 1e-12, 'RelTol', 1e-10);
verifyLessThan(tc, pt.c_d, 0.5*pt.si.JzJz_fluct);   % materially below the old upper bound
fprintf('eps_el fix: c_d=%.4g vs JzJz_fluct=%.4g (old formula overcounted %.3gx)\n', ...
    pt.c_d, pt.si.JzJz_fluct, pt.si.JzJz_fluct/max(abs(pt.c_d), 1e-300));
end
