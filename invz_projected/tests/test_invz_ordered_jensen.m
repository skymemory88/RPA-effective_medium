function tests = test_invz_ordered_jensen
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_jensen_mode_root_and_suppression(testCase)
% jensen mode at 2.85 T: ordered, converged, moment suppressed below bare, hmf below
% the bare fixed point, si.hz == hmf (fixed-field contract); at 3.30 T: PM early return.
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true, 'ordered_mode', 'jensen');
pt = invz_solve_point_ordered(ion, T, [2.85 0 0], Jnu, o);
verifyTrue(testCase, pt.is_ordered && pt.converged);
verifyEqual(testCase, pt.ordered_mode, 'jensen');
verifyEqual(testCase, pt.si.hz, pt.hmf, 'AbsTol', 1e-12);      % P0-1 contract on the final state
ptb = invz_solve_point_ordered(ion, T, [2.85 0 0], Jnu, rmfield(o, 'ordered_mode'));
verifyLessThan(testCase, pt.m0, ptb.m0);                       % fluctuation-suppressed moment
verifyLessThan(testCase, pt.hmf, ptb.si.hz);                   % root below the bare fixed point
verifyTrue(testCase, isfinite(pt.D_uni));                      % pole observable exposed
% P1-D (round-3 DISCRIMINATIVE form): the closure and D_uni identities alone cannot
% distinguish the elastic Gstat from the ordinary-Dyson value (both satisfy them), so:
% (a) pt.G(1) must equal the profile's closed Gstat_star;
% (b) pt.G(1) must equal an INDEPENDENT recomputation of Gstat from the final state;
% (c) pt.G(1) must DIFFER from the ordinary-Dyson static value at this finite moment.
verifyEqual(testCase, pt.G(1), pt.hmf_prof.Gstat_star, 'RelTol', 1e-8);
C = invz_const();  beta = 1/(C.kB*T);
c0e = invz_chi0z(pt.si, T, 0, struct('elastic', true));
c0i = invz_chi0z(pt.si, T, 0, struct('elastic', false));
X = real(c0e(:, :, 1));
fb = X(3,1) * (ion.Jxx0 / (1 - ion.Jxx0*X(1,1))) * X(1,3);     % legacy_x: 1x1 block (SS4a)
G0bare0 = -(X(3,3) + fb);  G0inel0 = -real(c0i(3,3,1));
Gs_ind = invz_gstat_ordered(pt.tl, pt.lambda(1:2), pt.K(1), pt.Sigma(1), beta, ...
                            G0inel0, G0bare0 - G0inel0);
verifyEqual(testCase, pt.G(1), Gs_ind, 'RelTol', 1e-8);        % independent recomputation
G_dyson = (-X(3,3)) / (1 + pt.Sigma(1) + pt.K(1)*(-X(3,3)));   % the OLD ordinary value
verifyGreaterThan(testCase, abs(pt.G(1) - G_dyson)/abs(G_dyson), 1e-6);  % must differ
verifyEqual(testCase, pt.D_uni, 1 + (6.4e-3 - pt.K(1))*pt.G(1), 'RelTol', 1e-10);
% round-5 P1-B: the EXPORTED tuple must satisfy the EMT lattice closure -- the formula
% consistency above proves what Gstat(K) is, this proves the exported K actually closes
Gq = pt.G(1) ./ (1 + (Jnu - pt.K(1)).*pt.G(1));
verifyEqual(testCase, mean(Gq), pt.G(1), 'RelTol', 1e-8);
verifyLessThan(testCase, pt.final_resid, 1e-8);                % outer-tol revalidation held
pth = invz_solve_point_ordered(ion, T, [3.30 0 0], Jnu, o);
verifyFalse(testCase, pth.is_ordered);                         % no ordered state above Bc_1z
end

function test_jensen_mode_guards(testCase)
ion = invz_ion();  Jnu = linspace(-2e-3, 6.0e-3, 24).';
base = struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true, 'ordered_mode', 'jensen');
o1 = base;  o1.forced_moment = true;
verifyError(testCase, @() invz_solve_point_ordered(ion, 0.31, [2 0 0.5], Jnu, o1), ...
            'invz:orderedMode');                               % forced_moment forbidden
verifyError(testCase, @() invz_solve_point_ordered(ion, 0.31, [2 0 0.5], Jnu, base), ...
            'invz:orderedMode');                               % longitudinal field forbidden (P1-6)
end
