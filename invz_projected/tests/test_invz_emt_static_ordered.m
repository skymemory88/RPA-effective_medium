function tests = test_invz_emt_static_ordered
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function tl = ordered_tl(hz)
ion = invz_ion();  T = 0.31;  Bx = [2.85 0 0];
tl = invz_twolevel_ordered(ion, T, Bx, hz, struct('Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x'));
end

function test_m_zero_reproduces_direct_solve(testCase)
% GATE C1 (exactness anchor, round-2 P0-A form): at m = 0 (elastic weight 0) the hybrid
% degenerates to the ordinary Dyson structure ON WHATEVER G0inel0 the caller supplies --
% invz_emt_scalar's direct solve is exact there. Check at BOTH a two-level-scale and a
% full-electronuclear-scale weight (the measured PM value is ~ -251 meV^-1): the closure
% must reproduce the direct solve's K(0) and G(0) for the SAME (G0, Sigma) inputs.
T = 0.31;  C = invz_const();  beta = 1/(C.kB*T);
tl = ordered_tl(1e-9);  tl.m = 0;
Sigma0 = 0.27;  lam = [0.012; 0.004];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
for G0inel = [-tl.M2*tl.g0, -251.4]
    med = invz_emt_scalar(G0inel, Sigma0, Jnu, struct());     % one-frequency direct solve
    [K0, Gs, out] = invz_emt_static_ordered(tl, lam, Sigma0, Jnu, 0, beta, 6.4e-3, ...
                                            G0inel, 0, struct());
    verifyTrue(testCase, out.converged);
    verifyEqual(testCase, K0, med.K(1), 'RelTol', 1e-9);
    verifyEqual(testCase, Gs, med.G(1), 'RelTol', 1e-9);
end
end

function test_closure_residual_is_enforced(testCase)
% GATE C2: at the returned K0 the EMT consistency holds -- the q-average of the
% elastic-inserted G(q,0) equals the site Gstat (the "R eq 8" convention) -- with
% genuinely ordered (m > 0) hybrid weights of full-electronuclear scale.
T = 0.31;  C = invz_const();  beta = 1/(C.kB*T);
tl = ordered_tl(0.15);                                    % genuinely ordered (m > 0)
Sigma0 = 0.2;  lam = [0.01; 0.003];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
G0inel = -230.0;  G0el = -12.0;                           % full-scale test weights
[K0, Gs, out] = invz_emt_static_ordered(tl, lam, Sigma0, Jnu, 0, beta, 6.4e-3, ...
                                        G0inel, G0el, struct());
verifyTrue(testCase, out.converged);
Gq = Gs ./ (1 + (Jnu - K0).*Gs);
verifyEqual(testCase, mean(Gq), Gs, 'RelTol', 1e-9);      % independent recomputation
verifyLessThan(testCase, out.resid, 1e-10);
verifyEqual(testCase, out.D_uni, 1 + (6.4e-3 - K0)*Gs, 'RelTol', 1e-14);
end

function test_bare_limit_closure(testCase)
% GATE C3: with Sigma0 = 0, lam = 0 the closure agrees with the direct solve at the
% m = 0 degenerate limit and stays finite/converged at finite m; opts threading works.
T = 0.31;  C = invz_const();  beta = 1/(C.kB*T);
Jnu = linspace(-2e-3, 6.0e-3, 24).';
tl = ordered_tl(0.15);
[gi, ge] = deal(-tl.M2*tl.g0, -tl.m^2*beta*(1 - tl.n01^2));
[K0m, Gsm, outm] = invz_emt_static_ordered(tl, [0;0], 0, Jnu, 0, beta, 6.4e-3, gi, ge, ...
                                           struct('resid_tol', 1e-12, 'maxit', 400));
verifyTrue(testCase, outm.converged && isfinite(K0m) && isfinite(Gsm));
verifyLessThan(testCase, outm.resid, 1e-11);              % tighter tol was honored
tl0 = ordered_tl(1e-9);  tl0.m = 0;
med = invz_emt_scalar(-tl0.M2*tl0.g0, 0, Jnu, struct());
[K00, ~, ~] = invz_emt_static_ordered(tl0, [0;0], 0, Jnu, 0, beta, 6.4e-3, ...
                                      -tl0.M2*tl0.g0, 0, struct());
verifyEqual(testCase, K00, med.K(1), 'RelTol', 1e-9);
end
