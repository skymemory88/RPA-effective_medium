function tests = test_invz_gstat_ordered
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function tl = ordered_tl(hz)
% Ordered two-level params at an imposed molecular field.
ion = invz_ion();  T = 0.31;  Bx = [2.85 0 0];
tl = invz_twolevel_ordered(ion, T, Bx, hz, struct('Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x'));
end

function [gi, ge] = twolevel_weights(tl, beta)
% Two-level static weights: the parametrization in which J 2.28 is EXACT (SS3 check a).
gi = -tl.M2*tl.g0;
ge = -tl.m^2 * beta*(1 - tl.n01^2);
end

function tl = twolevel_fix(h, D0, a, M0, beta)
% CLOSED analytic 2x2 fixture (round-3 P0-1): H = diag([D0/2, -D0/2]) - h*Jz2 with the
% FIXED traceless operator Jz2 = [a M0; M0 -a]. Unlike invz_twolevel_ordered -- whose
% doublet re-embeds in the full CF manifold, adding subspace-drift terms that break the
% J 2.31 identity by up to 31% (review-measured) -- this model is closed, so the
% identity is EXACT. Diagonals of V'*Jz2*V stay +/-m (traceless), so <Jz> = m*n01.
Jz2 = [a, M0; M0, -a];
H = [D0/2, 0; 0, -D0/2] - h*Jz2;
[V, E] = eig((H+H')/2, 'vector');  [E, ix] = sort(real(E));  V = V(:, ix);
p = exp(-beta*(E - E(1)));  p = p/sum(p);
Mz = V'*Jz2*V;
tl = struct('m', real(Mz(1,1)), 'M2', abs(Mz(1,2))^2, 'n01', p(1) - p(2), ...
            'Delta', E(2) - E(1), 'g0', 0);
tl.g0 = 2*tl.n01/tl.Delta;                              % g(0) = 2*n01*Delta/Delta^2
end

function test_bare_limit_and_fd_sign_anchor(testCase)
% GATE 1 (sign anchor, J 2.31; round-3 P0-1): the identity
%   -G0bare = M2*g0 + m^2*beta*(1-n01^2) = d(m*n01)/d(hmf)
% is exact ONLY for a closed two-level model -- built analytically here, NEVER from
% invz_twolevel_ordered (see twolevel_fix header).
T = 0.31;  C = invz_const();  beta = 1/(C.kB*T);
D0 = 0.2;  a = 2.0;  M0 = 3.0;  hz = 0.15;  d = 1e-6;
tlf = @(h) twolevel_fix(h, D0, a, M0, beta);
tl  = tlf(hz);  tlp = tlf(hz + d);  tlm = tlf(hz - d);
[gi, ge] = twolevel_weights(tl, beta);
[Gs, out] = invz_gstat_ordered(tl, [0; 0], 0, 0, beta, gi, ge);
verifyEqual(testCase, out.xi, 1, 'AbsTol', 1e-14);                 % tanh(0)=0, denom 1
verifyEqual(testCase, Gs, out.G0bare, 'RelTol', 1e-14);            % elastic + inelastic bare
fd = (tlp.m*tlp.n01 - tlm.m*tlm.n01) / (2*d);                      % d<Jz>_closed/d hmf
verifyEqual(testCase, -out.G0bare, fd, 'RelTol', 1e-7);            % EXACT for the closed model
verifyEqual(testCase, out.r, 1, 'RelTol', 1e-14);                  % bare integrand ratio
end

function test_m_zero_recovers_resummed_pm_static(testCase)
% GATE 2 (J 2.30 at w=0): at m = 0 the elastic weight vanishes and
% Gstat = G0inel/(1+Sigma0+K0*G0inel);  GATE 3 (K-cancellation, HTML 32):
% Gtil0 = Gstat/(1-K0*Gstat) = G0inel/(1+Sigma0);  hence r = 1+Sigma0.
% Run in BOTH parametrizations: two-level weights AND a generic full-sector weight,
% since the m->0 algebra must hold for any G0inel0 (the hybrid's whole point, SS3).
T = 0.31;  C = invz_const();  beta = 1/(C.kB*T);
tl = ordered_tl(1e-9);  tl.m = 0;                       % PM limit of the two-level params
Sigma0 = 0.31;  K0 = 0.018;  lam = [0.012; 0.004];      % generic nonzero test values
for G0inel = [-tl.M2*tl.g0, -251.4]                     % two-level and full-scale weights
    [Gs, out] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel, 0);
    verifyEqual(testCase, Gs, G0inel/(1 + Sigma0 + K0*G0inel), 'RelTol', 1e-13);
    verifyEqual(testCase, out.Gtil0, G0inel/(1 + Sigma0), 'RelTol', 1e-13);
    verifyEqual(testCase, out.r, 1 + Sigma0, 'RelTol', 1e-13);
end
end

function test_xi_formula_direct(testCase)
% GATE 4: xi transcription check against the closed formula (J 2.29) with
% hand-computed values -- independent arithmetic, not a copy of the source.
T = 0.31;  C = invz_const();  beta = 1/(C.kB*T);
tl = ordered_tl(0.15);
Sigma0 = 0.2;  K0 = 0.02;  lam = [0.01; 0.003];
[gi, ge] = twolevel_weights(tl, beta);
[~, out] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, gi, ge);
num = 1 + tanh(tl.m^2*tl.n01^2*beta*K0 - tl.M2*beta*lam(1));
den = 1 + (4*tl.n01^2*K0*tl.g0 + 2*lam(2) + tl.g0*lam(1))*tl.M2/tl.n01^2;
verifyEqual(testCase, out.xi, num/den, 'RelTol', 1e-14);
verifyGreaterThan(testCase, out.xi, 0);                 % bounded resummation stays positive
end
