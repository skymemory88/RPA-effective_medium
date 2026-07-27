function tests = test_invz_strict_identities
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function [ion, T, Bx, Jnu, o] = fx()
ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
           'static_medium', 'strict_1z_dyson_ref');
end

% G1a: dm/dh = -G0bare (J 2.31), panelwise on the profile.
function test_G1a_dm_dh_equals_minus_G0bare(testCase)
[ion, T, Bx, Jnu, o] = fx();
o.nH = 129;
[~, p] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyEqual(testCase, p.status, 'ok');
dm = diff(p.m)./diff(p.hgrid);
gb = 0.5*(p.G0bare(1:end-1) + p.G0bare(2:end));
ok = isfinite(dm) & isfinite(gb);
verifyLessThanOrEqual(testCase, max(abs(dm(ok) + gb(ok))./max(1, abs(gb(ok)))), 1e-6, ...
    'dm/dh = -G0bare must hold to panel order');
end

% G1b: Delta F/Delta h against the trapezoidal average of crit -- i.e. F' = crit (spec §A).
function test_G1b_F_prime_equals_crit(testCase)
[ion, T, Bx, Jnu, o] = fx();
o.nH = 129;
[~, p] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyEqual(testCase, p.status, 'ok');
dF = diff(p.F)./diff(p.hgrid);
cm = 0.5*(p.crit(1:end-1) + p.crit(2:end));
ok = isfinite(dF) & isfinite(cm);
verifyLessThanOrEqual(testCase, max(abs(dF(ok) - cm(ok))./max(1, abs(cm(ok)))), 1e-6);
end

% G1c: dF/dm = crit/chi_path, chi_path = -G0bare.
function test_G1c_dF_dm_equals_crit_over_chi_path(testCase)
[ion, T, Bx, Jnu, o] = fx();
o.nH = 129;
[~, p] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyEqual(testCase, p.status, 'ok');
dFdm = diff(p.F)./diff(p.m);
cm   = 0.5*(p.crit(1:end-1) + p.crit(2:end));
chp  = -0.5*(p.G0bare(1:end-1) + p.G0bare(2:end));
ok = isfinite(dFdm) & isfinite(cm) & chp > 0;
verifyLessThanOrEqual(testCase, ...
    max(abs(dFdm(ok) - cm(ok)./chp(ok))./max(1, abs(cm(ok)./chp(ok)))), 1e-6);
end

% G1d: second-order convergence under nH refinement (prereg §5).
function test_G1d_second_order_convergence(testCase)
[ion, T, Bx, Jnu, o] = fx();
e = nan(1, 3);  nHs = [33 65 129];
for k = 1:3
    ok_ = o;  ok_.nH = nHs(k);
    [~, p] = invz_hmf_ordered(ion, T, Bx, Jnu, ok_);
    verifyEqual(testCase, p.status, 'ok');
    dF = diff(p.F)./diff(p.hgrid);
    cm = 0.5*(p.crit(1:end-1) + p.crit(2:end));
    m_ = isfinite(dF) & isfinite(cm);
    e(k) = max(abs(dF(m_) - cm(m_)));
end
verifyGreaterThanOrEqual(testCase, e(1), e(2));
verifyGreaterThanOrEqual(testCase, e(2), e(3));
verifyLessThanOrEqual(testCase, e(3), 1e-6);
% A separate pure nested-uniform-grid test pins second-order trapezoid convergence. The
% adaptive production grid is not required to show a [3,5] ratio.
end

% G1d companion: second-order is pinned on a genuinely nested smooth quadrature fixture,
% separately from the non-nested adaptive production grids.
function test_G1d_nested_trapezoid_is_second_order(testCase)
n = [17 33 65];  e = nan(size(n));
Iexact = exp(1) - 1;
for k = 1:numel(n)
    x = linspace(0, 1, n(k));
    e(k) = abs(trapz(x, exp(x)) - Iexact);
end
ratio = e(1:end-1)./e(2:end);
verifyGreaterThan(testCase, min(ratio), 3.5);
verifyLessThan(testCase, max(ratio), 4.5);
end

% G2: the two legs coincide at m -> 0, WITHIN the frozen K tolerances (prereg §7.3 -- not
% bitwise: the callers reach Gref through different expressions).
function test_G2_onset_coincidence_at_m_zero(testCase)
[ion, T, Bx, Jnu, o] = fx();
Jscale = max(abs(Jnu));  K_atol = 1e-14;  K_rtol = 1e-12;
ptp = invz_solve_point(ion, T, Bx, Jnu, o);                     % PM leg, m = 0
[~, p] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyTrue(testCase, isfinite(p.slope0) && isfinite(ptp.crit));
% the dimensionless masses agree at onset
verifyEqual(testCase, p.slope0, ptp.crit, 'AbsTol', 1e-6);
% and so do the references and the media
G0pm = -ptp.chi0cc0;
GrefP = invz_static_medium_reference(G0pm, ptp.Sigma0, 'strict_1z_dyson_ref');
GrefO = invz_static_medium_reference(p.G0bare_pm0, p.Sigma0_pm0, 'strict_1z_dyson_ref');
verifyEqual(testCase, GrefO, GrefP, 'RelTol', 1e-6);
mom = invz_coupling_moments(Jnu);
KP = invz_medium_moment_closure(GrefP, mom, 'strict_1z_dyson_ref');
KO = invz_medium_moment_closure(GrefO, mom, 'strict_1z_dyson_ref');
Kgate = K_atol + K_rtol*max([abs(KO), abs(KP), Jscale]);
verifyLessThanOrEqual(testCase, abs(KO - KP), Kgate);
verifyLessThanOrEqual(testCase, abs(p.K0_pm0-ptp.K(1)), Kgate);
end

% G3: the pinned m = 0 identity r = 1 + Sigma0 survives under strict, for ANY K0.
function test_G3_r_equals_one_plus_sigma_under_strict(testCase)
beta = 1/(0.0862*0.31);
tl = struct('Delta', 0.02, 'M2', 0.8, 'm', 0, 'n01', tanh(0.02*beta/2), 'g0', 1);
tl.g0 = 2*tl.n01/tl.Delta;
for K0 = [0, 1e-3, 5e-3, 0.05]
    [~, out] = invz_gstat_ordered(tl, [0.01; 0.02], K0, 0.25, beta, -300, 0);
    verifyEqual(testCase, out.r, 1.25, 'RelTol', 1e-12, sprintf('K0 = %g', K0));
end
end
