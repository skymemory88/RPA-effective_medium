function tests = test_invz_transverse_mf_threading
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function fx = fixture()
fx = struct('Jnu', linspace(-2e-3, 6.0e-3, 24).', 'Jcc0', 6.4e-3);
end

function test_solve_point_legacy_default_identical(testCase)
fx = fixture();  ion = invz_ion();
p1 = invz_solve_point(ion, 0.31, [3 0 0], fx.Jnu, struct('J0eff', fx.Jcc0));
p2 = invz_solve_point(ion, 0.31, [3 0 0], fx.Jnu, ...
                      struct('J0eff', fx.Jcc0, 'transverse_mf', 'legacy_x'));
verifyTrue(testCase, isequaln(p1, p2));
end

function test_solve_point_vector_consistent(testCase)
% si and tl must carry the same MF model; hy is live in both.
fx = fixture();  ion = invz_ion();
pt = invz_solve_point(ion, 0.31, [3 1 0], fx.Jnu, ...
                      struct('J0eff', fx.Jcc0, 'transverse_mf', 'vector_ab'));
verifyEqual(testCase, pt.si.transverse_mf, 'vector_ab');
verifyEqual(testCase, pt.tl.transverse_mf, 'vector_ab');
verifyTrue(testCase, abs(pt.si.hy) > 1e-8);
end

function test_solve_auto_inplane_vector(testCase)
% In-plane vector field routes transversely (Bz = 0) and forwards the mode. At T=0.31K,
% Bx=[5 0.5 0] the transverse field is large enough to suppress the FM order (verified against
% the live solver: phase = 2, invz_solve_point's paramagnetic branch, is the actual outcome
% here -- phase is a numeric code (0 = no solution, 1 = ordered accepted, 2 = paramagnetic
% fallback accepted), not the string 'para').
fx = fixture();  ion = invz_ion();
[pt, phase] = invz_solve_auto(ion, 0.31, [5 0.5 0], fx.Jnu, ...
                              struct('J0eff', fx.Jcc0, 'transverse_mf', 'vector_ab'));
verifyEqual(testCase, phase, 2);
verifyEqual(testCase, pt.si.transverse_mf, 'vector_ab');
end

function test_solve_point_c4_axes(testCase)
% cc channel is C4-invariant: a/b-axis solves must agree in vector mode.
fx = fixture();  ion = invz_ion();
o = struct('J0eff', fx.Jcc0, 'transverse_mf', 'vector_ab');
px = invz_solve_point(ion, 0.31, [4 0 0], fx.Jnu, o);
py = invz_solve_point(ion, 0.31, [0 4 0], fx.Jnu, o);
verifyTrue(testCase, px.converged && py.converged);
verifyEqual(testCase, py.Sigma0, px.Sigma0, 'AbsTol', 1e-9);
verifyEqual(testCase, py.crit,   px.crit,   'AbsTol', 1e-9);
end

function test_chi_tensor_ref_vector_mode(testCase)
ion = invz_ion();  w = (0:0.01:0.6).';
o = struct('Jsel', 6.4e-3, 'Jaa0', 3.5e-3, 'eta', 0.02, 'transverse_mf', 'vector_ab');
R = invz_chi_tensor_ref(ion, 0.1, 3*[cosd(20) sind(20) 0], w, o);
verifyTrue(testCase, isfinite(R.eps_spec) && isfinite(R.Epeak_ten));
end

function test_twolevel_ordered_vector(testCase)
ion = invz_ion();
tl = invz_twolevel_ordered(ion, 0.31, [2 1 0], 0.01, ...
                           struct('transverse_mf', 'vector_ab'));
verifyEqual(testCase, tl.transverse_mf, 'vector_ab');
end
