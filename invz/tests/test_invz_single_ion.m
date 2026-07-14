function tests = test_invz_single_ion
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
end

function test_zero_field_electronic_levels(testCase)
% Rønnow 2007 Table II: ground doublet, first excited singlet at 11 K.
ion = invz_ion(); C = invz_const();
si = invz_single_ion(ion, 1.0, [0 0 0], struct('hyp', false));
verifyEqual(testCase, numel(si.E), 17);
verifyLessThan(testCase, si.E(2) - si.E(1), 1e-9);                  % degenerate doublet
gap_K = (si.E(3) - si.E(1)) / C.kB;
verifyGreaterThan(testCase, gap_K, 9);  verifyLessThan(testCase, gap_K, 13);
% Ising moment of the doublet: eigvals of the 2x2 Jz block = ±<Jz>, g|| = 2*gL*<Jz> = 13.78
blk = si.Mz(1:2,1:2);
mu  = sort(abs(eig(blk)));
verifyEqual(testCase, mu(2), 5.51, 'RelTol', 0.10);
end

function test_hyperfine_manifold(testCase)
% 136 states; lowest 16 (doublet x nuclear octet) spread ~ 2*A*I*<Jz> ≈ 0.13 meV (1.5 K)
ion = invz_ion();
si = invz_single_ion(ion, 0.5, [0 0 0], struct('hyp', true));
verifyEqual(testCase, numel(si.E), 136);
spread = si.E(16) - si.E(1);
verifyGreaterThan(testCase, spread, 0.10);  verifyLessThan(testCase, spread, 0.16);
verifyGreaterThan(testCase, si.E(17) - si.E(16), 0.3);   % gap to next manifold (11 K = 0.95 meV scale)
end

function test_transverse_field_splitting_and_para(testCase)
ion = invz_ion();
d = zeros(1,3); bx = [1 3 6];
for k = 1:3
    si = invz_single_ion(ion, 0.31, [bx(k) 0 0], struct('hyp', false));
    d(k) = si.E(2) - si.E(1);
    verifyLessThan(testCase, abs(si.Jexp(3)), 1e-8);      % paramagnetic: <Jz>=0
    verifyLessThan(testCase, abs(si.Mz(1,1)), 1e-6);      % m ≈ 0
end
verifyGreaterThan(testCase, d(2), d(1));  verifyGreaterThan(testCase, d(3), d(2));
verifyGreaterThan(testCase, d(3), 0.1);   verifyLessThan(testCase, d(3), 1.0);
% mean field converged: hx nonzero under transverse field
si = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false));
verifyGreaterThan(testCase, abs(si.hx), 1e-6);
end

function test_ordered_mean_field(testCase)
% Ordered (opts.order) mode adds the longitudinal ordering mean field and finds the
% spontaneous moment <Jz> = m0 below the boundary; it relaxes to ~0 above it, and the
% default (paramagnetic) path is unchanged.
ion = invz_ion();  T = 0.31;
% deep in the ordered phase: sizable spontaneous moment
si_lo = invz_single_ion(ion, T, [2 0 0], struct('hyp', false, 'order', true));
verifyGreaterThan(testCase, si_lo.Jexp(3), 1.0);
verifyGreaterThan(testCase, abs(si_lo.hz), 1e-4);              % longitudinal MF switched on
% well above the boundary: no spontaneous moment
si_hi = invz_single_ion(ion, T, [6 0 0], struct('hyp', false, 'order', true));
verifyLessThan(testCase, abs(si_hi.Jexp(3)), 0.1);
% paramagnetic path (order off) is unaffected: <Jz> stays zero even at low field
si_p = invz_single_ion(ion, T, [2 0 0], struct('hyp', false));
verifyLessThan(testCase, abs(si_p.Jexp(3)), 1e-8);
end

function test_hyp_operators_are_electronic(testCase)
% Tasks 3-4 depend on Mx/My/Mz being kron(J_electronic, eye(nuclear)) in the
% eigenbasis, i.e. the response operators act as electronic-J tensor
% nuclear-identity in the product (uncoupled) basis.
ion = invz_ion();
si = invz_single_ion(ion, 0.5, [2 0 0], struct('hyp', true));
oJ = stevens_ops(8);
Mz_prod = si.V * si.Mz * si.V';
verifyLessThan(testCase, max(abs(Mz_prod - kron(oJ.Jz, eye(8))), [], 'all'), 1e-10);
Mx_prod = si.V * si.Mx * si.V';
verifyLessThan(testCase, max(abs(Mx_prod - kron(oJ.Jx, eye(8))), [], 'all'), 1e-10);
end
