function tests = test_invz_ordered_phase
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_ordered_solve_and_soft_mode(testCase)
% SLOW. Ferromagnetic-phase 1/z solve on the real lattice: it converges to a stable (crit>0)
% gapped state below Bc, the spontaneous moment and the soft mode both grow deeper into the
% ordered phase, the mode never softens to zero, and above Bc the ordered solve reports no
% moment (fall back to the paramagnet). cf. Rønnow 2007 Fig 2 / Kovačević 2016 Fig 3d.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
cmd = "qVec_generator(ion.a, 'mode','grid','grid',[16 16 16],'range',[-0.5 0.5])";
[~, qvec] = evalc(cmd);
qvec = qvec(any(abs(qvec) > 1e-12, 2), :);
[Jnu, info] = invz_jq_modes(ion, qvec, struct('dpRng', 30, 'cache', true));
T = 0.31;  w = (0.005:0.005:0.7).';
o = struct('hyp', true, 'J0eff', info.Jcc0);

p2 = invz_solve_point_ordered(ion, T, 2.0, Jnu(:), o);   % deep ordered
p4 = invz_solve_point_ordered(ion, T, 4.0, Jnu(:), o);   % near the boundary
verifyTrue(testCase, p2.is_ordered && p2.converged);
verifyTrue(testCase, p4.is_ordered && p4.converged);
verifyGreaterThan(testCase, p2.crit, 0);                 % stable (gapped) ordered state
verifyGreaterThan(testCase, p4.crit, 0);
verifyGreaterThan(testCase, p2.crit, p4.crit);           % gap closes toward Bc
verifyGreaterThan(testCase, p2.m0, p4.m0);               % moment larger deeper in the FM phase

E2 = soft_peak(ion, T, 2.0, p2, w, info.Jcc0);
E4 = soft_peak(ion, T, 4.0, p4, w, info.Jcc0);
verifyGreaterThan(testCase, E2, E4);                     % mode hardens deeper into ordered phase
verifyGreaterThan(testCase, E4, 0.05);                   % gapped: never softens to zero
verifyLessThan(testCase, E2, 0.6);

pHi = invz_solve_point_ordered(ion, T, 6.0, Jnu(:), o);
verifyFalse(testCase, pHi.is_ordered);                   % paramagnet above the boundary
end

function E = soft_peak(ion, T, B, pt, w, Jsel)
out = invz_chi_realaxis(ion, T, B, pt, w, struct('Jsel', Jsel));
[~, ip] = max(imag(out.chi_cc_q(1, :)));
E = w(ip);
end
