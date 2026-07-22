function tests = test_invz_qcp_closure
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_fm_pm_poles_close_at_bc1z(testCase)
% Diagnosis regression 2, un-waived -- POLE-BASED (P1-5): the ordered static inverse
% response D_ord -> 0 from below and the PM mass crit -> 0 from above must vanish at
% the SAME field, the moment vanishing continuously; plus broadening/mesh refinement
% stability of the near-boundary spectrum. Minutes-scale -> INVZ_SLOW-gated.
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW=1 to run (jensen sweep).');
ion = invz_ion();  T = 0.31;
fields = linspace(2.70, 3.50, 17);                   % brackets Bc_1z ~ 3.0 (synthetic)
w = (0.005:0.005:0.6).';
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
S = invz_spectra_map(ion, T, fields, w, ...
        struct('Jnu', Jnu, 'info', info, 'verbose', false));

fm = S.phase_1z == 1;  pm = S.phase_1z == 2;
verifyGreaterThanOrEqual(testCase, nnz(fm), 3);
verifyGreaterThanOrEqual(testCase, nnz(pm), 3);
step = fields(2) - fields(1);
% (i) the label flip agrees with the PM mass zero-crossing within one grid step
cp = S.crit_pm;  okp = isfinite(cp) & pm;
Bc_crit = interp1(cp(okp), fields(okp), 0, 'linear', 'extrap');
verifyLessThan(testCase, abs(S.Bc_1z - Bc_crit), 2*step);
% (ii) POLE closure from BOTH sides: D_ord (ordered static inverse response) monotone
% decreasing toward the boundary and extrapolating to zero at the same field as crit
okd = isfinite(S.D_ord) & fm;
verifyGreaterThanOrEqual(testCase, nnz(okd), 3);
Dd = S.D_ord(okd);  Bd = fields(okd);
verifyTrue(testCase, all(Dd > 0) && all(diff(Dd) < 0));    % positive, closing from below
B_dn = interp1(Dd, Bd, 0, 'linear', 'extrap');
verifyLessThan(testCase, abs(B_dn - Bc_crit), 2*step);
verifyLessThan(testCase, abs(B_dn - S.Bc_1z), 2*step);
% (iii) the ordered moment vanishes continuously toward the boundary
mo = S.m_1z(fm);
verifyTrue(testCase, all(diff(mo) < 0) && mo(end) < 0.5*mo(1));
% (iv) broadening/mesh refinement via a LOCAL BRANCH TRACKER on the inverse response
% (round-3 P1-6): the soft branch is the LOWEST-frequency LOCAL minimum of
% |1 - J0eff*chi~0(w)|, chi~0 = chi0cc_w./(1+Sigma_w) -- a branch identity, not a
% global argmin (which can hop between electronuclear satellites under refinement).
% The refined minimum must lie in the SAME basin and shift only within tolerance.
% NOTE (documented non-identity): at finite m the real-axis pipeline excludes the
% elastic sector, so the w -> 0 limit of this inverse response equals D_uni only in
% the m -> 0 boundary limit; D_ord/crit above remain the load-bearing closure gates --
% this subtest checks refinement stability of the finite-frequency pole only.
Blast = fields(find(fm, 1, 'last'));
o = struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true, 'ordered_mode', 'jensen');
ptj = invz_solve_point_ordered(ion, T, [Blast 0 0], Jnu, o);
w2 = (0.0025:0.0025:0.6).';
c1 = struct('Jsel', 6.4e-3, 'eta', 5e-3,   'Jxx0', ion.Jxx0, 'Jshape', 0, 'hyp', true, ...
            'transverse_mf', 'legacy_x', 'si', ptj.si);
c2 = c1;  c2.eta = 2.5e-3;
o1 = invz_chi_realaxis(ion, T, [Blast 0 0], ptj, w,  c1);
o2 = invz_chi_realaxis(ion, T, [Blast 0 0], ptj, w2, c2);
dinv1 = abs(1 - 6.4e-3 * o1.chi0cc_w(:) ./ (1 + o1.Sigma_w(:)));
dinv2 = abs(1 - 6.4e-3 * o2.chi0cc_w(:) ./ (1 + o2.Sigma_w(:)));
i1 = find(islocalmin(dinv1), 1, 'first');               % lowest-frequency LOCAL minimum
i2 = find(islocalmin(dinv2), 1, 'first');
verifyTrue(testCase, ~isempty(i1) && ~isempty(i2));
% same-basin check: the refined minimum must fall inside the coarse minimum's basin
% (between the coarse grid's neighboring local maxima around i1)
lm = find(islocalmax(dinv1));
blo = max([w(lm(lm < i1)); 0]);  bhi = min([w(lm(lm > i1)); w(end)]);
verifyTrue(testCase, w2(i2) > blo && w2(i2) < bhi, 'refined minimum left the basin');
verifyLessThan(testCase, abs(w(i1) - w2(i2)), max(4*(w(2)-w(1)), 0.1*max(w(i1), w2(i2))));
% (v) sweep-direction invariance: reversed field order reproduces labels and moments
S2 = invz_spectra_map(ion, T, fliplr(fields), w, ...
        struct('Jnu', Jnu, 'info', info, 'verbose', false));
verifyEqual(testCase, fliplr(S2.phase_1z), S.phase_1z);
verifyEqual(testCase, fliplr(S2.m_1z), S.m_1z, 'RelTol', 1e-6, 'AbsTol', 1e-9);
end
