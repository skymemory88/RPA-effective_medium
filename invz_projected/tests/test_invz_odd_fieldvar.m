function tests = test_invz_odd_fieldvar
%TEST_INVZ_ODD_FIELDVAR T3.1 Tier-2 internal-field covariance C (E3).
% Structure (real symmetric PSD, tail control, two-basis cross-check at one
% generic q), C4 at Bx -> 0 (C_aa = C_bb), gated static approximation with the
% full-vs-static difference logged once (plan T3.1: never silent).
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));   % repo root: MF_dipole, exchange
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_fieldvar_structure(testCase)
ion = invz_ion();  T = 1.8;   % PM-guaranteed (> no-ODD Tc(0) = 1.743 K)
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 10, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
pt = invz_solve_point(ion, T, 0.5, [], o);
assumeTrue(testCase, pt.converged);
[C, info] = invz_odd_fieldvar(ion, pt, S, T, struct());
verifyEqual(testCase, C, C.', 'AbsTol', 1e-15*max(abs(C(:)) + 1e-30));
e = eig(C);
verifyGreaterThan(testCase, min(e), -1e-14*max(e));       % PSD
verifyGreaterThan(testCase, max(e), 0);                    % nonzero variance
verifyLessThan(testCase, info.tail_share, 0.05);           % w_n^-2 tail under control
fprintf('T3.1: C_aa = %.4g meV^2, heq = %.3f T (Dollberg ~0.4 T), tail %.2e\n', ...
    C(1,1), info.heq_T, info.tail_share);

% --- Two-basis cross-check (controller directive: one q, both bases): the
% mode-basis contraction the function uses (y_alpha = Vc_alpha' * u_nu column
% norms / overlaps, weighted by S_nu) must equal the brute-force SUBLATTICE-
% basis contraction tr[V_ac(q) * Scc(q) * V_cb(q)] with Scc(q) formed
% explicitly. Same converged state: dJ rebuilt from pt.odd.Xp is
% deterministic, so both sides diagonalize the identical Vcc + dJ. Checked at
% an on-axis q AND a generic low-symmetry q, AbsTol 1e-14 * scale. Measured
% fact (T3.1): the per-q integrand is REAL at every q of this mesh (inversion
% symmetry; max rel imag ~1e-16), so the +/-q cancellation the function
% asserts is trivially satisfied here and the off-diagonal (per-q rel ~0.25)
% cancels to ~0 in the q-sum.
[~, id] = invz_odd_fieldvar(ion, pt, S, T, struct('debug', true));
dJ = invz_odd_deltaJ(S.Vca, S.Vcb, pt.odd.Xp);
[~, wts, beta] = invz_matsubara(T, 40);
iqs = [1, find(all(qvec > 0, 2) & qvec(:,1) ~= qvec(:,2), 1)];   % on-axis + generic
for iq = iqs
    M = S.Vcc(:,:,iq) + dJ(:,:,iq);  M = (M + M')/2;
    [U, ev] = eig(M, 'vector');
    [ev, isrt] = sort(real(ev));  U = U(:, isrt);
    Gq = pt.G ./ (1 + (ev.' - pt.K).*pt.G);                % [nwn, 4] EMT lattice propagator
    Snu = -(wts.' * Gq)/beta;                              % equal-time S_nu(q) >= 0 (sum-rule sign)
    Scc = U * diag(Snu) * U';                              % mode -> sublattice basis
    Wa = S.Vca(:,:,iq);  Wb = S.Vcb(:,:,iq);
    Cq_bf = [trace(Wa'*Scc*Wa), trace(Wa'*Scc*Wb); ...
             trace(Wb'*Scc*Wa), trace(Wb'*Scc*Wb)];        % E3 integrand, V_ac = Vca' etc.
    verifyEqual(testCase, id.Cq(:,:,iq), Cq_bf, 'AbsTol', 1e-14*max(abs(Cq_bf(:))));
end
end

function test_fieldvar_aa_bb_at_zero_Bx(testCase)
% C_aa = C_bb at Bx -> 0 (C4). Small field + PM-guaranteed temperature.
ion = invz_ion();  T = 1.8;
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 10, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
pt = invz_solve_point(ion, T, 0.05, [], o);
assumeTrue(testCase, pt.converged);
C = invz_odd_fieldvar(ion, pt, S, T, struct());
verifyEqual(testCase, C(1,1), C(2,2), 'RelTol', 5e-2);
end

function test_fieldvar_static_approx_logged(testCase)
% The kB*T*chi(0) shortcut is gated and its error is measured once (plan T3.1).
ion = invz_ion();  T = 1.8;
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 10, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
pt = invz_solve_point(ion, T, 0.5, [], o);
assumeTrue(testCase, pt.converged);
[Cf, ~]  = invz_odd_fieldvar(ion, pt, S, T, struct());
[Cs, is] = invz_odd_fieldvar(ion, pt, S, T, struct('static_approx', true));
verifyTrue(testCase, is.static_approx);
fprintf('T3.1 static-approx rel diff: %.3f\n', norm(Cs - Cf)/norm(Cf));
verifyTrue(testCase, isfinite(norm(Cs - Cf)/norm(Cf)));
end
