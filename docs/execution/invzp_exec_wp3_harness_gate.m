function R = invzp_exec_wp3_harness_gate(save_path)
%INVZP_EXEC_WP3_HARNESS_GATE Gate for the exact-cluster coefficient harness (fix doc Sec. 12).
%
% Execution packet S8 of docs/execution/invzp_plan_execution_diary.md.
%
%   W1 EXTRACTOR SELF-TEST. invzf_taylor_extract on a function whose Taylor
%      coefficients are known in closed form, exp(a*J) with a_k = a^k/k!. Establishes
%      that the comparator's own error bars are real before any physics is graded.
%
%   W2 PARITY FIXTURE (h = 0). The exact two-site free energy must be EVEN in J,
%      because flipping X on one site is a local unitary that maps J -> -J at zero
%      source. Every odd extracted coefficient must vanish inside its own error bar.
%      This tests the extractor's ability to certify a coefficient as zero, which is
%      the harder half of a coefficient gate.
%
%   W3 MANIFEST vs EXACT. The local-data-only manifest (invzf_pair_local_manifest)
%      is graded against the exact cluster at both a zero-field and a finite-field
%      fixture. Orders 0-2 are derived at any field; order 3 is derived only at
%      h = 0; order 4 is not derived and must be reported ungraded, not passed.
%
%   W4 NEGATIVE CONTROL. A mean-field pair functional -- self-consistent
%      F_MF(J) = 2 f_loc(h + J m) + J m^2 -- is graded by the SAME extractor. It
%      reproduces orders 0 and 1 and the "C1C1|C2" part of order 2, and must FAIL at
%      order 2 by exactly the single-ring term -sum_n C_n^2/(2 beta). A harness that
%      cannot reject a known-incomplete candidate is worthless, so this failure is a
%      required outcome, and its SIZE is checked against the predicted ring term.
if nargin < 1, save_path = ''; end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);
addpath(fullfile(root, 'invz_functional'));
addpath(fullfile(root, 'invz_common'));

fprintf('=== WP3 exact-cluster coefficient harness gate ===\n\n');

% ---------------- W1 : extractor self-test -------------------------------------------------
aexp = 0.7;  K = 4;
sx = invzf_taylor_extract(@(J) exp(aexp*J), K, struct('r0', 0.05));
kfac = factorial((0:K).');
truth = aexp.^(0:K).' ./ kfac;
w1_rel = abs(sx.a - truth)./max(abs(truth), realmin);
% Two separate claims are tested, because they are separate claims:
%  (i) the orders the manifest actually grades (0-2) are extracted to <= 1e-10 relative;
% (ii) the REPORTED error bar is a usable bound at every order, i.e. the true error never
%      exceeds it by more than a small factor. The extraction error grows like r^-k, so
%      order 4 is genuinely only ~1e-6 accurate at r0 = 0.05; the requirement is that the
%      harness KNOWS that, not that it be accurate there.
w1_bar_ratio = abs(sx.a - truth)./max(sx.err, realmin);
w1_low_ok = all(w1_rel(1:3) < 1e-10);
w1_bar_ok = all(w1_bar_ratio < 10);
w1_pass = w1_low_ok && w1_bar_ok;
fprintf('W1 extractor self-test on exp(%.2f*J), K = %d:\n', aexp, K);
for k = 0:K
    fprintf('   a_%d: extracted %+.15g  exact %+.15g  rel %.2e  err_bar %.2e  |err|/bar %.2f\n', ...
        k, sx.a(k+1), truth(k+1), w1_rel(k+1), sx.err(k+1), w1_bar_ratio(k+1));
end
fprintf('   orders 0-2 within 1e-10: %s;  error bar honest (ratio < 10) at every order: %s\n', ...
    pf(w1_low_ok), pf(w1_bar_ok));
fprintf('   fit residual %.2e, worst Vandermonde cond %.2e, %d evaluations -> %s\n\n', ...
    sx.fit_resid, sx.cond, sx.n_evals, pf(w1_pass));

% ---------------- fixtures -----------------------------------------------------------------
fx0 = struct('Delta', 0.30, 'M', 1.00, 'h', 0.00, 'beta', 20.0);   % zero source
fxh = struct('Delta', 0.30, 'M', 1.00, 'h', 0.05, 'beta', 20.0);   % finite source
exo = struct('r0', 0.04, 'n_radii', 3, 'pad', 4);

% ---------------- W2 : parity at h = 0 -----------------------------------------------------
H0 = invzf_cluster_coeff_harness(fx0, struct('K', 5, 'extract', exo));
odd = 2:2:6;                                    % indices of orders 1,3,5
odd_a = H0.exact.a(odd);  odd_e = H0.exact.err(odd);
a2mag = abs(H0.exact.a(3));
w2_pass = all(abs(odd_a) <= max(odd_e, 1e-10*a2mag));
fprintf('W2 parity fixture (h = 0, Delta = %.2f, beta = %.1f): odd coefficients must vanish\n', ...
    fx0.Delta, fx0.beta);
for j = 1:numel(odd)
    fprintf('   a_%d = %+.6e   err_bar %.2e   |a| / |a_2| = %.2e\n', ...
        odd(j)-1, odd_a(j), odd_e(j), abs(odd_a(j))/a2mag);
end
fprintf('   even reference a_2 = %+.15g  -> %s\n\n', H0.exact.a(3), pf(w2_pass));

% ---------------- W3 : manifest vs exact ---------------------------------------------------
Hh = invzf_cluster_coeff_harness(fxh, struct('K', 4, 'extract', exo));
fprintf('W3 manifest (local data only) vs exact cluster\n');
w3_pass = true;
for tag = {'h = 0', 'h = 0.05'}
    if strcmp(tag{1}, 'h = 0'), Hx = H0; else, Hx = Hh; end
    fprintf('  --- %s ---   m = %+.10g   C_0 = %+.10g   sum C_n^2 = %+.10g (trunc err %.1e)\n', ...
        tag{1}, Hx.manifest.m, Hx.manifest.C0, Hx.manifest.sumC2, Hx.manifest.sumC2_trunc_err);
    for i = 1:min(5, numel(Hx.comparison))
        c = Hx.comparison(i);
        if strcmp(c.verdict, 'ungraded')
            fprintf('   order %d: %-9s exact %+.15g  (%s)\n', c.order, c.verdict, c.exact, c.reason);
        else
            fprintf('   order %d: %-9s exact %+.15g  manifest %+.15g  |diff| %.3e  band %.3e\n', ...
                c.order, c.verdict, c.exact, c.manifest, c.abs_diff, c.band);
        end
        if strcmp(c.verdict, 'MISMATCH'), w3_pass = false; end
    end
end
fprintf('  -> %s (a MISMATCH here would reject the manifest, not the cluster)\n\n', pf(w3_pass));

% ---------------- W4 : negative control, mean-field pair -----------------------------------
mf = @(J) mf_pair_free_energy(fxh.Delta, fxh.M, fxh.h, fxh.beta, J);
HC = invzf_cluster_coeff_harness(fxh, struct('K', 3, 'extract', exo, ...
    'candidate', mf, 'candidate_name', 'mean-field pair'));
ring_pred = -Hh.manifest.sumC2/(2*fxh.beta);
c2 = HC.candidate(3);
w4_orders_ok = strcmp(HC.candidate(1).verdict, 'match') && strcmp(HC.candidate(2).verdict, 'match');
w4_rejects   = strcmp(c2.verdict, 'MISMATCH');
ring_rel = abs((c2.candidate - c2.exact) - (-ring_pred))/abs(ring_pred);
w4_pass = w4_orders_ok && w4_rejects && ring_rel < 1e-6;
fprintf('W4 negative control: self-consistent mean-field pair at h = %.2f\n', fxh.h);
for i = 1:numel(HC.candidate)
    c = HC.candidate(i);
    fprintf('   order %d: %-9s exact %+.15g  candidate %+.15g  |diff| %.6e\n', ...
        c.order, c.verdict, c.exact, c.candidate, c.abs_diff);
end
fprintf('   predicted missing single-ring term  -sum C_n^2/(2 beta) = %+.15g\n', ring_pred);
fprintf('   measured order-2 shortfall                              = %+.15g\n', c2.candidate - c2.exact);
fprintf('   relative agreement of the shortfall with the ring term  = %.3e -> %s\n\n', ring_rel, pf(w4_pass));

allpass = w1_pass && w2_pass && w3_pass && w4_pass;
fprintf('OVERALL: %s   (W1 %s, W2 %s, W3 %s, W4 %s)\n', ...
    pf(allpass), pf(w1_pass), pf(w2_pass), pf(w3_pass), pf(w4_pass));

R = struct('W1', struct('pass', w1_pass, 'rel', w1_rel, 'extract', sx), ...
    'W2', struct('pass', w2_pass, 'odd_a', odd_a, 'odd_err', odd_e, 'harness', H0), ...
    'W3', struct('pass', w3_pass, 'harness_h0', H0, 'harness_h', Hh), ...
    'W4', struct('pass', w4_rejects && w4_orders_ok, 'ring_pred', ring_pred, ...
        'shortfall', c2.candidate - c2.exact, 'ring_rel', ring_rel, 'harness', HC), ...
    'all_pass', allpass);
if ~isempty(save_path), save(save_path, 'R'); fprintf('saved: %s\n', save_path); end
end

function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end

function F = mf_pair_free_energy(Delta, M, h, beta, J)
%MF_PAIR_FREE_ENERGY Self-consistent mean-field free energy of the two-site bond.
% Decoupling X1*X2 -> m*X2 + X1*m - m^2 gives H_MF = sum_i Hloc(h + J*m) + J*m^2 with
% m = <X> solved self-consistently. Deliberately INCOMPLETE: it retains the
% moment-moment and C1C1|C2 structures but no ring, which is what W4 detects.
m = local_moment(Delta, M, h, beta);
for it = 1:200
    mnew = local_moment(Delta, M, h + J*m, beta);
    if abs(mnew - m) < 1e-15, m = mnew; break; end
    m = m + 0.7*(mnew - m);
end
loc = invzf_cluster_exact(Delta, M, h + J*m, 0, beta, 0, 1);
F = 2*loc.F + J*m^2;
end

function m = local_moment(Delta, M, h, beta)
loc = invzf_cluster_exact(Delta, M, h, 0, beta, 0, 1);
m = loc.m_site(1);
end
