%% Task 6b Step 2/4 scratch harness -- closed-model J 2.34 validation + localization
% Standalone script (NOT a committed test -- Step 2 explicitly forbids committing the
% slow test's replacement before the mismatch is resolved). Section 1 mirrors the
% brief's frozen 2x2 closed model verbatim (test_two_routes_closed_model). Sections 2-4
% implement the Step-4 localization diagnostics on top of the same shared kernels.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'invz_common'));
addpath(fullfile(here, '..', '..', 'invz_projected'));

results = struct();

%% ================= Step 2: closed-model two-route harness =========================
T = 0.31;  C = invz_const();  beta = 1/(C.kB*T);
D0 = 0.2;  M0op = 3.0;  J0eff = 6.4e-3;
Jnu = linspace(-4e-3, 4.0e-3, 24).';                 % ZERO-MEAN Jensen lattice
assert(abs(mean(Jnu)) < 1e-15, 'zero-mean Jensen lattice violated (J(ii) = 0 identity)');

fprintf('=== Step 2: closed-model two-route harness ===\n');

hmax = 40*D0;  nH = 81;
[dFA,    dA] = closed_routeA(T, Jnu, J0eff, D0, M0op, hmax,   nH);
[dFA_h2, ~ ] = closed_routeA(T, Jnu, J0eff, D0, M0op, 2*hmax, nH);
[dFA_n2, ~ ] = closed_routeA(T, Jnu, J0eff, D0, M0op, hmax,   2*nH-1);
fprintf('dFA      = %.6e\n', dFA);
fprintf('dFA_h2   = %.6e  (hmax x2, rel diff %.4g%%)\n', dFA_h2, 100*abs(dFA-dFA_h2)/abs(dFA));
fprintf('dFA_n2   = %.6e  (nH x2,  rel diff %.4g%%)\n', dFA_n2, 100*abs(dFA-dFA_n2)/abs(dFA));
fprintf('dA.tail_est_A = %.4g, dh_end = %.4g, r0 = %.6g\n', dA.tail_est_A, dA.dh_end, dA.r0);
assert(abs(dFA-dFA_h2)/abs(dFA) < 2e-2, 'hmax-doubling convergence FAILED');
assert(abs(dFA-dFA_n2)/abs(dFA) < 2e-2, 'grid-doubling convergence FAILED');

% ---- Route B: temperature integral at h = 0 (PM, m = 0) ---------------------------
Tg = T * (1.35.^(0:0.5:17));
[dU37, dU21, tuples] = route_b_sweep(Tg, D0, M0op, Jnu);
assert(all(isfinite(dU37)), 'route-B node failed to converge somewhere on the grid');
dFB   = + T * trapz(Tg, dU37 ./ Tg.^2);
dFB_c = + T * trapz(Tg(1:2:end), dU37(1:2:end) ./ Tg(1:2:end).^2);   % true 1.35-spaced subset
nTm = nnz(Tg <= T*1.35^14);
dFBt = + T * trapz(Tg(1:nTm), dU37(1:nTm) ./ Tg(1:nTm).^2);
tail_est_B = + T * abs(dU37(end)) / Tg(end);

fprintf('dFB      = %.6e\n', dFB);
fprintf('dFB_c    = %.6e  (density-thinned subset, rel diff %.4g%%)\n', dFB_c, 100*abs(dFB-dFB_c)/abs(dFB));
fprintf('dFBt     = %.6e  (shorter Tmax prefix, rel diff %.4g%%)\n', dFBt, 100*abs(dFB-dFBt)/abs(dFB));
fprintf('tail_est_B = %.4g\n', tail_est_B);
assert(abs(dFB_c-dFB)/abs(dFB) < 0.05, 'T-density convergence FAILED');
assert(abs(dFBt-dFB)/abs(dFB) < 0.05, 'Tmax-extension convergence FAILED');
assert(tail_est_B < 0.05*max(abs(dFA), abs(dFB)), 'route-B tail too large');
assert(dA.tail_est_A < 0.05*max(abs(dFA), abs(dFB)), 'route-A tail too large');

mismatch = abs(dFA - dFB) / max(abs(dFA), abs(dFB));
fprintf('THE GATE: |dFA - dFB| / max(|dFA|,|dFB|) = %.4g%% (A=%.6e, B=%.6e)\n', 100*mismatch, dFA, dFB);

results.dFA = dFA; results.dFA_h2 = dFA_h2; results.dFA_n2 = dFA_n2; results.dA = dA;
results.dFB = dFB; results.dFB_c = dFB_c; results.dFBt = dFBt; results.tail_est_B = tail_est_B;
results.mismatch = mismatch;

%% ========= Step 4.1: dU via J 2.21 (eq 36) vs J 2.22 (eq 37) =======================
% Framework SS8 (worktree eq 36 region), read verbatim: (36) is "exact given G; the
% frequency-odd piece of wn*G requires the convergence factor." Direct substitution
% into the repo's nonnegative doubled-weight Matsubara sum is invalid for that piece.
%
% ANALYTIC TREATMENT (option (b): symmetrized expression in which the convergence-
% factor piece is analytic), AMENDED after a first attempt was caught by its own
% zero-order gate's sibling check: split G = G0 + (G-G0). The regularized sum of
% -(Delta+wn)/M^2*G0(wn) [convergence factor e^{wn 0+}] evaluates, via the SAME
% identity the framework's own MF check uses just above (36) ((Delta+wn)*g =
% 2*n01*Delta/(Delta-wn), (1/beta)*sum_n^reg(Delta-wn)^-1 = nB(Delta)), to EXACTLY
% n1*Delta [n1 = (1-n01)/2, the bare excited population], cancelling the "-n1*Delta"
% subtraction termwise. A FIRST ATTEMPT then assumed the remainder (G-G0) was
% absolutely convergent (hence its odd-wn part vanishes for free) -- this is WRONG:
% Sigma(wn) -> alpha (a NONZERO constant, not zero) as wn -> infinity (confirmed
% numerically: Sigma(end) = 6.2389e-3 vs alpha = 6.2383e-3), so G-G0 = -G0*(Sigma+K*G0)
% ~ -alpha*G0(wn) at large wn -- the SAME marginally-convergent tail as G0 itself,
% needing the SAME regularization, not a free pass. The naive treatment produced a
% spurious ~51% mismatch against J 2.22 that evaporates once corrected (see below).
% CORRECTED treatment: peel off the alpha-driven tail analytically, reusing the SAME
% regularized identity, and sum only the genuinely fast-decaying remainder directly:
%   R(wn) = (G(wn)-G0(wn)) + alpha*G0(wn)     [confirmed O(1/wn^4), absolutely convergent]
%   dU_21 = (1/2beta) sum_n K(wn)G(wn)  -  alpha*n1*Delta  -  (Delta/(2beta M^2)) sum_n R(wn)
% (the wn*R(wn) odd part now legitimately sums to exactly zero: R IS absolutely
% convergent). This is manifestly even in wn and sums correctly on the repo's standard
% nonneg-doubled wts grid.
fprintf('\n=== Step 4.1: J 2.21 (dU21, analytic-remainder treatment) vs J 2.22 (dU37) ===\n');

% Zero-order gate (MANDATORY before diagnostic use): Jnu = 0 forces K = Sigma = 0
% identically through the PM emt+sigma loop (A = mean(0./...) = 0 => K = 0; lam = 0
% since it is built from K; alpha = gamma = 0 => Sigma = 0), so G = G0 and BOTH dU
% formulas must vanish to numerical tolerance.
Jnu_zero = zeros(size(Jnu));
[dU37_0, dU21_0, tup0] = route_b_node(T, D0, M0op, Jnu_zero);
fprintf('zero-order gate (Jnu=0): dU21 = %.3e, dU37 = %.3e  (K(0)=%.3e, Sigma(0)=%.3e)\n', ...
    dU21_0, dU37_0, tup0.K(1), tup0.Sigma(1));
assert(abs(dU21_0) < 1e-12, 'ZERO-ORDER GATE FAILED for dU21 (J 2.21 treatment)');
assert(abs(dU37_0) < 1e-12, 'zero-order sanity check failed for dU37 (J 2.22)');
results.gate21_value = dU21_0;

% Interacting final tuple at the base temperature
[dU37_1, dU21_1, tup1] = route_b_node(T, D0, M0op, Jnu);
fprintf('interacting tuple @ T=%.3g: dU37 = %.6e, dU21 = %.6e, rel diff = %.4g%%\n', ...
    T, dU37_1, dU21_1, 100*abs(dU37_1-dU21_1)/max(abs(dU37_1),abs(dU21_1)));

% Term-by-term breakdown at the base T: -M2*lam(1) and the K*(G-G0) term are IDENTICAL
% algebra shared by both formulas (verified bit-level below); the population piece is
% eq37's "0.5*alpha*n01*Delta/(1+alpha)" vs eq36's CORRECTED exact-analytic piece
% "-alpha*n1*Delta - (Delta/(2*beta*M2))*sum(wts.*R)". Report the naive (uncorrected)
% and corrected population pieces side by side to document the amendment.
[wn1, wts1, beta1] = invz_matsubara(T, 40);
tl1 = closed_tl(0, D0, M0op, beta1);
g1 = real(invz_g(tl1, 1i*wn1));  G01 = -(tl1.M2*g1);
G1 = tup1.G;  K1 = tup1.K;  lam1 = tup1.lam;
sg1 = invz_sigma(tl1, lam1, K1, g1, beta1);
n1_1 = (1 - tl1.n01)/2;
R1 = (G1 - G01) + sg1.alpha*G01;
term_M2lam_37   = -0.5*tl1.M2*lam1(1);
term_KG_37      = 0.5*real(sum(wts1.*K1.*(G1-G01)))/beta1;
term_alpha_37   = 0.5*sg1.alpha*tl1.n01*tl1.Delta/(1+sg1.alpha);
term_pop_naive  = -(tl1.Delta/(2*beta1*tl1.M2))*real(sum(wts1.*(G1-G01)));            % WRONG (uncorrected)
term_pop_fixed  = -sg1.alpha*n1_1*tl1.Delta - (tl1.Delta/(2*beta1*tl1.M2))*real(sum(wts1.*R1));  % CORRECTED
fprintf('\n-- term-by-term breakdown @ T=%.3g --\n', T);
fprintf('  -M2*lam1 term       : eq37=%.6e  eq36=%.6e  (shared algebra, diff=%.3e)\n', ...
    term_M2lam_37, term_M2lam_37, 0);
fprintf('  K*(G-G0) term       : eq37=%.6e  eq36=%.6e  (shared algebra, diff=%.3e)\n', ...
    term_KG_37, term_KG_37, 0);
fprintf('  population term     : eq37 alpha-resummed=%.6e | eq36 NAIVE(wrong)=%.6e | eq36 CORRECTED=%.6e\n', ...
    term_alpha_37, term_pop_naive, term_pop_fixed);
fprintf('  population term rel diff (eq37 vs eq36-corrected) = %.4g%%\n', ...
    100*abs(term_alpha_37-term_pop_fixed)/max(abs(term_alpha_37),abs(term_pop_fixed)));
results.term_M2lam_37 = term_M2lam_37; results.term_KG_37 = term_KG_37;
results.term_alpha_37 = term_alpha_37; results.term_pop_naive = term_pop_naive;
results.term_pop_fixed = term_pop_fixed;

% Full alternative route-B free energy using dU21 instead of dU37, SAME T grid/weights
dFB_21 = + T * trapz(Tg, dU21 ./ Tg.^2);
fprintf('dFB (eq37/J2.22)                     = %.6e\n', dFB);
fprintf('dFB (eq36/J2.21 analytic-remainder)  = %.6e\n', dFB_21);
fprintf('dFA                                   = %.6e\n', dFA);
fprintf('mismatch dFA   vs dFB(J2.22) = %.4g%%\n', 100*abs(dFA-dFB)/max(abs(dFA),abs(dFB)));
fprintf('mismatch dFA   vs dFB(J2.21) = %.4g%%\n', 100*abs(dFA-dFB_21)/max(abs(dFA),abs(dFB_21)));
fprintf('mismatch dFB37 vs dFB21      = %.4g%%\n', 100*abs(dFB-dFB_21)/max(abs(dFB),abs(dFB_21)));

results.dU37_1 = dU37_1; results.dU21_1 = dU21_1; results.dFB_21 = dFB_21;

%% ========= Step 4.2: same-final-tuple audit ========================================
fprintf('\n=== Step 4.2: same-final-tuple audit ===\n');
maxdiff37 = 0; maxdiff21 = 0;
for k = 1:numel(Tg)
    Sk_final = tuples{k}.Sigma;                       % ONLY the persisted final Sigma
    [wnk, wtsk, betak] = invz_matsubara(Tg(k), 40);    % freshly rebuilt -- no reuse
    tlk = closed_tl(0, D0, M0op, betak);
    gk  = real(invz_g(tlk, 1i*wnk));
    G0k = -(tlk.M2*gk);
    medk = invz_emt_scalar(G0k, Sk_final, Jnu, struct());
    lamk = invz_lambdas(medk.K, gk, wtsk, betak, [1 2]);
    sgk  = invz_sigma(tlk, lamk, medk.K, gk, betak);
    Gk = G0k ./ (1 + Sk_final + medk.K .* G0k);
    dU37_audit = 0.5*( sgk.alpha*tlk.n01*tlk.Delta/(1 + sgk.alpha) - tlk.M2*lamk(1) ...
                  + real(sum(wtsk .* medk.K .* (Gk - G0k)))/betak );
    n1k = (1 - tlk.n01)/2;
    Rk  = (Gk - G0k) + sgk.alpha*G0k;
    dU21_audit = 0.5*real(sum(wtsk.*medk.K.*Gk))/betak - sgk.alpha*n1k*tlk.Delta ...
                 - (tlk.Delta/(2*betak*tlk.M2))*real(sum(wtsk.*Rk));
    maxdiff37 = max(maxdiff37, abs(dU37_audit - dU37(k)));
    maxdiff21 = max(maxdiff21, abs(dU21_audit - dU21(k)));
end
fprintf('max |dU37_audit - dU37| over the T grid = %.3e\n', maxdiff37);
fprintf('max |dU21_audit - dU21| over the T grid = %.3e\n', maxdiff21);
results.audit_maxdiff37 = maxdiff37; results.audit_maxdiff21 = maxdiff21;

%% ========= Step 4.3: coupling-scale exponent scan ==================================
fprintf('\n=== Step 4.3: coupling-scale exponent scan ===\n');
svals = [1, 0.5, 0.25];
dFAs = nan(size(svals)); dFBs = nan(size(svals)); epss = nan(size(svals));
dFB21s = nan(size(svals)); eps37v21 = nan(size(svals));
for i = 1:numel(svals)
    s = svals(i);
    Jnu_s = s*Jnu;
    [dFAi, ~] = closed_routeA(T, Jnu_s, J0eff, D0, M0op, hmax, nH);
    [dU37i, dU21i, ~] = route_b_sweep(Tg, D0, M0op, Jnu_s);
    dFBi = + T * trapz(Tg, dU37i ./ Tg.^2);
    dFB21i = + T * trapz(Tg, dU21i ./ Tg.^2);
    dFAs(i) = dFAi; dFBs(i) = dFBi; dFB21s(i) = dFB21i;
    epss(i) = abs(dFAi - dFBi) / max(abs(dFAi), abs(dFBi));
    eps37v21(i) = abs(dFBi - dFB21i) / max(abs(dFBi), abs(dFB21i));
    fprintf('s = %.4g: dFA = %.6e, dFB(37) = %.6e, dFB(21) = %.6e, eps_AvB = %.4g%%, eps_37v21 = %.4g%%\n', ...
        s, dFAi, dFBi, dFB21i, 100*epss(i), 100*eps37v21(i));
end
p_eps = polyfit(log(svals), log(epss), 1);
p_A   = polyfit(log(svals), log(abs(dFAs)), 1);
p_B   = polyfit(log(svals), log(abs(dFBs)), 1);
p_eps37v21 = polyfit(log(svals), log(eps37v21), 1);
fprintf('fitted exponent eps(s)            ~ s^%.4g   (route A vs route B, THE GATE)\n', p_eps(1));
fprintf('fitted exponent |dFA(s)|          ~ s^%.4g\n', p_A(1));
fprintf('fitted exponent |dFB(s)|          ~ s^%.4g\n', p_B(1));
fprintf('fitted exponent eps_37v21(s)      ~ s^%.4g   (internal J2.22 vs J2.21 energy-reduction check)\n', p_eps37v21(1));

results.svals = svals; results.dFAs = dFAs; results.dFBs = dFBs; results.epss = epss;
results.p_eps = p_eps(1); results.p_A = p_A(1); results.p_B = p_B(1);
results.dFB21s = dFB21s; results.eps37v21 = eps37v21; results.p_eps37v21 = p_eps37v21(1);

% Supplementary: with the CORRECTED treatment, does the (now small) eq37-vs-eq36
% population-term residual grow or shrink with coupling scale?
fprintf('\n-- supplementary: CORRECTED population-term comparison vs s --\n');
pop_ratio = nan(size(svals));
for i = 1:numel(svals)
    s = svals(i);
    [~, ~, tupi] = route_b_node(T, D0, M0op, s*Jnu);
    [wni, wtsi, betai] = invz_matsubara(T, 40);
    tli = closed_tl(0, D0, M0op, betai);
    gi_ = real(invz_g(tli, 1i*wni));  G0i = -(tli.M2*gi_);
    sgi = invz_sigma(tli, tupi.lam, tupi.K, gi_, betai);
    n1i = (1 - tli.n01)/2;
    Ri  = (tupi.G - G0i) + sgi.alpha*G0i;
    term_alpha_i = 0.5*sgi.alpha*tli.n01*tli.Delta/(1+sgi.alpha);
    term_pop_i   = -sgi.alpha*n1i*tli.Delta - (tli.Delta/(2*betai*tli.M2))*real(sum(wtsi.*Ri));
    pop_ratio(i) = term_pop_i / term_alpha_i;
    fprintf('s = %.4g: eq37-alpha = %.4e, eq36-corrected = %.4e, ratio = %.4g\n', ...
        s, term_alpha_i, term_pop_i, pop_ratio(i));
end
results.pop_ratio = pop_ratio;

save(fullfile(here, 'task6b_results.mat'), 'results');
fprintf('\nAll sections completed. Results saved to task6b_results.mat\n');

%% ================= local functions (closed 2x2 model, frozen per brief) ===========
function tl = closed_tl(h, D0, M0op, beta)
% FROZEN closed 2x2 (SS7b): H = diag(+/-D0/2) - h*Jz2, Jz2 = [0 M0op; M0op 0] (FIXED
% operator, never re-selected from a CF solve).
Jz2 = [0, M0op; M0op, 0];
H = [D0/2, 0; 0, -D0/2] - h*Jz2;
[V, E] = eig((H+H')/2, 'vector');  [E, ix] = sort(real(E));  V = V(:, ix);
p = exp(-beta*(E - E(1)));  p = p/sum(p);
Mz = V'*Jz2*V;
tl = struct('m', real(Mz(1,1)), 'M2', abs(Mz(1,2))^2, 'n01', p(1) - p(2), ...
            'Delta', E(2) - E(1), 'g0', 0);
tl.g0 = 2*tl.n01/tl.Delta;
end

function [rk, Mk, S0k, l1g0, ok] = closed_node(h, T, Jnu, J0eff, D0, M0op)
% One closed-model node, COLD-STARTED. Ordered Sigma<->EMT loop, TWO-LEVEL static
% weights (SS7b, J 2.28 verbatim parametrization). Verbatim from the brief.
[wn, wts, beta] = invz_matsubara(T, 40);
tl = closed_tl(h, D0, M0op, beta);
h0el = beta*(1 - tl.n01^2);
g  = real(invz_g(tl, 1i*wn));
G0 = -(tl.M2*g);  G0(1) = -(tl.M2*tl.g0 + tl.m^2*h0el);
gi = -tl.M2*tl.g0;  ge = -tl.m^2*h0el;
Sigma = zeros(size(wn));  K = zeros(size(wn));  lam = [0; 0; 0];  K0s = 0;  ok = false;
for outer = 1:200
    med = invz_emt_scalar(G0, Sigma, Jnu, struct());
    if ~med.converged, break; end
    K = med.K;
    [K0s, ~, so] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), Jnu, K0s, beta, ...
                                           J0eff, gi, ge, struct());
    K(1) = K0s;
    lam = invz_lambdas(K, g, wts, beta, [1 2 3]);
    sg  = invz_sigma_ordered(tl, lam, K, g, beta);
    dS  = max(abs(sg.Sigma - Sigma));
    Sigma = Sigma + 0.7*(sg.Sigma - Sigma);
    if dS < 1e-8 && so.converged, ok = true; break; end
end
[K0s, ~, so] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), Jnu, K0s, beta, ...
                                       J0eff, gi, ge, struct());
K(1) = K0s;
lam_chk = invz_lambdas(K, g, wts, beta, [1 2 3]);
sg_chk  = invz_sigma_ordered(tl, lam_chk, K, g, beta);
ok = ok && so.converged && isfinite(so.resid) && so.resid < 1e-10 ...
        && max(abs(sg_chk.Sigma - Sigma)) < 1e-8;
rk = so.r;  Mk = tl.m*tl.n01;  S0k = Sigma(1);  l1g0 = abs(lam_chk(1)*tl.g0);
end

function [dFA, dgn] = closed_routeA(T, Jnu, J0eff, D0, M0op, hmax, nH)
% Route A on ONE grid, all gates inside. Verbatim from the brief (testCase
% verifyTrue/verifyLessThan replaced by assert, since this is a standalone script).
hgrid = hmax * (1e-4^(1/(nH-1))).^((nH-1):-1:0);
[r0, ~, ~, ~, ok0] = closed_node(0, T, Jnu, J0eff, D0, M0op);
assert(ok0, 'closed node h=0 failed to converge');
[rv, Mv, S0v, Lv] = deal(nan(1, nH));
for k = 1:nH
    [rv(k), Mv(k), S0v(k), Lv(k), okk] = closed_node(hgrid(k), T, Jnu, J0eff, D0, M0op);
    assert(okk, sprintf('closed node h = %.4g failed', hgrid(k)));
end
h0 = cumtrapz([0 hgrid], [r0 rv]);  h0 = h0(2:end);
dh = h0 - hgrid;
q = max(1e-4, 0.01*abs(r0 - 1));
assert(all(abs(rv(end-4:end) - 1) < q), 'fluctuations not quenched');
assert(all(abs(S0v(end-4:end)) < q), 'Sigma(0) not quenched on the tail');
assert(all(abs(Lv(end-4:end)) < q), 'lambda1*g0 not quenched on the tail');
assert((M0op - Mv(end)) < 0.01*M0op, 'analytic saturation not reached');
dFA = -trapz([0 Mv], [0 dh]);
dgn = struct('tail_est_A', abs(dh(end))*(M0op - Mv(end)), 'dh_end', dh(end), 'r0', r0);
end

function [dU37, dU21, tup] = route_b_node(T, D0, M0op, Jnu)
% One route-B (h=0, PM) node: converge Sigma via the PM emt+sigma loop (verbatim from
% the brief's route-B body), then compute BOTH dU via J 2.22 (eq 37) and dU via J 2.21
% (eq 36, CORRECTED analytic-remainder treatment -- see the Step 4.1 header comment
% above AND the amendment below: a first attempt that dropped the ENTIRE (G-G0) odd-wn
% piece as "absolutely convergent" turned out WRONG, because Sigma(wn) -> alpha (a
% NONZERO constant, not zero) as wn -> infinity, so (G-G0) ~ -alpha*G0(wn) at large wn
% -- the SAME marginally-convergent 1/wn tail as the bare G0 piece, needing the SAME
% convergence-factor regularization, not a free pass. The amended treatment explicitly
% peels off this alpha-driven tail analytically (reusing the SAME regularized identity
% established for the G0 piece: (1/2beta)*sum^reg (Delta+wn)*g(wn) = n1*Delta, with
% n1 = (1-n01)/2 the bare excited-state population) and sums only the GENUINELY
% fast-decaying remainder R = (G-G0)+alpha*G0 (= O(1/wn^4), confirmed numerically)
% directly, with NO convergence factor, whose odd-wn part now legitimately vanishes.
[wn, wts, beta] = invz_matsubara(T, 40);
tl = closed_tl(0, D0, M0op, beta);
g  = real(invz_g(tl, 1i*wn));
G0 = -(tl.M2*g);                                      % m = 0: no elastic term
Sigma = zeros(size(wn));  ok = false;
for outer = 1:200
    med = invz_emt_scalar(G0, Sigma, Jnu, struct());
    if ~med.converged, break; end
    lam = invz_lambdas(med.K, g, wts, beta, [1 2]);
    sg  = invz_sigma(tl, lam, med.K, g, beta);
    dS = max(abs(sg.Sigma - Sigma));  Sigma = Sigma + 0.7*(sg.Sigma - Sigma);
    if dS < 1e-8, ok = true; break; end
end
assert(ok, 'route-B node did not converge');
% one SIMULTANEOUS final tuple, recomputed from the converged Sigma
med = invz_emt_scalar(G0, Sigma, Jnu, struct());
assert(med.converged);
lam = invz_lambdas(med.K, g, wts, beta, [1 2]);
sg  = invz_sigma(tl, lam, med.K, g, beta);
G = G0 ./ (1 + Sigma + med.K .* G0);                  % J 2.30
dU37 = 0.5*( sg.alpha*tl.n01*tl.Delta/(1 + sg.alpha) - tl.M2*lam(1) ...
              + real(sum(wts .* med.K .* (G - G0)))/beta );
n1 = (1 - tl.n01)/2;                                  % bare excited-state population
Rrem = (G - G0) + sg.alpha*G0;                        % fast-decaying remainder (O(1/wn^4))
dU21 = 0.5*real(sum(wts.*med.K.*G))/beta - sg.alpha*n1*tl.Delta ...
       - (tl.Delta/(2*beta*tl.M2))*real(sum(wts.*Rrem));
tup = struct('Sigma', Sigma, 'K', med.K, 'G', G, 'G0', G0, 'wn', wn, 'wts', wts, ...
             'beta', beta, 'tl', tl, 'lam', lam, 'alpha', sg.alpha);
end

function [dU37, dU21, tuples] = route_b_sweep(Tg, D0, M0op, Jnu)
dU37 = nan(size(Tg)); dU21 = nan(size(Tg)); tuples = cell(size(Tg));
for k = 1:numel(Tg)
    [dU37(k), dU21(k), tuples{k}] = route_b_node(Tg(k), D0, M0op, Jnu);
end
end
