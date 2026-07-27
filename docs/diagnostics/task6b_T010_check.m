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
T = 0.10;  C = invz_const();  beta = 1/(C.kB*T);
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

fprintf('T010 CHECK: dFA = %.6e, dFB = %.6e, mismatch = %.4f%%\n', dFA, dFB, 100*abs(dFA-dFB)/max(abs(dFA),abs(dFB)));
return
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
