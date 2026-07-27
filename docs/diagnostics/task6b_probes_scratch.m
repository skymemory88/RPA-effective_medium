%% Task 6b Probes (P1-P4): targeted localization of the closed-model 13.65% mismatch
% Standalone diagnostic script (NOT a test, no commits, no tolerance decisions). Adapts
% the Step-2/4 harness (task6b_closed_model_scratch.m) with route-A ablations that
% surgically modify the ordered elastic static sector (Gstat/Gtil0/r, invz_gstat_ordered
% + invz_emt_static_ordered, SS3/SS7b) while leaving route B (plain-PM dU, J 2.22, h=0)
% completely untouched. All variants share the closed 2x2 fixture, tl(h), and hgrid with
% the frozen baseline harness so numbers are directly comparable to task-6b-report.md.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'invz_common'));
addpath(fullfile(here, '..', '..', 'invz_projected'));

R = struct();

%% ================= Fixture (frozen, identical to task-6b) ==========================
T = 0.31;  C = invz_const();  beta = 1/(C.kB*T); %#ok<NASGU>
D0 = 0.2;  M0op = 3.0;  J0eff = 6.4e-3;
Jnu = linspace(-4e-3, 4.0e-3, 24).';
assert(abs(mean(Jnu)) < 1e-15, 'zero-mean Jensen lattice violated');
hmax = 40*D0;  nH = 81;

fprintf('=== Baseline reproduction (sanity gate before trusting any variant) ===\n');
[dFA_base, dgn_base] = closed_routeA_generic(T, Jnu, J0eff, D0, M0op, hmax, nH, 'baseline');
fprintf('dFA_base (this harness)   = %.6e\n', dFA_base);
fprintf('dFA      (task-6b-report) = -5.231699e-04  (expected match)\n');
relmis = abs(dFA_base - (-5.231699e-04))/abs(-5.231699e-04);
fprintf('rel diff vs task-6b-report baseline = %.4g%%\n', 100*relmis);
assert(relmis < 1e-4, 'BASELINE REIMPLEMENTATION DOES NOT MATCH task-6b-report -- STOP, fix the harness');

% Route B is untouched by every route-A ablation below -- compute ONCE and reuse.
Tg = T * (1.35.^(0:0.5:17));
dU37 = nan(size(Tg));
for k = 1:numel(Tg)
    [dU37(k), ~, ~] = route_b_node(Tg(k), D0, M0op, Jnu);
end
dFB = + T * trapz(Tg, dU37 ./ Tg.^2);
fprintf('dFB (route B, unaffected by route-A ablations) = %.6e\n', dFB);
mismatch_base = 100*abs(dFA_base - dFB)/max(abs(dFA_base), abs(dFB));
fprintf('BASELINE gate mismatch = %.4g%%  (report reference: 13.65%%)\n\n', mismatch_base);
R.dFA_base = dFA_base; R.dFB = dFB; R.mismatch_base = mismatch_base;

fmt = @(name, dFAv) fprintf('%-10s dFA=%.6e  dFB=%.6e  mismatch=%.4g%%  (Delta from baseline %.4g%% : %+.4g pp)\n', ...
    name, dFAv, dFB, 100*abs(dFAv-dFB)/max(abs(dFAv),abs(dFB)), mismatch_base, ...
    100*abs(dFAv-dFB)/max(abs(dFAv),abs(dFB)) - mismatch_base);

%% ================= P1: elastic-sector ablations of route A =========================
fprintf('=== P1: elastic-sector ablations ===\n');
[dFA_xi1,  dgn_xi1]  = closed_routeA_generic(T, Jnu, J0eff, D0, M0op, hmax, nH, 'xi1');
[dFA_noel, dgn_noel] = closed_routeA_generic(T, Jnu, J0eff, D0, M0op, hmax, nH, 'noel');
[dFA_pmc,  dgn_pmc]  = closed_routeA_generic(T, Jnu, J0eff, D0, M0op, hmax, nH, 'pmcancel');
fmt('P1a(xi=1)',  dFA_xi1);
fmt('P1b(noel)',  dFA_noel);
fmt('P1c(pmcxl)', dFA_pmc);
R.dFA_xi1 = dFA_xi1; R.dFA_noel = dFA_noel; R.dFA_pmc = dFA_pmc;
R.dgn_xi1 = dgn_xi1; R.dgn_noel = dgn_noel; R.dgn_pmc = dgn_pmc;

%% ================= P2: static closure vs direct-solve K0 ============================
fprintf('\n=== P2: K(1) from ordinary direct solve instead of the elastic closure ===\n');
[dFA_p2, dgn_p2] = closed_routeA_generic(T, Jnu, J0eff, D0, M0op, hmax, nH, 'p2');
fmt('P2', dFA_p2);
R.dFA_p2 = dFA_p2; R.dgn_p2 = dgn_p2;

%% ================= P3: ordered-Sigma static-approximation share =====================
fprintf('\n=== P3: PM Sigma(0) substituted in the static slot only ===\n');
[dFA_p3, dgn_p3] = closed_routeA_generic(T, Jnu, J0eff, D0, M0op, hmax, nH, 'p3');
fmt('P3', dFA_p3);
R.dFA_p3 = dFA_p3; R.dgn_p3 = dgn_p3;

fprintf('\n-- P3a: Sigma(0) difference profile vs h (ordered Sigma_ordered(0) vs PM Sigma(0), same K/lam/tl) --\n');
fprintf('%12s %16s %16s %16s\n', 'h', 'Sigma0_ordered', 'Sigma0_PM', 'diff');
idxshow = round(linspace(1, nH, 12));
for kk = idxshow
    fprintf('%12.4g %16.6e %16.6e %16.6e\n', dgn_base.hgrid(kk), dgn_base.S0v(kk), dgn_base.S0v_PM(kk), ...
        dgn_base.S0v(kk) - dgn_base.S0v_PM(kk));
end

%% ================= P4: m -> 0 consistency spot-check ================================
fprintf('\n=== P4: m->0 consistency spot-check (three smallest grid nodes vs PM anchor) ===\n');
[r0b, ~, S0_0, ~, K0_0, ok0] = closed_node_v2(0, T, Jnu, J0eff, D0, M0op, 'baseline');
[~, ~, tup_pm] = route_b_node(T, D0, M0op, Jnu);
r_pm = 1 + tup_pm.Sigma(1);
fprintf('h=0        (ordered) : r=%.8f  Sigma0=%.6e  K0=%.6e  (ok=%d)\n', r0b, S0_0, K0_0, ok0);
fprintf('h=0        (PM route B) : r=%.8f  Sigma0=%.6e  K0=%.6e\n', r_pm, tup_pm.Sigma(1), tup_pm.K(1));
fprintf('  |Delta r|=%.3e  |Delta Sigma0|=%.3e  |Delta K0|=%.3e\n', ...
    abs(r0b-r_pm), abs(S0_0-tup_pm.Sigma(1)), abs(K0_0-tup_pm.K(1)));
for kk = 1:3
    h = dgn_base.hgrid(kk);
    fprintf('h=%.4e (ordered, node %d) : r=%.8f  Sigma0=%.6e  K0=%.6e\n', ...
        h, kk, dgn_base.rv(kk), dgn_base.S0v(kk), dgn_base.K0v(kk));
end
R.P4 = struct('r0_ordered', r0b, 'r_pm', r_pm, 'S0_0_ordered', S0_0, 'S0_pm', tup_pm.Sigma(1), ...
              'K0_0_ordered', K0_0, 'K0_pm', tup_pm.K(1), ...
              'h3', dgn_base.hgrid(1:3), 'r3', dgn_base.rv(1:3), 'S03', dgn_base.S0v(1:3), 'K03', dgn_base.K0v(1:3));

%% ================= Summary table =====================================================
fprintf('\n=== SUMMARY ===\n');
fprintf('%-12s %14s %14s %10s %10s\n', 'variant', 'dFA', 'dFB', 'mismatch%', 'Delta_pp');
rows = {'baseline', dFA_base; 'P1a_xi1', dFA_xi1; 'P1b_noel', dFA_noel; 'P1c_pmcxl', dFA_pmc; ...
        'P2', dFA_p2; 'P3', dFA_p3};
for i = 1:size(rows,1)
    name = rows{i,1}; dFAv = rows{i,2};
    mism = 100*abs(dFAv-dFB)/max(abs(dFAv),abs(dFB));
    fprintf('%-12s %14.6e %14.6e %10.4g %+10.4g\n', name, dFAv, dFB, mism, mism-mismatch_base);
end

save(fullfile(here, 'task6b_probes_results.mat'), 'R');
fprintf('\nAll probes completed. Results saved to task6b_probes_results.mat\n');

%% ================= local functions ===================================================
function tl = closed_tl(h, D0, M0op, beta)
Jz2 = [0, M0op; M0op, 0];
H = [D0/2, 0; 0, -D0/2] - h*Jz2;
[V, E] = eig((H+H')/2, 'vector');  [E, ix] = sort(real(E));  V = V(:, ix);
p = exp(-beta*(E - E(1)));  p = p/sum(p);
Mz = V'*Jz2*V;
tl = struct('m', real(Mz(1,1)), 'M2', abs(Mz(1,2))^2, 'n01', p(1) - p(2), ...
            'Delta', E(2) - E(1), 'g0', 0);
tl.g0 = 2*tl.n01/tl.Delta;
end

function [Gstat, out] = gstat_core(tl, lam, K0, Sigma0, beta, G0inel0, G0el0, xiMode)
% Mirrors invz_gstat_ordered.m exactly when xiMode = 'auto'; xiMode = scalar overrides xi.
m = tl.m; M2 = tl.M2; n01 = tl.n01; g0 = tl.g0;
if ischar(xiMode) && strcmp(xiMode, 'auto')
    xi = (1 + tanh(m^2*n01^2*beta*K0 - M2*beta*lam(1))) / ...
         (1 + (4*n01^2*K0*g0 + 2*lam(2) + g0*lam(1))*M2/n01^2);
else
    xi = xiMode;
end
Gstat  = G0inel0/(1 + Sigma0 + K0*G0inel0) + xi*G0el0;
G0bare = G0inel0 + G0el0;
Gtil0  = Gstat/(1 - K0*Gstat);
out = struct('xi', xi, 'G0bare', G0bare, 'Gtil0', Gtil0, 'r', G0bare/Gtil0);
end

function [K0, Gstat, out] = emt_static_core(tl, lam, Sigma0, Jnu_flat, K0_seed, beta, J0eff, G0inel0, G0el0, xiMode)
% Mirrors invz_emt_static_ordered.m exactly when xiMode = 'auto' and G0el0 is the true
% elastic weight; xiMode / G0el0 may be overridden by the caller to implement P1a/P1b.
rtol = 1e-10; maxit = 200; mix = 0.5;
Jf = Jnu_flat(:);
K0 = K0_seed;
for it = 1:maxit
    Gs = gstat_core(tl, lam, K0, Sigma0, beta, G0inel0, G0el0, xiMode);
    Gq = Gs ./ (1 + (Jf - K0).*Gs);
    Gbar = mean(Gq);
    if abs(Gbar - Gs) < rtol, break; end
    K0_new = mean(Jf .* Gq) / Gbar;
    dK = abs(K0_new - K0);
    if dK < 4*eps(abs(K0)), break; end
    K0 = K0 + mix*(K0_new - K0);
end
[Gstat, go] = gstat_core(tl, lam, K0, Sigma0, beta, G0inel0, G0el0, xiMode);
Gq = Gstat ./ (1 + (Jf - K0).*Gstat);
out = go;
out.D_uni = 1 + (J0eff - K0)*Gstat;
out.resid = abs(mean(Gq) - Gstat);
out.iters = it;
out.converged = out.resid < rtol;
end

function [rk, Mk, S0k, l1g0, K0k, ok, extra] = closed_node_v2(h, T, Jnu, J0eff, D0, M0op, kind)
%CLOSED_NODE_V2 extra returns the final converged tuple (K, lam, tl, g, beta) so callers
% (e.g. the P3a Sigma(0) profile) can recompute alternative Sigma formulas post hoc
% without perturbing the loop that produced them.
% One closed-model node under a named ablation ('baseline'|'xi1'|'noel'|'pmcancel'|'p2'|'p3').
[wn, wts, beta] = invz_matsubara(T, 40);
tl = closed_tl(h, D0, M0op, beta);
h0el = beta*(1 - tl.n01^2);
g  = real(invz_g(tl, 1i*wn));
G0 = -(tl.M2*g);  G0(1) = -(tl.M2*tl.g0 + tl.m^2*h0el);
gi = -tl.M2*tl.g0;  ge = -tl.m^2*h0el;
G0bare_true = gi + ge;

ge_arg = ge;   if strcmp(kind, 'noel'), ge_arg = 0; end
xiMode = 'auto'; if strcmp(kind, 'xi1'), xiMode = 1; end

Sigma = zeros(size(wn));  K = zeros(size(wn));  lam = [0; 0; 0];  K0s = 0;  ok = false;  so = struct('converged', false);
for outer = 1:200
    med = invz_emt_scalar(G0, Sigma, Jnu, struct());
    if ~med.converged, break; end
    K = med.K;
    if strcmp(kind, 'p2')
        % P2: keep K(1) = med.K(1), the ORDINARY direct solve on G0(1)=G0bare_true --
        % skip the elastic static closure entirely.
        Gt = med.G(1)/(1 - K(1)*med.G(1));
        so = struct('converged', true, 'resid', 0, 'Gtil0', Gt);
    else
        if strcmp(kind, 'p3')
            sig_pm = invz_sigma(tl, lam(1:2), K, g, beta);
            Sigma0_arg = sig_pm.Sigma(1);
        else
            Sigma0_arg = Sigma(1);
        end
        [K0s, ~, so] = emt_static_core(tl, lam(1:2), Sigma0_arg, Jnu, K0s, beta, J0eff, gi, ge_arg, xiMode);
        K(1) = K0s;
    end
    lam = invz_lambdas(K, g, wts, beta, [1 2 3]);
    sg  = invz_sigma_ordered(tl, lam, K, g, beta);
    dS  = max(abs(sg.Sigma - Sigma));
    Sigma = Sigma + 0.7*(sg.Sigma - Sigma);
    if dS < 1e-8 && so.converged, ok = true; break; end
end
% final refresh (same simultaneous-final-tuple contract as the baseline harness)
if strcmp(kind, 'p2')
    med = invz_emt_scalar(G0, Sigma, Jnu, struct());
    K = med.K;
    Gt = med.G(1)/(1 - K(1)*med.G(1));
    so = struct('converged', med.converged, 'resid', 0, 'Gtil0', Gt);
else
    if strcmp(kind, 'p3')
        sig_pm = invz_sigma(tl, lam(1:2), K, g, beta);
        Sigma0_arg = sig_pm.Sigma(1);
    else
        Sigma0_arg = Sigma(1);
    end
    [K0s, ~, so] = emt_static_core(tl, lam(1:2), Sigma0_arg, Jnu, K0s, beta, J0eff, gi, ge_arg, xiMode);
    K(1) = K0s;
end
lam_chk = invz_lambdas(K, g, wts, beta, [1 2 3]);
sg_chk  = invz_sigma_ordered(tl, lam_chk, K, g, beta);
ok = ok && so.converged && max(abs(sg_chk.Sigma - Sigma)) < 1e-8;

if strcmp(kind, 'pmcancel')
    rk = 1 + Sigma(1);                       % J 2.32's K-cancelled PM form, r = 1+Sigma(0)
else
    rk = G0bare_true / so.Gtil0;
end
Mk = tl.m*tl.n01; S0k = Sigma(1); K0k = K(1); l1g0 = abs(lam_chk(1)*tl.g0);
extra = struct('K', K, 'lam', lam_chk, 'tl', tl, 'g', g, 'beta', beta);
end

function [dFA, dgn] = closed_routeA_generic(T, Jnu, J0eff, D0, M0op, hmax, nH, kind)
hgrid = hmax * (1e-4^(1/(nH-1))).^((nH-1):-1:0);
[r0, ~, ~, ~, ~, ok0] = closed_node_v2(0, T, Jnu, J0eff, D0, M0op, kind);
if ~ok0, warning('task6b_probes:node0', '[%s] closed node h=0 failed to converge', kind); end
[rv, Mv, S0v, Lv, K0v] = deal(nan(1, nH));
S0v_PM = nan(1, nH);
nfail = 0;
for k = 1:nH
    [rv(k), Mv(k), S0v(k), Lv(k), K0v(k), okk, extrak] = closed_node_v2(hgrid(k), T, Jnu, J0eff, D0, M0op, kind);
    if ~okk, nfail = nfail + 1; end
    % P3a diagnostic: PM Sigma(0) at the SAME final (K, lam, tl) as this node's converged
    % tuple -- computed for every kind (cheap), read off the baseline call in the main script.
    sig_pm_k = invz_sigma(extrak.tl, extrak.lam(1:2), extrak.K, extrak.g, extrak.beta);
    S0v_PM(k) = sig_pm_k.Sigma(1);
end
if nfail > 0
    warning('task6b_probes:nodeFail', '[%s] %d/%d nodes failed to converge', kind, nfail, nH);
end
h0 = cumtrapz([0 hgrid], [r0 rv]);  h0 = h0(2:end);
dh = h0 - hgrid;
q = max(1e-4, 0.01*abs(r0 - 1));
quench_r   = all(abs(rv(end-4:end) - 1) < q);
quench_S0  = all(abs(S0v(end-4:end)) < q);
quench_L   = all(abs(Lv(end-4:end)) < q);
sat_ok     = (M0op - Mv(end)) < 0.01*M0op;
if ~(quench_r && quench_S0 && quench_L && sat_ok)
    warning('task6b_probes:quench', ...
        '[%s] quench/saturation criteria not all satisfied (r:%d S0:%d L:%d sat:%d) -- reporting anyway (investigation, not a gate)', ...
        kind, quench_r, quench_S0, quench_L, sat_ok);
end
dFA = -trapz([0 Mv], [0 dh]);
dgn = struct('tail_est_A', abs(dh(end))*(M0op - Mv(end)), 'dh_end', dh(end), 'r0', r0, ...
             'hgrid', hgrid, 'rv', rv, 'Mv', Mv, 'S0v', S0v, 'S0v_PM', S0v_PM, 'Lv', Lv, 'K0v', K0v, ...
             'quench_r', quench_r, 'quench_S0', quench_S0, 'quench_L', quench_L, 'sat_ok', sat_ok, 'nfail', nfail);
end

function [dU37, dU21, tup] = route_b_node(T, D0, M0op, Jnu)
[wn, wts, beta] = invz_matsubara(T, 40);
tl = closed_tl(0, D0, M0op, beta);
g  = real(invz_g(tl, 1i*wn));
G0 = -(tl.M2*g);
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
med = invz_emt_scalar(G0, Sigma, Jnu, struct());
lam = invz_lambdas(med.K, g, wts, beta, [1 2]);
sg  = invz_sigma(tl, lam, med.K, g, beta);
G = G0 ./ (1 + Sigma + med.K .* G0);
dU37 = 0.5*( sg.alpha*tl.n01*tl.Delta/(1 + sg.alpha) - tl.M2*lam(1) ...
              + real(sum(wts .* med.K .* (G - G0)))/beta );
dU21 = NaN;   % not needed for the probes (P1-P4 do not touch dU21/eq36)
tup = struct('Sigma', Sigma, 'K', med.K, 'G', G, 'G0', G0, 'lam', lam, 'alpha', sg.alpha, 'tl', tl);
end
