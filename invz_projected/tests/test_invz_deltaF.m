function tests = test_invz_deltaF
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_bare_limit_deltaF_zero(testCase)
% Route A with force_bare: dh = h0 - hmf = 0 identically, so BOTH dF and the
% non-claiming endpoint diagnostic are exactly zero (P0-3: the round-1 draft asserted
% tail > 0, which contradicts dh == 0). Field names per the dF_partial contract
% (stage-2 task 6b, step 1): out.dF_partial, out.endpoint_dh.
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
[dF, out] = invz_deltaF_ordered(ion, T, [2.85 0 0], Jnu, ...
    struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true, 'force_bare', true));
verifyEqual(testCase, dF, 0, 'AbsTol', 1e-12);
verifyEqual(testCase, out.dF_partial, 0, 'AbsTol', 1e-12);
verifyEqual(testCase, out.endpoint_dh, 0, 'AbsTol', 1e-15);
end

function test_two_routes_closed_model(testCase)
% J 2.34 two-route identity in the CLOSED 2x2 model (plan SS7b, stage-2 task 6b step 6 --
% amends the step-2 harness test of the same name). One moment M = m*n01 serves J 2.31
% (Gate 1), J 2.33, the work term, and J 2.34's measure -- the conjugacy chain the
% rejected SS7a draft broke. Saturation is VERIFIED against the ANALYTIC M_sat = M0op.
% FAST: seconds of 2x2 work.
%
% Blocker-review resolution (invzp_task6b_blocker_review_Codex.md, adjudicated against
% the published Jensen paper): the closed-model two-route mismatch (~13.7%, essentially
% scale-free in the coupling) is a same-retained-order static-elastic approximation
% residual -- NOT a demonstrated implementation defect (Jensen's own full-machinery
% published check achieves only 2-3% on physical HoF3 at low T, PRB 49, 11833). Per that
% resolution, NO percentage-equality gate (10% or 15%) may be committed here as a
% physics gate -- the old verifyEqual(dFA, dFB, 'RelTol', 0.10) gate is REMOVED. What IS
% derivable and IS enforced below: every exact-limit, convergence, and quench gate from
% the original harness (unchanged), plus (1) the ORDER-CONSISTENCY gate
% (leading-order-in-coupling scaling) and (2) the approximation-fingerprint regression,
% both in place of the retired equality gate. See invzp_QCP_diagnosis.md and
% invz_deltaF_ordered.m for the full scope-limitation record.
T = 0.31;
D0 = 0.2;  M0op = 3.0;  J0eff = 6.4e-3;
Jnu = linspace(-4e-3, 4.0e-3, 24).';                  % ZERO-MEAN: enforces Jensen's
verifyLessThan(testCase, abs(mean(Jnu)), 1e-15);      % no-self-site identity J(ii) = 0
% (load-bearing, second SS7 review P0-1 -- a nonzero mean is NOT a valid Jensen lattice:
%  measured consequence was opposite-sign routes)

% ---- Route A with EXECUTED convergence sweeps (third-review P1) --------------------
hmax = 40*D0;  nH = 81;
[dFA,    dA] = closed_routeA(T, Jnu, J0eff, D0, M0op, hmax,   nH,     testCase);
[dFA_h2, ~ ] = closed_routeA(T, Jnu, J0eff, D0, M0op, 2*hmax, nH,     testCase);
[dFA_n2, ~ ] = closed_routeA(T, Jnu, J0eff, D0, M0op, hmax,   2*nH-1, testCase);
verifyEqual(testCase, dFA, dFA_h2, 'RelTol', 2e-2);   % hmax-doubling convergence
verifyEqual(testCase, dFA, dFA_n2, 'RelTol', 2e-2);   % grid-doubling convergence

% ---- Route B: temperature integral at h = 0 (PM, m = 0) ----------------------------
% DENSITY-CONVERGED grid (ratio sqrt(1.35)); EVERY node must converge -- silently
% deleting failed temperatures corrupts both convergence checks and can false-pass
% (third-review P1). The DENSE value enters the gate.
Tg = T * (1.35.^(0:0.5:17));
dU = route_b_sweep(testCase, Tg, D0, M0op, Jnu);
verifyTrue(testCase, all(isfinite(dU)));              % NO skipped nodes in a blocking test
dFB   = + T * trapz(Tg, dU ./ Tg.^2);                 % the DENSE value: THE gate input
dFB_c = + T * trapz(Tg(1:2:end), dU(1:2:end) ./ Tg(1:2:end).^2);   % true 1.35-spaced subset
verifyEqual(testCase, dFB_c, dFB, 'RelTol', 0.05);    % temperature-DENSITY convergence
nTm = nnz(Tg <= T*1.35^14);                           % shorter-ceiling prefix (intact grid)
dFBt = + T * trapz(Tg(1:nTm), dU(1:nTm) ./ Tg(1:nTm).^2);
verifyEqual(testCase, dFBt, dFB, 'RelTol', 0.05);     % Tmax-extension convergence
tail_est_B = + T * abs(dU(end)) / Tg(end);
verifyLessThan(testCase, tail_est_B, 0.05*max(abs(dFA), abs(dFB)));
verifyLessThan(testCase, dA.tail_est_A, 0.05*max(abs(dFA), abs(dFB)));

% ---- ORDER-CONSISTENCY GATE (blocker-review resolution -- replaces the invalid
% 10%/15% equality gate; plan SS7b "the engineering gate that IS derivable now") -------
% Run at s = 1, 1/2, 1/4 under Jnu -> s*Jnu; fit the |dFA|, |dFB|, |dFA-dFB| exponents
% and the relative-mismatch exponent p_eps. The bands below are ENGINEERING TEST
% TOLERANCES encoding the OBSERVED leading-order structure (measured: p_A ~ 2.00,
% p_B ~ 2.01, p_eps ~ 0.06 -- both routes are clean leading-order-s^2 free-energy
% corrections and their mismatch is scale-free) -- they are explicitly NOT analytic
% error bounds.
svals = [1, 0.5, 0.25];
dFAs = nan(size(svals));  dFBs = nan(size(svals));
dFAs(1) = dFA;  dFBs(1) = dFB;                        % reuse the base (s=1) computation
for i = 2:numel(svals)
    s = svals(i);
    [dFAs(i), ~] = closed_routeA(T, s*Jnu, J0eff, D0, M0op, hmax, nH, testCase);
    dUi = route_b_sweep(testCase, Tg, D0, M0op, s*Jnu);
    verifyTrue(testCase, all(isfinite(dUi)));
    dFBs(i) = + T * trapz(Tg, dUi ./ Tg.^2);
end
epss = abs(dFAs - dFBs) ./ max(abs(dFAs), abs(dFBs));
p_A   = polyfit(log(svals), log(abs(dFAs)), 1);
p_B   = polyfit(log(svals), log(abs(dFBs)), 1);
p_D   = polyfit(log(svals), log(abs(dFAs - dFBs)), 1);
p_eps = polyfit(log(svals), log(epss), 1);
verifyGreaterThanOrEqual(testCase, p_A(1), 1.9);
verifyLessThanOrEqual(testCase, p_A(1), 2.1);
verifyGreaterThanOrEqual(testCase, p_B(1), 1.9);
verifyLessThanOrEqual(testCase, p_B(1), 2.1);
verifyGreaterThanOrEqual(testCase, p_D(1), 1.9);
verifyLessThanOrEqual(testCase, p_D(1), 2.1);
verifyLessThanOrEqual(testCase, abs(p_eps(1)), 0.15);

% ---- APPROXIMATION-FINGERPRINT REGRESSION (non-gating diagnostic; plan SS7b stage-2
% task 6b step 6 item 2) -------------------------------------------------------------
% Pins dFA/dFB at T = 0.31 K (the s=1 base-fixture values already computed above -- NO
% recomputation) and the T = 0.10 K cross-check pair to their recorded values
% (task-6b-report.md; task6b_T010_check.m -- reproduced twice, controller + reviewer,
% digit-identical). This detects CODE DRIFT in the closed-model kernels -- it is NOT a
% two-route thermodynamic validation. The ~13.7% route residual is the documented
% SAME-RETAINED-ORDER static-elastic approximation fingerprint: Jensen's own
% full-machinery published check (PRB 49, 11833) achieves only 2-3% on the physical HoF3
% lattice at low T (with his stated high-T elastic caveat); see plan SS7b and
% invzp_QCP_diagnosis.md.
verifyEqual(testCase, dFA, -5.2317e-4, 'RelTol', 0.02);
verifyEqual(testCase, dFB, -6.0585e-4, 'RelTol', 0.02);

[dFA10, ~] = closed_routeA(0.10, Jnu, J0eff, D0, M0op, hmax, nH, testCase);
Tg10 = 0.10 * (1.35.^(0:0.5:17));
dU10 = route_b_sweep(testCase, Tg10, D0, M0op, Jnu);
verifyTrue(testCase, all(isfinite(dU10)));
dFB10 = + 0.10 * trapz(Tg10, dU10 ./ Tg10.^2);
verifyEqual(testCase, dFA10, -5.136449e-4, 'RelTol', 0.02);
verifyEqual(testCase, dFB10, -6.006005e-4, 'RelTol', 0.02);
end

%% ================= local functions (closed 2x2 model, frozen per plan SS7b) =========
function tl = closed_tl(h, D0, M0op, beta)
% FROZEN closed 2x2 (SS7b): H = diag(+/-D0/2) - h*Jz2, Jz2 = [0 M0op; M0op 0] (FIXED
% operator, never re-selected from a CF solve). Purely off-diagonal Jz2: m(0) = 0 (true
% PM anchor) and M_sat = M0op ANALYTICALLY (the Jz2 eigenvalue).
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
% One closed-model node, COLD-STARTED (third-review P2: node independence explicit --
% no warm-start seed). The ordered Sigma<->EMT loop with TWO-LEVEL static weights: the
% exact SS3-check-(a) parametrization in which J 2.28 is verbatim. Mirrors the
% production round-5 final-tuple contract (third-review P1): every dynamic solve
% requires med.converged; the final static refresh keeps its K0 in K(1); lambdas and
% the ordered Sigma target are RECOMPUTED from that exported K and must close to the
% outer tolerance; the returned diagnostics come from the revalidated final tuple.
[wn, wts, beta] = invz_matsubara(T, 40);
tl = closed_tl(h, D0, M0op, beta);
h0el = beta*(1 - tl.n01^2);
g  = real(invz_g(tl, 1i*wn));
G0 = -(tl.M2*g);  G0(1) = -(tl.M2*tl.g0 + tl.m^2*h0el);      % elastic in the static slot
gi = -tl.M2*tl.g0;  ge = -tl.m^2*h0el;
Sigma = zeros(size(wn));  K = zeros(size(wn));  lam = [0; 0; 0];  K0s = 0;  ok = false;
for outer = 1:200
    med = invz_emt_scalar(G0, Sigma, Jnu, struct());
    if ~med.converged, break; end
    K = med.K;
    [K0s, ~, so] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), Jnu, K0s, beta, ...
                                           J0eff, gi, ge, struct('warn', false));
    K(1) = K0s;
    lam = invz_lambdas(K, g, wts, beta, [1 2 3]);
    sg  = invz_sigma_ordered(tl, lam, K, g, beta);
    dS  = max(abs(sg.Sigma - Sigma));
    Sigma = Sigma + 0.7*(sg.Sigma - Sigma);
    if dS < 1e-8 && so.converged, ok = true; break; end
end
% Simultaneous FINAL tuple: refresh the static slot at the final Sigma, keep ITS K0,
% then REVALIDATE lambdas and the ordered Sigma target against the exported K.
[K0s, ~, so] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), Jnu, K0s, beta, ...
                                       J0eff, gi, ge, struct('warn', false));
K(1) = K0s;
lam_chk = invz_lambdas(K, g, wts, beta, [1 2 3]);
sg_chk  = invz_sigma_ordered(tl, lam_chk, K, g, beta);
ok = ok && so.converged && isfinite(so.resid) && so.resid < 1e-10 ...
        && max(abs(sg_chk.Sigma - Sigma)) < 1e-8;
rk = so.r;  Mk = tl.m*tl.n01;  S0k = Sigma(1);  l1g0 = abs(lam_chk(1)*tl.g0);
end

function [dFA, dgn] = closed_routeA(T, Jnu, J0eff, D0, M0op, hmax, nH, testCase)
% Route A on ONE grid, all gates inside -- called three times so the SS7b convergence
% sweeps are EXECUTED assertions, not comments (third-review P1).
hgrid = hmax * (1e-4^(1/(nH-1))).^((nH-1):-1:0);
[r0, ~, ~, ~, ok0] = closed_node(0, T, Jnu, J0eff, D0, M0op);
verifyTrue(testCase, ok0);
[rv, Mv, S0v, Lv] = deal(nan(1, nH));
for k = 1:nH
    [rv(k), Mv(k), S0v(k), Lv(k), okk] = closed_node(hgrid(k), T, Jnu, J0eff, D0, M0op);
    verifyTrue(testCase, okk, sprintf('closed node h = %.4g failed', hgrid(k)));
end
h0 = cumtrapz([0 hgrid], [r0 rv]);  h0 = h0(2:end);
dh = h0 - hgrid;
q = max(1e-4, 0.01*abs(r0 - 1));
verifyTrue(testCase, all(abs(rv(end-4:end) - 1) < q), 'fluctuations not quenched');
verifyTrue(testCase, all(abs(S0v(end-4:end)) < q), 'Sigma(0) not quenched on the tail');
verifyTrue(testCase, all(abs(Lv(end-4:end)) < q), 'lambda1*g0 not quenched on the tail');
verifyLessThan(testCase, M0op - Mv(end), 0.01*M0op);          % analytic saturation reached
dFA = -trapz([0 Mv], [0 dh]);                                  % one moment throughout
dgn = struct('tail_est_A', abs(dh(end))*(M0op - Mv(end)), 'dh_end', dh(end), 'r0', r0);
end

function dU = route_b_sweep(testCase, Tg, D0, M0op, Jnu)
% Route-B sweep (stage-2 task 6b step 6 item 1 -- ported as its own localfunction so the
% order-consistency scan and the approximation-fingerprint pins can reuse it at
% different coupling scales and temperatures without duplicating the PM loop).
dU = nan(size(Tg));
for k = 1:numel(Tg)
    dU(k) = route_b_node(testCase, Tg(k), D0, M0op, Jnu);
end
end

function dU = route_b_node(testCase, T, D0, M0op, Jnu)
% One route-B (h=0, PM) node: dU per J 2.22 (eq 37), G reconstructed via J 2.30, on a
% Matsubara grid rebuilt fresh per temperature (round-2 P0-B).
[wn, wts, beta] = invz_matsubara(T, 40);
tl = closed_tl(0, D0, M0op, beta);
g  = real(invz_g(tl, 1i*wn));
G0 = -(tl.M2*g);                                      % m = 0: no elastic term
Sigma = zeros(size(wn));  ok = false;
for outer = 1:200                                     % PM loop: emt + invz_sigma (lam [1 2])
    med = invz_emt_scalar(G0, Sigma, Jnu, struct());
    if ~med.converged, break; end
    lam = invz_lambdas(med.K, g, wts, beta, [1 2]);
    sg  = invz_sigma(tl, lam, med.K, g, beta);
    dS = max(abs(sg.Sigma - Sigma));  Sigma = Sigma + 0.7*(sg.Sigma - Sigma);
    if dS < 1e-8, ok = true; break; end
end
verifyTrue(testCase, ok, sprintf('route-B node T = %.4g K did not converge', T));
% one SIMULTANEOUS final tuple, recomputed from the converged Sigma (third-review P1)
med = invz_emt_scalar(G0, Sigma, Jnu, struct());
verifyTrue(testCase, med.converged);
lam = invz_lambdas(med.K, g, wts, beta, [1 2]);
sg  = invz_sigma(tl, lam, med.K, g, beta);
G = G0 ./ (1 + Sigma + med.K .* G0);                  % J 2.30
dU = 0.5*( sg.alpha*tl.n01*tl.Delta/(1 + sg.alpha) - tl.M2*lam(1) ...
              + real(sum(wts .* med.K .* (G - G0)))/beta );
end
