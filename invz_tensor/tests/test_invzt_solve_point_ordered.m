function tests = test_invzt_solve_point_ordered
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
ion = invz_ion();
g   = invzt_qgrid(16, 'halfopen');
tc.TestData.ion = ion;
tc.TestData.lat = invzt_jq_tensor(ion, g, struct('dpRng', 30, 'cache', true));
end

function test_deep_ordered_point(tc)
% Deep FM point: spontaneous moment, converged Sigma, stable (crit >= 0 within tol).
pt = invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0], tc.TestData.lat, struct());
verifyTrue(tc, pt.is_ordered);
verifyTrue(tc, pt.converged);
verifyEqual(tc, pt.moment_branch, 'spontaneous');
verifyGreaterThan(tc, abs(pt.m0), 1.0);              % LiHoF4 FM moment is O(5) here
verifyTrue(tc, isfinite(pt.Sigma0) && all(isfinite(pt.Sigma)));
verifyTrue(tc, isfinite(pt.alpha_m));
verifyGreaterThan(tc, pt.crit, -1e-3);               % broken-symmetry state locally stable
verifyEqual(tc, pt.J0z_used, tc.TestData.lat.info.Jcc0);
verifyEqual(tc, pt.Jxx0_used, tc.TestData.lat.info.Jaa0);
verifyLessThan(tc, pt.sumrule_rel, 0.1);
fprintf('ordered a1 @ 3T: m0=%.4f Sigma0=%.5f crit=%.5f alpha_m=%.4g sumrule=%.3g\n', ...
    pt.m0, pt.Sigma0, pt.crit, pt.alpha_m, pt.sumrule_rel);
end

function test_converged_state_is_self_consistent(tc)
% Review P2-2: on a converged exit, the returned (Sigma, K, lambda, alpha_m) must
% describe the SAME state -- re-evaluating one medium+self-energy pass at pt.Sigma
% must reproduce pt.Sigma to the outer tolerance (check-before-mix loop ordering).
ion = tc.TestData.ion;  lat = tc.TestData.lat;
pt = invzt_solve_point_ordered(ion, 0.1, [3.0 0 0], lat, struct());
verifyTrue(tc, pt.converged);
sg2 = local_one_pass(ion, 0.1, [3.0 0 0], lat, pt);   % helper below: one g(Sigma) map step
verifyLessThan(tc, max(abs(sg2.Sigma - pt.Sigma)), 1e-6);   % >= tol_outer scale, not machine
end

function sg = local_one_pass(ion, T, B, lat, pt)
% One medium + moment-form self-energy evaluation AT pt.Sigma (no mixing): the
% fixed-point map g(Sigma) evaluated at the returned state. WHOLE-CC (2026-07-20
% amendment): no dominant/rest split -- mirrors the solver's own medium step.
[wn, wts, beta] = invz_matsubara(T, 40);
c0 = invz_chi0z(pt.si, T, 1i*wn, struct('elastic', true));
c0_cc = real(squeeze(c0(3,3,:)));
g = real(invz_g(pt.tl, 1i*wn));
ctil = c0 ./ reshape(1 + pt.Sigma, 1, 1, numel(wn));
Gcc = invzt_gcc_lattice(ctil, lat);
Gloc = -Gcc(:);
G0til = -(c0_cc ./ (1 + pt.Sigma));
K = 1 ./ Gloc - 1 ./ G0til;
lam = invz_lambdas(K, g, wts, beta, [1 2 3]);
sg = invz_sigma_ordered(pt.tl, lam, K, g, beta);
end

function test_pm_early_return(tc)
% WELL above the bare-MF boundary (~5.0 T -- NOT the QPT at 4.65-4.70; measured
% m0 = 1.17 at 4.8 T, so 4.8 would NOT early-return) the order-mode MF relaxes
% to ~0: paramagnetic early return.
pt = invzt_solve_point_ordered(tc.TestData.ion, 0.1, [5.5 0 0], tc.TestData.lat, struct());
verifyFalse(tc, pt.is_ordered);
verifyFalse(tc, pt.converged);
verifyEqual(tc, pt.moment_branch, 'none');
verifyTrue(tc, isnan(pt.Sigma0) && isnan(pt.crit));
verifyTrue(tc, isempty(pt.tl));
verifyLessThan(tc, abs(pt.m0), 1e-2);
end

function test_mf_moment_persists_past_QPT(tc)
% Locks the P0-1 finding into the suite: the bare-MF moment is STILL nonzero at
% the 4.8 T PM anchor (measured m0 = 1.1717) -- which is exactly why phase
% selection must NOT be ordered-first (see invzt_solve_auto, Task 5).
pt = invzt_solve_point_ordered(tc.TestData.ion, 0.1, [4.8 0 0], tc.TestData.lat, struct());
verifyTrue(tc, abs(pt.m0) > 1.0);
fprintf('P0-1 lock: bare-MF m0(4.8 T) = %.4f (QPT is at 4.65-4.70 T)\n', pt.m0);
end

function test_longitudinal_rejected(tc)
verifyError(tc, @() invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0.1], ...
    tc.TestData.lat, struct()), 'invzt:orderedLongitudinal');
end

function test_full_dress_a3_rejected(tc)
% Full-dress 'a3' is PERMANENTLY rejected (136-state vertex budget-refused).
% This assertion stays valid after Task 7D adds 'a3d' to the mode gate.
verifyError(tc, @() invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0], ...
    tc.TestData.lat, struct('mode', 'a3')), 'invzt:orderedMode');
end

function test_nlevels_std_only(tc)
verifyError(tc, @() invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0], ...
    tc.TestData.lat, struct('nlevels', 'three')), 'invzt:orderedNlevels');
end

function test_split_knobs_rejected(tc)
% 2026-07-20 amendment: the ordered medium is WHOLE-CC (no dominant/rest split);
% the PM solver's 'Esplit'/'chi_rest' knobs have no meaning here and must fail loud.
verifyError(tc, @() invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0], ...
    tc.TestData.lat, struct('Esplit', 0.4)), 'invzt:orderedSplitKnobs');
verifyError(tc, @() invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0], ...
    tc.TestData.lat, struct('chi_rest', false)), 'invzt:orderedSplitKnobs');
end

% =====================================================================================
%  Task 7D: a3d ordered solve -- full-response fixed-rank dominant-vertex hybrid.
%
%  AFFORDABLE ANCHOR (2026-07-20). The 7D brief proposed T = 1.0 K, B = 2.0 T. MEASURED
%  during 7D: 1.0 K / 2.0 T is ORDERED (m0 = 4.67) but the a1 ordered solver's plain
%  damped-mixing outer loop does NOT converge there -- it oscillates (crit drifts
%  0.63 -> 0.47 -> 0.23 as mix goes 0.7 -> 0.5 -> 0.3, iters hit the cap). The a1
%  ordered solver has no Anderson acceleration (unlike invzt_solve_point), and 2.0 T
%  sits on the near-singular ordered-side RPA where plain Picard fails (same failure
%  mode documented in invzt_solve_point lines 271-277). B = 2.5 T likewise fails.
%  B = 3.0 T CONVERGES cleanly for BOTH a1 (23 iters, crit = 0.224, sumrule = 0.068)
%  and a3d, still ordered (m0 = 3.24), still affordable (T = 1.0 K, nwn = 20). The
%  gates below are UNWEAKENED; only the field moved 2.0 -> 3.0 T so both maps converge.
%
%  COST: the FIRST a3d solve builds the compact cc;cc Gamma4 (~2.9 min at this anchor);
%  it is session-cached (gamma4cc_cached), so the second/third a3d tests reuse it.
% =====================================================================================

function test_ordered_a3d_point(tc)
% a3d Matsubara solve on the affordable anchor: converged physical objects,
% honest surface (no fabricated Jensen fields), consistent criticality.
ion = tc.TestData.ion;
g8  = invzt_qgrid(8, 'halfopen');
lat8 = invzt_jq_tensor(ion, g8, struct('dpRng', 15, 'cache', true));
pt = invzt_solve_point_ordered(ion, 1.0, [3.0 0 0], lat8, ...
    struct('mode', 'a3d', 'Ecut', 10));
if ~pt.is_ordered, fprintf('anchor relaxed PM -- pick a lower B\n'); end
verifyTrue(tc, pt.is_ordered && pt.converged);
verifyTrue(tc, all(isfinite(pt.Vcc)) && all(isfinite(pt.chi_til(:))));
verifyTrue(tc, isnan(pt.alpha) && isnan(pt.alpha_m));        % NOT fabricated
verifyEqual(tc, pt.vb.n_vertex, 16);
verifyTrue(tc, isfinite(pt.eps_el) && isfinite(pt.c_d));
% criticality consistency (7E gate): crit is recomputed from the RETURNED chi_til
cr = invzt_crit_static(real((pt.chi_til(:,:,1) + pt.chi_til(:,:,1)')/2), lat8.JtGamma);
verifyEqual(tc, pt.crit, cr, 'AbsTol', 1e-10);
% real-axis stays rejected (pt.mode 'a3d' ~= 'a1' triggers the A1-only gate)
verifyError(tc, @() invzt_chi_realaxis(ion, 1.0, [3.0 0 0], pt, ...
    linspace(0, 0.6, 11).', struct()), 'invzt:realaxisMode');
fprintf('a3d anchor: m0=%.4f crit=%.5f chi_share=%.4f eps_el=%.3g\n', ...
    pt.m0, pt.crit, pt.vb.chi_share, pt.eps_el);
end

function test_a3d_vs_a1_approximation_control(tc)
% 7E approximation-control gate: a3d vs ordered a1 at the same affordable
% anchor -- REPORT the spread (the beyond-Jensen content + basis truncation),
% and require same-sign, same-order crit (a wildly different crit means a bug,
% not physics). Sigma0 is NaN for a3d (honest surface); the DIAGNOSTIC
% self-energy is pt.Sigma_cc_equiv, compared against a1's pt.Sigma0.
ion = tc.TestData.ion;
g8  = invzt_qgrid(8, 'halfopen');
lat8 = invzt_jq_tensor(ion, g8, struct('dpRng', 15, 'cache', true));
p1 = invzt_solve_point_ordered(ion, 1.0, [3.0 0 0], lat8, struct('Ecut', 10));
p3 = invzt_solve_point_ordered(ion, 1.0, [3.0 0 0], lat8, ...
    struct('mode', 'a3d', 'Ecut', 10));
verifyTrue(tc, p1.converged && p3.converged);
verifyEqual(tc, sign(p3.crit), sign(p1.crit));
verifyLessThan(tc, abs(p3.crit - p1.crit), max(0.5*abs(p1.crit), 0.05));
fprintf('a3d vs a1: crit %.5f / %.5f, dSigma_cc(0) = %.4g\n', ...
    p3.crit, p1.crit, real(p3.Sigma_cc_equiv(1)) - p1.Sigma0);
end

function test_a3d_self_generation(tc)
% 7D-rework HARD GATE 1: the returned Vcc must be GENERATED BY the returned Kmat.
% The a3d fixed point iterates V on the COMPLETE hybrid map, so at convergence
% Vcc(n) == (1/2beta) sum_l Gamma16_cc;cc(n,l) * Kmat_cc(l) EXACTLY -- recompute the
% compact cc;cc contraction in-test from pt.a3d.G4cc + pt.Kmat and match pt.Vcc.
% (Against the pre-rework isolated-map + post-hoc-EMT build this FAILS by O(1): the
% returned Kmat is NOT the medium that generated the returned Vcc.)
ion = tc.TestData.ion;
g8  = invzt_qgrid(8, 'halfopen');
lat8 = invzt_jq_tensor(ion, g8, struct('dpRng', 15, 'cache', true));
pt = invzt_solve_point_ordered(ion, 1.0, [3.0 0 0], lat8, ...
    struct('mode', 'a3d', 'Ecut', 10));
[~, ~, beta] = invz_matsubara(1.0, 10);
vgen = a3d_contract_cc(pt.a3d.G4cc, pt.Kmat, pt.a3d.Lmax, beta);
verifyLessThan(tc, max(abs(vgen - pt.Vcc)), 1e-6);
fprintf('self-gen: max|contract(Gamma16,Kmat) - Vcc| = %.3e\n', max(abs(vgen - pt.Vcc)));
end

function test_a3d_fixed_point_consistency(tc)
% 7D-rework HARD GATE 2: one re-evaluation of the COMPLETE hybrid map at the returned
% state reproduces Vcc/Kmat/chi_til within tolerance. opts.reeval = pt seeds the vertex
% Vmat fixed point with the returned converged Vmat and runs a single outer pass of the
% FULL map (dress + chi_base + EMT + contraction), returning the one-pass images.
ion = tc.TestData.ion;
g8  = invzt_qgrid(8, 'halfopen');
lat8 = invzt_jq_tensor(ion, g8, struct('dpRng', 15, 'cache', true));
pt = invzt_solve_point_ordered(ion, 1.0, [3.0 0 0], lat8, ...
    struct('mode', 'a3d', 'Ecut', 10));
re = invzt_solve_point_ordered(ion, 1.0, [3.0 0 0], lat8, ...
    struct('mode', 'a3d', 'Ecut', 10, 'reeval', pt));
verifyLessThan(tc, max(abs(re.Vcc - pt.Vcc)), 1e-6);
verifyLessThan(tc, max(abs(re.Kmat(:) - pt.Kmat(:))), 1e-6);
verifyLessThan(tc, max(abs(re.chi_til(:) - pt.chi_til(:))), 1e-6);
end

function test_a3d_reduction_identity(tc)
% 7D-rework HARD GATE 3 (+ legacy regression for the invzt_sigma_tensor chi_base
% extension): with chi_base = 0 the COMPLETE-map machinery reduces EXACTLY to the
% pre-rework isolated 16-state solve. Call the vertex solve twice -- chi_base absent
% vs chi_base = zeros(3,3,nwn) -- and require bit-scale agreement (Vmat/chi_til/Kmat).
ion = tc.TestData.ion;
g8  = invzt_qgrid(8, 'halfopen');
lat8 = invzt_jq_tensor(ion, g8, struct('dpRng', 15, 'cache', true));
[si_vb, wn, beta, ~, vb] = a3d_ingredients(ion, 1.0, 3.0, 10, lat8);
verifyEqual(tc, vb.n_vertex, 16);
o0 = struct('dress', 'dominant', 'dom_basis', vb, 'rank_tol', 1e-12);
st_none = invzt_sigma_tensor(si_vb, 1.0, lat8, wn, beta, o0);
oz = o0;  oz.chi_base = complex(zeros(3, 3, numel(wn)));
st_zero = invzt_sigma_tensor(si_vb, 1.0, lat8, wn, beta, oz);
verifyLessThan(tc, max(abs(st_zero.Vmat(:)    - st_none.Vmat(:))),    1e-12);
verifyLessThan(tc, max(abs(st_zero.chi_til(:) - st_none.chi_til(:))), 1e-12);
verifyLessThan(tc, max(abs(st_zero.Kmat(:)    - st_none.Kmat(:))),    1e-12);
end

function test_a3d_map_choice_report(tc)
% 7D-rework GATE 4 (REPORT, no threshold): the measured effect of iterating V on the
% COMPLETE hybrid map vs the pre-rework isolated map (+ post-hoc EMT), at the anchor --
% dVcc(0) relative and dcrit. This records what the map choice changes (the spectator
% feedback the rework restores). Both solves must converge.
ion = tc.TestData.ion;
g8  = invzt_qgrid(8, 'halfopen');
lat8 = invzt_jq_tensor(ion, g8, struct('dpRng', 15, 'cache', true));
[si_vb, wn, beta, chi_base, vb, ~, si] = a3d_ingredients(ion, 1.0, 3.0, 10, lat8);
chi_full = invz_chi0z(si, 1.0, 1i*wn, struct('elastic', true));
chi_dom  = invz_chi0z(si_vb, 1.0, 1i*wn, struct('elastic', true));
o0 = struct('dress', 'dominant', 'dom_basis', vb, 'rank_tol', 1e-12);
st_iso = invzt_sigma_tensor(si_vb, 1.0, lat8, wn, beta, o0);     % as-built isolated map
oh = o0;  oh.chi_base = chi_base;
st_hyb = invzt_sigma_tensor(si_vb, 1.0, lat8, wn, beta, oh);     % complete hybrid map
verifyTrue(tc, st_iso.converged && st_hyb.converged);
chi_ab  = chi_full + (st_iso.chi_til - chi_dom);                 % as-built returned chi_til
crit_ab = invzt_crit_static(real((chi_ab(:,:,1)+chi_ab(:,:,1)')/2), lat8.JtGamma);
crit_h  = invzt_crit_static(real((st_hyb.chi_til(:,:,1)+st_hyb.chi_til(:,:,1)')/2), lat8.JtGamma);
Vab = squeeze(st_iso.Vmat(3,3,:));  Vh = squeeze(st_hyb.Vmat(3,3,:));
fprintf('map choice (isolated -> complete): dVcc(0)_rel = %.4g, dcrit = %.5f (crit %.5f -> %.5f)\n', ...
    abs(Vh(1)-Vab(1))/max(abs(Vab(1)),1e-30), crit_h-crit_ab, crit_ab, crit_h);
end

function test_a3d_production_point(tc)
% PRODUCTION a3d point, INVZ_SLOW-gated (CORE Incomplete allowlist member).
% production build ~= 7.4 h -- launch via nohup outside an interactive session;
% see ODD-LOG SS A-ordered. NOT run in an interactive/CI harness (the build dies
% silently mid-run). The 7C compact vertex + budget guards keep it feasible offline.
% assumeTrue (not early-return) so the skip REGISTERS as Incomplete -- repo slow-gate
% convention (test_invzt_a4_ladder, test_invzt_critical_T/_parity).
assumeTrue(tc, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = tc.TestData.ion;
g16 = invzt_qgrid(16, 'halfopen');
lat16 = invzt_jq_tensor(ion, g16, struct('dpRng', 30, 'cache', true));
pt = invzt_solve_point_ordered(ion, 0.1, [3.0 0 0], lat16, struct('mode', 'a3d'));
verifyTrue(tc, pt.is_ordered && pt.converged);
verifyEqual(tc, pt.vb.n_vertex, 16);
verifyTrue(tc, all(isfinite(pt.chi_til(:))) && isfinite(pt.crit));
fprintf('a3d PRODUCTION @ 0.1 K/3 T: m0=%.4f crit=%.5f chi_share=%.4f Sigma_cc(0)=%.5f\n', ...
    pt.m0, pt.crit, pt.vb.chi_share, real(pt.Sigma_cc_equiv(1)));
end

% ------------------------------- a3d local helpers ----------------------------------
function v = a3d_contract_cc(G4cc, Kmat, Lmax, beta)
% Recompute the compact cc;cc vertex contraction (the SAME formula
% invzt_sigma_tensor's contract_vertex_cc uses): V_cc(n) = (1/2beta) sum_l Kcc(l)
% G4cc(n,l), Kcc(l) = Kmat(3,3,|l|+1) (locked transpose relation, identity for cc).
nl = 2*Lmax + 1;  Kcc = complex(zeros(1, nl));
for li = 1:nl, l = li - Lmax - 1;  Kcc(li) = Kmat(3, 3, abs(l) + 1); end
v = sum(G4cc .* Kcc, 2) / (2*beta);   % [nwn,1]
end

function [si_vb, wn, beta, chi_base, vb, lat_eff, si] = a3d_ingredients(ion, T, B, Ecut, lat)
% Reconstruct the a3d branch's vertex-solve ingredients (front half + hybrid offset),
% so the gate tests can call invzt_sigma_tensor directly on the 16-state basis.
J0z = lat.info.Jcc0;  Jxx0 = lat.info.Jaa0;
[wn, ~, beta] = invz_matsubara(T, Ecut);
si = invz_single_ion(ion, T, [B 0 0], struct('hyp', true, 'order', true, ...
    'J0z', J0z, 'Jxx0', Jxx0, 'transverse_mf', 'legacy_x'));
vb = invzt_ordered_vertex_basis(ion, T, si, struct());
si_vb = vb.si_vertex;
chi_full = invz_chi0z(si,    T, 1i*wn, struct('elastic', true));
chi_dom  = invz_chi0z(si_vb, T, 1i*wn, struct('elastic', true));
chi_base = chi_full - chi_dom;
lat_eff  = lat;
end
