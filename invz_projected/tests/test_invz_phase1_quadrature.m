function tests = test_invz_phase1_quadrature
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function ion = make_ion()
ion = invz_ion();
end

function J = frozen_J_ref()
J = 0.006424435656;
end

function c = dummy_couplings(nrows)
% Cheap stand-in for invz_phase1_couplings' output, for tests that exercise items 1/3/4 (which
% depend only on the GRID g, never on c) without paying a real invz_jq_modes evaluation.
c.Jnu_unflat = zeros(nrows, 4);
c.Jnu_flat   = zeros(4*nrows, 1);
c.w_flat     = ones(4*nrows, 1) / (4*nrows);
c.info       = struct();
c.J0eff  = 1;  c.Jcc0 = 1;  c.maxJnu = 1;
end

% ===========================================================================================
% invz_phase1_offsets: structural sanity (deterministic, instant).
% ===========================================================================================
function test_offsets_structure(testCase)
offs = invz_phase1_offsets();
verifyEqual(testCase, numel(offs), 8);
verifyEqual(testCase, {offs.tag}, {'000','h','k','l','hk','hl','kl','hkl'});
verifyEqual(testCase, offs(1).flags, logical([0 0 0]));
verifyEqual(testCase, offs(8).flags, logical([1 1 1]));
verifyEqual(testCase, offs(2).flags, logical([1 0 0]));
verifyEqual(testCase, offs(4).flags, logical([0 0 1]));
% all 8 combinations are distinct
allflags = cat(1, offs.flags);
verifyEqual(testCase, size(unique(allflags, 'rows'), 1), 8);
end

% ===========================================================================================
% Item 1 (uniqueness): the brief's own required numeric coverage -- half-open N=16 exactly
% 4096 distinct; legacy inclusive N=16 exactly 3375 distinct (the documented, expected-FAIL
% baseline). Grid-only (invz_phase1_qgrid never calls invz_jq_modes), so real N=16 is cheap here;
% a lightweight dummy `c` avoids paying for a real coupling evaluation this test does not need.
% ===========================================================================================
function test_halfopen_n16_distinct_4096_every_offset(testCase)
ion = make_ion();
offs = invz_phase1_offsets();
for k = 1:8
    g = invz_phase1_qgrid(ion, 16, offs(k).flags, 'halfopen', 'P_complete');
    res = invz_phase1_checks(ion, g, dummy_couplings(g.nominal), frozen_J_ref());
    verifyEqual(testCase, res.item1.nominal, 4096);
    verifyEqual(testCase, res.item1.distinct, 4096, sprintf('halfopen offset %s', offs(k).tag));
    verifyTrue(testCase, res.item1.pass, sprintf('halfopen offset %s must PASS item 1', offs(k).tag));
end
end

function test_legacy_n16_distinct_3375_documented_baseline(testCase)
% Prereg: "legacy inclusive grid is expected to FAIL (report N^3 nominal vs (N-1)^3 distinct --
% for N=16, 4096 vs 3375); that failure is the documented baseline, not a stop." This builder
% (invz_phase1_qgrid.m) additionally documents that EVERY offset, not just baseline, shows this
% under its ax0/ax1-constant-shift construction (see that file's header) -- asserted here too.
ion = make_ion();
offs = invz_phase1_offsets();
for k = 1:8
    g = invz_phase1_qgrid(ion, 16, offs(k).flags, 'legacy_inclusive', 'P_complete');
    res = invz_phase1_checks(ion, g, dummy_couplings(g.nominal), frozen_J_ref());
    verifyEqual(testCase, res.item1.nominal, 4096);
    verifyEqual(testCase, res.item1.distinct, 3375, sprintf('legacy offset %s', offs(k).tag));
    verifyFalse(testCase, res.item1.pass, sprintf('legacy offset %s is expected to FAIL item 1', offs(k).tag));
end
end

function test_tol_uniq_boundary(testCase)
% Item 1's own tolerance (tol_uniq = 1e-12) classifies a synthetic near-duplicate pair correctly:
% differing by < 1e-12 collapses to 1 distinct point; differing by > 1e-12 stays 2 distinct.
ion = make_ion();
g_dup = struct('qvec', [0.1 0.2 0.3; 0.1 0.2 0.3+5e-13], 'w', [0.5;0.5], 'n_gamma', 0, ...
    'nominal', 2, 'N', NaN, 'offsetFlags', logical([0 0 0]), 'convention', 'halfopen', 'gammaPolicy', 'P_complete');
res_dup = invz_phase1_checks(ion, g_dup, dummy_couplings(2), frozen_J_ref());
verifyEqual(testCase, res_dup.item1.distinct, 1, 'points 5e-13 apart (< tol_uniq) must collapse');

g_apart = struct('qvec', [0.1 0.2 0.3; 0.1 0.2 0.3+5e-12], 'w', [0.5;0.5], 'n_gamma', 0, ...
    'nominal', 2, 'N', NaN, 'offsetFlags', logical([0 0 0]), 'convention', 'halfopen', 'gammaPolicy', 'P_complete');
res_apart = invz_phase1_checks(ion, g_apart, dummy_couplings(2), frozen_J_ref());
verifyEqual(testCase, res_apart.item1.distinct, 2, 'points 5e-12 apart (> tol_uniq) must stay distinct');
end

% ===========================================================================================
% Item 4 (weight normalization): sum(w)==1 within 1e-12, both Gamma policies, several configs.
% Grid-only (cheap even at N=16).
% ===========================================================================================
function test_weights_sum_to_one_both_gamma_policies(testCase)
ion = make_ion();
offs = invz_phase1_offsets();
convs = {'halfopen', 'legacy_inclusive'};
Ns = [12 16];
for ci = 1:numel(convs)
    for ni = 1:numel(Ns)
        for k = [1 4 8]   % a few representative offsets, not all 8x2x2 (still thorough, stays fast)
            gC = invz_phase1_qgrid(ion, Ns(ni), offs(k).flags, convs{ci}, 'P_complete');
            resC = invz_phase1_checks(ion, gC, dummy_couplings(size(gC.qvec,1)), frozen_J_ref());
            verifyTrue(testCase, resC.item4.pass, sprintf('%s N=%d offset %s P_complete', convs{ci}, Ns(ni), offs(k).tag));
            verifyLessThanOrEqual(testCase, resC.item4.abs_err, 1e-12);

            gD = invz_phase1_qgrid(ion, Ns(ni), offs(k).flags, convs{ci}, 'P_drop');
            resD = invz_phase1_checks(ion, gD, dummy_couplings(size(gD.qvec,1)), frozen_J_ref());
            verifyTrue(testCase, resD.item4.pass, sprintf('%s N=%d offset %s P_drop', convs{ci}, Ns(ni), offs(k).tag));
            verifyLessThanOrEqual(testCase, resD.item4.abs_err, 1e-12);
        end
    end
end
end

% ===========================================================================================
% Item 3 (cardinality + Gamma): P_drop removes exactly the P_complete Gamma count and
% renormalizes; post-drop Gamma count is re-verified at 0. Grid-only.
% ===========================================================================================
function test_p_drop_removes_exactly_gamma_count_and_renormalizes(testCase)
ion = make_ion();
offs = invz_phase1_offsets();
for k = 1:8
    gC = invz_phase1_qgrid(ion, 16, offs(k).flags, 'halfopen', 'P_complete');
    gD = invz_phase1_qgrid(ion, 16, offs(k).flags, 'halfopen', 'P_drop');
    verifyEqual(testCase, gC.n_gamma, gD.n_gamma, sprintf('offset %s', offs(k).tag));
    verifyEqual(testCase, size(gD.qvec,1), gC.nominal - gC.n_gamma, sprintf('offset %s', offs(k).tag));
    resD = invz_phase1_checks(ion, gD, dummy_couplings(size(gD.qvec,1)), frozen_J_ref());
    verifyTrue(testCase, resD.item3.pass, sprintf('offset %s P_drop item3', offs(k).tag));
    verifyEqual(testCase, resD.item3.n_gamma_after_drop, 0);
    verifyEqual(testCase, abs(sum(gD.w) - 1) <= 1e-12, true);
end
% halfopen: Gamma present ONLY in the [0 0 0] baseline offset (prereg item 3 expectation).
n_gamma_by_offset = arrayfun(@(o) invz_phase1_qgrid(ion, 16, o.flags, 'halfopen', 'P_complete').n_gamma, offs);
verifyEqual(testCase, n_gamma_by_offset, [1 0 0 0 0 0 0 0]);
end

% ===========================================================================================
% Item 2 (periodicity): sorted-branch spectrum at q+G matches q within tolerance, on a small
% sample. Uses a small dpRng-bearing invz_jq_modes call directly (cheap: 5 q's x 6 G's = 35 pts).
% ===========================================================================================
function test_periodicity_sample_passes(testCase)
ion = make_ion();
res = invz_phase1_check_periodicity(ion, 30);
verifyGreaterThan(testCase, res.n_pairs, 0);   % non-vacuous
verifyTrue(testCase, res.pass);
verifyLessThanOrEqual(testCase, res.max_violation_margin, 0);
verifyEqual(testCase, res.AbsTol_J, 1e-10);
verifyEqual(testCase, res.RelTol_J, 1e-8);
end

% ===========================================================================================
% invz_phase1_gate / invz_phase1_refinement_gate: the "gate is non-vacuous" requirement -- a
% synthetic pair that IS within tolerance and one that IS NOT are both classified correctly, for
% both the raw pairwise gate and the 3-rung refinement gate (finest-comparison + spread
% non-increasing).
% ===========================================================================================
function test_gate_shape_and_energy_classify_synthetic_pairs(testCase)
[p_close, ~, ~] = invz_phase1_gate('shape', 1.000000, 1.0000005, []);
verifyTrue(testCase, p_close);
[p_far, ~, ~] = invz_phase1_gate('shape', 1.0, 2.0, []);
verifyFalse(testCase, p_far);

Jref = frozen_J_ref();
[p_close_e, ~, ~] = invz_phase1_gate('energy', Jref, Jref + 1e-11, Jref);
verifyTrue(testCase, p_close_e);
[p_far_e, ~, ~] = invz_phase1_gate('energy', Jref, Jref*2, Jref);
verifyFalse(testCase, p_far_e);

verifyError(testCase, @() invz_phase1_gate('bogus', 1, 2, Jref), 'invz:phase1Config');
verifyError(testCase, @() invz_phase1_gate('energy', 1, 2, -1), 'invz:phase1Config');
end

function test_refinement_gate_non_vacuous(testCase)
% Converging sequence (finest step smaller than preceding, and within tolerance) -> PASS.
res_pass = invz_phase1_refinement_gate('shape', 1.0000, 1.00001, 1.000011, []);
verifyTrue(testCase, res_pass.finest_pass);
verifyTrue(testCase, res_pass.spread_nonincreasing);
verifyTrue(testCase, res_pass.pass);

% Finest comparison outside tolerance -> finest_pass false -> overall FAIL.
res_fail_tol = invz_phase1_refinement_gate('shape', 1.0, 1.1, 1.5, []);
verifyFalse(testCase, res_fail_tol.finest_pass);
verifyFalse(testCase, res_fail_tol.pass);

% Finest step within tolerance BUT growing relative to the preceding step -> spread_nonincreasing
% false -> overall FAIL even though finest_pass alone would read PASS.
res_fail_spread = invz_phase1_refinement_gate('shape', 1.000000, 1.0000001, 1.000005, []);
verifyTrue(testCase, res_fail_spread.finest_pass, 'finest step alone is within the shape tolerance');
verifyFalse(testCase, res_fail_spread.spread_nonincreasing, 'but the step grew vs the preceding one');
verifyFalse(testCase, res_fail_spread.pass);
end

function test_offset_spread_pairwise_gate_non_vacuous(testCase)
agree8 = 1.0 + 1e-9*(1:8);         % all mutually within the shape tolerance
res_agree = invz_phase1_offset_spread('shape', agree8, []);
verifyTrue(testCase, res_agree.pairwise_pass);

outlier8 = [1.0 1.0 1.0 1.0 1.0 1.0 1.0 5.0];   % one offset far outside tolerance
res_outlier = invz_phase1_offset_spread('shape', outlier8, []);
verifyFalse(testCase, res_outlier.pairwise_pass);
verifyEqual(testCase, res_outlier.spread, 4.0, 'AbsTol', 1e-12);
end

% ===========================================================================================
% Determinism: identical inputs, no RNG/Date anywhere, must yield byte-identical point sets and
% weights across repeated calls (prereg "Execution order"; brief's own explicit test bullet).
% ===========================================================================================
function test_determinism(testCase)
ion = make_ion();
for gp = {'P_complete','P_drop'}
    g1 = invz_phase1_qgrid(ion, 16, [1 0 1], 'legacy_inclusive', gp{1});
    g2 = invz_phase1_qgrid(ion, 16, [1 0 1], 'legacy_inclusive', gp{1});
    verifyEqual(testCase, g1.qvec, g2.qvec);
    verifyEqual(testCase, g1.w, g2.w);
    verifyEqual(testCase, g1.n_gamma, g2.n_gamma);
end
offs1 = invz_phase1_offsets();
offs2 = invz_phase1_offsets();
verifyEqual(testCase, {offs1.tag}, {offs2.tag});
verifyEqual(testCase, cat(1,offs1.flags), cat(1,offs2.flags));
end

% ===========================================================================================
% Offset [0 0 0] reproduces a DIRECT qVec_generator call bit-for-bit for halfopen (the
% construction's central correctness claim, invz_phase1_qgrid.m header) -- NOT for legacy_
% inclusive, where the uniform BZ-wrap deliberately maps the raw +0.5 face to -0.5 (that
% difference IS what lets item 1 detect the historical duplicate face; asserted here too, so this
% documented asymmetry stays regression-tested rather than silently assumed).
% ===========================================================================================
function test_offset_000_matches_direct_call(testCase)
ion = make_ion();
N = 8;
g = invz_phase1_qgrid(ion, N, [0 0 0], 'halfopen', 'P_complete');
qdirect = qVec_generator(ion.a, 'mode','grid','grid',[N N N],'range',[-0.5 0.5], ...
    'endpoint', false, 'verbose', false);
verifyEqual(testCase, g.qvec, qdirect, 'AbsTol', 0, ...
    'halfopen offset 000 must be bit-identical (values AND row order) to a direct qVec_generator call');

gL = invz_phase1_qgrid(ion, N, [0 0 0], 'legacy_inclusive', 'P_complete');
qdirectL = qVec_generator(ion.a, 'mode','grid','grid',[N N N],'range',[-0.5 0.5], ...
    'endpoint', true, 'verbose', false);
verifyNotEqual(testCase, gL.qvec, qdirectL);
verifyEqual(testCase, max(abs(gL.qvec(:) - qdirectL(:))), 1.0, 'AbsTol', 1e-12, ...
    'the only difference from the raw legacy grid must be the wrap of the +0.5 face down to -0.5');
end

% ===========================================================================================
% invz_phase1_couplings: end-to-end sanity on a cheap real config (N=6, dpRng=30) -- Jcc0/J0eff
% agree, maxJnu equals Jcc0 under P_complete (Gamma's own branch is the global max) and is
% strictly smaller under P_drop (Gamma removed). Real invz_jq_modes call, kept small for speed.
% ===========================================================================================
function test_couplings_evaluator_energy_scalars(testCase)
ion = make_ion();
gC = invz_phase1_qgrid(ion, 6, [0 0 0], 'halfopen', 'P_complete');
cC = invz_phase1_couplings(ion, gC, 30);
verifyEqual(testCase, cC.J0eff, cC.Jcc0);
verifyEqual(testCase, cC.maxJnu, cC.Jcc0, 'AbsTol', 1e-9, 'Gamma kept: the grid max must equal Jcc0');
verifyEqual(testCase, numel(cC.Jnu_flat), 4*size(gC.qvec,1));
verifyEqual(testCase, cC.w_flat, repmat(gC.w, 4, 1));

gD = invz_phase1_qgrid(ion, 6, [0 0 0], 'halfopen', 'P_drop');
cD = invz_phase1_couplings(ion, gD, 30);
verifyLessThan(testCase, cD.maxJnu, cD.Jcc0, 'Gamma dropped: the grid max must be strictly below Jcc0');

res = invz_phase1_checks(ion, gC, cC, frozen_J_ref());
verifyEqual(testCase, res.item5.energy.Jcc0, cC.Jcc0);
verifyEqual(testCase, res.item5.norm.mean, res.item5.raw.mean / frozen_J_ref(), 'AbsTol', 1e-15);
verifyEqual(testCase, numel(res.item5.raw.q), 5);
end

% ===========================================================================================
% invz_phase1_qgrid input validation (fast, no grid built).
% ===========================================================================================
function test_qgrid_rejects_bad_inputs(testCase)
ion = make_ion();
verifyError(testCase, @() invz_phase1_qgrid(ion, 16, [0 0 0], 'bogus', 'P_complete'), 'invz:phase1Config');
verifyError(testCase, @() invz_phase1_qgrid(ion, 16, [0 0 0], 'halfopen', 'bogus'), 'invz:phase1Config');
verifyError(testCase, @() invz_phase1_qgrid(ion, 16, [0 0], 'halfopen', 'P_complete'), 'invz:phase1Config');
verifyError(testCase, @() invz_phase1_qgrid(ion, 1.5, [0 0 0], 'halfopen', 'P_complete'), 'invz:phase1Config');
end
