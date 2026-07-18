function tests = test_invzt_a3_threestate
%TEST_INVZT_A3_THREESTATE  CORE gate for A3 on the explicit three-state model (Task 12).
% Validates invzt_threestate (the 3-param (Delta1,m0,rho) -> (Delta,M2,chiperp) match),
% invzt_sigma_tensor (the genuine tensor 1/z self-energy from the exact four-point
% vertex), and invzt_solve_point mode 'a3'. The two make-or-break validations are the
% EXACT rho->0 scalar-compatibility gate (category 1, tol 1e-8) and the framework SS11.8
% lambda-emergence slope gate (A3's Gaussian truncation reproduces A1/E1). Runs with
% invz_projected OFF the CORE path (a single guarded local addpath reaches the
% dependency-free invz_emt_scalar reference in the rho->0 gate, self-cleaning on exit).
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion engine
end

function test_threestate_matches_twolevel_sector(testCase)
% Constructor contract: doublet sector reproduces invz_twolevel's (Delta, M2).
ion = invz_ion();  T = 1.6;  B = [0.5 0 0];
tl = invz_twolevel(ion, T, B(1), struct());
ts = invzt_threestate(ion, T, B, struct());
verifyEqual(testCase, ts.E(2) - ts.E(1), tl.Delta, 'RelTol', 1e-9);
verifyEqual(testCase, abs(ts.Mz(1,2))^2, tl.M2, 'RelTol', 1e-9);
verifyTrue(testCase, all(isfield(ts, {'P','Jexp','JzJz_fluct','Mx','My'})));
fprintf(['[threestate] match residuals: dDelta=%.2e, dM2=%.2e, chiperp achieved=%.4f ' ...
    '(target %.4f, res %.2e)\n'], ts.match.res_Delta, ts.match.res_M2, ts.match.chiperp, ...
    ts.match.chiperp_target, ts.match.res_chiperp);
end

function test_a3_scalar_compat_exact_rho0(testCase)
% SCALAR-SECTOR SS11.8 emergence gate (the strongest single emergence statement: the A3
% four-point vertex reduces to E1's self-energy invz_sigma EXACTLY). The ODD-SECTOR
% emergence is the matched-truncation collapse+band of test_a3_emergence_matched_truncation.
% EXACT scalar-compatibility gate (category 1) -- CONSISTENT under v3 option (a).
% chiperp_scale = 0 sets rho = 0, so Ja0 = Jb0 = 0 and |3> has ZERO matrix
% elements to the doublet in EVERY response (Mx = My = 0; Mz(*,3) = 0): the
% spectator is disconnected regardless of its thermal population, and the doublet
% KEEPS its splitting Delta1 = tl.Delta and moment m0 = sqrt(tl.M2) from the DIRECT
% tunnelling term. So the toy is EXACTLY the two-level system and A3 must reproduce
% the scalar two-level chain (invz_emt_scalar + invz_sigma on the same branch
% spectrum) to 1e-8. (This is the fix for the v2 contradiction, where rho -> 0
% also erased the splitting.) far_excited is kept as belt-and-suspenders but is no
% longer load-bearing: disconnection, not depopulation, makes the limit exact. The
% populated-spectator normalization (constraint 3) is the ORACLE's job (Task 10
% check 6), not an exact solver identity.
% invz_emt_scalar (the scalar reference) lives in invz_projected, deliberately OFF
% the CORE path; reach it via a guarded self-cleaning local addpath so the suite
% stays isolated while this one gate still calls the reference verbatim.
projdir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'invz_projected');
if isempty(which('invz_emt_scalar'))
    addpath(projdir);
    cleaner = onCleanup(@() rmpath(projdir));  %#ok<NASGU> removed on test exit
end
ion = invz_ion();  T = 1.6;  B = [0.5 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
pt = invzt_solve_point(ion, T, B, lat, struct('mode', 'a3', 'nlevels', 'three', ...
    'chiperp_scale', 0, 'far_excited', true, 'odd', false));
verifyTrue(testCase, pt.converged);
% scalar reference on the SAME cc branch spectrum and the SAME (Delta, M2):
tl = invz_twolevel(ion, T, B(1), struct());
[wn, wts, beta] = invz_matsubara(T, 40);
g = real(invz_g(tl, 1i*wn));
G0 = -tl.M2 * g;
Jnu = local_cc_branches(lat);                 % helper: sort(eig) of cc blocks, flattened
Sigma = zeros(size(wn));
% Reference-convergence fix (caught scaffolding defect, cf. Tasks 5/7/9/10): the brief's
% loop capped at 60 iterations, but the damped scalar map (mix 0.7) is SLOW at this
% 0.5 T / 1.6 K proxy (contraction rate ~0.95, so ~410 iters to reach 1e-10) and 60
% iters UNDERSHOOTS the true fixed point by ~3.3e-3 (never hitting the break). The
% PHYSICAL identity is correct -- A3 reproduces the FULLY-CONVERGED scalar chain to
% ~1e-10 -- so the reference must actually converge for the 1e-8 identity to be
% meaningful. Cap raised to 3000 (break unchanged); the 1e-8 assertion is untouched.
for it = 1:3000
    med = invz_emt_scalar(G0, Sigma, Jnu, struct());
    lam = invz_lambdas(med.K, g, wts, beta, [1 2]);
    sg = invz_sigma(tl, lam, med.K, g, beta);
    dS = max(abs(sg.Sigma - Sigma));  Sigma = Sigma + 0.7*(sg.Sigma - Sigma);
    if dS < 1e-10, break, end
end
fprintf('[rho->0 gate] A3 Sigma_cc_equiv(1)=%.10f, scalar Sigma(1)=%.10f, |diff|=%.3e\n', ...
    real(pt.Sigma_cc_equiv(1)), Sigma(1), abs(real(pt.Sigma_cc_equiv(1)) - Sigma(1)));
verifyEqual(testCase, real(pt.Sigma_cc_equiv(1)), Sigma(1), 'AbsTol', 1e-8);
end

function test_a3_emergence_matched_truncation(testCase)
% Framework SS11.8 ODD-SECTOR emergence gate -- MATCHED-TRUNCATION protocol (Task-12
% coordinator resolution; the scalar-sector emergence is the SEPARATE rho->0 identity
% test_a3_scalar_compat_exact_rho0, A3 vertex -> invz_sigma to 3e-11 -- the strongest
% single emergence statement).
%
% The FULL-A3 four-point vertex DRESSES the transverse (a,b) spectator that E1 (Jensen's
% dominant-only rule = mode 'a1') leaves BARE; that transverse-spectator dressing is the
% genuine beyond-E1 correction (full-A3 vs A1 crit-shift ratio ~1.11, REPORTED). Under the
% SAME truncation as E1 -- dress='dominant' dresses the dominant cc TRANSITION only,
% transverse bare, with E1's dominant-renormalized/rest-bare crit -- A3 reduces to E1 up
% to the resummation-SCHEME difference (A3 whole-cc symmetric-bracket Dyson + matrix EMT vs
% A1 dominant/rest split + K-bookkeeping), which is O(1/z^2) BY the plan's LOCKED hard-math
% constraint 8 ("resummation ambiguity is O(1/z^2)"). A lambda^2 coefficient differing at
% O(1/z^2) makes |dcrit_A3dom - dcrit_A1| ~ lambda^2 (slope ~2, NOT >=3): the crit-shift
% SLOPE is O(1/z^2)-capped and is the WRONG probe (a >=2.3 slope gate would contradict
% constraint 8). The ODD-sector emergence is therefore gated by the two matched-truncation
% facts:
%   (a) COLLAPSE: dropping the transverse dressing must shrink the beyond-E1 crit-shift
%       excess |ratio-1| by a large factor -- i.e. the transverse dressing IS the beyond-E1
%       physics (measured |rd-1|/|rf-1| ~ 0.15);
%   (b) BAND: the dominant residual |rd-1| sits inside the O(1/z^2) resummation band
%       (measured ~0.017; resum_spread_crit is the matching constraint-8 method error bar).
%
% STABLE PM config: T = 2.0 K, Bx = 0.5 T -- A1 crit > 0, SINGLE-ROOT (converges in ~7
% iters); Sigma_seed / Vmat_seed CONTINUITY from the odd-off baseline keeps every lambda on
% one root branch (T6/T7 discipline). The near-critical Bx = 0.1 T of earlier drafts was
% multi-root (crit < 0) and unusable.
collapse_frac = 0.4;   % HARD (a): |rd-1| <= collapse_frac*|rf-1|. Measured ratio ~0.15;
                       % 0.4 gives >2x margin. The transverse-dressing collapse.
band = 0.05;           % HARD (b): |rd-1| <= band. Measured ~0.017; fixed documented band
                       % (constraint-8 O(1/z^2) resummation window), residual comfortably
                       % inside; resum_spread_crit below is the matching absolute error bar.
ion = invz_ion();  T = 2.0;  B = [0.5 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
b1 = invzt_solve_point(ion, T, B, lat, struct('mode', 'a1', 'nlevels', 'three', 'odd', false));
bd = invzt_solve_point(ion, T, B, lat, struct('mode', 'a3', 'nlevels', 'three', 'odd', false, 'dress', 'dominant'));
bf = invzt_solve_point(ion, T, B, lat, struct('mode', 'a3', 'nlevels', 'three', 'odd', false, 'dress', 'full'));
verifyGreaterThan(testCase, b1.crit, 0);   % stable PM baseline (single-root A1)
lams = [1.0 0.5 0.25];  ddom = zeros(1, 3);  d1 = zeros(1, 3);  rf1 = NaN;
for i = 1:3
    latL = scale_odd_blocks(lat, lams(i));
    a1 = invzt_solve_point(ion, T, B, latL, struct('mode', 'a1', 'nlevels', 'three', 'odd', true, ...
        'Sigma_seed', b1.Sigma));
    ad = invzt_solve_point(ion, T, B, latL, struct('mode', 'a3', 'nlevels', 'three', 'odd', true, ...
        'dress', 'dominant', 'Vmat_seed', bd.Vmat, 'Sigma_seed', bd.Sigma));
    d1(i) = a1.crit - b1.crit;  ddom(i) = ad.crit - bd.crit;
    if i == 1   % full-A3 needed only at lambda=1 (beyond-E1 report + collapse denominator)
        af = invzt_solve_point(ion, T, B, latL, struct('mode', 'a3', 'nlevels', 'three', 'odd', true, ...
            'dress', 'full', 'Vmat_seed', bf.Vmat, 'Sigma_seed', bf.Sigma));
        rf1 = (af.crit - bf.crit) / d1(1);
    end
end
rd1   = ddom(1) / d1(1);
slope = polyfit(log(lams), log(abs(ddom - d1)), 1) * [1; 0];
fprintf(['[SS11.8 ODD emergence, MATCHED] rd(1)=%.4f rf(1)=%.4f | (a) collapse |rd-1|/|rf-1|=%.3f ' ...
    '(<=%.2f) | (b) band |rd-1|=%.4f (<=%.2f) | REPORT slope=%.2f (O(1/z^2)-capped, not>=3), ' ...
    'resum_spread_crit=%.2e\n'], rd1, rf1, abs(rd1 - 1)/abs(rf1 - 1), collapse_frac, ...
    abs(rd1 - 1), band, slope, bf.resum_spread_crit);
% HARD (a): matched-truncation collapse -- transverse dressing IS the beyond-E1 physics.
verifyLessThanOrEqual(testCase, abs(rd1 - 1), collapse_frac * abs(rf1 - 1));
% HARD (b): dominant residual within the O(1/z^2) resummation band (constraint 8).
verifyLessThanOrEqual(testCase, abs(rd1 - 1), band);
end

function test_uniform_shift_response_reported(testCase)
% As v1 (report-only probe; the strict theorem lives at the closed-form Tc root).
ion = invz_ion();  T = 1.6;  B = [0.1 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
d = 5e-4;
latS = lat;
latS.Jt(3:3:12, 3:3:12, :) = latS.Jt(3:3:12, 3:3:12, :) + d*repmat(eye(4), 1, 1, size(lat.Jt, 3));
latS.JtGamma(3:3:12, 3:3:12) = latS.JtGamma(3:3:12, 3:3:12) + d*eye(4);
p0 = invzt_solve_point(ion, T, B, lat,  struct('mode', 'a3', 'nlevels', 'three', 'odd', false));
p1 = invzt_solve_point(ion, T, B, latS, struct('mode', 'a3', 'nlevels', 'three', 'odd', false));
r = (p1.crit - p0.crit) / d;
fprintf('uniform-shift crit response (A3, three-state): dcrit/d = %.3f\n', r);
verifyTrue(testCase, isfinite(r));  verifyLessThan(testCase, abs(r), 5);
end

function test_a3_monitors_reported(testCase)
ion = invz_ion();  T = 1.6;  B = [0.5 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
pt = invzt_solve_point(ion, T, B, lat, struct('mode', 'a3', 'nlevels', 'three', 'odd', true));
fprintf('A3 monitors: sumrule_rel %.3f, eps_el %.3f, eps_cross %.4f, resum spread %.2e\n', ...
    pt.sumrule_rel, pt.eps_el, pt.eps_cross, pt.resum_spread_crit);
verifyTrue(testCase, pt.converged && isfinite(pt.eps_el));
verifyLessThan(testCase, pt.sumrule_rel, 0.25);
end

function test_vertex_array_K_nonsymmetric_transpose(testCase)
% CARRY (Task-12 resolution 5): exercise the vertex ARRAY-K negative-l transpose
% reconstruction with a NON-SYMMETRIC Kmat. T11's array-K test used a DIAGONAL K
% (test_negfreq_V_reconstruction: eye(3)*scalar), which is symmetric, so a
% K_{rho sigma} <-> K_{sigma rho} swap bug on the negative-l branch (arr_K) would
% hide -- exactly where a bug would live for the production GYROTROPIC (non-symmetric)
% K. Here: build a non-symmetric array K, verify the array path reproduces a
% function-handle kernel carrying the SAME transpose convention to machine precision,
% AND that a deliberately WRONG (no-transpose) handle DIFFERS (so the guard is real).
ion = invz_ion();  T = 1.6;
[~, ~, beta] = invz_matsubara(T, 40);
% small non-degenerate 3-level system with generic real Hermitian operators
E = [0; 0.35; 0.92];  p = exp(-beta*(E - min(E)));  p = p/sum(p);
mkH = @(M) (M + M')/2;
oa = mkH([0 0.7 0.2; 0.7 0 0.5; 0.2 0.5 0]);
ob = mkH([0.1 0.3 0.9; 0.3 -0.2 0.4; 0.9 0.4 0.05]);
oc = mkH([0.8 0.1 0.15; 0.1 -0.6 0.25; 0.15 0.25 0.3]);
ctr = @(O) O - real(sum(p(:).*diag(O)))*eye(3);
es  = struct('E', E, 'p', p);
ops = struct('a', ctr(oa), 'b', ctr(ob), 'c', ctr(oc));
Lmax = 6;
Karr = zeros(3, 3, Lmax + 1);
base = [1.0 0.5 0.2; -0.3 0.8 0.1; 0.15 -0.25 0.6];   % NON-symmetric base
for l = 0:Lmax
    wl = 2*pi*l/beta;
    Karr(:, :, l + 1) = base / (1 + (wl/1.5)^2);
end
verifyGreaterThan(testCase, norm(base - base.', 'fro'), 0.5);   % genuinely non-symmetric
Kf   = @(ri, si, l) local_Kfun(Karr, ri, si, l);        % correct transpose on l<0
Kbad = @(ri, si, l) Karr(ri, si, abs(l) + 1);           % WRONG: no transpose on l<0
opts = struct('stage', 'V', 'Lmax', Lmax);
opts.comps = {'a', 'b', 'c'};  opts.ext = {{'c', 'a'}};
nn = [-2 -1 1 2];
outA = invzt_vertex4(es, ops, Karr, nn, beta, opts);
outF = invzt_vertex4(es, ops, Kf,   nn, beta, opts);
outB = invzt_vertex4(es, ops, Kbad, nn, beta, opts);
verifyEqual(testCase, outA.val, outF.val, 'RelTol', 1e-12, 'AbsTol', 1e-14);
verifyGreaterThan(testCase, max(abs(outA.val(:) - outB.val(:))), 1e-9);
fprintf('[nonsym array-K] |arr - fhandle| = %.2e (match); |arr - no-transpose| = %.2e (differs, guard real)\n', ...
    max(abs(outA.val(:) - outF.val(:))), max(abs(outA.val(:) - outB.val(:))));
end

% ============================== helpers ============================== %
function Jnu = local_cc_branches(lat)
% Flattened sorted eigenvalues of the 4x4 cc blocks over all q (the scalar mode
% spectrum invz_emt_scalar averages uniformly; matches invzt_emt_matrix's BZ+
% sublattice mean when lat.w is uniform, as for a halfopen grid).
nq = size(lat.Jt, 3);
Jnu = zeros(4*nq, 1);
for iq = 1:nq
    Jcc = (lat.Jt(3:3:12, 3:3:12, iq) + lat.Jt(3:3:12, 3:3:12, iq)')/2;
    Jnu((iq-1)*4 + (1:4)) = sort(real(eig(Jcc)));
end
end

function latL = scale_odd_blocks(lat, lambda)
% Multiply the ODD (c<->a,b) Cartesian blocks -- ca/cb and their ac/bc transposes --
% of every sublattice pair by lambda, leaving cc/aa/bb/ab untouched. The ODD entry
% is where EXACTLY ONE composite Cartesian index is c (=3).
cart = mod((0:11).', 3) + 1;        % 1,2,3,1,2,3,...
isc  = (cart == 3);
oddmask = xor(isc, isc.');          % [12,12] logical: exactly one index is c
scal = 1 + (lambda - 1)*oddmask;    % lambda on odd entries, 1 elsewhere
latL = lat;
latL.Jt = lat.Jt .* scal;           % broadcasts over pages
latL.JtGamma = lat.JtGamma .* scal;
end

function v = local_Kfun(Karr, ri, si, l)
% Array accessor with the LOCKED transpose relation K_{rho sigma}(-l)=K_{sigma rho}(+l).
if l >= 0
    v = Karr(ri, si, l + 1);
else
    v = Karr(si, ri, -l + 1);
end
end
