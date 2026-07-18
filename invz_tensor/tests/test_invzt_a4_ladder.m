function tests = test_invzt_a4_ladder
%TEST_INVZT_A4_LADDER  CORE gate for A4 -- the basis-defined state-space ladder (Task 13).
% Validates invzt_rung_basis (multiplet-complete electronic CF subspaces + the xI8
% nuclear product), invzt_run_ladder (the DATA-ONLY, budget-refusing driver), and the
% invzt_gamma4 == invzt_vertex4 stage-'Gamma' cross-channel equivalence the ladder relies
% on at larger N (the cached Gamma4 precompute is quoted as physics). Runs with
% invz_projected OFF the CORE path. The production slow ladder (e6 at ~10 h) is gated
% INVZ_SLOW and left to the user per the Task-13 budget policy.
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion engine
end

function test_ladder_fast_two_rungs(testCase)
ion = invz_ion();
out = invzt_run_ladder(ion, struct('rungs', {{'three', 'e3'}}));
verifyEqual(testCase, numel(out.rungs), 2);
verifyTrue(testCase, all(isfinite(out.crit_shift_odd)));
verifyTrue(testCase, all(out.converged));
verifyTrue(testCase, all(out.chi0_virtual_deficit >= 0));
fprintf('ladder fast: crit shifts %s\n', mat2str(out.crit_shift_odd(:).', 4));
end

function test_rung_basis_multiplet_complete(testCase)
% v3 (review Other 6): assert COMPLETE-MULTIPLET inclusion and record the ACTUAL
% dimension, rather than hard-coding the nominal label size. At Bx = 0 the CF
% ground doublet is degenerate; the e3 basis must hold BOTH doublet states (no
% split multiplet at the cut) and expose dim_actual + multiplet_complete. The
% nuclear product multiplies by the full I=7/2 space (8).
ion = invz_ion();
rb = invzt_rung_basis(ion, 'e3');
verifyTrue(testCase, isfield(rb, 'dim_actual') && rb.dim_actual == size(rb.projector, 2));
verifyGreaterThanOrEqual(testCase, rb.dim_actual, 3);                 % nominal 3; larger only if a multiplet completes
verifyLessThan(testCase, abs(rb.E_basis(2) - rb.E_basis(1)), 1e-6);   % ground doublet intact (not split by the cut)
verifyTrue(testCase, rb.multiplet_complete);                         % no partial multiplet at the top edge
rb8 = invzt_rung_basis(ion, 'e3xI8');
verifyEqual(testCase, rb8.dim_actual, 8*rb.dim_actual);              % complete nuclear product
end

function test_ladder_production_slow(testCase)
% HEADLINE (report): rung table at production settings + small-B-proxy Tc for
% the largest affordable rung. Cross-validation comparators (REPORT, never
% tune): projected Tier-1+2 = 1.509 K, DeltaTc = 0.2341 K; A3-beyond-Gaussian
% share vs the projected Tier-2 share (~2.8%).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
out = invzt_run_ladder(ion, struct('rungs', {{'three','e3','e6','e17','e3xI8'}}, ...
    'production', true, 'tc', true, 'budget_hours', 12));
verifyTrue(testCase, isstruct(out) && ~isempty(out.rungs));
end

function test_sigma_tensor_assembles_dim6(testCase)
% REGRESSION GUARD (coordinator review): the A3 assembly must be N-adaptive in the
% STATE space (numel(si.E)), NOT the toy's 3. A hardcoded eye(3) operator-centering in
% invzt_sigma_tensor crashed the full-A3 solve for e6 (6x6 si.Mx) with
% MATLAB:sizeDimensionsMustMatch -- a bug that hid behind INVZ_SLOW because T12 only
% exercised N=3. Feed a SYNTHETIC dim-6 si through invzt_sigma_tensor (full dress) on a
% TINY Matsubara grid and assert it ASSEMBLES (Vmat [3,3,nwn], finite) -- a cheap
% structural guard (no full/production solve), the regression that should have existed.
ion = invz_ion();  T = 1.6;  C = invz_const();  beta = 1/(C.kB*T);
g = invzt_qgrid(4, 'halfopen');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 8, 'cache', false));
N = 6;                                           % > 3: the dimension that used to crash
E = [0; 0.02; 0.5; 0.9; 1.2; 1.6];               % distinct energies (ground near-pair + spread)
p = exp(-beta*(E - E(1)));  p = p/sum(p);
rng(20260718);
si = struct('E', E, 'P', p, 'Mx', local_herm(N), 'My', local_herm(N), 'Mz', local_herm(N));
si.Jexp = [real(sum(p.*diag(si.Mx))); real(sum(p.*diag(si.My))); real(sum(p.*diag(si.Mz)))];
si.JzJz_fluct = real(sum(p.*diag(si.Mz*si.Mz))) - si.Jexp(3)^2;
wn = 2*pi*(0:1).'/beta;                          % nwn = 2 -> Lmax = 1, cheap (n,l) grid
out = invzt_sigma_tensor(si, T, lat, wn, beta, struct('max_outer', 3, 'dress', 'full'));
verifyEqual(testCase, size(out.Vmat), [3 3 numel(wn)]);   % external Cartesian tensor stays 3x3
verifyTrue(testCase, all(isfinite(out.Vmat(:))));
fprintf('[dim-6 guard] invzt_sigma_tensor assembled N=%d: Vmat %s, max|V|=%.3g\n', ...
    N, mat2str(size(out.Vmat)), max(abs(out.Vmat(:))));
end

function test_gamma4_equals_vertex4_cross_channels(testCase)
% Task-12 review Minor #1 (brief resolution 8): the rho->0 scalar gate exercised only
% the cc channel, but the ladder uses invzt_gamma4 (the cached Gamma4 precompute) at
% larger N with the FULL-A3 rf quoted as physics. Assert gamma4 reproduces the
% invzt_vertex4 stage-'Gamma' connected cumulant over the a,b (aa,ab,ac,bc,cc) external
% pairs and ALL internal (rho,sigma) channels on a small dense (n,l) grid to ~1e-11 --
% the equivalence invzt_sigma_tensor's precompute-once-contract-cheap loop rests on.
T = 1.6;  [~, ~, beta] = invz_matsubara(T, 40);
% small non-degenerate 3-level system with generic real Hermitian, centred operators
E = [0; 0.41; 0.95];  p = exp(-beta*(E - min(E)));  p = p/sum(p);
mkH = @(M) (M + M')/2;
ctr = @(O) O - real(sum(p(:).*diag(O)))*eye(3);
ops = struct('a', ctr(mkH([0 0.7 0.2; 0.7 0 0.5; 0.2 0.5 0])), ...
             'b', ctr(mkH([0.1 0.3 0.9; 0.3 -0.2 0.4; 0.9 0.4 0.05])), ...
             'c', ctr(mkH([0.8 0.1 0.15; 0.1 -0.6 0.25; 0.15 0.25 0.3])));
es = struct('E', E, 'p', p);
comps = {'a', 'b', 'c'};
ext = {{'a','a'}, {'a','b'}, {'a','c'}, {'b','c'}, {'c','c'}};
next = (0:1).';                       % n = 0, 1
Lmax = 1;  lvals = (-Lmax:Lmax).';    % l = -1, 0, 1
G4 = invzt_gamma4(es, ops, ext, comps, next, lvals, beta);
nlgrid = [reshape(repmat(next.', numel(lvals), 1), [], 1), repmat(lvals, numel(next), 1)];
worst = 0;
for ip = 1:numel(ext)
    for ri = 1:3
        for si = 1:3
            vo = invzt_vertex4(es, ops, [], [], beta, struct('stage', 'Gamma', ...
                'quad', {{ext{ip}{1}, ext{ip}{2}, comps{ri}, comps{si}}}, 'nl', nlgrid));
            % nlgrid rows are (n-major): [n0 l-1; n0 l0; n0 l1; n1 l-1; ...]; G4 is [.,.,.,nwn,nl]
            g4v = zeros(size(nlgrid, 1), 1);
            for k = 1:size(nlgrid, 1)
                in = find(next  == nlgrid(k, 1), 1);
                il = find(lvals == nlgrid(k, 2), 1);
                g4v(k) = G4(ip, ri, si, in, il);
            end
            worst = max(worst, max(abs(g4v - vo.val)));
        end
    end
end
fprintf('[gamma4==vertex4] worst |Gamma4 - vertex4 Gamma| over aa/ab/ac/bc/cc x {a,b,c}^2 = %.3e\n', worst);
verifyLessThan(testCase, worst, 1e-11);
end

% ============================== helpers ============================== %
function H = local_herm(N)
% Generic N x N Hermitian operator (state-space dimension N, for the dim-6 guard).
A = randn(N) + 1i*randn(N);
H = (A + A')/2;
end
