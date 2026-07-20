function tests = test_invzt_vertex
%TEST_INVZT_VERTEX  CORE gate for the A3 tensor vertex engine (Task 11).
% Validates invzt_kernels / invzt_vertex4 / invzt_vertex3 against the committed
% mpmath oracle invz_tensor/tests/fixtures/vertex_oracle.json (Task 10, factored_ok
% = false -> dense-only).  Runs with invz_projected OFF the path (CORE isolation).
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion engine
testCase.TestData.fx = jsondecode(fileread(fullfile(here, 'fixtures', 'vertex_oracle.json')));
end

% ===================================================================== %
% 1. KERNEL ROWS  (I2/I3 incl. repeated-node Hermite + large-beta*Delta rows)
% ===================================================================== %
function test_kernel_rows_vs_oracle(testCase)
% Whole-value per constraint 4/5: exprel-stable phi + Hermite limits + large-beta.
k  = invzt_kernels();
kr = testCase.TestData.fx.kernel_rows;
if ~iscell(kr), kr = num2cell(kr); end
worst = 0;
for i = 1:numel(kr)
    r = kr{i};
    args = r.args_re(:) + 1i * r.args_im(:);
    beta = r.beta;
    switch r.kernel
        case 'I2', v = k.I2(args(1), args(2), beta);
        case 'I3', v = k.I3(args(1), args(2), args(3), beta);
        otherwise, error('unknown kernel %s', r.kernel);
    end
    ref = r.re + 1i * r.im;
    e = relerr(v, ref);
    worst = max(worst, e);
    verifyLessThan(testCase, e, 1e-11, ...
        sprintf('%s(beta=%g) args=%s', r.kernel, beta, mat2str(r.args_re(:).', 3)));
end
fprintf('[kernels] %d rows, worst rel residual = %.3e\n', numel(kr), worst);
end

% ===================================================================== %
% 2. DENSE four-point F and connected Gamma4 vs oracle -- PRIMARY GATE
% ===================================================================== %
function test_dense_F_Gamma_vs_oracle(testCase)
fx = testCase.TestData.fx;
worstF = 0; worstG = 0; nF = 0; nG = 0;
for i = 1:numel(fx.vertex_rows)
    row = fx.vertex_rows(i);
    parts = strsplit(row.tags, ';');
    stage = parts{1};
    if ~ismember(stage, {'F', 'Gamma'}), continue; end
    munu = parts{2}; rs = parts{3}; system = parts{4};
    d  = sysdata(fx, system);
    es = struct('E', d.E, 'p', d.p);
    opts = struct('stage', stage);
    opts.quad = {munu(1), munu(2), rs(1), rs(2)};
    opts.nl   = [row.n, row.l];
    out = invzt_vertex4(es, d.ops, [], [], d.beta, opts);
    v   = out.val(1);
    ref = row.value_re + 1i * row.value_im;
    e   = relerr(v, ref);
    verifyLessThan(testCase, e, 1e-9, sprintf('%s system=%s (n,l)=(%d,%d)', ...
        row.tags, system, row.n, row.l));
    if strcmp(stage, 'F'), worstF = max(worstF, e); nF = nF + 1;
    else,                  worstG = max(worstG, e); nG = nG + 1; end
end
fprintf('[dense F]     %d rows, worst rel residual = %.3e\n', nF, worstF);
fprintf('[dense Gamma] %d rows, worst rel residual = %.3e\n', nG, worstG);
end

% ===================================================================== %
% 3. DENSE contracted V vs oracle (two-level Jensen system rows)
% ===================================================================== %
function test_dense_V_jensen_vs_oracle(testCase)
fx = testCase.TestData.fx;
S  = fx.systems.jensen_2lvl;
beta = S.beta; Ecc = S.E(:);
pcc  = boltz(Ecc, beta);
X    = centered(opmat(S.ops.c), pcc);
es   = struct('E', Ecc, 'p', pcc);
ops  = struct('c', X);
W    = 1.3;
Kf   = @(ri, si, l) kjen(ri, si, l, beta, W);
opts = struct('stage', 'V', 'Lmax', 400);
opts.comps = {'c'};
opts.ext   = {{'c', 'c'}};
worst = 0; nrow = 0;
for i = 1:numel(fx.vertex_rows)
    row = fx.vertex_rows(i);
    % NOTE: the exact-match prefix (trailing ';') is required -- without it this
    % filter would ALSO swallow the 'V;cc;jensen2lvlord_b14/b3;...' ordered rows
    % (round-2 P1-6: an over-broad filter must not silently mix systems).
    if ~startsWith(row.tags, 'V;cc;jensen2lvl;'), continue; end
    out = invzt_vertex4(es, ops, Kf, row.n, beta, opts);
    v   = out.val(1, 1);
    ref = row.value_re + 1i * row.value_im;
    e   = relerr(v, ref);
    verifyLessThan(testCase, e, 1e-9, sprintf('V;cc;jensen n=%d', row.n));
    worst = max(worst, e); nrow = nrow + 1;
end
verifyGreaterThan(testCase, nrow, 0, 'no V;cc;jensen2lvl rows matched (filter regression)');
fprintf('[dense V jensen] %d rows, worst rel residual = %.3e\n', nrow, worst);
end

% ===================================================================== %
% 3b. DENSE contracted V vs oracle (ordered two-level Jensen system rows)
% ===================================================================== %
function test_dense_V_jensen_ordered_vs_oracle(testCase)
% Ordered two-level oracle rows (split doublet, diagonal +-m, off-diagonal M; a3d
% prereq, J 2.26-2.29). TWO scalar-beta systems (_b14, _b3; round-2 P1-6 -- sysdata
% expects scalar S.beta). Same dense four-point contraction path as the PM Jensen
% test above; this is a pure dense-vs-dense oracle gate (NOT the moment-form
% approximation diagnostic -- that is test_jensen_ordered_static_diagnostic below).
fx = testCase.TestData.fx;
W = 1.3;
for tagc = {'b14', 'b3'}
    tag = tagc{1};
    S    = fx.systems.(['jensen_2lvl_ordered_' tag]);
    beta = S.beta; Ecc = S.E(:);
    pcc  = boltz(Ecc, beta);
    X    = centered([S.m sqrt(S.M2); sqrt(S.M2) -S.m], pcc);
    es   = struct('E', Ecc, 'p', pcc);
    ops  = struct('c', X);
    Kf   = @(ri, si, l) kjen(ri, si, l, beta, W);
    opts = struct('stage', 'V', 'Lmax', 400);
    opts.comps = {'c'};
    opts.ext   = {{'c', 'c'}};
    prefix = ['V;cc;jensen2lvlord_' tag ';'];
    worst = 0; nrow = 0;
    for i = 1:numel(fx.vertex_rows)
        row = fx.vertex_rows(i);
        if ~startsWith(row.tags, prefix), continue; end
        out = invzt_vertex4(es, ops, Kf, row.n, beta, opts);
        v   = out.val(1, 1);
        ref = row.value_re + 1i * row.value_im;
        e   = relerr(v, ref);
        verifyLessThan(testCase, e, 1e-9, sprintf('%s n=%d', row.tags, row.n));
        worst = max(worst, e); nrow = nrow + 1;
    end
    verifyGreaterThan(testCase, nrow, 0, sprintf('no ordered V rows matched tag %s (filter regression)', tag));
    fprintf('[dense V jensen ordered %s] %d rows, worst rel residual = %.3e\n', tag, nrow, worst);
end
end

% ===================================================================== %
% 4. KMS large-beta*Delta whole contracted value (constraint 4; betaDelta<=200)
% ===================================================================== %
function test_kms_large_betaDelta(testCase)
fx = testCase.TestData.fx;
Delta = 1.0; M = 0.9;
opts = struct('stage', 'V', 'Lmax', 30);
opts.comps = {'c'};
opts.ext   = {{'c', 'c'}};
worst = 0;
for i = 1:numel(fx.kms_rows)
    row  = fx.kms_rows(i);
    bD   = row.beta_Delta; n = row.n;
    beta = bD / Delta;
    Ecc  = [0; Delta]; pcc = boltz(Ecc, beta);
    X    = centered([0 M; M 0], pcc);
    es   = struct('E', Ecc, 'p', pcc);
    ops  = struct('c', X);
    Kf   = @(ri, si, l) kkms(ri, si, l, beta);
    out  = invzt_vertex4(es, ops, Kf, n, beta, opts);
    v    = out.val(1, 1);
    ref  = row.V_re + 1i * row.V_im;
    e    = relerr(v, ref);
    verifyLessThan(testCase, e, 1e-11, sprintf('KMS V_cc(%d) betaDelta=%g', n, bD));
    worst = max(worst, e);
    verifyTrue(testCase, all(isfinite([real(v) imag(v)])));
end
fprintf('[KMS betaDelta<=200] %d rows, worst rel residual = %.3e\n', numel(fx.kms_rows), worst);
end

% ===================================================================== %
% 5. Two-level Jensen bridge:  V == G0 .* Sigma  (invz_twolevel/invz_sigma)
% ===================================================================== %
function test_jensen_bridge(testCase)
% invz_common code UNCHANGED: build tl via invz_twolevel, self-energy via
% invz_lambdas/invz_sigma, and check the dense four-point contraction reproduces
% G0.*Sigma to RelTol 1e-6.
ion = invz_ion(); T = 1.6; Bx = 0.5;
tl  = invz_twolevel(ion, T, Bx, struct());
C   = invz_const(); beta = 1 / (C.kB * T);
Delta = tl.Delta; M2 = tl.M2;
W = 1.3;
% internal Matsubara sum for the lambdas (k = 0..Lmax, doubling weights):
Lmax = 400; kk = (0:Lmax).'; wts = [1; 2 * ones(Lmax, 1)];
wln  = 2 * pi * kk / beta;
Kint = 0.37 ./ (1 + (wln / W).^2);
gint = invz_g(tl, 1i * wln);
lam  = invz_lambdas(Kint, gint, wts, beta, [1 2]);
% external frequencies n = 0,1,2:
next = [0; 1; 2]; wext = 2 * pi * next / beta;
Kext = 0.37 ./ (1 + (wext / W).^2);
gext = invz_g(tl, 1i * wext);
sig  = invz_sigma(tl, lam, Kext, gext, beta);
G0   = -M2 * gext;
Vjen = G0 .* sig.Sigma;
% independent: exact four-point contraction of the cc channel
Ecc = [0; Delta]; pcc = boltz(Ecc, beta);
X   = centered([0 sqrt(M2); sqrt(M2) 0], pcc);
es  = struct('E', Ecc, 'p', pcc); ops = struct('c', X);
Kf  = @(ri, si, l) kjen(ri, si, l, beta, W);
opts = struct('stage', 'V', 'Lmax', 400);
opts.comps = {'c'}; opts.ext = {{'c', 'c'}};
out = invzt_vertex4(es, ops, Kf, next.', beta, opts);
Vme = out.val(1, :).';
verifyEqual(testCase, Vme, Vjen, 'RelTol', 1e-6);
fprintf('[Jensen bridge] Delta=%.4f meV, max rel diff V vs G0*Sigma = %.3e\n', ...
    Delta, max(abs(Vme - Vjen) ./ max(1, abs(Vjen))));
end

% ===================================================================== %
% 5b. Ordered two-level Jensen DIAGNOSTIC: STATIC gate (n=0) + finite-freq REPORT
% ===================================================================== %
function test_jensen_ordered_static_diagnostic(testCase)
% Ordered moment-form comparison (J 2.26-2.29) -- STATIC GATE + FINITE-FREQUENCY
% REPORT (round-2 P0-2: the moment form is n=0-exact/finite-n-approximate BY
% DERIVATION -- g(wm+-wn) -> g(wm) -- so equality at n>0 must NOT be gated).
% This is an approximation diagnostic for a3d, NOT a correctness gate on the
% dense vertex (that gate is the oracle-row test above).
Delta = 1.0;  M2 = 0.81;  m = 0.6;  W = 1.3;
for ib = 1:2
    betas = [14, 3];  beta = betas(ib);
    n01 = tanh(beta*Delta/2);
    tl = struct('Delta', Delta, 'M2', M2, 'm', m, 'n01', n01, 'g0', 2*n01/Delta);
    Lmax = 400; kk = (0:Lmax).'; wts = [1; 2*ones(Lmax, 1)];
    wln = 2*pi*kk/beta;  Kint = 0.37 ./ (1 + (wln/W).^2);
    gint = invz_g(tl, 1i*wln);
    lam = invz_lambdas(Kint, gint, wts, beta, [1 2 3]);
    next = [0; 1; 2];  wext = 2*pi*next/beta;
    Kext = 0.37 ./ (1 + (wext/W).^2);
    gext = invz_g(tl, 1i*wext);
    sig = invz_sigma_ordered(tl, lam, Kext, gext, beta);
    h0  = beta*(1 - n01^2);                      % elastic term (J 2.8), n = 0 slot only
    G0  = -M2*gext;  G0(1) = G0(1) - m^2*h0;     % ordered bare propagator (J 2.7)
    Vjen = G0 .* sig.Sigma;
    % independent: exact four-point contraction, m != 0 (centered operator)
    Ecc = [0; Delta];  pcc = boltz(Ecc, beta);
    X   = centered([m sqrt(M2); sqrt(M2) -m], pcc);
    es  = struct('E', Ecc, 'p', pcc);  ops = struct('c', X);
    Kf  = @(ri, si, l) kjen(ri, si, l, beta, W);
    opts = struct('stage', 'V', 'Lmax', 400);
    opts.comps = {'c'};  opts.ext = {{'c', 'c'}};
    out = invzt_vertex4(es, ops, Kf, next.', beta, opts);
    Vme = out.val(1, :).';
    res = abs(Vme - Vjen) ./ max(abs(Vjen), 1e-300);
    if ib == 1
        % n = 0 STATIC GATE, low-elastic-weight regime. Tolerance derived from
        % the omitted elastic/static-resummation scale (30x margin over the
        % omitted-term ratio; floor 1e-4). Measured naive residual: 3.13e-5.
        tol0 = max(1e-4, 30 * m^2*h0 / (M2*abs(gext(1))));
        verifyLessThan(testCase, res(1), tol0);
        % derivation-shape check: the n=0 identity is MUCH better than the
        % n>0 static-value approximation (their measured ratio is ~1e3).
        verifyLessThan(testCase, res(1), 0.1*res(2));
        fprintf('[ordered static gate bD=14] res(n=0)=%.3e (tol %.3e); REPORT n=1: %.3e, n=2: %.3e\n', ...
            res(1), tol0, res(2), res(3));
    else
        fprintf(['[ordered diagnostic bD=3, REPORT ONLY] res = [%.3e %.3e %.3e] -- ' ...
            'elastic-sector resummation ambiguity (constraint 7), a3d input\n'], res);
    end
end
end

% ===================================================================== %
% 6. Degenerate-doublet regularity: Jc = diag(m,-m,0), exact E1=E2 -> no NaN
% ===================================================================== %
function test_degenerate_doublet_no_nan(testCase)
% Bx=0 doublet: E1 = E2 exactly.  Jc = diag(m,-m,0) is ALREADY centered at equal
% populations (no bogus centering line).  The path sum uses initial-state weights
% p_r (not D_pq/(eps_q-eps_p) ratios), so exact degeneracy stays finite.
d2 = 0.9; m = 0.85; beta = 1.7;
E  = [0; 0; d2]; p = boltz(E, beta);
Jc = diag([m, -m, 0]);
verifyEqual(testCase, real(sum(p .* diag(Jc))), 0, 'AbsTol', 1e-14);  % already centered
es = struct('E', E, 'p', p); ops = struct('c', Jc);
% un-contracted connected Gamma4 at repeated-node (n,l):
optsG = struct('stage', 'Gamma');
optsG.quad = {'c', 'c', 'c', 'c'};
optsG.nl   = [0 0; 0 1; 1 -1; 2 1; 1 1];
outG = invzt_vertex4(es, ops, [], [], beta, optsG);
verifyTrue(testCase, all(isfinite(outG.val)), 'Gamma NaN at exact degeneracy');
% contracted V at the doublet:
Kf = @(ri, si, l) kkms(ri, si, l, beta);
optsV = struct('stage', 'V', 'Lmax', 20);
optsV.comps = {'c'}; optsV.ext = {{'c', 'c'}};
outV = invzt_vertex4(es, ops, Kf, [0 1 2], beta, optsV);
verifyTrue(testCase, all(isfinite(outV.val(:))), 'V NaN at exact degeneracy');
fprintf('[degeneracy] Gamma4 & V finite at exact E1=E2 doublet (max|V|=%.3e)\n', max(abs(outV.val(:))));
end

% ===================================================================== %
% 7. Negative-frequency reconstruction (constraint 9): C rows + transpose
% ===================================================================== %
function test_negfreq_C_and_transpose(testCase)
fx = testCase.TestData.fx;
d  = sysdata(fx, 'fieldon_complex_3lvl');
es = struct('E', d.E, 'p', d.p); beta = d.beta;
lab = struct('a', 'a', 'b', 'b', 'c', 'c');
worst = 0; nrow = 0;
for i = 1:numel(fx.negfreq_rows)
    row = fx.negfreq_rows(i);
    parts = strsplit(row.tags, ';');
    if ~strcmp(parts{1}, 'C'), continue; end
    pair = parts{2};
    opts = struct('stage', 'C');
    opts.quad = {pair(1), pair(2)};
    opts.nl   = row.n;
    out = invzt_vertex4(es, d.ops, [], [], beta, opts);
    v   = out.val(1);
    ref = row.value_re + 1i * row.value_im;
    e   = relerr(v, ref);
    verifyLessThan(testCase, e, 1e-9, sprintf('%s n=%d', row.tags, row.n));
    worst = max(worst, e); nrow = nrow + 1;
end
% LOCKED transpose relation C_{rho sigma}(-l) = C_{sigma rho}(+l) (complex, off-diag)
labs = {'a', 'b', 'c'};
tmax = 0;
for ii = 1:3
    for jj = 1:3
        for l = 1:3
            on  = struct('stage', 'C'); on.quad = {labs{ii}, labs{jj}}; on.nl = -l;
            op  = struct('stage', 'C'); op.quad = {labs{jj}, labs{ii}}; op.nl =  l;
            cn  = invzt_vertex4(es, d.ops, [], [], beta, on);
            cp  = invzt_vertex4(es, d.ops, [], [], beta, op);
            tmax = max(tmax, abs(cn.val(1) - cp.val(1)));
        end
    end
end
verifyLessThan(testCase, tmax, 1e-11, 'C transpose relation (constraint 9)');
fprintf('[negfreq C] %d rows, worst rel residual = %.3e; transpose |C(-l)-C^T(l)| max = %.3e\n', ...
    nrow, worst, tmax);
end

% ===================================================================== %
% 7b. Negative-n contracted V (array Kmat -> transpose reconstruction of -l)
% ===================================================================== %
function test_negfreq_V_reconstruction(testCase)
fx = testCase.TestData.fx;
d  = sysdata(fx, 'fieldon_complex_3lvl');
es = struct('E', d.E, 'p', d.p); beta = d.beta;
Lmax = 80;
labs = {'a', 'b', 'c'};
idx  = struct('a', 1, 'b', 2, 'c', 3);
% build the internal kernel as an ARRAY K[ri,si,l], l = 0..Lmax; the engine
% reconstructs negative l via the transpose relation (constraint 9).
Karr = zeros(3, 3, Lmax + 1);
for l = 0:Lmax
    wl = 2 * pi * l / beta;
    Karr(:, :, l + 1) = eye(3) * (0.3 / (1 + (wl / 1.5)^2));
end
opts = struct('stage', 'V', 'Lmax', Lmax);
opts.comps = {'a', 'b', 'c'};
worst = 0; nrow = 0;
for i = 1:numel(fx.negfreq_rows)
    row = fx.negfreq_rows(i);
    parts = strsplit(row.tags, ';');
    if ~strcmp(parts{1}, 'V'), continue; end
    pair = parts{2};
    opts.ext = {{pair(1), pair(2)}};
    out = invzt_vertex4(es, d.ops, Karr, row.n, beta, opts);
    v   = out.val(1, 1);
    ref = row.value_re + 1i * row.value_im;
    e   = relerr(v, ref);
    verifyLessThan(testCase, e, 1e-8, sprintf('%s n=%d', row.tags, row.n));
    worst = max(worst, e); nrow = nrow + 1;
    % transpose relation on the contracted vertex: V_{mu nu}(-1) = V_{nu mu}(1)
    opts.ext = {{pair(2), pair(1)}};
    outp = invzt_vertex4(es, d.ops, Karr, -row.n, beta, opts);
    verifyLessThan(testCase, abs(v - outp.val(1, 1)), 1e-8, ...
        sprintf('V transpose %s', row.tags));
end
fprintf('[negfreq V] %d rows, worst rel residual = %.3e (array-K transpose reconstruction)\n', ...
    nrow, worst);
end

% ===================================================================== %
% 8. factored path disabled (factored_ok = false -> invzt:factoredUnproven)
% ===================================================================== %
function test_factored_impl_errors(testCase)
fx = testCase.TestData.fx;
verifyFalse(testCase, logical(fx.factored.factored_ok), ...
    'fixture factored_ok must be false for this branch of the plan');
d  = sysdata(fx, 'random_real_3lvl');
es = struct('E', d.E, 'p', d.p);
opts = struct('impl', 'factored', 'stage', 'F');
opts.quad = {'a', 'a', 'a', 'a'}; opts.nl = [0 0];
verifyError(testCase, @() invzt_vertex4(es, d.ops, [], [], d.beta, opts), ...
    'invzt:factoredUnproven');
fprintf('[factored] disabled -> invzt:factoredUnproven (factored_ok=%d)\n', fx.factored.factored_ok);
end

% ===================================================================== %
% 9. three-point vertex vs independent Gauss-Legendre simplex quadrature
% ===================================================================== %
function test_vertex3_vs_quadrature(testCase)
fx = testCase.TestData.fx;
d  = sysdata(fx, 'random_real_3lvl');
es = struct('E', d.E, 'p', d.p); beta = d.beta;
quad = {'a', 'b', 'c'};                 % external a, internal (b,c)
lvals = [0; 1; 2];
opts = struct(); opts.quad = quad;
[dm, T3] = invzt_vertex3(es, d.ops, lvals, beta, opts);
verifyTrue(testCase, isfinite(dm));
worst = 0;
for il = 1:numel(lvals)
    ref = ta_brute(es, d.ops, quad, lvals(il), beta, 90);
    e   = relerr(T3(il), ref);
    verifyLessThan(testCase, e, 1e-7, sprintf('T3(l=%d) vs quadrature', lvals(il)));
    worst = max(worst, e);
end
fprintf('[vertex3] worst T3-vs-quadrature rel residual = %.3e ; dm=%.4f\n', worst, real(dm));
end

function test_vertex3_degenerate_no_nan(testCase)
E = [0; 0; 0.9]; beta = 1.7; p = boltz(E, beta);
ops = struct('c', diag([0.85, -0.85, 0]));
[dm, T3] = invzt_vertex3(struct('E', E, 'p', p), ops, [0; 1; 2], beta, ...
    struct('quad', {{'c', 'c', 'c'}}));
verifyTrue(testCase, isfinite(dm) && all(isfinite(T3)), 'vertex3 NaN at degeneracy');
end

% ============================== helpers ============================== %
function d = sysdata(fx, name)
% Rebuild (E, p, ops) for a fixture system, reconstructing the operators the oracle
% did NOT serialise (populated_spectator's cc channel; the v2 scalar reference).
if strcmp(name, 'v2_scalar_ref')
    Xm = [0.3 0.8 0.45; 0.8 -0.5 0.6; 0.45 0.6 0.2];
    E = [0; 0.6; 2.9]; beta = 1.7; p = boltz(E, beta);
    d = struct('E', E, 'p', p, 'beta', beta);
    d.ops = struct('c', centered(Xm, p));
    return
end
S = fx.systems.(name);
beta = S.beta; E = S.E(:); p = boltz(E, beta);
d = struct('E', E, 'p', p, 'beta', beta);
d.ops = struct();
if isfield(S, 'ops')
    fns = fieldnames(S.ops);
    for i = 1:numel(fns), d.ops.(fns{i}) = opmat(S.ops.(fns{i})); end
end
if strcmp(name, 'populated_spectator_3lvl')
    Jp = zeros(3); Jp(1, 2) = 1; Jp(2, 1) = 1;   % pure 0<->1 channel, M=1
    d.ops.c = centered(Jp, p);
end
end

function O = opmat(s)
O = s.re + 1i * s.im;
end

function p = boltz(E, beta)
E = E(:); w = exp(-beta * (E - min(E))); p = w / sum(w);
end

function O = centered(O, p)
av = real(sum(p(:) .* diag(O)));
O = O - av * eye(size(O, 1));
end

function e = relerr(a, b)
e = abs(a - b) / max([1, abs(a), abs(b)]);
end

function K = kjen(ri, si, l, beta, W)
% Jensen cc-channel Lorentzian kernel (single internal channel).
if ri == 1 && si == 1
    wl = 2 * pi * l / beta;
    K = 0.37 / (1 + (wl / W)^2);
else
    K = 0;
end
end

function K = kkms(ri, si, l, beta)
% KMS-test cc-channel Lorentzian kernel.
if ri == 1 && si == 1
    wl = 2 * pi * l / beta;
    K = 0.3 / (1 + (wl / 2.0)^2);
else
    K = 0;
end
end

function v = ta_brute(es, ops, quad, l, beta, ng)
% Independent three-point value by Gauss-Legendre quadrature over the ordered
% simplex t1 > t2 > 0 (both internal-leg orderings), mirroring invzt_vertex3.
E = es.E(:); p = es.p(:); N = numel(E);
Om = ops.(quad{1}); Or = ops.(quad{2}); Os = ops.(quad{3});
zl = 1i * 2 * pi * l / beta;
[gx, gw] = gauleg(ng);
a = 0.5 * (gx + 1); wa = 0.5 * gw;
[A1, A2] = ndgrid(a, a); [Wa1, Wa2] = ndgrid(wa, wa);
T1 = beta * A1; T2 = T1 .* A2;
Wj = Wa1 .* Wa2 * beta^2 .* A1;                 % simplex jacobian
val = complex(zeros(size(T1)));
for r = 1:N
    for s = 1:N
        Ors = Or(r, s);
        if Ors == 0, continue; end
        for u = 1:N
            w = p(r) * Ors * Os(s, u) * Om(u, r);
            if w == 0, continue; end
            val = val + w * exp(E(r) * T1 - E(s) * (T1 - T2) - E(u) * T2);
        end
    end
end
v = sum(sum(Wj .* val .* exp(zl * (T1 - T2)))) ...
  + sum(sum(Wj .* val .* exp(zl * (T2 - T1))));
end

function [x, w] = gauleg(n)
% Gauss-Legendre nodes/weights on [-1,1] via the Golub-Welsch eigenproblem.
i = 1:n-1; b = i ./ sqrt(4 * i.^2 - 1);
J = diag(b, 1) + diag(b, -1);
[V, D] = eig(J);
x = diag(D); [x, idx] = sort(x);
w = 2 * (V(1, idx).^2); w = w(:); x = x(:);
end
