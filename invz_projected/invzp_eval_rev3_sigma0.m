function invzp_eval_rev3_sigma0()
%INVZP_EVAL_REV3_SIGMA0 Self-contained test of rev3's claims on regenerated node data.
% Regenerates the mix_outer=0.3 predictor+33-node HMF sweeps at Bx = 1 and 3 T
% (exact replica of invz_solve_point_ordered's node kernel and warm-start chain;
% deterministic, reproduces the 2026-07-29 states), then per node evaluates:
%   pred    = K0*g0*(M2 - 3*m^2)   -- rev3's closed form for Sigma0 (n01->1, const K)
%   3m2/M2  -- rev3's sign-crossing variable (Sigma0 < 0 iff > 1)
%   J0|G0b| -- rev3's h-dagger admissibility ratio J0eff*|G0bare(0;h)| (single-ion only)
%   domshare= M2*g0/|G0bare|       -- dominant-sector static share (whole-cc remainder)
%   nneg    -- # mesh points with 1+(Jq-K0)Gstat < 0 at the final state (class-A marker)
% Output: console table + diag_rev3_check.mat. Throwaway measurement script.
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(here, '..')); addpath(fullfile(here, '..', 'invz_common'));
ion = invz_ion();
co = struct('grid', [16 16 16], 'cache', false, 'dipole', 'ewald', ...
    'gridConvention', 'halfopen', 'gridOffset', [0 0 0], ...
    'gammaPolicy', 'P_drop', ...
    'ewald', struct('alpha', 0.3, 'r_cut', 16, 'g_cut', 3, ...
                    'boundary', 'conducting_k0_omitted'));
[J, info] = invz_bz_couplings(ion, co);
Jf = J(:);  J0eff = info.Jcc0;  Jxx0 = info.Jaa0;
T = 0.1;  hyp = true;
[wn, wts, beta] = invz_matsubara(T, 40);
ctx = struct('wn', wn, 'wts', wts, 'beta', beta, 'hyp', hyp, 'Jxx0', Jxx0, ...
    'J0eff', J0eff, 'mixo', 0.3, 'tolo', 1e-8, 'maxo', 500);
out = struct('Bx', {}, 'nodes', {});
for Bx = [1.0 3.0]
    Bvec = [Bx 0 0];
    sib = invz_single_ion(ion, T, Bvec, struct('hyp', hyp, 'Jxx0', Jxx0, ...
        'order', true, 'J0z', J0eff));
    hmax = 1.25*abs(sib.hz);
    nH = 33;  ratio = (1e-3)^(1/(nH-1));
    hgrid = [0, hmax * ratio.^((nH-1):-1:0)];
    fprintf('\n=== Bx=%g T | mix03 (regenerated) ===\n', Bx);
    fprintf('idx  h            conv Sigma0     K0         pred       3m2/M2   J0|G0b|   domshare  nneg\n');
    Sigma = [];  K0s = 0;  nodes = {};
    for k = 1:numel(hgrid)
        [nd, Sigma, K0s] = eval_node_lean(ion, T, Bx, hgrid(k), Jf, ctx, Sigma, K0s);
        tl = nd.tl;
        nd.pred = nd.K0 * tl.g0 * (tl.M2 - 3*tl.m^2);
        nd.r3   = 3*tl.m^2 / tl.M2;
        nd.jr   = J0eff * abs(nd.G0bare0);
        nd.ds   = tl.M2*tl.g0 / abs(nd.G0bare0);
        nodes{end+1} = nd; %#ok<AGROW>
        fprintf('%3d  %-12.5g %d    %-10.4g %-10.5g %-10.4g %-8.4g %-9.4g %-9.3g %d\n', ...
            k-1, nd.h, nd.conv, nd.Sigma0, nd.K0, nd.pred, nd.r3, nd.jr, nd.ds, nd.nneg);
    end
    out(end+1) = struct('Bx', Bx, 'nodes', {nodes}); %#ok<AGROW>
end
save(fullfile(fileparts(here), 'diag_rev3_check.mat'), 'out', 'info', 'T');
fprintf('\nSaved diag_rev3_check.mat\n');
end

function [nd, Sigma, K0s] = eval_node_lean(ion, T, Bx, hp, Jf, ctx, Sigma, K0s)
% Exact replica of production eval_node (invz_solve_point_ordered), no traces.
Bvec = [Bx 0 0];
si = invz_single_ion(ion, T, Bvec, struct('hyp', ctx.hyp, 'Jxx0', ctx.Jxx0, 'hz_fixed', hp));
tl = invz_twolevel_ordered(ion, T, Bx, hp, struct('Jxx0', ctx.Jxx0));
c0 = invz_chi0z(si, T, 1i*ctx.wn, struct('elastic', true));
G0 = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si, T, 1i*ctx.wn(1), struct('elastic', false));
G0i = -real(c0i(3,3,1));
X = real(c0(:, :, 1));
fb = X(3,1) * (ctx.Jxx0 / (1 - ctx.Jxx0*X(1,1))) * X(1,3);
G0b = -(X(3,3) + fb);
G0e = G0b - G0i;
g = real(invz_g(tl, 1i*ctx.wn));
if isempty(Sigma), Sigma = zeros(size(ctx.wn)); end
K = zeros(size(ctx.wn));  lam = [0; 0; 0];  ok = false;
for outer = 1:ctx.maxo
    med = invz_emt_scalar(G0, Sigma, Jf, struct());
    K = med.K;
    [K0s, ~, sout] = emt_static_replica(tl, lam(1:2), Sigma(1), Jf, K0s, ctx.beta, ctx.J0eff, G0i, G0e);
    K(1) = K0s;
    lam = invz_lambdas(K, g, ctx.wts, ctx.beta, [1 2 3]);
    sg = invz_sigma_ordered(tl, lam, K, g, ctx.beta);
    dS = max(abs(sg.Sigma - Sigma));
    Sigma = Sigma + ctx.mixo*(sg.Sigma - Sigma);
    if dS < ctx.tolo && sout.converged, ok = true; break; end
end
[K0s, Gsk, so] = emt_static_replica(tl, lam(1:2), Sigma(1), Jf, K0s, ctx.beta, ctx.J0eff, G0i, G0e);
ok = ok && so.converged && isfinite(so.resid) && so.resid < 1e-10;
den = 1 + (Jf - K0s).*Gsk;
nd = struct('h', hp, 'conv', ok, 'iters', outer, 'K0', K0s, 'Sigma0', Sigma(1), ...
    'Gstat', Gsk, 'G0inel0', G0i, 'G0el0', G0e, 'G0bare0', G0b, 'm_si', si.Jexp(3), ...
    'tl', tl, 'lam', lam, 'nneg', sum(den < 0), 'den_min', min(abs(den)));
end

function [K0, Gstat, out] = emt_static_replica(tl, lam2, Sigma0, Jf, K0_seed, beta, J0eff, G0i, G0e)
% Verbatim semantics of the production local invz_emt_static_ordered.
rtol = 1e-10;  tol = 0;  maxit = 200;  mix = 0.5;
K0 = K0_seed;
for it = 1:maxit
    Gs = invz_gstat_ordered(tl, lam2, K0, Sigma0, beta, G0i, G0e);
    Gq = Gs ./ (1 + (Jf - K0).*Gs);
    Gbar = mean(Gq);
    if abs(Gbar - Gs) < rtol, break; end
    K0_new = mean(Jf .* Gq) / Gbar;
    dK = abs(K0_new - K0);
    if dK < max(tol, 4*eps(abs(K0))), break; end
    K0 = K0 + mix*(K0_new - K0);
end
[Gstat, go] = invz_gstat_ordered(tl, lam2, K0, Sigma0, beta, G0i, G0e);
Gq = Gstat ./ (1 + (Jf - K0).*Gstat);
out = go;
out.resid = abs(mean(Gq) - Gstat);
out.iters = it;
out.converged = out.resid < rtol;
end
