function P = invzp_exec_s10_outer_eigenvalue(Bx, alphas, save_prefix)
%INVZP_EXEC_S10_OUTER_EIGENVALUE Measure the outer-map Jacobian eigenvalue from damping.
%
% Execution packet S10 of docs/execution/invzp_plan_execution_diary.md.
% MEASUREMENT ONLY. Decides, for one field, whether the outer Sigma<->K map's
% failures are curable by damping AT ALL, rather than searching for a lucky rung.
%
% THE IDENTITY BEING USED. Damped Picard on the outer map is
%     u_{k+1} = (1 - a) u_k + a G(u_k),
% so with e_k = u_k - u* and J = dG/du at the fixed point,
%     e_{k+1} = [(1-a) I + a J] e_k =: M e_k.
% The recorded outer residual is the RAW map step r_k = G(u_k) - u_k ~ (J - I) e_k.
% J - I and M are both polynomials in J, hence commute, so r obeys the same
% recursion: r_{k+1} = M r_k. The tail ratio of |r| therefore measures
%     rho(a) = |mu|,  mu = 1 - a + a*lambda,
% with lambda the dominant eigenvalue of J. For REAL lambda this is a straight
% line in a:
%
%   lambda > 1  ->  rho(a) = 1 + a(lambda - 1):  intercept 1, slope POSITIVE,
%                   rho > 1 for EVERY a in (0,1]. The fixed point is linearly
%                   unstable and NO damping value can converge onto it. Damping
%                   and iteration budget are the wrong tools; the fix is Newton
%                   (which is indifferent to the sign of lambda - 1) or
%                   continuation through the fold.
%   lambda < 1  ->  rho(a) = |1 - a(1 - lambda)|: slope NEGATIVE in a. The map is
%                   merely slow/non-contractive at large a and SOME a converges.
%
% So the sign of d rho / d a discriminates "damping can never work here" from
% "damping needs the right alpha", and lambda = 1 + slope reads off the fit.
%
% WHY THIS IS A MEASUREMENT AND NOT MACHINERY. It adds no solver, changes no
% production path, and reuses invzp_exec_s1_failure_census unchanged. Its only
% output is a number, lambda, and the qualitative verdict that follows from its
% position relative to 1.
%
% CAVEATS, STATED BECAUSE THEY BOUND THE CONCLUSION.
%  - rho is a LOCAL linearization diagnostic. If the iterate is wandering across
%    the Gstat pole rather than sitting near a fixed point, the linear model does
%    not apply; that case is detected here as a large tail-ratio IQR and a failed
%    straight-line fit, and is REPORTED as such rather than being fitted anyway.
%  - a complex-conjugate dominant pair gives an oscillating ratio; also flagged
%    by the IQR, and by the alternation count.
%  - lambda is the dominant eigenvalue only. A subdominant unstable direction
%    would not be seen.
%
%   Bx           field in T (default 1.0, the deep-ordered failure found by S9)
%   alphas       mix_outer values to probe (default [0.15 0.25 0.35 0.50])
%   save_prefix  prefix for the per-alpha trace .mat files (required: the traces
%                are large and are what the ratio is computed from)
if nargin < 1 || isempty(Bx), Bx = 1.0; end
if nargin < 2 || isempty(alphas), alphas = [0.15 0.25 0.35 0.50]; end
if nargin < 3 || isempty(save_prefix)
    save_prefix = fullfile('docs','execution','out', sprintf('s10_B%s', strrep(num2str(Bx),'.','p')));
end

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); addpath(fullfile(root, 'docs', 'execution'));

fprintf('=== S10 outer-map eigenvalue probe, Bx = %.4f T ===\n', Bx);
fprintf('alphas = [%s], max_outer = 1000\n\n', num2str(alphas));

runs = struct('alpha', {}, 'mat', {}, 'n_failed', {}, 'failed_ids', {});
for k = 1:numel(alphas)
    mp = sprintf('%s_a%s.mat', save_prefix, strrep(num2str(alphas(k)),'.','p'));
    so = struct('mix_outer', alphas(k), 'max_outer', 1000);
    txt = evalc('Rk = invzp_exec_s1_failure_census(Bx, mp, so);'); %#ok<NASGU>
    fid = find(~Rk.tab.accepted).';
    runs(end+1) = struct('alpha', alphas(k), 'mat', mp, ...
        'n_failed', Rk.meta.n_failed, 'failed_ids', fid); %#ok<AGROW>
    fprintf('alpha = %.2f -> %s, failed %d/%d, ids [%s]\n', alphas(k), ...
        Rk.meta.hmf_status, Rk.meta.n_failed, Rk.meta.n_nodes, num2str(fid));
end

% Nodes that fail at EVERY alpha -- the ones whose lambda is worth measuring.
core = runs(1).failed_ids;
for k = 2:numel(runs), core = intersect(core, runs(k).failed_ids); end
fprintf('\nnodes failing at every alpha: [%s]\n\n', num2str(core));
if isempty(core)
    fprintf('no common failing node; lambda is not measurable on a shared node.\n');
    P = struct('Bx', Bx, 'runs', runs, 'core', core, 'rows', [], 'verdict', 'no-common-node');
    return
end

rows = struct('node', {}, 'alpha', {}, 'rho_med', {}, 'rho_iqr', {}, ...
    'n_iter', {}, 'gdenom_min', {}, 'gdenom_max', {}, 'sign_flips', {});
for c = core(:).'
    for k = 1:numel(runs)
        S = load(runs(k).mat, 'trace');
        it = S.trace.iters;
        sel = it([it.node_id] == c);
        if numel(sel) < 50, continue; end
        r = abs([sel.resid_map].');
        gd = [sel.gstat_local_denom].';
        % Tail only: the linearization is about the LATE iterate, and the first
        % few hundred iterations include the transient from the seed.
        tail = max(2, numel(r)-299):numel(r);
        ratio = r(tail)./max(r(tail-1), realmin);
        q = quantile(ratio, [0.25 0.5 0.75]);
        dK = diff([sel.K0].');
        rows(end+1) = struct('node', c, 'alpha', runs(k).alpha, ...
            'rho_med', q(2), 'rho_iqr', q(3)-q(1), 'n_iter', numel(r), ...
            'gdenom_min', min(gd(tail)), 'gdenom_max', max(gd(tail)), ...
            'sign_flips', nnz(diff(sign(dK)) ~= 0)); %#ok<AGROW>
    end
end

fprintf('%-6s %-7s %-11s %-11s %-8s %-10s %-10s %s\n', ...
    'node', 'alpha', 'rho_median', 'rho_IQR', 'n_iter', 'gdenom_lo', 'gdenom_hi', 'signflips');
for i = 1:numel(rows)
    fprintf('%-6d %-7.2f %-11.6f %-11.3e %-8d %-10.3g %-10.3g %d\n', ...
        rows(i).node, rows(i).alpha, rows(i).rho_med, rows(i).rho_iqr, ...
        rows(i).n_iter, rows(i).gdenom_min, rows(i).gdenom_max, rows(i).sign_flips);
end

% ---- straight-line fit rho(a) = 1 + a(lambda - 1), per node -------------------
fprintf('\n--- fit rho(a) = c + a*s   (real lambda => c ~ 1, lambda = 1 + s) ---\n');
fits = struct('node', {}, 'slope', {}, 'intercept', {}, 'lambda', {}, ...
    'resid', {}, 'max_iqr', {}, 'verdict', {});
for c = unique([rows.node])
    m = [rows.node] == c;
    a = [rows(m).alpha].'; rh = [rows(m).rho_med].';
    if numel(a) < 2, continue; end
    A = [ones(numel(a),1), a];
    p = A\rh;
    res = norm(A*p - rh, inf);
    maxiqr = max([rows(m).rho_iqr]);
    lam = 1 + p(2);
    % A large tail-ratio IQR means the ratio is not settling to a single |mu|,
    % so the linear model does not describe this node and lambda is NOT reported
    % as meaningful -- the pole-crossing / complex-pair case.
    if maxiqr > 0.25
        vd = 'linearization-invalid (ratio not settling; suspect pole crossing or complex pair)';
    elseif abs(p(1) - 1) > 0.15
        vd = 'intercept far from 1 (linear model questionable)';
    elseif p(2) > 0.02
        vd = 'UNSTABLE: lambda > 1, NO damping value converges';
    elseif p(2) < -0.02
        vd = 'stable-but-slow: some alpha converges';
    else
        vd = 'marginal: lambda ~ 1';
    end
    fits(end+1) = struct('node', c, 'slope', p(2), 'intercept', p(1), ...
        'lambda', lam, 'resid', res, 'max_iqr', maxiqr, 'verdict', vd); %#ok<AGROW>
    fprintf('node %3d: intercept %+.4f  slope %+.4f  => lambda = %+.4f  (fit resid %.2e, max IQR %.2e)\n', ...
        c, p(1), p(2), lam, res, maxiqr);
    fprintf('          %s\n', vd);
end

P = struct('Bx', Bx, 'alphas', alphas, 'runs', runs, 'core', core, ...
    'rows', rows, 'fits', fits);
sp = [save_prefix '_summary.mat'];
save(sp, 'P');
fprintf('\nsaved: %s\n', sp);
end
