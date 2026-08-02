function out = invzt_run_ladder(ion, opts)
%INVZT_RUN_LADDER  A4 basis-defined state-space ladder driver (DATA ONLY, budget-refusing).
%
%   out = INVZT_RUN_LADDER(ion, opts) climbs the historical Task-13 state-space ladder --
%   'three' (toy) then the
%   basis-content rungs 'e3'/'e6'/'e17' (lowest 3/6/17 CF states, multiplet-complete) and
%   their 'xI8' nuclear products -- and reports how the tensor 1/z ODD physics evolves with
%   BASIS CONTENT, cross-validated against the projected Tier-1+2 numbers. It RETURNS A
%   STRUCT and WRITES NO FILES; INVZT_REPORT_LADDER serializes a completed run for the
%   controller to paste wherever the run is being recorded (this driver writes nothing).
%
%   BUDGET REFUSAL (LOCKED, consumes Task-11's measured dense-vertex scaling). Each rung's
%   PROJECTED one-solve cost is estimated (conservative N^4 anchored at the Task-11 e6
%   boundary, 9.86 h; the canonical rungs quote T11's recorded projection). Any rung whose
%   projection exceeds opts.budget_hours (12) is REFUSED -- recorded in out.skipped_rungs,
%   NOT run. Under the LOCKED T11 verdict this affords {three, e3, e6} and refuses {e17
%   (~196 h), every xI8, e17xI8 = 136 (~5.9e5 h)}. Reaching the budget is not a failure; it
%   bounds which rungs the ladder reports.
%
%   PER-RUNG MEASUREMENT (at the validation point, seed-continuous). At one (T, B) point
%   [default the LOCKED validation point 1.6 K, (0.1 0 0) T, 6^3 / dpRng 10] the driver runs
%   the A1 bridge, A3 dominant-dress (E1-matched truncation), and A3 full-dress, each with
%   the ODD (c<->a,b) blocks OFF then ON. Seeds carry the odd-off fixed point into the
%   odd-on solve (single-root continuity, T6/T7 discipline) and across rungs. It records,
%   per rung:
%     crit_shift_odd  = crit(a3-full, odd on) - crit(a3-full, odd off)   [+ = ODD lowers Tc,
%                       the projected DeltaTc>0 direction: at fixed T below Tc0 the ODD system
%                       is less deep in the ordered phase, so its criticality margin is higher]
%     rf, rd          = full-A3 / A1 and dominant-dress-A3 / A1 odd-shift ratios; rf is the
%                       beyond-E1 transverse-spectator dressing, rd the matched E1 truncation
%     matched_collapse= |rd-1|/|rf-1| (dropping the transverse dressing shrinks the beyond-E1
%                       excess -- the T12 collapse, reported here at the ladder point)
%     eps_el, eps_cross, sumrule_rel, resum_spread   (A3 monitors, constraints 7/8/9/6)
%     chi0_virtual_deficit = ||chi0_full136 - chi0_rung||_F / ||chi0_full136||_F, static (T,B):
%                       a DIAGNOSTIC of the virtual intermediate states the rung's basis is
%                       missing -- NOT a bound on the missing four-point vertex.
%     dim_actual, multiplet_complete, N, npaths (=6 N^4 vertex path count), t_sec.
%   opts.tc true additionally computes the SMALL-Bx-PROXY Tc (0.05 T) per run rung via
%   INVZT_TC_PM_EXTRAP (tensor A3 true-B=0 Tc is DEFERRED, v3 zero-field policy -- EVERY
%   rung's tc is the proxy, never a true-zero-field number).
%
%   OPTIONS (getf defaults):
%     rungs        {'three','e3'}  cellstr of rungs to attempt (budget may refuse some).
%     T            1.6             validation temperature (K).
%     B            [0.1 0 0]       validation field (T); transverse component must be > 1e-6.
%     ngrid        6 (prod 8)      q-grid points per axis (halfopen).
%     conv         'halfopen'      q-grid convention.
%     dipole       'ewald'         'ewald' | 'bruteforce'.
%     ewald        certified       exact Ewald controls (Ewald backend only).
%     dpRng        10 (prod 20)    brute-force lattice-sum range (diagnostic only).
%     production   false           true bumps ngrid/brute dpRng to production sizes.
%     Ecut         40              Matsubara cutoff (meV).
%     budget_hours 12              refusal threshold (projected one-solve hours).
%     tc           false           add the small-Bx-proxy Tc per rung (SLOW; production).
%     tc_bx        0.05            transverse proxy field for the Tc extrapolation (T).
%     tc_grid      1.30:1/30:1.80  ascending T-grid (K) for the proxy Tc extrapolation.
%     Esplit       0.4653          dominant/rest split energy (meV) for dress='dominant'.
%     cache        true            cache the lattice tensor.
%     verbose      true            progress printing (stdout only; no file writes).
%
%   See also INVZT_RUNG_BASIS, INVZT_SOLVE_POINT, INVZT_REPORT_LADDER,
%   INVZT_TC_PM_EXTRAP, INVZT_SIGMA_TENSOR.
if nargin < 2 || isempty(opts), opts = struct(); end
production = getf(opts, 'production', false);
rungs   = getf(opts, 'rungs',   {'three', 'e3'});
if ischar(rungs) || isstring(rungs), rungs = cellstr(rungs); end
T       = getf(opts, 'T', 1.6);
B       = invz_field_vec(getf(opts, 'B', [0.1 0 0]));
ngrid   = getf(opts, 'ngrid', tern(production, 8, 6));
conv    = getf(opts, 'conv',  'halfopen');
dpRng   = getf(opts, 'dpRng', tern(production, 20, 10));
dipole  = getf(opts, 'dipole', 'ewald');
eopts   = getf(opts, 'ewald', invzt_ewald_defaults(ion));
Ecut    = getf(opts, 'Ecut', 40);
budget  = getf(opts, 'budget_hours', 12);
do_tc   = getf(opts, 'tc', false);
tc_bx   = getf(opts, 'tc_bx', 0.05);
tc_grid = getf(opts, 'tc_grid', 1.30:(1/30):1.80);
Esplit  = getf(opts, 'Esplit', 0.4653);
cache   = getf(opts, 'cache', true);
verbose = getf(opts, 'verbose', true);

% --- lattice (built once; shared by every rung) ----------------------------------------
g   = invzt_qgrid(ngrid, conv);
latopts = lattice_opts(dipole, eopts, dpRng, cache);
lat = invzt_jq_tensor(ion, g, latopts);

% --- full-136 static chi0 at (T,B): the virtual-completeness diagnostic reference -------
si_full = invz_single_ion(ion, T, B, struct('hyp', true, 'transverse_mf', 'legacy_x', ...
    'Jxx0', lat.info.Jaa0));
c0_full = invz_chi0z(si_full, T, 0, struct('elastic', true));
nrm_full = norm(c0_full, 'fro');

% --- attempt each rung: budget-refuse or run+measure -----------------------------------
skipped = struct('name', {}, 'dim', {}, 'N', {}, 'projected_hours', {}, 'reason', {});
R = struct([]);                                          % per-run-rung result accumulator
seeds = struct('a1', [], 'dom_V', [], 'dom_S', [], 'full_V', [], 'full_S', []);
for ir = 1:numel(rungs)
    nl = char(rungs{ir});
    [N, dim_actual, mcomplete] = rung_dim(ion, nl);
    [phours, hbasis] = rung_cost_hours(nl, N);
    if phours > budget
        skipped(end+1) = struct('name', nl, 'dim', dim_actual, 'N', N, ...
            'projected_hours', phours, 'reason', sprintf('projected %.4g h > budget %g h (%s)', ...
            phours, budget, hbasis));  %#ok<AGROW>
        if verbose
            fprintf('[ladder] REFUSE %-8s N=%-3d projected %.4g h > budget %g h -- skipped.\n', ...
                nl, N, phours, budget);
        end
        continue;
    end
    if verbose
        fprintf('[ladder] RUN    %-8s N=%-3d dim_actual=%-3d (projected %.3g h <= budget) ...\n', ...
            nl, N, dim_actual, phours);
    end
    t0 = tic;
    [res, seeds] = run_one_rung(ion, T, B, lat, nl, Esplit, Ecut, seeds, verbose);
    res.t_sec = toc(t0);
    res.name = nl;  res.N = N;  res.dim_actual = dim_actual;  res.multiplet_complete = mcomplete;
    res.projected_hours = phours;  res.hours_basis = hbasis;  res.npaths = 6 * N^4;
    % --- virtual-completeness diagnostic: rung static chi0 vs full-136 ------------------
    c0_rung = invz_chi0z(res.si, T, 0, struct('elastic', true));
    res.chi0_virtual_deficit = norm(c0_full - c0_rung, 'fro') / max(nrm_full, eps);
    res.chi0_rung_ccdiag = real(c0_rung(3,3));
    % --- small-Bx-proxy Tc (optional; DEFERRED true-B=0, proxy only) --------------------
    res.tc_proxy = NaN;  res.tc_used = [NaN NaN];
    if do_tc
        [res.tc_proxy, res.tc_used] = rung_proxy_tc(ion, lat, nl, tc_bx, tc_grid, Ecut, Esplit, verbose);
    end
    if isempty(R), R = res; else, R(end+1) = res; end   %#ok<AGROW>
    if verbose
        fprintf(['[ladder]   %-8s crit_shift_odd=%+.4e conv=%d | rf=%.3f rd=%.3f collapse=%.3f | ' ...
            'eps_el=%.3f eps_cross=%.4f sumrule=%.3f resum=%.2e | vdef=%.3f (%.1fs)\n'], ...
            nl, res.crit_shift_odd, res.converged, res.rf, res.rd, res.matched_collapse, ...
            res.eps_el, res.eps_cross, res.sumrule_rel, res.resum_spread, res.chi0_virtual_deficit, res.t_sec);
    end
end

% --- assemble the scalar-struct/array-field output -------------------------------------
% Empty-safe (Minor #1): when EVERY requested rung is budget-refused (e.g. rungs={'e17'}),
% R is an empty struct array -- {R.name}/[R.field] would error, so build empty column
% fields directly and still return a clean struct carrying skipped_rungs + meta.
nr = numel(R);
out = struct();
if nr > 0
    out.rungs = {R.name}.';
    out.dim_actual = colv([R.dim_actual]);
    out.multiplet_complete = colv(logical([R.multiplet_complete]));
    out.N = colv([R.N]);
    out.npaths = colv([R.npaths]);
    out.projected_hours = colv([R.projected_hours]);
    out.crit_shift_odd = colv([R.crit_shift_odd]);
    out.crit_oddoff = colv([R.crit_oddoff]);
    out.crit_oddon  = colv([R.crit_oddon]);
    out.converged = colv(logical([R.converged]));
    out.rf = colv([R.rf]);
    out.rd = colv([R.rd]);
    out.matched_collapse = colv([R.matched_collapse]);
    out.a1_odd_shift = colv([R.a1_odd_shift]);
    out.eps_el = colv([R.eps_el]);
    out.eps_cross = colv([R.eps_cross]);
    out.sumrule_rel = colv([R.sumrule_rel]);
    out.resum_spread = colv([R.resum_spread]);
    out.chi0_virtual_deficit = colv([R.chi0_virtual_deficit]);
    out.chi0_rung_ccdiag = colv([R.chi0_rung_ccdiag]);
    out.t_sec = colv([R.t_sec]);
    out.tc_proxy = colv([R.tc_proxy]);
    out.tc_used = reshape([R.tc_used], 2, nr).';
else
    out.rungs = cell(0, 1);
    for f = {'dim_actual','N','npaths','projected_hours','crit_shift_odd','crit_oddoff', ...
             'crit_oddon','rf','rd','matched_collapse','a1_odd_shift','eps_el','eps_cross', ...
             'sumrule_rel','resum_spread','chi0_virtual_deficit','chi0_rung_ccdiag','t_sec','tc_proxy'}
        out.(f{1}) = zeros(0, 1);
    end
    out.multiplet_complete = false(0, 1);
    out.converged = false(0, 1);
    out.tc_used = zeros(0, 2);
end
out.skipped_rungs = skipped;
% provenance (no file writes): the report serializer + controller consume this.
out.meta = struct('T', T, 'B', B(:).', 'ngrid', ngrid, 'conv', conv, ...
    'dipole', lat.info.dipole, 'dpRng', lat.info.dpRng, ...
    'Ecut', Ecut, 'production', production, 'budget_hours', budget, 'tc', do_tc, ...
    'tc_bx', tc_bx, 'Esplit', Esplit, 'qhash', g.hash, ...
    'Jcc0', lat.info.Jcc0, 'Jaa0', lat.info.Jaa0, ...
    'chi0_full136_ccdiag', real(c0_full(3,3)), 'chi0_full136_perpdiag', real(c0_full(1,1)), ...
    'nrm_chi0_full136', nrm_full, 'date', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    'git', local_git_hash(), 'nrun', nr, 'nskip', numel(skipped));
end

function opts = lattice_opts(backend, eopts, dpRng, cache)
opts = struct('dipole', backend, 'cache', cache);
if strcmp(backend, 'ewald')
    opts.ewald = eopts;
elseif strcmp(backend, 'bruteforce')
    opts.dpRng = dpRng;
else
    error('invzt:ladderDipoleBackend', ...
        'opts.dipole must be ''ewald'' or ''bruteforce'' (got ''%s'').', backend);
end
end

% ======================================================================================= %
%  One rung: 6 seed-continuous solves (a1/a3dom/a3full x odd off/on) at (T,B).
% ======================================================================================= %
function [res, seeds] = run_one_rung(ion, T, B, lat, nl, Esplit, Ecut, seeds, verbose)
base = struct('nlevels', nl, 'Ecut', Ecut, 'Esplit', Esplit);
o = @(m, od, extra) mergestruct(base, mergestruct(struct('mode', m, 'odd', od), extra));

% --- odd OFF baselines (warm-started across rungs where a seed exists) ------------------
b1 = invzt_solve_point(ion, T, B, lat, o('a1', false, seedopt('Sigma_seed', seeds.a1)));
bd = invzt_solve_point(ion, T, B, lat, o('a3', false, ...
    mergestruct(struct('dress', 'dominant'), seedopt2('Vmat_seed', seeds.dom_V, 'Sigma_seed', seeds.dom_S))));
bf = invzt_solve_point(ion, T, B, lat, o('a3', false, ...
    mergestruct(struct('dress', 'full'), seedopt2('Vmat_seed', seeds.full_V, 'Sigma_seed', seeds.full_S))));

% --- odd ON (warm-started from THIS rung's odd-off fixed point: single-root continuity) -
a1 = invzt_solve_point(ion, T, B, lat, o('a1', true, struct('Sigma_seed', b1.Sigma)));
ad = invzt_solve_point(ion, T, B, lat, o('a3', true, ...
    struct('dress', 'dominant', 'Vmat_seed', bd.Vmat, 'Sigma_seed', bd.Sigma)));
af = invzt_solve_point(ion, T, B, lat, o('a3', true, ...
    struct('dress', 'full', 'Vmat_seed', bf.Vmat, 'Sigma_seed', bf.Sigma)));

% --- carry seeds to the next rung (odd-off fixed points) --------------------------------
seeds.a1 = b1.Sigma;  seeds.dom_V = bd.Vmat;  seeds.dom_S = bd.Sigma;
seeds.full_V = bf.Vmat;  seeds.full_S = bf.Sigma;

% --- odd shifts and ratios (rf beyond-E1, rd matched E1) -------------------------------
d1 = a1.crit - b1.crit;                                  % A1 odd shift
ddom = ad.crit - bd.crit;                                % A3 dominant-dress odd shift
dfull = af.crit - bf.crit;                               % A3 full-dress odd shift
res = struct();
res.crit_oddoff = bf.crit;  res.crit_oddon = af.crit;
res.crit_shift_odd = dfull;                              % + = ODD lowers Tc (projected DeltaTc direction)
res.a1_odd_shift = d1;
if abs(d1) > 0
    res.rf = dfull / d1;  res.rd = ddom / d1;
else
    res.rf = NaN;  res.rd = NaN;
end
res.matched_collapse = abs(res.rd - 1) / max(abs(res.rf - 1), eps);
% converged = the PRIMARY a3-full physics (off & on); the ratio inputs are NaN-tolerant.
res.converged = bf.converged && af.converged;
res.aux_converged = b1.converged && a1.converged && bd.converged && ad.converged;
% monitors from the full-A3 odd-on solve (constraints 7/8/9/6)
res.eps_el = af.eps_el;  res.eps_cross = af.eps_cross;
res.sumrule_rel = af.sumrule_rel;  res.resum_spread = af.resum_spread_crit;
res.si = bf.si;                                          % for the virtual-deficit diagnostic
if verbose && ~res.aux_converged
    fprintf('[ladder]   NOTE %-8s an auxiliary (a1/a3dom) solve did not converge; rf/rd reported as-is.\n', nl);
end
end

% ======================================================================================= %
%  Rung Hilbert dimension (N drives the vertex cost) + reported dim_actual.
% ======================================================================================= %
function [N, dim_actual, mcomplete] = rung_dim(ion, nl)
if strcmp(nl, 'three')
    N = 3;  dim_actual = 3;  mcomplete = true;           % synthetic toy (not a truncation)
else
    rb = invzt_rung_basis(ion, nl);
    N = rb.dim_actual;  dim_actual = rb.dim_actual;  mcomplete = rb.multiplet_complete;
end
end

% ======================================================================================= %
%  Projected one-solve cost (h). Consumes Task-11's measured DENSE-vertex scaling: the
%  canonical rungs quote T11's recorded projection; any other uses the conservative N^4
%  anchored at the T11 e6 boundary (9.86 h at N=6), which is >= the recorded piecewise
%  model everywhere and so reproduces the LOCKED verdict at any budget.
% ======================================================================================= %
function [hours, basis] = rung_cost_hours(nl, N)
switch nl
    case 'e6',     hours = 9.86;    basis = 'T11-recorded';
    case 'e17',    hours = 196;     basis = 'T11-recorded';
    case 'e17xI8', hours = 5.9e5;   basis = 'T11-recorded';
    otherwise,     hours = 9.86 * (N/6)^4;  basis = 'N4-anchored-e6(9.86h)';
end
end

% ======================================================================================= %
%  Small-Bx-proxy Tc for one rung (a3-full, ODD on). True-B=0 Tc is DEFERRED (v3): this is
%  the 0.05 T proxy for EVERY rung. Warm-starts each T from the previous grid point.
% ======================================================================================= %
function [tc, used] = rung_proxy_tc(ion, lat, nl, bx, Tg, Ecut, Esplit, verbose)
Bp = [bx 0 0];
persistentseed = struct('V', [], 'S', []);
    function [crit, ok] = critfun(Tt)
        try
            ex = struct('dress', 'full', 'Ecut', Ecut, 'Esplit', Esplit);
            if ~isempty(persistentseed.V), ex.Vmat_seed = persistentseed.V; ex.Sigma_seed = persistentseed.S; end
            pt = invzt_solve_point(ion, Tt, Bp, lat, mergestruct(struct('mode', 'a3', ...
                'nlevels', nl, 'odd', true), ex));
            crit = pt.crit;  ok = pt.converged && isfinite(crit);
            if ok, persistentseed.V = pt.Vmat;  persistentseed.S = pt.Sigma; end
        catch ME
            crit = NaN;  ok = false;
            if verbose, fprintf('[ladder]   tc(%s) T=%.3f: %s\n', nl, Tt, ME.message); end
        end
    end
try
    [tc, used] = invzt_tc_pm_extrap(@critfun, Tg);
catch ME
    tc = NaN;  used = [NaN NaN];
    if verbose, fprintf('[ladder]   tc(%s): no PM window (%s)\n', nl, ME.message); end
end
end

% ------------------------------------- small helpers ----------------------------------- %
function v = colv(x)
v = x(:);
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function s = mergestruct(s, t)
% Shallow overlay of struct t onto struct s (t wins).
f = fieldnames(t);
for i = 1:numel(f), s.(f{i}) = t.(f{i}); end
end

function s = seedopt(field, val)
% {field: val} if val nonempty, else empty struct (no seed -> cold start).
if isempty(val), s = struct(); else, s = struct(field, val); end
end

function s = seedopt2(f1, v1, f2, v2)
s = struct();
if ~isempty(v1), s.(f1) = v1; end
if ~isempty(v2), s.(f2) = v2; end
end

function h = local_git_hash()
try
    here = fileparts(mfilename('fullpath'));
    [st, o] = system(sprintf('git -C "%s" rev-parse --short HEAD 2>/dev/null', here));
    if st == 0, h = strtrim(o); else, h = 'unknown'; end
catch
    h = 'unknown';
end
end
