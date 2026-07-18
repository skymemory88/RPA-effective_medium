function out = perf_vertex_scaling(opts)
%PERF_VERTEX_SCALING  A3 dense-vertex PERFORMANCE GATE (Task 11, INVZ_SLOW).
%
%   out = PERF_VERTEX_SCALING()  measures the wall time of one DENSE Vmat
%   assembly (all 6 independent external Cartesian pairs, internal channels
%   {a,b,c} against a synthetic K) with invzt_vertex4 at N = 3, 6, 12, 24
%   (synthetic random-Hermitian spectra), fits the dense scaling exponent
%   (expected ~ N^4 per time-ordering), and PROJECTS the cost of one production
%   A3 solve (~30 outer iterations) at each ladder rung -- three, e3, e6, e17 and
%   the xI8 products up to 136 (e3xI8=24, e6xI8=48, e17xI8=136).
%
%   STOP RULE (two independent triggers, brief resolution 6):
%     (i)  FACTORIZATION NOT ESTABLISHED -- Task-10 check 10 recorded
%          factored_ok = false, so the O(N^3) 'factored' path is disabled and the
%          dense O(N^4) engine is the sole costed path.  Large rungs may then be
%          budget-refused.  (Recorded here; no factored projection is emitted.)
%     (ii) BUDGET -- any rung whose PROJECTED one-solve time exceeds
%          opts.budget_hours (12) is flagged for Task-13 REFUSAL.  Reaching the
%          budget is NOT a failure; it bounds which rungs Task 13 reports.
%
%   This is a DEV-TIME driver (never run by a test suite): it returns a data
%   struct and prints the scaling table + per-rung projection + verdict for the
%   controller to log to docs/ODD-LOG.md (Section A3).  No file writes.
%
%   opts (all optional; defaults are production-representative):
%     .N_list      [3 6 12 24]   measured Hilbert dimensions
%     .beta        7.73          inverse T (meV^-1) ~ 1.5 K
%     .nwn_meas    1             external Matsubara freqs used in the measurement
%     .Lmax_meas   1             internal l-cutoff used in the measurement (l=-1..1)
%     .nwn_prod    51            production external freqs (Ecut=40, T~1.5K, invz_matsubara)
%     .Lmax_prod   50            production internal l-cutoff
%     .outer_iters 30            outer self-consistent iterations per A3 solve
%     .budget_hours 12           refusal threshold
%     .seed        20260718
%
%   See also INVZT_VERTEX4, INVZT_KERNELS.

if nargin < 1 || isempty(opts), opts = struct(); end
gf = @(f, d) def(opts, f, d);
N_list      = gf('N_list',      [3 6 12 24]);
beta        = gf('beta',        7.73);
nwn_meas    = gf('nwn_meas',    1);
Lmax_meas   = gf('Lmax_meas',   1);
nwn_prod    = gf('nwn_prod',    51);
Lmax_prod   = gf('Lmax_prod',   50);
outer_iters = gf('outer_iters', 30);
budget_hours= gf('budget_hours',12);
seed        = gf('seed',        20260718);

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));                 % invz_tensor (invzt_vertex4/kernels)
fxpath = fullfile(here, '..', 'fixtures', 'vertex_oracle.json');
fx = jsondecode(fileread(fxpath));
factored_ok = logical(fx.factored.factored_ok);

extpairs = {{'a','a'},{'a','b'},{'a','c'},{'b','b'},{'b','c'},{'c','c'}};   % 6 independent
comps    = {'a','b','c'};
wn_meas  = 1:nwn_meas;

% ---- warm-up (JIT / first-call overhead off the timed path) ----
run_assembly(3, beta, extpairs, comps, wn_meas, Lmax_meas, seed);

% ---- measure ----
t_meas = zeros(numel(N_list), 1);
fprintf('\n=== A3 DENSE VERTEX PERFORMANCE GATE (Task 11) ===\n');
fprintf('measured workload per assembly: %d external pairs, %d internal channels (%d pairs), ', ...
    numel(extpairs), numel(comps), numel(comps)^2);
fprintf('nwn=%d, l in [-%d,%d]\n\n', nwn_meas, Lmax_meas, Lmax_meas);
fprintf('  %4s   %12s\n', 'N', 't_assy [s]');
for i = 1:numel(N_list)
    N = N_list(i);
    tt = tic;
    run_assembly(N, beta, extpairs, comps, wn_meas, Lmax_meas, seed);
    t_meas(i) = toc(tt);
    fprintf('  %4d   %12.4f\n', N, t_meas(i));
end

% ---- fit log-log scaling ----
cf   = polyfit(log(N_list(:)), log(t_meas(:)), 1);
pexp = cf(1);
Apre = exp(cf(2));
slopes = diff(log(t_meas(:))) ./ diff(log(N_list(:)));      % pairwise local slopes
asym = slopes(end);
fprintf('\nfit  t_assy(N) ~ A * N^p :  p = %.3f  (A = %.3e s)\n', pexp, Apre);
fprintf('pairwise local slopes %s -> asymptotic %.3f\n', mat2str(slopes.', 3), asym);
fprintf(['NOTE: measured wall-time exponent (~%.1f) is BELOW the N^4 arithmetic complexity ', ...
    'because\n  the innermost eigenstate index is vectorised (SIMD amortises per-call ', ...
    'overhead);\n  the local slope is CLIMBING toward 4, so the budget gate uses the ', ...
    'CONSERVATIVE\n  theoretical N^4 anchored at the largest measured N (N=%d).\n'], asym, N_list(end));

% ---- project one full A3 solve at each rung ----
% production one-solve time = t_assy(N) * (nwn_prod/nwn_meas)
%                             * ((2 Lmax_prod+1)/(2 Lmax_meas+1)) * outer_iters
% Budget basis (conservative): anchor at the largest measured N and extrapolate at
% the theoretical N^4.  Empirical fit shown alongside for transparency.
freq_scale  = (nwn_prod / nwn_meas) * ((2*Lmax_prod + 1) / (2*Lmax_meas + 1));
solve_mult  = freq_scale * outer_iters;                       % t_assy -> one full solve
t_solve_anchor = t_meas(end) * solve_mult;                    % one solve at N = N_list(end)
Nanchor = N_list(end);
rungs = {'three',3; 'e3',3; 'e6',6; 'e17',17; 'e3xI8',24; 'e6xI8',48; 'e17xI8',136};
nr = size(rungs, 1);
rung_name = cell(nr,1); rung_N = zeros(nr,1); rung_hours = zeros(nr,1);
rung_hours_emp = zeros(nr,1); fits = false(nr,1);
fprintf('\nprojection basis: nwn_prod=%d, l in [-%d,%d], outer_iters=%d, budget=%d h\n', ...
    nwn_prod, Lmax_prod, Lmax_prod, outer_iters, budget_hours);
% Budget gate = MAX of the two models per rung (conservative both directions:
% empirical fit dominates BELOW the anchor, N^4 dominates ABOVE it).
fprintf('\n  %-8s %5s   %16s   %14s   %12s   %s\n', 'rung', 'N', 't_solve N^4 [h]', ...
    'emp-fit [h]', 'gate [h]', 'budget');
for i = 1:nr
    name = rungs{i,1}; N = rungs{i,2};
    hrs     = (t_solve_anchor * (N/Nanchor)^4) / 3600;        % conservative N^4
    hrs_emp = (Apre * N^pexp * solve_mult) / 3600;            % empirical fit
    hrs_gate = max(hrs, hrs_emp);
    ok = hrs_gate <= budget_hours;
    rung_name{i} = name; rung_N(i) = N; rung_hours(i) = hrs; rung_hours_emp(i) = hrs_emp; fits(i) = ok;
    fprintf('  %-8s %5d   %16.4g   %14.4g   %12.4g   %s\n', name, N, hrs, hrs_emp, hrs_gate, ...
        tern(ok, 'FIT', 'REFUSE (>12h -> Task 13)'));
end

% ---- verdict ----
fprintf('\n--- VERDICT ---\n');
fprintf('STOP TRIGGER (i)  FACTORIZATION: %s (factored_ok=%d) -> dense-only engine.\n', ...
    tern(factored_ok, 'ESTABLISHED', 'NOT ESTABLISHED'), factored_ok);
refused = rung_name(~fits);
if isempty(refused)
    fprintf('STOP TRIGGER (ii) BUDGET: no rung exceeds %d h.\n', budget_hours);
else
    fprintf('STOP TRIGGER (ii) BUDGET: rungs refused (>%d h) for Task 13: %s\n', ...
        budget_hours, strjoin(refused, ', '));
end
fprintf('affordable rungs: %s\n', strjoin(rung_name(fits), ', '));
fprintf('optimisation backlog (Task-13 next steps): prove/use the factored identity, ');
fprintf('transition/Liouville-space contraction, time-simplex matrix quadrature, ');
fprintf('symmetry blocking by electronuclear quantum numbers, omega-grid compression, ');
fprintf('tail subtraction, cached transition sums.\n');

out = struct();
out.N_list = N_list(:); out.t_meas = t_meas(:);
out.exponent = pexp; out.prefactor_s = Apre; out.local_slopes = slopes(:);
out.asymptotic_slope = asym;
out.rungs = rung_name; out.rung_N = rung_N;
out.rung_solve_hours = rung_hours;          % conservative N^4 (budget basis)
out.rung_solve_hours_empfit = rung_hours_emp;
out.rung_fits_budget = fits; out.refused_rungs = refused;
out.affordable_rungs = rung_name(fits);
out.factored_ok = factored_ok;
out.params = struct('beta', beta, 'nwn_meas', nwn_meas, 'Lmax_meas', Lmax_meas, ...
    'nwn_prod', nwn_prod, 'Lmax_prod', Lmax_prod, 'outer_iters', outer_iters, ...
    'budget_hours', budget_hours, 'freq_scale', freq_scale, 'seed', seed);
end

% ---------------------------------------------------------------------- %
function run_assembly(N, beta, extpairs, comps, wn, Lmax, seed)
% One synthetic dense Vmat assembly at Hilbert dimension N.
rng(seed + N);
E = sort(5 * rand(N, 1));                       % spread ~5 meV, no accidental degeneracy
E = E - min(E);
w = exp(-beta * E); p = w / sum(w);
ops = struct();
for c = 1:numel(comps)
    A = randn(N) + 1i * randn(N); A = (A + A') / 2;   % random Hermitian
    av = real(sum(p .* diag(A))); A = A - av * eye(N);
    ops.(comps{c}) = A;
end
nc = numel(comps);
Karr = zeros(nc, nc, Lmax + 1);
for l = 0:Lmax
    wl = 2 * pi * l / beta;
    Karr(:, :, l + 1) = eye(nc) * (0.3 / (1 + (wl / 2.0)^2));
end
es = struct('E', E, 'p', p);
opts = struct('stage', 'V', 'Lmax', Lmax);
opts.comps = comps; opts.ext = extpairs;
invzt_vertex4(es, ops, Karr, wn, beta, opts);
end

function v = def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
