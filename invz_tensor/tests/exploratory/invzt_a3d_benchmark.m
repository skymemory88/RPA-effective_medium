function out = invzt_a3d_benchmark(opts)
%INVZT_A3D_BENCHMARK  7C Step-3 GATE: compact cc;cc dominant-vertex timing + N^4*nwn*nl
%   scaling fit + extrapolation to the 0.1 K production grid vs the 12 h budget.
%
%   out = INVZT_A3D_BENCHMARK()  builds the ORDERED electronuclear single ion at the deep
%   FM anchor (0.1 K, 3.0 T), and for N = 8, 10, 13, 16 takes the fixed-rank field-adapted
%   vertex basis (INVZT_ORDERED_VERTEX_BASIS, n_vertex = N), centers its Jz by the consumer
%   rule, and TIMES one compact cc;cc four-point build INVZT_GAMMA4({c},{{c,c}}) on a REDUCED
%   Matsubara grid (Ecut = 10). This is the dominant cost of one A3-dominant solve: the
%   vertex is built ONCE and cached (INVZT_SIGMA_TENSOR gamma4cc_cached); only the cheap
%   O(nwn*nl) K-contraction repeats each outer iteration. The tiling only partitions an
%   l-independent walk, so it does not change the arithmetic measured here.
%
%   It fits t_build(N) ~ A*N^p at the fixed measurement grid, then extrapolates to the
%   production grid (nwn = 740, nl = 1479) two ways:
%     (i)  DIRECT: production rank is 16 (= the largest measured N), so
%          t_prod = t_build(16) * (nwn_prod*nl_prod)/(nwn_meas*nl_meas)  -- grid-only scaling;
%     (ii) CONSERVATIVE N^4: anchor at the largest measured N and scale (N/Nanchor)^4.
%   HARD GATE (brief Step 3): the extrapolated one-build wall-clock at rank 16 must be
%   <= budget_hours (12). If it fails -> STOP (the follow-on options -- on-the-fly
%   contraction, Matsubara compression/tail treatment, transition/time-domain factorization
%   -- become the next plan); do NOT promise the production run.
%
%   DEV-TIME driver (never run by a test suite): returns a data struct and prints the table
%   for the controller to log to the ODD-LOG. No file writes. NOTE gpuDeviceCount('all') = 0
%   on this Apple-silicon machine -- the CPU path is the only verified path (round-2 P1-2).
%
%   opts (all optional):
%     .T           0.1          anchor temperature (K)
%     .Bx          3.0          anchor transverse field (T)
%     .N_list      [8 10 13 16] measured vertex ranks (n_vertex)
%     .Ecut_meas   10           reduced measurement Matsubara cutoff (meV)
%     .Lmax_meas   []           internal l cutoff for the measurement (default nwn_meas-1);
%                               set smaller to speed the measurement (extrapolation rescales)
%     .nwn_prod    740          production external Matsubara count
%     .nl_prod     1479         production internal l-grid size (= 2*Lmax_prod+1)
%     .Nprod       16           production vertex rank (the gate rank)
%     .budget_hours 12          HARD gate threshold
%     .grid_label  '16^3 halfopen, dpRng 30'   lattice provenance for the report
%
%   See also INVZT_SIGMA_TENSOR, INVZT_ORDERED_VERTEX_BASIS, INVZT_GAMMA4,
%   PERF_VERTEX_SCALING.
if nargin < 1 || isempty(opts), opts = struct(); end
gf = @(f, d) def(opts, f, d);
T           = gf('T',            0.1);
Bx          = gf('Bx',           3.0);
N_list      = gf('N_list',       [8 10 13 16]);
Ecut_meas   = gf('Ecut_meas',    10);
Lmax_meas   = gf('Lmax_meas',    []);
nwn_prod    = gf('nwn_prod',     740);
nl_prod     = gf('nl_prod',      1479);
Nprod       = gf('Nprod',        16);
budget_hours= gf('budget_hours', 12);
grid_label  = gf('grid_label',   '16^3 halfopen, dpRng 30');

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));                 % invz_tensor
addpath(fullfile(here, '..', '..', '..'));           % repo root
addpath(fullfile(here, '..', '..', '..', 'invz_common'));

% --- ordered single ion at the deep FM anchor (same si_full the a3d consumes) --------
ion = invz_ion();
g   = invzt_qgrid(16, 'halfopen');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 30, 'cache', true));   % loads from cache
J0z = lat.info.Jcc0;  Jxx0 = lat.info.Jaa0;
siopts  = struct('hyp', true, 'order', true, 'J0z', J0z, 'Jxx0', Jxx0, 'transverse_mf', 'legacy_x');
si_full = invz_single_ion(ion, T, [Bx 0 0], siopts);

% --- reduced measurement grid (Ecut_meas) --------------------------------------------
[wn_meas, ~, beta] = invz_matsubara(T, Ecut_meas);
nwn_meas = numel(wn_meas);
if isempty(Lmax_meas), Lmax_meas = nwn_meas - 1; end
nl_meas  = 2*Lmax_meas + 1;
next  = (0:nwn_meas-1).';
lvals = (-Lmax_meas:Lmax_meas).';

fprintf('\n=== A3d COMPACT DOMINANT-VERTEX BENCHMARK (7C Step 3) ===\n');
fprintf('anchor: T = %.3g K, Bx = %.3g T, lattice %s (n_full = %d)\n', ...
    T, Bx, grid_label, numel(si_full.E));
fprintf('measurement grid: Ecut = %g meV -> nwn = %d, l in [-%d,%d] (nl = %d)\n', ...
    Ecut_meas, nwn_meas, Lmax_meas, Lmax_meas, nl_meas);
fprintf('production grid:  nwn = %d, nl = %d ; production rank Nprod = %d ; budget = %g h\n\n', ...
    nwn_prod, nl_prod, Nprod, budget_hours);

% --- warm-up (JIT off the timed path) -------------------------------------------------
vb0 = invzt_ordered_vertex_basis(ion, T, si_full, struct('n_vertex', min(N_list)));
[es0, ops0] = vertex_es_ops(vb0);
invzt_gamma4(es0, ops0, {{'c','c'}}, {'c'}, next, lvals, beta);

% --- measure one compact cc;cc build per rank -----------------------------------------
t_meas    = zeros(numel(N_list), 1);
chi_share = zeros(numel(N_list), 1);
var_share = zeros(numel(N_list), 1);
fprintf('  %4s   %12s   %10s   %10s\n', 'N', 't_build [s]', 'chi_share', 'var_share');
for i = 1:numel(N_list)
    N = N_list(i);
    vb = invzt_ordered_vertex_basis(ion, T, si_full, struct('n_vertex', N));
    [es, ops] = vertex_es_ops(vb);
    tt = tic;
    invzt_gamma4(es, ops, {{'c','c'}}, {'c'}, next, lvals, beta);
    t_meas(i) = toc(tt);
    chi_share(i) = vb.chi_share;  var_share(i) = vb.var_share;
    fprintf('  %4d   %12.4f   %10.5f   %10.5f\n', N, t_meas(i), chi_share(i), var_share(i));
end

% --- fit log-log scaling --------------------------------------------------------------
cf   = polyfit(log(N_list(:)), log(t_meas(:)), 1);
pexp = cf(1);  Apre = exp(cf(2));
slopes = diff(log(t_meas(:))) ./ diff(log(N_list(:)));
fprintf('\nfit  t_build(N) ~ A*N^p :  p = %.3f  (A = %.3e s)   [expected ~4; vectorised', pexp, Apre);
fprintf(' innermost index depresses the measured slope]\n');
fprintf('pairwise local slopes: %s\n', mat2str(slopes.', 3));

% --- extrapolate to the production grid -----------------------------------------------
grid_ratio = (nwn_prod * nl_prod) / (nwn_meas * nl_meas);
Nanchor    = N_list(end);
[~, iprod] = min(abs(N_list - Nprod));
if N_list(iprod) == Nprod
    t_prod_direct = t_meas(iprod) * grid_ratio;         % measured rank, grid-only scaling
else
    t_prod_direct = Apre * Nprod^pexp * grid_ratio;     % fit if Nprod not measured
end
t_prod_N4  = t_meas(end) * grid_ratio * (Nprod/Nanchor)^4;
h_direct   = t_prod_direct / 3600;
h_N4       = t_prod_N4     / 3600;
h_gate     = max(h_direct, h_N4);
fits       = h_gate <= budget_hours;

fprintf('\nextrapolation to production (one cached vertex build; contraction O(nwn*nl) negligible):\n');
fprintf('  grid_ratio = (%d*%d)/(%d*%d) = %.4g\n', nwn_prod, nl_prod, nwn_meas, nl_meas, grid_ratio);
fprintf('  rank %d: DIRECT %.3g h | N^4-anchored %.3g h | GATE(max) %.3g h  -> %s (<= %g h)\n', ...
    Nprod, h_direct, h_N4, h_gate, tern(fits, 'FIT', 'STOP'), budget_hours);
% larger ranks for context (convergence-study candidates)
fprintf('  context (N^4-anchored one-build hours): ');
for N = [16 24 32]
    fprintf('N=%d: %.3g h   ', N, t_meas(end)*grid_ratio*(N/Nanchor)^4/3600);
end
fprintf('\n');

fprintf('\n--- VERDICT (7C Step 3 HARD GATE) ---\n');
if fits
    fprintf('PASS: extrapolated one-build wall-clock at rank %d = %.3g h <= %g h budget.\n', ...
        Nprod, h_gate, budget_hours);
else
    fprintf(['STOP: extrapolated one-build wall-clock at rank %d = %.3g h EXCEEDS %g h. Do NOT ' ...
        'promise the production run; escalate the follow-on plan (on-the-fly contraction, ' ...
        'Matsubara compression/tail treatment, transition/time-domain factorization).\n'], ...
        Nprod, h_gate, budget_hours);
end

out = struct();
out.N_list = N_list(:);  out.t_build_s = t_meas(:);
out.chi_share = chi_share(:);  out.var_share = var_share(:);
out.exponent = pexp;  out.prefactor_s = Apre;  out.local_slopes = slopes(:);
out.grid_ratio = grid_ratio;
out.hours_direct = h_direct;  out.hours_N4 = h_N4;  out.hours_gate = h_gate;
out.fits_budget = fits;  out.budget_hours = budget_hours;
out.params = struct('T', T, 'Bx', Bx, 'Ecut_meas', Ecut_meas, 'nwn_meas', nwn_meas, ...
    'Lmax_meas', Lmax_meas, 'nl_meas', nl_meas, 'nwn_prod', nwn_prod, 'nl_prod', nl_prod, ...
    'Nprod', Nprod, 'grid_label', grid_label);
end

% ---------------------------------------------------------------------- %
function [es, ops] = vertex_es_ops(vb)
% Reduced (es, ops) for the compact cc;cc build: ground-shifted E, normalized p, and the
% CONSUMER-CENTERED Jz (7B rule Mz - <Jz>_p * I), exactly as invzt_sigma_tensor forms them.
Ed  = vb.E(:) - vb.E(1);
pd  = vb.p(:);  pd = pd / sum(pd);
opc = vb.Mz - real(sum(pd .* diag(vb.Mz)))*eye(numel(Ed));
es  = struct('E', Ed, 'p', pd);
ops = struct('c', opc);
end

function v = def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
