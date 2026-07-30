%INVZT_RUN_PHASE_DIAGRAM  Full-tensor 1/z PM-side phase boundary: Bc(T) + Tc(B).
%
% TWO-REGIME SEARCH, mirroring the projected invz_run_phase_diagram: fixed-T
% FIELD CUTS (invzt_critical, low-T branch, one Bc per Ts entry) plus fixed-B
% TEMPERATURE CUTS (invzt_critical_T, the near-vertical branch approaching
% Tc(0) where a field cut crosses the boundary at a glancing angle and is
% ill-conditioned -- one Tc per Bs entry). Both kinds run in one flat parfor.
%
% SCOPE: PM side only -- there is still no ordered-phase tensor solve, so
% nothing below the boundary is computed. Unlike the projected T-cut (which
% cannot bracket with ODD on), the tensor A1 solver converges
% metastable PM points inside the ordered phase, so the T-cut brackets with
% odd on/off alike.
%
% OPTION NAMESPACES: solve_opts carries solver physics knobs shared by every
% point solve; the two finders get their OWN control knobs (Btol in tesla vs
% Ttol/Twidth/Tgridstep in kelvin -- both finders name their option 'tol',
% in different units, so they are merged deliberately at the call boundary,
% never shared).
%
% ERROR POLICY: the sweep absorbs ONLY per-point invzt:bracket (a genuine
% no-crossing outcome once Ts/Bs/Twin/Brange are preflighted below); every
% other identifier that ESCAPES the finders rethrows. The finders' own
% internal sampler additionally classifies shared-engine physics signals as
% invalid/ordered votes (their committed policy, documented in their headers).

addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_common'));

ion = invz_ion();

% ---- knobs ------------------------------------------------------------------
% Ts = linspace(0.4, 1.4, 11);       % fixed-T FIELD-CUT grid (K); low-T branch.
%                                        % [] disables field cuts.
% Bs = [];
Brange = [0.1 6];                   % field-cut [Blo Bhi] bracket (T). Blo > 0:
                                       % mode 'a1' forbids exact zero transverse field.
Btol   = 0.02;                         % field-cut bracket tolerance (TESLA) ->
                                       % invzt_critical opts.tol

Ts = [];                                       
Bs = linspace(0.01, 5.51, 56);       % fixed-B TEMPERATURE-CUT fields (T); the
                                       % near-vertical branch. [] disables T-cuts.
Twin   = [];                           % [] -> adaptive T-cut window anchored at the
                                       % small-Bx proxy Tc0 (computed below); or an
                                       % explicit HARD [Tlo Thi] bound (K) forwarded
                                       % to invzt_critical_T (no sliding; skips the
                                       % proxy anchor).
Ttol   = 0.005;                        % T-cut refinement tolerance (KELVIN) ->
                                       % invzt_critical_T opts.tol
Twidth = 1.0;                          % T-cut adaptive-window width (K)
Tgridstep = 1/30;                      % T-cut coarse-grid step (K)
gridN  = 16;  gridConv = 'halfopen';   % invzt_qgrid(gridN, gridConv)
dpRng  = 50;                           % invzt_jq_tensor coupling-sum range
useParallel = true;                    % false -> force serial
solve_opts  = struct('mode', 'a1', 'odd', true, 'nlevels', 'std', 'dress', 'full');
                                       % sigma_floor may be added here too; defaults to
                                       % invzt_critical's own -0.5 (single-sourced)

show_proxy_anchor     = true;          % plot the small-Bx Tc(0) proxy marker. The proxy
                                       % is COMPUTED regardless whenever Bs is nonempty
                                       % and Twin is [] (it anchors the adaptive T-cut
                                       % window); this knob only controls the marker.
Bproxy                = 0.05;          % transverse field (T) for that proxy
Ts_proxy              = 1.40:1/30:2.00;% proxy extrapolation grid (K); must contain >= 2
                                       % converged PM points at Bproxy (i.e. reach above
                                       % Tc0 ~ 1.56 K at production grids). NB proxy
                                       % samples are COLD (not branch-tracked): keep this
                                       % grid reaching WELL above Tc0 -- a too-low anchor
                                       % pushes the T-cuts onto the adaptive slide-up
                                       % path (final review; execution finding E2).
show_projected_anchor = false;         % OPT-IN cross-model comparator: the PROJECTED
                                       % closed-form Tc(0) (invz_odd_zero_field).
                                       % Off by default because it (a) puts
                                       % invz_projected on the path, (b) requires
                                       % ion.demag == 0, (c) solves 7 variants per
                                       % grid, and (d) is a DIFFERENT ODD treatment
                                       % (projected Tier-1 chi_perp-mediated dJ) from
                                       % this driver's structural tensor `odd` flag --
                                       % a comparator, NOT the same quantity.
% -----------------------------------------------------------------------------

% Preflights (invzt:bracket doubles as the finders' arg-validation id; absorbing
% it per point without these checks would turn a typo into a silent all-NaN
% sweep). invzt_str, not mat2str: mat2str throws on structs while BUILDING the
% intended error message.
assert(isnumeric(Brange) && isreal(Brange) && numel(Brange) == 2 && ...
    all(isfinite(Brange)) && Brange(2) > Brange(1) && Brange(1) > 0, ...
    'invzt_run_phase_diagram:Brange', ...
    'Brange must be finite [Blo Bhi] with 0 < Blo < Bhi (mode ''a1'' forbids B = 0); got %s.', ...
    invzt_str(Brange));
assert(isempty(Ts) || (isnumeric(Ts) && isreal(Ts) && isvector(Ts) && ...
    all(isfinite(Ts)) && all(Ts > 0)), 'invzt_run_phase_diagram:Ts', ...
    'Ts must be empty or a finite positive real vector; got %s.', invzt_str(Ts));
assert(isempty(Bs) || (isnumeric(Bs) && isreal(Bs) && isvector(Bs) && ...
    all(isfinite(Bs)) && all(Bs > 0)), 'invzt_run_phase_diagram:Bs', ...
    'Bs must be empty or a finite positive real vector (mode ''a1'' forbids B = 0); got %s.', ...
    invzt_str(Bs));
assert(isempty(Twin) || (isnumeric(Twin) && isreal(Twin) && numel(Twin) == 2 && ...
    all(isfinite(Twin)) && Twin(2) > Twin(1) && Twin(1) > 0 && Twin(2) > 0.02), ...
    'invzt_run_phase_diagram:Twin', ...
    ['Twin must be [] (adaptive) or a finite HARD [Tlo Thi] bound with 0 < Tlo < Thi ' ...
     'reaching above the 0.02 K solve floor; got %s.'], invzt_str(Twin));
% Finder control knobs: validate BEFORE the lattice build and proxy (the
% finders would reject or misuse them only inside a sweep job, after the
% expensive setup -- second review R4). invzt_critical does NOT validate its
% opts.tol at all, so a bad Btol would silently corrupt the bisection.
posk = @(x) isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x > 0;
assert(posk(Btol), 'invzt_run_phase_diagram:Btol', ...
    'Btol must be a finite positive real scalar (TESLA); got %s.', invzt_str(Btol));
assert(posk(Ttol), 'invzt_run_phase_diagram:Ttol', ...
    'Ttol must be a finite positive real scalar (KELVIN); got %s.', invzt_str(Ttol));
assert(posk(Twidth), 'invzt_run_phase_diagram:Twidth', ...
    'Twidth must be a finite positive real scalar (KELVIN); got %s.', invzt_str(Twidth));
assert(posk(Tgridstep), 'invzt_run_phase_diagram:Tgridstep', ...
    'Tgridstep must be a finite positive real scalar (KELVIN); got %s.', invzt_str(Tgridstep));

if show_projected_anchor
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_projected'));
    assert(ion.demag == 0, 'invzt:anchorDemag', ...
        ['show_projected_anchor requires ion.demag == 0 (invz_odd_blocks is ' ...
         'intrinsic-only). Set ion.demag = 0 or disable the anchor.']);
end

g   = invzt_qgrid(gridN, gridConv);
lat = invzt_jq_tensor(ion, g, struct('dpRng', dpRng, 'cache', true));

% ---- small-Bx proxy Tc(0): plot marker AND the adaptive T-cut window anchor --
Tc0_proxy = NaN;
need_anchor = ~isempty(Bs) && isempty(Twin);
if show_proxy_anchor || need_anchor
    critfun = @(T) invzt_crit_at(ion, T, [Bproxy 0 0], lat, solve_opts);
    try
        Tc0_proxy = invzt_tc_pm_extrap(critfun, Ts_proxy);
    catch err
        if ~strcmp(err.identifier, 'invzt:tcNoWindow'), rethrow(err); end
        fprintf('  Tc(0) proxy: no PM extrapolation window on Ts_proxy (skipped)\n');
    end
end
if need_anchor
    assert(isfinite(Tc0_proxy), 'invzt_run_phase_diagram:tcAnchor', ...
        ['The T-cut jobs need an adaptive-window anchor but the small-Bx proxy did ' ...
         'not resolve on Ts_proxy. Extend Ts_proxy above Tc0, or set an explicit Twin.']);
end

% Per-finder opts, merged deliberately at the call boundary (both finders name
% their tolerance 'tol' -- tesla for the field cut, kelvin for the T-cut).
bopts = solve_opts;  bopts.tol = Btol;
topts = solve_opts;  topts.tol = Ttol;  topts.width = Twidth;  topts.gridstep = Tgridstep;
if ~isempty(Twin), topts.window = Twin; else, topts.Tc0 = Tc0_proxy; end

% ---- one flat parfor over both cut kinds -------------------------------------
nWorkers = 0;
if useParallel && ~isempty(ver('parallel')), nWorkers = Inf; end

nT = numel(Ts);  nB = numel(Bs);
jobs = [Ts(:).' Bs(:).'];              % one independent 1-D root find per job
kind = [ones(1, nT) 2*ones(1, nB)];    % 1 = Bc(T) at fixed T, 2 = Tc(B) at fixed B
res  = nan(1, nT + nB);
parfor (k = 1:nT+nB, nWorkers)
    t0 = tic;  v = jobs(k);  val = NaN;
    % val-then-assign: keeps the sliced output unconditionally assigned.
    try
        if kind(k) == 1
            val = invzt_critical(ion, v, lat, Brange, bopts);
        else
            val = invzt_critical_T(ion, v, lat, topts);
        end
    catch err
        % invzt:bracket is a genuine per-point outcome (inputs preflighted
        % above, so this window/bracket simply lacks a crossing). Anything
        % else is a MISCONFIGURATION and rethrows.
        if ~strcmp(err.identifier, 'invzt:bracket'), rethrow(err); end
        if kind(k) == 1
            fprintf('  T = %.2f K: no bracket in [%.2f %.2f] T (left NaN)\n', ...
                    v, Brange(1), Brange(2));
        else
            fprintf('  B = %.2f T: no valid T-crossing (left NaN)\n', v);
        end
    end
    if kind(k) == 1
        fprintf('  [%2d/%2d] T = %.2f K  ->  Bc = %.3f T   (%.0f s)\n', ...
                k, nT+nB, v, val, toc(t0));
    else
        fprintf('  [%2d/%2d] B = %.2f T  ->  Tc = %.3f K   (%.0f s)\n', ...
                k, nT+nB, v, val, toc(t0));
    end
    res(k) = val;
end
Bc  = res(1:nT);                       % low-T branch:  Bc at each Ts
TcB = res(nT+1:end);                   % high-T branch: Tc at each Bs

% ---- optional projected closed-form Tc(0) comparator --------------------------
Tc0_closed = NaN;
if show_projected_anchor
    % Track this driver's own odd setting -- comparing a tensor ODD-on curve
    % against a projected ODD-off anchor (or vice versa) would be misleading.
    anchor_mode = 'off';
    if ~isfield(solve_opts, 'odd') || ~isequal(solve_opts.odd, false), anchor_mode = 'full'; end
    % Same nominal N and dipole cutoff as this driver's lattice, but on a
    % different model with a different ODD treatment: a cross-model
    % COMPARATOR, never the same quantity.
    Tc0_closed = invz_odd_zero_field(ion, struct('mode', anchor_mode, ...
        'grids', {{gridN}}, 'dpRng', dpRng, 'cache', true));
end

figure; hold on;
plot(Bc, Ts, 'o-', 'DisplayName', 'tensor A1: field-cut B_c(T)');
if nB > 0
    plot(Bs, TcB, 's-', 'DisplayName', 'tensor A1: T-cut T_c(B)');
end
if show_proxy_anchor && isfinite(Tc0_proxy)
    plot(0, Tc0_proxy, 'k^', 'MarkerFaceColor', 'y', ...
         'DisplayName', sprintf('tensor small-B_x proxy T_c(0), %.2f T', Bproxy));
end
if isfinite(Tc0_closed)
    plot(0, Tc0_closed, 'ks', 'MarkerFaceColor', 'c', ...
         'DisplayName', 'projected closed-form T_c(0) (comparator; legacy-inclusive mesh)');
end
xlabel('B_c (T)'); ylabel('T (K)');
title('LiHoF_4 full-tensor 1/z phase boundary (PM side: field cuts + T cuts)');
legend('Location', 'southwest');

% Merged boundary, T-sorted, finite points only. Columns [T(K) B(T)]. Near the
% regime join both branches can contribute a point, so the curve is not
% strictly single-valued in T there (same note as the projected driver).
phase_boundary = sortrows([Ts(:) Bc(:); TcB(:) Bs(:)], 1);
phase_boundary = phase_boundary(all(isfinite(phase_boundary), 2), :);
