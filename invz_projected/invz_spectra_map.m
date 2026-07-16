function S = invz_spectra_map(ion, T, fields, w, opts)
%INVZ_SPECTRA_MAP chi''_cc(omega, B) at the uniform (q=0) mode over a field-magnitude sweep.
%   S = invz_spectra_map(ion, T, fields, w) sweeps the field magnitude |B| over `fields`
%   (Tesla, nonnegative) along opts.field_dir at fixed temperature T (K) and returns
%   chi''_cc(omega, B) across BOTH phases: at each field it tries the ordered 1/z solve
%   first, falling back to the paramagnetic solve, giving a single soft-mode map continuous
%   across the transition.
%
%   Returned maps (columns = fields):
%     S.chiz   [nw x nB]  1/z-renormalized chi''_cc  (moment-form below Bc, PM above)
%     S.chirpa [nw x nB]  bare-RPA chi''_cc          (Sigma = 0, matching phase)
%     With ion.demag ~= 0 both are the demag-corrected MEASURED observable (via
%     info.Jshape_cc, saturating instead of diverging); demag = 0 gives the intrinsic response.
%   Per-field diagnostics: S.phase (1 = moment-form (spontaneous FM below Bc, or field-induced
%   under a longitudinal tilt -- a rounded crossover, no sharp Bc), 2 = strict paramagnet,
%   0 = masked), S.Sigma0, and S.Epeak/S.Epeak_rpa (censored, parabolic-refined peak energy;
%   NaN at a non-positive or boundary maximum, via invz_peak_energy, shared with
%   invz_spectra_qpath). Fields with no solution at all (e.g. the degenerate doublet at
%   Bx -> 0) are left NaN and masked out.
%   S.field_dir [1x3]  normalized field direction actually used.
%   S.Bvec      [nB x 3] field vectors actually used (fields(:) * field_dir, dead band applied).
%
%   opts fields (all optional):
%     .hyp       (true)        include the nuclear hyperfine manifold in chi0
%     .grid      ([16 16 16])  q-grid for the effective-medium lattice sum
%     .dpRng     (30)          real-space dipole cutoff for invz_jq_modes
%     .eta       (5e-3)        real-axis broadening (meV) passed to invz_chi_realaxis
%     .parallel  (false)       distribute the field points over a parallel pool (parfor)
%     .peak_wmin (0)           meV; excludes a known low-frequency line (e.g. hyperfine
%                              pole) from the peak search; default is no exclusion
%     .Jnu, .info              precomputed coupling branches / info struct (skips the
%                              lattice sum; used by the tests)
%     .verbose   (true)        print one progress line per field
%     .field_dir ([1 0 0])     nonzero finite real 3-vector, normalized internally; sets the
%                              sweep direction of `fields` (|B|). Error invz:fieldDir if
%                              invalid. A nonzero z-component routes through the longitudinal
%                              (field-induced moment) solve once |Bz| clears .bz_tol.
%                              Validated envelope: ac-plane directions [cos(theta) 0 sin(theta)] with theta_c <= 5 deg (see docs/SESSION-2026-07-16-field-angle.md); By ~= 0 runs under the legacy x-only transverse MF and demag ~= 0 with tilt is unvalidated.
%     .bz_tol    (1e-9)        T; dead band on Bz -- resolved ONCE, applied to the field table
%                              BEFORE any solve, and forwarded to invz_solve_auto/one_field.
%     .solve_opts (struct())   merged into the per-field invz_solve_auto opts; fields
%                              J0eff/Jxx0/hyp are reserved (driver-owned) -> error invz:solveOpts.
%                              transverse_mf ('legacy_x' | 'none' | 'vector_ab', default
%                              'legacy_x') is a legal solve_opts field forwarded to the
%                              solvers. Under 'legacy_x' (x-only mean field) a nonzero
%                              b-axis (y) field component is C4-inconsistent and errors
%                              invz:transverseMF; use 'vector_ab' for genuine in-plane
%                              rotation.
%
%   Returns S.transverse_mf: the resolved MF mode string (echoes opts.solve_opts.transverse_mf,
%   default 'legacy_x').
%
%   Cost is one 1/z solve per field (~15-25 min for a 61-point sweep on a 16^3 grid, single
%   core). Compute S once, then replot freely (invz_plot_spectra_map / invz_run_spectra).

if nargin < 5, opts = struct(); end
hyp      = getf(opts, 'hyp', true);
grid     = getf(opts, 'grid', [16 16 16]);
dpRng    = getf(opts, 'dpRng', 30);
eta      = getf(opts, 'eta', 5e-3);
verbose  = getf(opts, 'verbose', true);
parallel = getf(opts, 'parallel', false);
wmin     = getf(opts, 'peak_wmin', 0);

fdir  = getf(opts, 'field_dir', [1 0 0]);
bztol = getf(opts, 'bz_tol', 1e-9);
sxtra = getf(opts, 'solve_opts', struct());
if any(isfield(sxtra, {'J0eff', 'Jxx0', 'hyp'}))
    error('invz:solveOpts', 'solve_opts fields J0eff/Jxx0/hyp are reserved (driver-owned).');
end
if ~isnumeric(fdir) || ~isreal(fdir) || numel(fdir) ~= 3 || ~all(isfinite(fdir)) || norm(fdir(:)) == 0
    error('invz:fieldDir', 'field_dir must be a nonzero finite real 3-vector.');
end
fdir = reshape(fdir, 1, 3) / norm(fdir);
if any(fields < 0)
    error('invz:fields', 'fields are sweep magnitudes |B| and must be nonnegative.');
end

fields = fields(:).';   w = w(:);
nB = numel(fields);     nw = numel(w);

BvecM = fields(:) * fdir;                        % [nB x 3] actual solve fields
BvecM(abs(BvecM(:, 3)) <= bztol, 3) = 0;         % dead band: identical rule to invz_solve_auto

tmf = getf(sxtra, 'transverse_mf', 'legacy_x');
if strcmp(tmf, 'legacy_x') && any(abs(BvecM(:, 2)) > 0)
    error('invz:transverseMF', ['field has a b-axis (y) component but transverse_mf is ' ...
        '''legacy_x'' (x-only mean field; C4-inconsistent, 17 ueV a/b asymmetry at 4 T). ' ...
        'Set opts.solve_opts.transverse_mf = ''vector_ab'' (or ''none'' for bare diagnostics).']);
end

if isfield(opts, 'Jnu') && isfield(opts, 'info')
    Jnu = opts.Jnu(:);   info = opts.info;
else
    [qc, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', grid, 'range', [-0.5 0.5], 'verbose', false);
    qc = qc(any(abs(qc) > 1e-12, 2), :);
    [Jnu, info] = invz_jq_modes(ion, qc, struct('dpRng', dpRng, 'cache', true));
    Jnu = Jnu(:);
end
Jcc0 = info.Jcc0;
Jaa0   = ion.Jxx0;  if isfield(info, 'Jaa0'),      Jaa0   = info.Jaa0;      end
Jshape = 0;         if isfield(info, 'Jshape_cc'), Jshape = info.Jshape_cc; end

% Sliced plain-array outputs (parfor-safe); packed into S after the sweep.
chizM   = nan(nw, nB);   chirpaM = nan(nw, nB);
Sig0    = nan(1, nB);    phaseC  = zeros(1, nB);

% nWorkers = 0 forces serial execution even inside a parfor, and works without the
% Parallel Computing Toolbox; Inf lets parfor use (and auto-create) the pool.
nWorkers = 0;
if parallel && ~isempty(ver('parallel')), nWorkers = Inf; end

sopts = sxtra;
sopts.hyp = hyp;  sopts.J0eff = Jcc0;  sopts.Jxx0 = Jaa0;  sopts.bz_tol = bztol;

parfor (k = 1:nB, nWorkers)
    [chizM(:, k), chirpaM(:, k), Sig0(k), phaseC(k)] = ...
        one_field(ion, T, BvecM(k, :), Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts, bztol);
    if verbose
        ph = {'masked (no converged solve)', 'moment-form (FM or field-induced)', 'paramagnet'};
        fprintf('  |B| = %5.2f T : %-34s Sigma0 = %s\n', fields(k), ph{phaseC(k)+1}, num2str(Sig0(k)));
    end
end

S = struct();
S.fields = fields;  S.w = w;  S.T = T;  S.info = info;
S.field_dir = fdir;  S.Bvec = BvecM;  S.transverse_mf = tmf;
S.chiz = chizM;  S.chirpa = chirpaM;
S.Sigma0 = Sig0;  S.phase = phaseC;
S.Epeak     = invz_peak_energy(chizM,   w, wmin);
S.Epeak_rpa = invz_peak_energy(chirpaM, w, wmin);
end

% -------------------------------------------------------------------------------------------
function [chiz, chirpa, Sigma0, phase] = one_field(ion, T, B, Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts, bztol)
%ONE_FIELD chi''_cc(omega) at one field via the shared ordered-first solve
% (invz_solve_auto); phase = 1 (moment-form: spontaneous FM or field-induced),
% 2 (strict PM), 0 (no accepted solution -> masked 1/z column).
% Jsel = Jcc0 is the strict-uniform observable, so the demag correction Jshape applies.
nw = numel(w);
chiz = nan(nw, 1);  chirpa = nan(nw, 1);  Sigma0 = NaN;  phase = 0;
copts = struct('Jsel', Jcc0, 'eta', eta, 'Jxx0', Jaa0, 'Jshape', Jshape, 'hyp', hyp);

[pt, phase] = invz_solve_auto(ion, T, B, Jnu, sopts);

if phase == 1                                     % --- moment-form branch (FM or induced) ---
    o  = invz_chi_realaxis(ion, T, B, pt, w, copts);   % reuses pt.si (moment-form eigenstates)
    chiz = imag(o.chi_cc_q(1, :)).';
    pt0 = struct('alpha', 0, 'alpha_m', 0, 'lambda', [0; 0; 0], 'tl', pt.tl, ...
                 'K', [], 'is_ordered', true, 'si', pt.si);
    c0opts = copts;  c0opts.npass = 1;  c0opts.chi0cc_w = o.chi0cc_w;   % share bare cc
    o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
    chirpa = imag(o0.chi_cc_q(1, :)).';
    Sigma0 = pt.Sigma0;
    return;
end

if abs(B(3)) > bztol
    % --- longitudinal failure: NEVER the strict-paramagnet overlay (its m = 0 gate
    % would raise invz:orderedPhase and abort the parfor -- review finding 5). If the
    % failed moment-branch pt still carries valid si/tl, compute the RPA-only overlay
    % from the ordered-style pt0; otherwise leave the whole column masked.
    if ~isempty(pt) && ~isempty(pt.si) && isfield(pt, 'tl') && ~isempty(pt.tl)
        pt0 = struct('alpha', 0, 'alpha_m', 0, 'lambda', [0; 0; 0], 'tl', pt.tl, ...
                     'K', [], 'is_ordered', true, 'si', pt.si);
        c0opts = copts;  c0opts.npass = 1;
        try
            o0 = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
            chirpa = imag(o0.chi_cc_q(1, :)).';
        catch err
            if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
        end
    end
    if ~isempty(pt) && isfield(pt, 'Sigma0'), Sigma0 = pt.Sigma0; end
    return;
end

% --- transverse paramagnetic side: unchanged historical logic --------------------------
if phase == 2 && ~isempty(pt), tl0 = pt.tl;  si0 = pt.si;
else, tl0 = invz_twolevel(ion, T, B, struct('Jxx0', Jaa0));  si0 = []; end
chi0cc = [];
try
    pt0 = struct('alpha', 0, 'lambda', [0; 0], 'tl', tl0, 'K', []);
    c0opts = copts;  c0opts.npass = 1;  c0opts.si = si0;
    o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
    chirpa = imag(o0.chi_cc_q(1, :)).';
    chi0cc = o0.chi0cc_w;                          % share the bare cc with the 1/z call
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
end
if ~isempty(pt) && isfield(pt, 'Sigma0'), Sigma0 = pt.Sigma0; end
if phase == 2                                     % --- converged paramagnetic 1/z ---
    copts1 = copts;
    if ~isempty(chi0cc), copts1.chi0cc_w = chi0cc; end
    o = invz_chi_realaxis(ion, T, B, pt, w, copts1);
    chiz = imag(o.chi_cc_q(1, :)).';
end
end
