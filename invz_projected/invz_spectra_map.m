function S = invz_spectra_map(ion, T, fields, w, opts)
%INVZ_SPECTRA_MAP chi''_cc(omega, B) at the uniform (q=0) mode over a field-magnitude sweep.
%   S = invz_spectra_map(ion, T, fields, w) sweeps the field magnitude |B| over `fields`
%   (Tesla, nonnegative) along opts.field_dir at fixed temperature T (K) and returns
%   chi''_cc(omega, B) across BOTH phases: at each field it tries the ordered 1/z solve
%   first, falling back to the paramagnetic solve, giving a single soft-mode map continuous
%   across the transition.
%
%   Returned maps (columns = fields):
%     S.chiz   [nw x nB]  1/z-renormalized chi''_cc at the 1/z theory's OWN phase: the
%                         strict-PM state wherever it is stable (crit > 0, i.e. above
%                         Bc_1z), the bare-MF moment-form state below Bc_1z (DIAGNOSTIC
%                         there -- the complete ordered 1/z state is Stage 2 of
%                         invzp_QCP_diagnosis.md and does not exist yet)
%     S.chirpa [nw x nB]  bare-RPA chi''_cc          (Sigma = 0, matching phase)
%     With ion.demag ~= 0 both are the demag-corrected MEASURED observable (via
%     info.Jshape_cc, saturating instead of diverging); demag = 0 gives the intrinsic response.
%   Per-field diagnostics: S.phase (1 = moment-form (spontaneous FM below Bc, or field-induced
%   under a longitudinal tilt -- a rounded crossover, no sharp Bc), 2 = strict paramagnet,
%   0 = masked). RPA and 1/z are DIFFERENT theories with different critical fields
%   (invzp_QCP_diagnosis.md). S.phase stays the AUTO dispatch (ordered-first bare-MF); the
%   Sigma = 0 overlay built from it approximates the RPA state only where the ordered 1/z
%   EMT converged to the bare boundary -- an RPA-independent dispatcher is a scheduled
%   follow-up (diagnosis regression 4 stays OPEN), hence the boundary output is named
%   Bc_auto, never Bc_rpa. S.phase_1z is the 1/z leg's own label (2 = stable PM 1/z state
%   used for S.chiz, ALWAYS gated on crit > 0; 1 = bare-MF ordered diagnostic below
%   Bc_1z; 0 = no valid 1/z label: masked column OR a spurious converged PM point with
%   crit <= 0; mirrors S.phase under a longitudinal tilt -- no strict-PM phase there).
%   S.suspect flags the spurious PM columns (S.phase = 2 with finite crit <= 0). Masking
%   before peak extraction: chiz is NaN where phase_1z = 0; chirpa is NaN where there is
%   NO accepted auto state (S.phase = 0 -- the legacy Sigma-zero fallback is computed but
%   no longer surfaced), S.suspect, or a non-finite PM mass (valid_auto_pm requires
%   S.phase = 2 AND finite crit_pm > 0, so a phase-2 column with non-finite crit_pm is
%   masked the same as suspect). S.crit_pm is the per-field PM 1/z mass
%   1 + Sigma0 - J0eff*chi0cc0 (NaN
%   where no PM solve returned a mass -- a failed auto column may still record a finite
%   value here, diagnostic only); S.Bc_auto (anchored only on valid, non-suspect PM
%   points) / S.Bc_1z are sweep-midpoint boundary
%   estimates (NaN when the sweep does not bracket the flip; precision = half the gap
%   between the bracketing labelled fields, so masked or suspect columns between them
%   widen it -- refine with invz_critical for a solver-grade Bc). S.Sigma0 is the Sigma0
%   of the state used for
%   S.chiz. S.Epeak/S.Epeak_rpa (censored, parabolic-refined peak energy;
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
%                              Validated envelope: ac-plane directions [cos(theta) 0 sin(theta)] with theta_c <= 5 deg (see docs/SESSION-2026-07-16-field-angle.md); By ~= 0 now errors under the legacy x-only transverse MF (invz:transverseMF) and is validated under vector_ab (peak observables, all tested in-plane angles; see docs/SESSION-2026-07-16-inplane-rotation.md); demag ~= 0 with tilt is unvalidated.
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
invz_check_solve_opts(sxtra);
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

tmf = invz_check_transverse_mf(sxtra, BvecM(:, 2));

if isfield(opts, 'Jnu') && isfield(opts, 'info')
    Jnu = opts.Jnu(:);   info = opts.info;
else
    [Jnu, info, ~] = invz_bz_couplings(ion, struct('grid', grid, 'dpRng', dpRng));
end
Jcc0 = info.Jcc0;
Jaa0   = ion.Jxx0;  if isfield(info, 'Jaa0'),      Jaa0   = info.Jaa0;      end
Jshape = 0;         if isfield(info, 'Jshape_cc'), Jshape = info.Jshape_cc; end

% Sliced plain-array outputs (parfor-safe); packed into S after the sweep.
chizM   = nan(nw, nB);   chirpaM = nan(nw, nB);
Sig0    = nan(1, nB);    phaseC  = zeros(1, nB);
ph1z    = zeros(1, nB);  critPM  = nan(1, nB);

% nWorkers = 0 forces serial execution even inside a parfor, and works without the
% Parallel Computing Toolbox; Inf lets parfor use (and auto-create) the pool.
nWorkers = 0;
if parallel && ~isempty(ver('parallel')), nWorkers = Inf; end

sopts = sxtra;
sopts.hyp = hyp;  sopts.J0eff = Jcc0;  sopts.Jxx0 = Jaa0;  sopts.bz_tol = bztol;

parfor (k = 1:nB, nWorkers)
    [chizM(:, k), chirpaM(:, k), Sig0(k), phaseC(k), ph1z(k), critPM(k)] = ...
        one_field(ion, T, BvecM(k, :), Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts, bztol);
    if verbose
        ph = {'masked', 'moment-form', 'paramagnet'};
        fprintf('  |B| = %5.2f T : auto-state %-11s 1z-state %-11s Sigma0 = %s\n', ...
                fields(k), ph{phaseC(k)+1}, ph{ph1z(k)+1}, num2str(Sig0(k)));
    end
end

S = struct();
S.fields = fields;  S.w = w;  S.T = T;  S.info = info;
S.field_dir = fdir;  S.Bvec = BvecM;  S.transverse_mf = tmf;
S.Sigma0 = Sig0;  S.phase = phaseC;   % auto (ordered-first bare-MF) dispatch -- see header
S.phase_1z = ph1z;  S.crit_pm = critPM;
S.suspect = (phaseC == 2) & isfinite(critPM) & (critPM <= 0);   % spurious below-Bc auto-PM
% Validity masks BEFORE the spectra are packed and the peaks extracted (re-review
% findings 2, R3-1; R4-3 hardening): chiz is masked where there is no valid 1/z label;
% the auto overlay is masked wherever there is NO accepted auto state (phase 0 --
% including the legacy Sigma-zero fallbacks, still computed inside one_field but no
% longer surfaced) or the PM mass is not confirmed positive (suspect OR non-finite).
% The plotter then hides them (NaN -> transparent) and invz_peak_energy never sees them.
valid_auto_pm = (phaseC == 2) & isfinite(critPM) & (critPM > 0);
invalid_auto  = (phaseC == 0) | ((phaseC == 2) & ~valid_auto_pm);
chizM(:, ph1z == 0)      = NaN;
chirpaM(:, invalid_auto) = NaN;
S.chiz = chizM;  S.chirpa = chirpaM;
S.Bc_auto = invz_boundary_field(fields, phaseC == 1, valid_auto_pm);   % suspect never anchors
S.Bc_1z   = invz_boundary_field(fields, ph1z  == 1, ph1z  == 2);
S.Epeak     = invz_peak_energy(chizM,   w, wmin);
S.Epeak_rpa = invz_peak_energy(chirpaM, w, wmin);
end

% -------------------------------------------------------------------------------------------
function [chiz, chirpa, Sigma0, phase, phase_1z, crit_pm] = one_field(ion, T, B, Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts, bztol)
%ONE_FIELD chi''_cc(omega) at one field -- TWO independently phased legs (QCP Stage 1,
% docs/superpowers/plans/2026-07-21-invzp-qcp-stage1-split-overlays.md):
%   auto/overlay leg (chirpa): the ordered-first invz_solve_auto selection (bare-MF
%     moment), UNCHANGED -- phase = 1 moment-form (FM or field-induced), 2 strict PM,
%     0 masked, historical semantics. The Sigma = 0 overlay built from it approximates
%     the RPA state only where the ordered EMT converged to the bare boundary; a
%     1/z-independent RPA dispatcher is a scheduled follow-up (diagnosis regression 4
%     stays open), hence the driver reports Bc_auto, never Bc_rpa.
%   1/z leg (chiz): the 1/z theory's own stability. Wherever the strict-PM solve converges
%     with crit = 1 + Sigma0 - J0eff*chi0cc0 > 0, the PM state is the consistent 1/z state
%     (phase_1z = 2) even though the bare-MF leg still orders; below Bc_1z the bare-MF
%     ordered 1/z curve is kept as a DIAGNOSTIC (phase_1z = 1) until the full ordered 1/z
%     state exists (diagnosis Stage 2). EVERY phase_1z = 2 label is gated on crit > 0;
%     phase_1z = 0 means no valid 1/z label (masked column OR a spurious converged PM
%     point with crit <= 0). Longitudinal tilts have no strict-PM phase
%     (rounded crossover): phase_1z mirrors phase there.
% crit_pm: PM 1/z mass at this field (NaN where no PM solve was attempted or returned).
% Jsel = Jcc0 is the strict-uniform observable, so the demag correction Jshape applies.
nw = numel(w);
chiz = nan(nw, 1);  chirpa = nan(nw, 1);  Sigma0 = NaN;  phase = 0;
phase_1z = 0;  crit_pm = NaN;
tmf = getf(sopts, 'transverse_mf', 'legacy_x');
% transverse_mf here so any fallback single-ion rebuild inside invz_chi_realaxis (when a
% branch below can't supply si) uses the SAME MF model as the solve, not the 'legacy_x'
% default (review finding I1; the C1 companion bug in invz_spectra_qpath.m).
copts = struct('Jsel', Jcc0, 'eta', eta, 'Jxx0', Jaa0, 'Jshape', Jshape, 'hyp', hyp, ...
               'transverse_mf', tmf);

[pt, phase] = invz_solve_auto(ion, T, B, Jnu, sopts);

if phase == 1                                     % --- moment-form branch (FM or induced) ---
    % 1/z leg: strict-PM stability probe, TRANSVERSE ONLY (the strict-PM solver's m = 0
    % gate would raise invz:orderedPhase under a longitudinal tilt -- review finding 5).
    ptp = [];
    if B(3) == 0
        try
            ptp = invz_solve_point(ion, T, B, Jnu, sopts);
        catch err
            if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
        end
    end
    if ~isempty(ptp) && isfield(ptp, 'crit') && isfinite(ptp.crit), crit_pm = ptp.crit; end
    % crit_pm doubles as the field-guarded crit: pm_valid never touches ptp.crit directly
    % (an early-return ptp may lack the field).
    pm_valid = ~isempty(ptp) && ptp.converged && isfinite(ptp.Sigma0) && ...
               isfinite(crit_pm) && crit_pm > 0;

    if pm_valid                         % --- Bc_1z < |B| < bare boundary: 1/z theory is PM ---
        phase_1z = 2;  Sigma0 = ptp.Sigma0;
        copts1 = copts;  copts1.si = ptp.si;      % PM chi0 differs from the ordered leg's:
        o1 = invz_chi_realaxis(ion, T, B, ptp, w, copts1);   % no chi0cc_w sharing possible
        chiz = imag(o1.chi_cc_q(1, :)).';
        c0opts = copts;  c0opts.npass = 1;        % RPA overlay stays on the bare-MF state
        o0 = invz_chi_realaxis(ion, T, B, invz_zero_sigma_overlay(pt), w, c0opts);
        chirpa = imag(o0.chi_cc_q(1, :)).';
    else                                     % --- below Bc_1z: ordered 1/z, DIAGNOSTIC ---
        phase_1z = 1;  Sigma0 = pt.Sigma0;
        o  = invz_chi_realaxis(ion, T, B, pt, w, copts);   % reuses pt.si (moment-form eigenstates)
        chiz = imag(o.chi_cc_q(1, :)).';
        pt0 = invz_zero_sigma_overlay(pt);
        c0opts = copts;  c0opts.npass = 1;  c0opts.chi0cc_w = o.chi0cc_w;   % share bare cc
        o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
        chirpa = imag(o0.chi_cc_q(1, :)).';
    end
    return;
end

if abs(B(3)) > bztol
    % --- longitudinal failure: NEVER the strict-paramagnet overlay (its m = 0 gate
    % would raise invz:orderedPhase and abort the parfor -- review finding 5). If the
    % failed moment-branch pt still carries valid si/tl, compute the RPA-only overlay
    % from the ordered-style pt0; otherwise leave the whole column masked.
    if ~isempty(pt) && ~isempty(pt.si) && isfield(pt, 'tl') && ~isempty(pt.tl)
        pt0 = invz_zero_sigma_overlay(pt);
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
else, tl0 = invz_twolevel(ion, T, B, struct('Jxx0', Jaa0, 'transverse_mf', tmf));  si0 = []; end
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
if ~isempty(pt) && isfield(pt, 'crit') && isfinite(pt.crit), crit_pm = pt.crit; end
if phase == 2                                     % --- converged paramagnetic 1/z ---
    if isfinite(crit_pm) && crit_pm > 0
        phase_1z = 2;                             % the auto PM state IS the 1/z leg here
    end                                           % crit <= 0: spurious below-Bc PM point
                                                  % (documented tensor-side failure mode) --
                                                  % phase_1z stays 0, column flagged suspect
    copts1 = copts;  copts1.si = pt.si;           % reuse the solve's si (C1 companion: pt is
    if ~isempty(chi0cc), copts1.chi0cc_w = chi0cc; end   % not is_ordered, so without si this
    o = invz_chi_realaxis(ion, T, B, pt, w, copts1);     % would silently rebuild at 'legacy_x'
    chiz = imag(o.chi_cc_q(1, :)).';
end
end
