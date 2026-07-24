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
%                         Bc_1z), and below Bc_1z -- under opts.ordered_1z = 'jensen'
%                         (default) -- the Jensen-consistent ordered 1/z state (J 2.28-2.33,
%                         diagnosis Stage 2 delivered on the projected path), or under the
%                         'bare' escape hatch the retired Stage-1 bare-MF diagnostic
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
%   used for S.chiz, ALWAYS gated on crit > 0; 1 = under 'jensen' the CONSISTENT ordered
%   1/z state (root-existence gated, not a diagnostic) -- under 'bare' the retired
%   Stage-1 bare-MF diagnostic; 0 = no valid 1/z label: masked column (including a jensen
%   solve whose H_MF root does not exist, or whose static closure failed) OR a spurious
%   converged PM point with crit <= 0 -- jensen failure NEVER silently falls back to
%   'bare' at a field where jensen actually applies; mirrors S.phase under a
%   longitudinal tilt -- no strict-PM phase there, and jensen (transverse/spontaneous
%   only by construction) is never attempted, so a tilted field-induced moment always
%   uses the bare diagnostic below Bc_1z regardless of opts.ordered_1z).
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
%   S.chiz. S.m_1z [1 x nB] is the ordered moment of the jensen 1/z state (NaN under
%   'bare' and on all non-ordered columns). S.D_ord [1 x nB] is the ordered static
%   inverse response 1 + (J0eff - K(1))*G(1) AT THE FINAL jensen STATE -- the pole
%   observable, which vanishes at Bc_1z from below so the FM and PM 1/z branches close
%   at the same field (NaN under 'bare' and on all non-ordered columns). S.Epeak/S.Epeak_rpa (censored, parabolic-refined peak energy;
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
%                              lattice sum; used by the tests). Must be supplied TOGETHER
%                              (error invz:spectraPrecomputedPartial if only one is present).
%     .dipole, .ewald          opt-in Ewald dipolar backend, forwarded BY PRESENCE into
%                              invz_bz_couplings when the couplings are computed (absent =>
%                              unchanged brute-force default). .dipole is 'bruteforce' |
%                              'ewald'; .ewald is a complete {alpha,r_cut,g_cut,boundary}
%                              struct, required only with .dipole = 'ewald'
%                              (docs/invzp_ewald_integration_map.md Sec.5A).
%     .gridConvention, .gridOffset, .gammaPolicy   opt-in BZ quadrature policy, forwarded BY
%                              PRESENCE into invz_bz_couplings when computed (any one present
%                              there switches to the invz_phase1_qgrid route; all absent here
%                              leaves the legacy grid unchanged and info.grid unset).
%     .cache     (true)        invz_jq_modes file cache; forwarded BY PRESENCE when computed.
%     .verbose   (true)        print one progress line per field
%     .field_dir ([1 0 0])     nonzero finite real 3-vector, normalized internally; sets the
%                              sweep direction of `fields` (|B|). Error invz:fieldDir if
%                              invalid. A nonzero z-component routes through the longitudinal
%                              (field-induced moment) solve once |Bz| clears .bz_tol.
%                              Validated envelope: ac-plane directions [cos(theta) 0 sin(theta)] with theta_c <= 5 deg (see docs/SESSION-2026-07-16-field-angle.md); By ~= 0 now errors under the legacy x-only transverse MF (invz:transverseMF) and is validated under vector_ab (peak observables, all tested in-plane angles; see docs/SESSION-2026-07-16-inplane-rotation.md); demag ~= 0 with tilt is unvalidated.
%     .bz_tol    (1e-9)        T; dead band on Bz -- resolved ONCE, applied to the field table
%                              BEFORE any solve, and forwarded to invz_solve_auto/one_field.
%     .ordered_1z ('jensen')   the 1/z leg's ordered-side solve below Bc_1z: 'jensen' (default)
%                              is the Stage-2 consistent ordered 1/z state (invz_solve_point_ordered
%                              opts.ordered_mode = 'jensen'); 'bare' reproduces the retired
%                              Stage-1 bare-MF diagnostic verbatim. Error invz:ordered1z if
%                              neither. This is driver-owned: solve_opts.ordered_mode is
%                              reserved and errors invz:solveOpts (P1-6) so the auto/overlay
%                              leg can never be flipped by a caller.
%     .solve_opts (struct())   merged into the per-field invz_solve_auto opts; fields
%                              J0eff/Jxx0/hyp are reserved (driver-owned) -> error invz:solveOpts.
%                              transverse_mf ('legacy_x' | 'none' | 'vector_ab', default
%                              'legacy_x') is a legal solve_opts field forwarded to the
%                              solvers. Under 'legacy_x' (x-only mean field) a nonzero
%                              b-axis (y) field component is C4-inconsistent and errors
%                              invz:transverseMF; use 'vector_ab' for genuine in-plane
%                              rotation.
%
%   Precomputed .Jnu/.info + backend/grid strictness (Ewald Step-5 Task 7): with NO explicit
%   .dipole/.ewald/.gridConvention/.gridOffset/.gammaPolicy request, a complete precomputed
%   info.dipole is simply carried through in S.info, and a provenance-less legacy synthetic
%   .info (e.g. only .Jcc0, as built by existing tests) is accepted unchanged and resolves to
%   brute force -- unaffected backward-compatible behavior. With an explicit request, the
%   precomputed .info must carry the matching COMPLETE provenance or the call errors:
%   invz:spectraBackendProvenanceMissing / invz:spectraGridProvenanceMissing (provenance
%   absent), invz:spectraBackendConflict / invz:spectraGridConflict (provenance present but
%   disagrees with the request). An invalid explicit .dipole/.ewald is rejected with the same
%   invz:jqModes* identifiers invz_jq_modes itself raises, even when precomputed .Jnu/.info
%   bypass the lattice sum entirely.
%
%   Returns S.transverse_mf: the resolved MF mode string (echoes opts.solve_opts.transverse_mf,
%   default 'legacy_x'). Returns S.ordered_1z: echoes opts.ordered_1z ('jensen' or 'bare').
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
o1z = getf(opts, 'ordered_1z', 'jensen');
if ~any(strcmp(o1z, {'jensen', 'bare'}))
    error('invz:ordered1z', 'ordered_1z must be ''jensen'' or ''bare''.');
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

tmf = invz_check_transverse_mf(sxtra, BvecM(:, 2));

% Ewald Step-5 Task 7 (docs/invzp_ewald_integration_map.md Sec.5A): backend/grid options are
% forwarded BY PRESENCE into invz_bz_couplings on the compute branch, and a precomputed
% opts.Jnu/opts.info pair is validated against any EXPLICIT backend/grid-policy request rather
% than trusted blindly -- see the local_* helpers below for the exact conflict rules.
hasBackendReq = isfield(opts, 'dipole') || isfield(opts, 'ewald');
hasGridReq    = isfield(opts, 'gridConvention') || isfield(opts, 'gridOffset') || isfield(opts, 'gammaPolicy');
if hasBackendReq
    % Validate opts.dipole/opts.ewald THEMSELVES even though the precomputed branch below may
    % never call invz_jq_modes: a malformed request must not escape checking just because no
    % lattice sum runs.
    [backendReq, eoptsReq] = local_validate_dipole_opts(opts);
end

hasJnuOpt  = isfield(opts, 'Jnu');
hasInfoOpt = isfield(opts, 'info');
if hasJnuOpt ~= hasInfoOpt
    error('invz:spectraPrecomputedPartial', ...
          'opts.Jnu and opts.info must be supplied together (exactly one was present).');
end

if hasJnuOpt && hasInfoOpt
    Jnu = opts.Jnu(:);   info = opts.info;
    if hasBackendReq, local_check_backend_provenance(info, backendReq, eoptsReq); end
    if hasGridReq,    local_check_grid_provenance(info, opts, grid);              end
else
    bzOpts = struct('grid', grid, 'dpRng', dpRng);
    if isfield(opts, 'dipole'),         bzOpts.dipole         = opts.dipole;         end
    if isfield(opts, 'ewald'),          bzOpts.ewald          = opts.ewald;          end
    if isfield(opts, 'gridConvention'), bzOpts.gridConvention = opts.gridConvention; end
    if isfield(opts, 'gridOffset'),     bzOpts.gridOffset     = opts.gridOffset;     end
    if isfield(opts, 'gammaPolicy'),    bzOpts.gammaPolicy    = opts.gammaPolicy;    end
    if isfield(opts, 'cache'),          bzOpts.cache          = opts.cache;          end
    [Jnu, info, ~] = invz_bz_couplings(ion, bzOpts);
end
Jcc0 = info.Jcc0;
Jaa0   = ion.Jxx0;  if isfield(info, 'Jaa0'),      Jaa0   = info.Jaa0;      end
Jshape = 0;         if isfield(info, 'Jshape_cc'), Jshape = info.Jshape_cc; end

% Sliced plain-array outputs (parfor-safe); packed into S after the sweep.
chizM   = nan(nw, nB);   chirpaM = nan(nw, nB);
Sig0    = nan(1, nB);    phaseC  = zeros(1, nB);
ph1z    = zeros(1, nB);  critPM  = nan(1, nB);
m1zM    = nan(1, nB);    DordM   = nan(1, nB);

% nWorkers = 0 forces serial execution even inside a parfor, and works without the
% Parallel Computing Toolbox; Inf lets parfor use (and auto-create) the pool.
nWorkers = 0;
if parallel && ~isempty(ver('parallel')), nWorkers = Inf; end

sopts = sxtra;
sopts.hyp = hyp;  sopts.J0eff = Jcc0;  sopts.Jxx0 = Jaa0;  sopts.bz_tol = bztol;

parfor (k = 1:nB, nWorkers)
    [chizM(:, k), chirpaM(:, k), Sig0(k), phaseC(k), ph1z(k), critPM(k), m1zM(k), DordM(k)] = ...
        one_field(ion, T, BvecM(k, :), Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts, bztol, o1z);
    if verbose
        ph = {'masked', 'moment-form', 'paramagnet'};
        fprintf('  |B| = %5.2f T : auto-state %-11s 1z-state %-11s Sigma0 = %s\n', ...
                fields(k), ph{phaseC(k)+1}, ph{ph1z(k)+1}, num2str(Sig0(k)));
    end
end

S = struct();
S.fields = fields;  S.w = w;  S.T = T;  S.info = info;
S.field_dir = fdir;  S.Bvec = BvecM;  S.transverse_mf = tmf;  S.ordered_1z = o1z;
S.Sigma0 = Sig0;  S.phase = phaseC;   % auto (ordered-first bare-MF) dispatch -- see header
S.phase_1z = ph1z;  S.crit_pm = critPM;  S.m_1z = m1zM;  S.D_ord = DordM;
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
function [chiz, chirpa, Sigma0, phase, phase_1z, crit_pm, m_1z, D_ord] = one_field(ion, T, B, Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts, bztol, o1z)
%ONE_FIELD chi''_cc(omega) at one field -- TWO independently phased legs (QCP Stage 1,
% docs/superpowers/plans/2026-07-21-invzp-qcp-stage1-split-overlays.md; ordered 1/z leg
% below Bc_1z upgraded to the Stage-2 jensen solve, docs/superpowers/plans/
% 2026-07-22-invzp-stage2-ordered-thermodynamics.md):
%   auto/overlay leg (chirpa): the ordered-first invz_solve_auto selection (bare-MF
%     moment), UNCHANGED -- phase = 1 moment-form (FM or field-induced), 2 strict PM,
%     0 masked, historical semantics. The Sigma = 0 overlay built from it approximates
%     the RPA state only where the ordered EMT converged to the bare boundary; a
%     1/z-independent RPA dispatcher is a scheduled follow-up (diagnosis regression 4
%     stays open), hence the driver reports Bc_auto, never Bc_rpa.
%   1/z leg (chiz): the 1/z theory's own stability. Wherever the strict-PM solve converges
%     with crit = 1 + Sigma0 - J0eff*chi0cc0 > 0, the PM state is the consistent 1/z state
%     (phase_1z = 2) even though the bare-MF leg still orders; below Bc_1z, TRANSVERSE ONLY
%     (B(3) == 0) -- under o1z = 'jensen' (default) -- the CONSISTENT ordered 1/z state is
%     solved (invz_solve_point_ordered opts.ordered_mode = 'jensen', J 2.28-2.33; diagnosis
%     Stage 2 delivered), gated on H_MF root existence and static-closure convergence,
%     NEVER a silent fallback to the bare diagnostic on failure (phase_1z stays 0 and the
%     column is masked instead); under o1z = 'bare' the retired Stage-1 bare-MF diagnostic
%     is kept verbatim. EVERY phase_1z = 2 label is gated on crit > 0; phase_1z = 0 means
%     no valid 1/z label (masked column OR a spurious converged PM point with crit <= 0).
%     Longitudinal tilts have no strict-PM phase (rounded crossover): phase_1z mirrors
%     phase there, and ALWAYS uses the bare diagnostic below Bc_1z regardless of o1z --
%     jensen is transverse/spontaneous only by invz_solve_point_ordered's own scope
%     contract (invz:orderedMode on any |B(3)| > bz_tol), so a tilted field-induced
%     moment is never even offered to it (this mirrors the strict-PM probe above, which
%     likewise only runs at B(3) == 0).
% crit_pm: PM 1/z mass at this field (NaN where no PM solve was attempted or returned).
% m_1z/D_ord: jensen ordered moment / ordered static inverse response (pole observable) at
%     this field -- NaN except on a converged jensen-ordered column (phase_1z = 1, o1z =
%     'jensen'); always NaN under o1z = 'bare'.
% Jsel = Jcc0 is the strict-uniform observable, so the demag correction Jshape applies.
nw = numel(w);
chiz = nan(nw, 1);  chirpa = nan(nw, 1);  Sigma0 = NaN;  phase = 0;
phase_1z = 0;  crit_pm = NaN;  m_1z = NaN;  D_ord = NaN;
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
    else                                     % --- below Bc_1z: ordered 1/z leg ---
        % jensen is TRANSVERSE/SPONTANEOUS only (invz_solve_point_ordered's own scope
        % guard, invz:orderedMode, on any |B(3)| > bz_tol): a longitudinal tilt's
        % field-induced moment falls through to the bare diagnostic below, matching
        % the auto/overlay leg's own B(3) == 0 gate on the strict-PM probe above --
        % NOT a "jensen failed" masked column, since jensen is never attempted there.
        if strcmp(o1z, 'jensen') && B(3) == 0   % Stage 2: consistent ordered state (H_MF)
            so2 = sopts;  so2.ordered_mode = 'jensen';
            ptj = [];
            try
                ptj = invz_solve_point_ordered(ion, T, B, Jnu, so2);
            catch err
                if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
            end
            if ~isempty(ptj) && ptj.is_ordered && ptj.converged && isfinite(ptj.Sigma0)
                phase_1z = 1;  Sigma0 = ptj.Sigma0;  m_1z = ptj.m0;  D_ord = ptj.D_uni;
                oj = invz_chi_realaxis(ion, T, B, ptj, w, copts);   % jensen si differs from
                chiz = imag(oj.chi_cc_q(1, :)).';                   % the auto pt's -- no sharing
            end                              % else phase_1z stays 0 -> chiz column masked
            pt0 = invz_zero_sigma_overlay(pt);                      % overlay: UNCHANGED auto state
            c0opts = copts;  c0opts.npass = 1;
            o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
            chirpa = imag(o0.chi_cc_q(1, :)).';
        else                                 % 'bare', or a longitudinal tilt under 'jensen':
                                              % Stage-1 diagnostic escape hatch
            phase_1z = 1;  Sigma0 = pt.Sigma0;
            o  = invz_chi_realaxis(ion, T, B, pt, w, copts);
            chiz = imag(o.chi_cc_q(1, :)).';
            pt0 = invz_zero_sigma_overlay(pt);
            c0opts = copts;  c0opts.npass = 1;  c0opts.chi0cc_w = o.chi0cc_w;
            o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
            chirpa = imag(o0.chi_cc_q(1, :)).';
        end
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

% -------------------------------------------------------------------------------------------
% Ewald Step-5 Task 7: precomputed-coupling provenance/conflict validation (shared logic,
% duplicated byte-for-byte in invz_spectra_map.m and invz_spectra_qpath.m -- no new shared file
% is introduced, per the Task-7 staging contract; docs/invzp_ewald_integration_map.md Sec.5A).
% local_validate_dipole_opts intentionally replicates invz_jq_modes.m's own
% local_resolve_dipole_backend checks AND error identifiers, so a malformed opts.dipole/
% opts.ewald request is rejected identically whether couplings are computed for real or a
% precomputed opts.Jnu/opts.info bypasses invz_jq_modes entirely.
% -------------------------------------------------------------------------------------------
function [backend, eopts] = local_validate_dipole_opts(opts)
if ~isfield(opts, 'dipole') || isempty(opts.dipole)
    backend = 'bruteforce';
else
    raw = opts.dipole;
    if isstring(raw) && isscalar(raw), raw = char(raw); end
    if ~(ischar(raw) && isrow(raw))
        error('invz:jqModesBackend', ...
            ['opts.dipole must be a scalar string/char naming a backend ' ...
             '(''bruteforce''|''ewald''); got class %s.'], class(opts.dipole));
    end
    if ~(strcmp(raw, 'bruteforce') || strcmp(raw, 'ewald'))
        error('invz:jqModesBackend', ...
            'unknown opts.dipole backend ''%s''; supported backends are ''bruteforce'' and ''ewald''.', raw);
    end
    backend = raw;
end

hasEwaldOpts = isfield(opts, 'ewald') && ~isempty(opts.ewald);
if hasEwaldOpts && ~strcmp(backend, 'ewald')
    error('invz:jqModesEwaldOptsUnexpected', ...
        ['opts.ewald was supplied but the resolved opts.dipole backend is ''%s'', not ''ewald''; ' ...
         'Ewald controls are only accepted with an explicit opts.dipole=''ewald'' request.'], backend);
end

eopts = [];
if strcmp(backend, 'ewald')
    if ~hasEwaldOpts || ~isstruct(opts.ewald) || ~isscalar(opts.ewald)
        error('invz:jqModesEwaldOptsFields', ...
            ['opts.dipole=''ewald'' requires a complete scalar struct opts.ewald with EXACTLY the ' ...
             'fields {alpha, r_cut, g_cut, boundary}; jq_modes does not synthesize frozen defaults.']);
    end
    want = sort({'alpha', 'r_cut', 'g_cut', 'boundary'});
    have = sort(reshape(fieldnames(opts.ewald), 1, []));
    if ~isequal(have, want)
        error('invz:jqModesEwaldOptsFields', ...
            'opts.ewald must have EXACTLY the fields {alpha, r_cut, g_cut, boundary}; got {%s}.', ...
            strjoin(reshape(fieldnames(opts.ewald), 1, []), ', '));
    end
    b = opts.ewald.boundary;
    if isstring(b) && isscalar(b), b = char(b); end
    if ~(ischar(b) && isrow(b) && strcmp(b, 'conducting_k0_omitted'))
        error('invz:jqModesEwaldBoundary', ...
            'opts.ewald.boundary must be ''conducting_k0_omitted''; got %s.', class(opts.ewald.boundary));
    end
    eopts = opts.ewald;
end
end

function tf = local_has_complete_dipole_provenance(info)
%LOCAL_HAS_COMPLETE_DIPOLE_PROVENANCE True iff info.dipole is a complete backend-provenance
% struct (invz_jq_modes' info.dipole contract: backend/ewald/q_reduction/primitive_schema).
tf = isfield(info, 'dipole') && isstruct(info.dipole) && isscalar(info.dipole) && ...
     all(isfield(info.dipole, {'backend', 'ewald', 'q_reduction', 'primitive_schema'}));
end

function tf = local_has_complete_grid_provenance(info)
%LOCAL_HAS_COMPLETE_GRID_PROVENANCE True iff info.grid is a complete BZ grid-policy provenance
% struct (invz_bz_couplings' info.grid contract).
tf = isfield(info, 'grid') && isstruct(info.grid) && isscalar(info.grid) && ...
     all(isfield(info.grid, {'schema', 'convention', 'offset', 'gammaPolicy', 'requested', ...
                             'nominal', 'retained', 'n_gamma', 'qhash'}));
end

function local_check_backend_provenance(info, backendReq, eoptsReq)
%LOCAL_CHECK_BACKEND_PROVENANCE An explicit opts.dipole/opts.ewald request was made alongside
% precomputed opts.Jnu/opts.info: require complete info.dipole and compare EXACTLY.
if ~local_has_complete_dipole_provenance(info)
    error('invz:spectraBackendProvenanceMissing', ['an explicit opts.dipole/opts.ewald request ' ...
        'was given together with precomputed opts.Jnu/opts.info, but opts.info lacks complete ' ...
        'info.dipole provenance (a struct with backend/ewald/q_reduction/primitive_schema); cannot ' ...
        'verify the precomputed couplings actually used the requested backend.']);
end
conflict = ~strcmp(info.dipole.backend, backendReq);
if ~conflict && strcmp(backendReq, 'ewald')
    conflict = ~isequaln(info.dipole.ewald, eoptsReq);
end
if conflict
    error('invz:spectraBackendConflict', ['the requested opts.dipole/opts.ewald backend/controls ' ...
        'conflict with the precomputed opts.info.dipole provenance (requested ''%s'', precomputed ' ...
        'info carries ''%s'').'], backendReq, info.dipole.backend);
end
end

function local_check_grid_provenance(info, opts, grid)
%LOCAL_CHECK_GRID_PROVENANCE An explicit gridConvention/gridOffset/gammaPolicy request was made
% alongside precomputed opts.Jnu/opts.info: require complete info.grid, resolve any omitted
% request member to invz_bz_couplings' own defaults, and compare EXACTLY (convention, offset,
% gammaPolicy, and the requested grid dimensions).
if ~local_has_complete_grid_provenance(info)
    error('invz:spectraGridProvenanceMissing', ['an explicit gridConvention/gridOffset/gammaPolicy ' ...
        'request was given together with precomputed opts.Jnu/opts.info, but opts.info lacks ' ...
        'complete info.grid provenance; cannot verify the precomputed couplings actually used the ' ...
        'requested grid policy.']);
end
convention = getf(opts, 'gridConvention', 'legacy_inclusive');
if isstring(convention) && isscalar(convention), convention = char(convention); end
gridOffset  = logical(reshape(getf(opts, 'gridOffset', [0 0 0]), 1, 3));
gammaPolicy = getf(opts, 'gammaPolicy', 'P_drop');
if isstring(gammaPolicy) && isscalar(gammaPolicy), gammaPolicy = char(gammaPolicy); end
gr = info.grid;
conflict = ~isequal(gr.convention, convention) || ~isequal(gr.offset, gridOffset) || ...
    ~isequal(gr.gammaPolicy, gammaPolicy) || ~isequal(gr.requested, reshape(grid, 1, 3));
if conflict
    error('invz:spectraGridConflict', ['the requested gridConvention/gridOffset/gammaPolicy/grid ' ...
        'dimensions conflict with the precomputed opts.info.grid provenance.']);
end
end
