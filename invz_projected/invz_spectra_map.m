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
%   (invzp_convg_diagnosis.md, especially Sections 2.4 and 8). S.phase stays the AUTO dispatch
%   (ordered-first bare-MF); the
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
%   S.chiz. S.m_1z [1 x nB] is the ordered moment of the jensen 1/z state and S.D_ord [1 x nB]
%   the ordered static inverse response 1 + (J0eff - K(1))*G(1) AT THE FINAL jensen STATE -- the
%   pole observable, which vanishes at Bc_1z from below so the FM and PM 1/z branches close
%   at the same field. BOTH ARE ALWAYS NaN UNDER 'bare'; otherwise they are recorded wherever
%   the jensen leg RETURNED A CLOSURE-CONSISTENT POINT, which is NOT the same column set as
%   phase_1z = 1 (review M3, the sibling F8 missed). Under a STRICT scheme they are also finite
%   on a column MASKED as 'pm_probe_unknown', where the jensen solve ran for DIAGNOSTICS only
%   (review F2), and they are RETAINED on a 'response_failed' column, where S.Sigma0 is reset to
%   NaN but these two are not: unlike Sigma0 they are not defined as properties of the state
%   used for S.chiz, so nulling them would discard a genuine measurement rather than enforce a
%   definition. Gate on S.phase_1z == 1 if you want labelled ordered columns only.
%   S.Epeak/S.Epeak_rpa (censored, parabolic-refined peak energy;
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
%     .field_continuation ('none')  'none' preserves independent per-field solves.
%                              'qcp_down' is an EXPERIMENTAL preview: solve the supplied
%                              increasing field table in descending execution order and pass
%                              the last accepted Jensen predictor's full Sigma/K0 state to the
%                              next lower field. This forces serial execution, retains the
%                              ordinary cold retry and every residual/stability gate, and does
%                              not certify branch identity. S records seed provenance.
%     .peak_wmin (0)           meV; excludes a known low-frequency line (e.g. hyperfine
%                              pole) from the peak search; default is no exclusion
%     .Jnu, .info              precomputed coupling branches / info struct (skips the
%                              lattice sum; used by the tests). Must be supplied TOGETHER
%                              (error invz:spectraPrecomputedPartial if only one is present).
%     .dipole, .ewald          opt-in Ewald dipolar backend, forwarded BY PRESENCE into
%                              invz_bz_couplings when the couplings are computed (absent =>
%                              unchanged brute-force default). .dipole is 'bruteforce' |
%                              'ewald'; .ewald is a complete {alpha,r_cut,g_cut,boundary}
%                              struct, required only with .dipole = 'ewald'.
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
%                              Validated envelope: ac-plane directions [cos(theta) 0 sin(theta)]
%                              with theta_c <= 5 deg; By ~= 0 now errors under the legacy x-only
%                              transverse MF (invz:transverseMF) and is validated under
%                              vector_ab (peak observables, all tested in-plane angles);
%                              demag ~= 0 with tilt is unvalidated.
%     .bz_tol    (1e-9)        T; dead band on Bz -- resolved ONCE, applied to the field table
%                              BEFORE any solve, and forwarded to invz_solve_auto/one_field.
%     .ordered_1z ('jensen')   the 1/z leg's ordered-side solve below Bc_1z: 'jensen' (default)
%                              is the Stage-2 consistent ordered 1/z state (invz_solve_point_ordered
%                              opts.ordered_mode = 'jensen'); 'bare' reproduces the retired
%                              Stage-1 bare-MF diagnostic verbatim. Error invz:ordered1z if
%                              neither. This is driver-owned: solve_opts.ordered_mode is
%                              reserved and errors invz:solveOpts (P1-6) so the auto/overlay
%                              leg can never be flipped by a caller.
%     .static_medium ('resummed')  omega_n = 0 static-medium scheme (spec SS4.2), resolved
%                              ONCE here by invz_check_static_medium and threaded into every
%                              per-field solve: 'resummed' (DEFAULT, legacy and numerically
%                              bit-identical) | 'strict_1z_dyson_ref' | 'strict_1z_bare_ref'.
%                              .ref_margin (1e-6) is the reference-denominator floor. BOTH are
%                              driver-owned: either one inside solve_opts errors invz:solveOpts.
%                              A STRICT scheme also switches the 1/z leg's phase dispatch to
%                              the three-way rule described under S.phase_1z_reason below;
%                              under 'resummed' the historical dispatch is preserved exactly.
%     .solve_opts (struct())   merged into the per-field invz_solve_auto opts; fields
%                              J0eff/Jxx0/hyp/ordered_mode/static_medium/ref_margin are
%                              reserved (driver-owned) -> error invz:solveOpts.
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
%   PHASE PROVENANCE (spec SS4.4). S.static_medium is the scalar REQUESTED scheme and is
%   always present. Per column: S.static_medium_used (the same string, except
%   'n/a_bare_escape' on a bare escape-hatch column, which is the retired Stage-1 bare-MF
%   diagnostic and NOT a result of the requested scheme); S.pm_probe_status
%   ('not_attempted' | 'nonconverged' | 'nonfinite' | 'boundary_band' | 'stable' |
%   'unstable' | 'recoverable_error'); S.pm_probe_error_id (the EXACT absorbed identifier,
%   never collapsed into 'nonconverged'); S.stability_1z; and S.phase_1z_reason -- a reason
%   for EVERY column, so no phase_1z = 0 is ever unexplained:
%     'pm' 'ordered' 'bare_escape_hatch'                    (labelled columns)
%     'unstable_endpoint' 'medium_out_of_domain' 'degenerate_doublet' 'solver_failed'
%     'pm_probe_unknown' 'boundary_indeterminate' 'not_attempted_longitudinal'
%     'bare_not_ordered' 'response_failed'                  (masked columns)
%   VOCABULARY EXTENSION (review M7). The historical design spec in Git history enumerated
%   phase_1z_reason as a CLOSED list of ELEVEN values; 'response_failed' (added by
%   review F3, so a chiz response failure under a strict scheme masks the column instead of
%   letting an all-NaN spectrum anchor Bc_1z while wearing a 'pm'/'ordered' label) makes TWELVE.
%   A Stage-4 G12/G16 consumer switching on this enum must handle it. It is unreachable today
%   (see S.response_error_id below), and its ratification is with the plan owner.
%   S.stability_1z is NOT a statement about the column's phase, and it does NOT have a single
%   source. On a 'pm' column it is hard-set TRUE by the PM STABILITY GATE that labelled the
%   column -- crit_pm > crit_tol under a strict scheme (invz_pm_verdict), crit_pm > 0 on the
%   'resummed' path -- with NO ordered solve involved at all, so it is true on EVERY 'pm'
%   column, on all three schemes and under o1z = 'bare' (measured; review I1 corrects the
%   earlier claim that this flag is always pt.stable_1z and is false under 'bare'). It is the
%   ordered 1/z endpoint's OWN flag (invz_solve_point_ordered's pt.stable_1z) only on
%   ORDERED/DIAGNOSTIC columns -- recorded whenever the jensen leg returned it, so it CAN be
%   true on a MASKED column, e.g. a strict 'unknown' column where the jensen solve ran for
%   diagnostics only, or 'unstable_endpoint'. It is false only where NEITHER source fired: a
%   column with no PM label and no ordered solve (bare escape hatch, longitudinal tilt, masked
%   auto state, and o1z = 'bare' where pt.stable_1z is absent from the returned point
%   altogether), or an ordered point whose own flag came back false. Only S.phase_1z labels a
%   column.
%   S.ordered_diag_reason is the ordered (jensen) solve's OWN reduced verdict, recorded even
%   where the phase rule does not adopt it: 'not_attempted' (the ordered leg never ran) |
%   'accepted' (a closure-consistent ordered point was returned) | otherwise the same
%   'degenerate_doublet' / 'medium_out_of_domain' / 'solver_failed' vocabulary. This is the
%   field that carries the KNOWABLE cause on a strict column masked as 'pm_probe_unknown',
%   where the PM probe deliberately dominates the phase reason (a failed PM probe must never
%   vote 'ordered'), and it is what the sweep-summary counters attribute on -- so the summary
%   can never report zero degenerate doublets for a column that is one.
%   S.response_error_id is the FIRST recoverable identifier absorbed by a response-function
%   (invz_chi_realaxis) boundary in that column, '' when none fired. It is expected to be ''
%   everywhere, because no whitelisted identifier is reachable through invz_chi_realaxis today
%   (task-15 report SS6); it exists so that if one ever becomes reachable the cause is recorded
%   rather than inferred. Under a STRICT scheme a failed chiz response call also MASKS the
%   column ('response_failed', Sigma0 -> NaN), so an all-NaN spectrum can never anchor Bc_1z
%   while wearing a 'pm'/'ordered' label. A failed RPA-OVERLAY call is recorded but never
%   changes S.phase: the auto label is a solve outcome, not an overlay-evaluation outcome.
%   S.Sigma0 IS SCHEME-DEPENDENT ON MASKED COLUMNS -- do not compare it there across schemes.
%   Under a strict scheme it is NaN wherever phase_1z = 0, because no state was adopted for
%   chiz and there is therefore no Sigma0 to report; 'resummed' keeps its historical behaviour
%   of reporting the auto point's pt.Sigma0 on the longitudinal-failure and transverse tails
%   even when the column is masked. Any scheme-vs-scheme comparison of S.Sigma0 must restrict
%   itself to labelled columns.
%   Under a STRICT scheme the 1/z leg dispatches on invz_pm_verdict's three-way rule: a
%   failed or in-band PM probe is 'unknown' and may run the jensen solver for DIAGNOSTICS
%   only -- it can never emit phase_1z = 1, because solver availability is not a phase
%   criterion. Under 'resummed' the historical dispatch is byte-preserved and these fields
%   are pure diagnostics. S.Bc_1z_interval/[S.Bc_1z_status] are the anchor-bracketed boundary
%   (invz_boundary_interval: 'valid' | 'widened' | 'unbracketed' | 'invalid'); they FEED
%   S.Bc_1z only under a strict scheme, while 'resummed' keeps invz_boundary_field's
%   historical scalar reduction. FROM THIS DRIVER 'invalid' can ONLY mean that the sweep's
%   `fields` are empty, non-finite, non-numeric, non-real or not strictly increasing, so the
%   reduction was never run at all. The helper's other documented 'invalid' cause -- overlapping
%   masks -- is UNREACHABLE here by construction, because the three masks are built by strcmp
%   against the single reason string each column carries and are therefore mutually exclusive.
%   Read invz_boundary_interval's docstring with that restriction in mind. ONE FURTHER
%   RESTRICTION ON 'valid' (review M6): only 'pm_probe_unknown' / 'boundary_indeterminate' feed
%   `unk`, so an INTERVENING column masked for a KNOWABLE cause does not widen the bracket --
%   'valid' therefore means "no INDETERMINATE column between the anchors", not "every column
%   between the anchors is labelled". A 'response_failed' column between the anchors would leave
%   the status 'valid' while the reported interval spans an unlabelled column. That matches how
%   'medium_out_of_domain' / 'degenerate_doublet' / 'solver_failed' are already treated (a
%   knowable cause is a reportable defect, not boundary uncertainty), and is unreachable today;
%   the sweep-summary warning is what surfaces those columns instead.
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
fcont    = getf(opts, 'field_continuation', 'none');
wmin     = getf(opts, 'peak_wmin', 0);

fdir  = getf(opts, 'field_dir', [1 0 0]);
bztol = getf(opts, 'bz_tol', 1e-9);
sxtra = getf(opts, 'solve_opts', struct());
invz_check_solve_opts(sxtra);
% Static-medium scheme resolved ONCE at this public entry by the sole authority (spec SS4.2)
% and threaded into every per-field solve below, so a sweep can never mix truncation orders
% across its columns. static_medium/ref_margin are DRIVER-OWNED at this level --
% invz_check_solve_opts above rejects either inside solve_opts. Absent => 'resummed' =>
% numerically identical to the pre-strict path, dispatch included.
sm = invz_check_static_medium(opts);
o1z = getf(opts, 'ordered_1z', 'jensen');
if ~any(strcmp(o1z, {'jensen', 'bare'}))
    error('invz:ordered1z', 'ordered_1z must be ''jensen'' or ''bare''.');
end
if ~(ischar(fcont) && isrow(fcont)) || ~any(strcmp(fcont, {'none', 'qcp_down'}))
    error('invz:fieldContinuation', ...
        'field_continuation must be ''none'' or ''qcp_down''.');
end
if strcmp(fcont, 'qcp_down') && ~strcmp(o1z, 'jensen')
    error('invz:fieldContinuation', ...
        'field_continuation ''qcp_down'' requires ordered_1z = ''jensen''.');
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

% Ewald Step-5 Task 7: backend/grid options are
% forwarded BY PRESENCE into invz_bz_couplings on the compute branch, and a precomputed
% opts.Jnu/opts.info pair is validated against any EXPLICIT backend/grid-policy request rather
% than trusted blindly -- see invz_check_coupling_opts.m (shared with invz_spectra_qpath.m,
% task-7 review dedup fix) for the exact conflict rules.
chk = invz_check_coupling_opts();
hasBackendReq = isfield(opts, 'dipole') || isfield(opts, 'ewald');
hasGridReq    = isfield(opts, 'gridConvention') || isfield(opts, 'gridOffset') || isfield(opts, 'gammaPolicy');
if hasBackendReq
    % Validate opts.dipole/opts.ewald THEMSELVES even though the precomputed branch below may
    % never call invz_jq_modes: a malformed request must not escape checking just because no
    % lattice sum runs.
    [backendReq, eoptsReq] = chk.validate_dipole_opts(opts);
end

hasJnuOpt  = isfield(opts, 'Jnu');
hasInfoOpt = isfield(opts, 'info');
if hasJnuOpt ~= hasInfoOpt
    error('invz:spectraPrecomputedPartial', ...
          'opts.Jnu and opts.info must be supplied together (exactly one was present).');
end

if hasJnuOpt && hasInfoOpt
    Jnu = opts.Jnu(:);   info = opts.info;
    if hasBackendReq, chk.check_backend_provenance(info, backendReq, eoptsReq); end
    if hasGridReq,    chk.check_grid_provenance(info, opts, grid);              end
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
% Per-column dispatch provenance (spec SS4.4). Cells slice under parfor exactly like the
% numeric arrays; one_field assigns every one of them on every return path.
smUsed  = cell(1, nB);   pmStat  = cell(1, nB);   pmEid = cell(1, nB);
reasonC = cell(1, nB);   stab1z  = false(1, nB);
% ordDiag: the ordered leg's OWN verdict, kept even where the PM-probe-dominant phase rule
% does not adopt it (review F2 -- paying for a diagnostic solve and discarding its verdict left
% the sweep counters reporting falsehoods). resEid: the identifier absorbed by a response-call
% boundary, so a NaN spectrum on a labelled column can never be causeless (review F3).
ordDiag = cell(1, nB);   resEid  = cell(1, nB);
seedOut = cell(1, nB);   seedSupplied = false(1, nB);
seedSourceField = nan(1, nB);

% nWorkers = 0 forces serial execution even inside a parfor, and works without the
% Parallel Computing Toolbox; Inf lets parfor use (and auto-create) the pool.
nWorkers = 0;
if parallel && ~isempty(ver('parallel')), nWorkers = Inf; end

sopts = sxtra;
sopts.hyp = hyp;  sopts.J0eff = Jcc0;  sopts.Jxx0 = Jaa0;  sopts.bz_tol = bztol;
% static_medium/ref_margin are driver-owned at this level: stamped here so BOTH point solvers
% resolve the SAME scheme. 'resummed' is their own default, so this is inert by construction
% on the protected path.
sopts.static_medium = sm.scheme;  sopts.ref_margin = sm.ref_margin;

if strcmp(fcont, 'qcp_down')
    if parallel
        warning('invz:fieldContinuationSerial', ...
            'field_continuation ''qcp_down'' is sequential; opts.parallel is ignored.');
    end
    lastSeed = [];  lastSeedField = NaN;
    for k = nB:-1:1
        sk = sopts;
        if ~isempty(lastSeed)
            sk.hmf_seed = lastSeed;
            seedSupplied(k) = true;
            seedSourceField(k) = lastSeedField;
        end
        [chizM(:, k), chirpaM(:, k), Sig0(k), phaseC(k), ph1z(k), critPM(k), ...
         m1zM(k), DordM(k), smUsed{k}, pmStat{k}, pmEid{k}, stab1z(k), ...
         reasonC{k}, ordDiag{k}, resEid{k}, seedOut{k}] = ...
            one_field(ion, T, BvecM(k, :), Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, ...
                      sk, bztol, o1z, sm);
        % Only a fully accepted ordered column advances physical-field branch history.
        % A failed field cannot poison the next seed, even if its predictor happened to close.
        if ph1z(k) == 1 && ~isempty(seedOut{k})
            lastSeed = seedOut{k};
            lastSeedField = fields(k);
        end
        if verbose
            ph = {'masked', 'moment-form', 'paramagnet'};
            fprintf('  |B| = %5.2f T : auto-state %-11s 1z-state %-11s Sigma0 = %s seed=%s\n', ...
                    fields(k), ph{phaseC(k)+1}, ph{ph1z(k)+1}, num2str(Sig0(k)), ...
                    local_seed_label(seedSupplied(k), seedSourceField(k)));
        end
    end
else
    parfor (k = 1:nB, nWorkers)
        [chizM(:, k), chirpaM(:, k), Sig0(k), phaseC(k), ph1z(k), critPM(k), ...
         m1zM(k), DordM(k), smUsed{k}, pmStat{k}, pmEid{k}, stab1z(k), ...
         reasonC{k}, ordDiag{k}, resEid{k}, seedOut{k}] = ...
            one_field(ion, T, BvecM(k, :), Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, ...
                      sopts, bztol, o1z, sm);
        if verbose
            ph = {'masked', 'moment-form', 'paramagnet'};
            fprintf('  |B| = %5.2f T : auto-state %-11s 1z-state %-11s Sigma0 = %s\n', ...
                    fields(k), ph{phaseC(k)+1}, ph{ph1z(k)+1}, num2str(Sig0(k)));
        end
    end
end

S = struct();
S.fields = fields;  S.w = w;  S.T = T;  S.info = info;
S.field_dir = fdir;  S.Bvec = BvecM;  S.transverse_mf = tmf;  S.ordered_1z = o1z;
S.field_continuation = fcont;
S.hmf_seed_supplied = seedSupplied;
S.hmf_seed_source_field = seedSourceField;
S.Sigma0 = Sig0;  S.phase = phaseC;   % auto (ordered-first bare-MF) dispatch -- see header
S.phase_1z = ph1z;  S.crit_pm = critPM;  S.m_1z = m1zM;  S.D_ord = DordM;
% Dispatch provenance (spec SS4.4): the REQUESTED scheme is scalar and mandatory; everything
% else is per column, so a masked column's cause is recorded rather than inferred.
S.static_medium = sm.scheme;      S.static_medium_used = smUsed;
S.pm_probe_status = pmStat;       S.pm_probe_error_id = pmEid;
S.stability_1z = stab1z;          S.phase_1z_reason = reasonC;
S.ordered_diag_reason = ordDiag;  S.response_error_id = resEid;
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
% --- 1/z boundary (spec SS4.4) --------------------------------------------------------------
% Anchor-bracketed reduction: only unknown/indeterminate columns BETWEEN the last ordered and
% the first stable-PM anchor widen the bracket. 'ordered' deliberately excludes
% bare_escape_hatch -- the retired Stage-1 diagnostic must never anchor a 1/z boundary.
unk = strcmp(S.phase_1z_reason, 'pm_probe_unknown') | ...
      strcmp(S.phase_1z_reason, 'boundary_indeterminate');
ord = strcmp(S.phase_1z_reason, 'ordered');       % excludes bare_escape_hatch
pm  = strcmp(S.phase_1z_reason, 'pm');
% invz_boundary_interval owns a HARD input contract on `fields`: NUMERIC, REAL, NONEMPTY,
% FINITE, STRICTLY INCREASING (five predicates -- see its line 33/38 guards). The driver has
% never constrained the sweep order OR required a nonempty sweep, so the precondition is
% screened HERE rather than allowed to become a new throw on the protected 'resummed' path.
%
% WHY THE SCREEN MUST BE A SUPERSET OF THE HELPER'S PRECONDITION (review F1, completed by M1).
% The call is pure and its outputs land in fields that did not exist before, so the ONLY way it
% can change default-path behaviour is by THROWING. G9 therefore holds exactly when this screen
% rejects everything the helper rejects. Mirror all five, in the helper's own order:
%   isnumeric  -- MEASURED GAP: a char `fields` sweeps fine (char*double is double) and reaches
%                 here, where all(isfinite) and diff(...) > 0 both pass, so it hit the helper and
%                 raised invz:boundaryInterval on the DEFAULT path, where b9082fd returned.
%   isreal     -- defensive: a genuinely complex `fields` dies earlier at invz_field_vec
%                 (invz:fieldVec, unchanged since b9082fd), and complex-with-zero-imaginary-part
%                 narrows to real in fields(:).', so no case is known to reach here. Mirrored
%                 anyway, because completeness of the mirror is what makes the G9 argument hold
%                 without a reachability case analysis.
%   ~isempty   -- FIRST and explicitly: all([]) is TRUE and 0 < 2 is TRUE, so an empty sweep
%                 otherwise slips past both later clauses (the original F1 defect).
%   isfinite / strictly increasing -- as before.
% The helper's precondition is right and is NOT loosened; only the screen was incomplete.
if isnumeric(fields) && isreal(fields) && ~isempty(fields) && all(isfinite(fields)) && ...
        (nB < 2 || all(diff(fields) > 0))
    [S.Bc_1z_interval, S.Bc_1z_status, Bc1z_strict] = invz_boundary_interval(fields, ord, pm, unk);
else
    S.Bc_1z_interval = [NaN NaN];  S.Bc_1z_status = 'invalid';  Bc1z_strict = NaN;
end
if sm.is_strict
    S.Bc_1z = Bc1z_strict;
    % ONE sweep-level summary, never a per-node warning: a 61-point sweep must not emit 61
    % near-identical messages, and the counts are what a reader needs to judge the map.
    %
    % COUNTER ATTRIBUTION (review F2). phase_1z_reason is deliberately PM-probe-dominant: a
    % failed PM probe must never vote 'ordered', so a column whose PM probe died reads
    % 'pm_probe_unknown' even when the DIAGNOSTIC ordered solve reduced a knowable cause. That
    % is the right per-column phase vocabulary, but it must not make the SUMMARY false -- at
    % B = 0 it reported "0 degenerate-doublet, 1 unknown-PM-probe" for a column that IS a
    % degenerate doublet, and Stage-4 G16 reads these counters. So a MASKED column whose
    % ordered_diag_reason names a knowable cause is attributed to that cause instead. One
    % cause per column, so the reported total is still a column count and never double counts.
    %
    % SCOPED TO THE FAILED-PROBE CASE ONLY (review M2). Keying on ph1z == 0 also re-attributed
    % 'boundary_indeterminate' columns, where the PM probe CONVERGED and merely sits inside the
    % +/-crit_tol band: there is no failed-probe gap to fill there, and re-labelling such a
    % column by its diagnostic jensen verdict would hide that it is the band column. Only
    % 'pm_probe_unknown' -- where the phase reason is PM-probe-dominant precisely BECAUSE the
    % probe told us nothing -- may borrow the ordered leg's verdict. `unk` itself is untouched,
    % so S.Bc_1z and S.Bc_1z_interval are unaffected; this is summary attribution only.
    diag_cause = strcmp(S.phase_1z_reason, 'pm_probe_unknown') & ...
                 ismember(S.ordered_diag_reason, {'medium_out_of_domain', 'degenerate_doublet'});
    cause = S.phase_1z_reason;
    cause(diag_cause) = S.ordered_diag_reason(diag_cause);
    n_dom = nnz(strcmp(cause, 'medium_out_of_domain'));
    n_deg = nnz(strcmp(cause, 'degenerate_doublet'));
    n_unk = nnz(unk & ~diag_cause);
    % A LOST BOUNDARY MUST NEVER BE SILENT (review F7). Under ordered_1z = 'bare' every ordered
    % column reports 'bare_escape_hatch', which by design cannot anchor a 1/z boundary, so `ord`
    % is empty and Bc_1z comes back NaN -- where 'resummed' returns a finite one from the same
    % run. With no masked columns the count-based trigger below never fired, so the difference
    % was invisible. Trigger on the lost boundary too whenever the sweep DID produce columns
    % that the historical reduction would have anchored on (bare escape hatches) or masked ones.
    n_bare = nnz(strcmp(S.phase_1z_reason, 'bare_escape_hatch'));
    n_msk  = nnz(ph1z == 0);
    % NO MASKED COLUMN MAY BE SILENT OR UNDERCOUNTED (review I2). n_dom/n_deg/n_unk between them
    % cover only 4 of the 9 masked reason strings: 'unstable_endpoint', 'solver_failed',
    % 'not_attempted_longitudinal', 'bare_not_ordered' and 'response_failed' contribute ZERO to
    % all three. So (a) the HEADLINE is n_msk, the true unlabelled-column count, with the three
    % attributed causes as a breakdown plus an explicit 'other' remainder that closes the sum,
    % and (b) n_msk itself arms the trigger. boundary_lost alone could not: it requires
    % ~isfinite(Bc_1z), so a strict sweep with a FINITE Bc_1z plus one 'solver_failed' column
    % emitted NO warning at all. The three counts are mutually exclusive and count only masked
    % columns, so n_oth >= 0 by construction.
    n_oth = n_msk - (n_dom + n_deg + n_unk);
    boundary_lost = ~isfinite(S.Bc_1z) && (n_bare > 0 || n_msk > 0);
    if n_msk > 0 || boundary_lost
        warning('invz:spectraMapMasked', ...
            ['1/z dispatch left %d of %d columns unlabelled: %d medium-out-of-domain, ' ...
             '%d degenerate-doublet, %d unknown-PM-probe, %d other; %d column(s) are ' ...
             'bare-escape-hatch (never a 1/z anchor). Bc_1z = %g (status ''%s'').'], ...
            n_msk, nB, n_dom, n_deg, n_unk, n_oth, n_bare, S.Bc_1z, S.Bc_1z_status);
    end
else
    S.Bc_1z = invz_boundary_field(fields, ph1z  == 1, ph1z  == 2);   % historical scalar (G9)
end
S.Epeak     = invz_peak_energy(chizM,   w, wmin);
S.Epeak_rpa = invz_peak_energy(chirpaM, w, wmin);
end

% -------------------------------------------------------------------------------------------
function [chiz, chirpa, Sigma0, phase, phase_1z, crit_pm, m_1z, D_ord, static_medium_used, ...
          pm_probe_status, pm_probe_error_id, stability_1z, phase_1z_reason, ...
          ordered_diag_reason, response_error_id, hmf_seed_out] = ...
    one_field(ion, T, B, Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts, bztol, o1z, sm)
%ONE_FIELD chi''_cc(omega) at one field -- TWO independently phased legs (QCP Stage 1;
% ordered 1/z leg below Bc_1z upgraded to the Stage-2 Jensen solve; see
% invz_projected/README.html and invzp_convg_diagnosis.md):
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
% sm is the resolved static-medium record. sm.is_strict switches the 1/z leg to the three-way
%     dispatch (spec SS4.4) at the TOP of this function and returns before the historical body
%     is reached, so the protected 'resummed' path keeps its exact order and arithmetic and the
%     new phase rule cannot leak into it. phase_1z_reason and the pm_probe_* / stability_1z /
%     static_medium_used provenance are populated on BOTH paths -- on 'resummed' they are pure
%     bookkeeping written alongside the existing decisions, never used to make one.
% Jsel = Jcc0 is the strict-uniform observable, so the demag correction Jshape applies.
nw = numel(w);
chiz = nan(nw, 1);  chirpa = nan(nw, 1);  Sigma0 = NaN;  phase = 0;
phase_1z = 0;  crit_pm = NaN;  m_1z = NaN;  D_ord = NaN;
pm_probe_status = 'not_attempted';  pm_probe_error_id = '';
stability_1z = false;  phase_1z_reason = 'solver_failed';
static_medium_used = sm.scheme;
% ordered_diag_reason records the ordered (jensen) leg's OWN verdict wherever that leg runs,
% including the strict 'unknown' arm where it is diagnostic only -- the phase reason there is
% PM-probe-dominant by design, so without this field the knowable cause is paid for and thrown
% away. 'not_attempted' is the honest default: the leg did not run on this column.
% response_error_id stays '' unless a response-function (invz_chi_realaxis) boundary absorbs a
% recoverable identifier, which no reachable path does today.
ordered_diag_reason = 'not_attempted';  response_error_id = '';
hmf_seed_out = [];
crit_tol = getf(sopts, 'crit_tol', 1e-6);   % the frozen band, shared with the ordered checker
tmf = getf(sopts, 'transverse_mf', 'legacy_x');
% transverse_mf here so any fallback single-ion rebuild inside invz_chi_realaxis (when a
% branch below can't supply si) uses the SAME MF model as the solve, not the 'legacy_x'
% default (review finding I1; the C1 companion bug in invz_spectra_qpath.m).
copts = struct('Jsel', Jcc0, 'eta', eta, 'Jxx0', Jaa0, 'Jshape', Jshape, 'hyp', hyp, ...
               'transverse_mf', tmf);

if sm.is_strict
    % =========================================================================================
    % STRICT ONE-SHOT CLOSURE: three-way phase dispatch (spec SS4.4). Gated, and returning
    % before the historical body, because an UNGATED three-way rule would classify the
    % below-Bc_1z PM probe -- which fails for the same pole reason under 'resummed' -- as
    % 'unknown' and mask MORE columns than today, breaking G9.
    % =========================================================================================
    [pt, phase] = invz_solve_auto(ion, T, B, Jnu, sopts);

    % (0) DEFERRED ROBUSTNESS ITEM, stated rather than hidden: this stage does not un-nest the
    % 1/z leg from auto availability. With no accepted auto state the strict column stays
    % masked. Stage-4 G16 must expose that as a limitation, never report it as a strict solve.
    if phase == 0
        if abs(B(3)) > bztol, phase_1z_reason = 'not_attempted_longitudinal';
        else,                 phase_1z_reason = 'solver_failed';  end
        return;
    end
    longit = abs(B(3)) > bztol;

    % (1) PM probe, TRANSVERSE ONLY (the strict-PM solver's m = 0 gate raises
    % invz:orderedPhase under a longitudinal tilt). ONE probe serves BOTH successful auto
    % states: a converged strict auto-PM point IS the probe, so it is neither skipped on an
    % auto-PM column (the previous draft's hole: phase == 2 with a negative PM mass never
    % reached the ordered leg) nor solved a second time.
    ptp = [];  pm_completed = true;
    if ~longit
        if phase == 2 && ~isempty(pt) && isfield(pt, 'converged') && pt.converged && ...
                isfinite(pt.Sigma0) && strcmp(getf(pt, 'static_medium', ''), sm.scheme)
            ptp = pt;
        else
            [ptp, pm_completed, pm_probe_error_id] = invz_try_solver_call( ...
                @() invz_solve_point(ion, T, B, Jnu, sopts));
        end
    end
    if ~isempty(ptp) && isfield(ptp, 'crit') && isfinite(ptp.crit), crit_pm = ptp.crit; end
    if ~pm_completed
        % The exact identifier is already in pm_probe_error_id and must NEVER be collapsed
        % into 'nonconverged' -- a recoverable domain signal and a failed iteration are
        % different causes with different follow-ups.
        pm_probe_status = 'recoverable_error';  pm_consistent = false;
    else
        pm_probe_status = local_pm_status(ptp, crit_pm, crit_tol);
        pm_consistent = ~isempty(ptp) && isfield(ptp, 'converged') && ptp.converged && ...
                        isfinite(ptp.Sigma0);
    end
    verdict = invz_pm_verdict(crit_pm, pm_consistent, crit_tol);

    % (2) three-way dispatch -------------------------------------------------------------
    if strcmp(verdict, 'pm')
        phase_1z = 2;  phase_1z_reason = 'pm';  stability_1z = true;  Sigma0 = ptp.Sigma0;
        c1 = copts;  c1.si = ptp.si;             % PM chi0 differs from the ordered leg's
        [o1, ok1, e1] = invz_try_solver_call(@() invz_chi_realaxis(ion, T, B, ptp, w, c1));
        if ok1
            chiz = imag(o1.chi_cc_q(1, :)).';
        else
            % REVIEW F3. Discarding this identifier reproduced, inside the new code, exactly
            % the defect this task exists to remove: phase_1z was already 2 and the reason
            % already 'pm', so the column would ANCHOR Bc_1z while carrying an all-NaN spectrum
            % and no recorded cause. Mask it, name the cause, keep the identifier. Sigma0 goes
            % back to NaN because it is documented as the Sigma0 of the state used for chiz and
            % there is now no such state. Unreachable today -- no whitelisted identifier is
            % reachable through invz_chi_realaxis (task-15 report SS6) -- and the strict
            % provenance test pins that by asserting response_error_id is empty, so if it ever
            % becomes reachable the pin fails instead of the honesty thesis.
            phase_1z = 0;  phase_1z_reason = 'response_failed';
            Sigma0 = NaN;  response_error_id = e1;
        end
    elseif longit || strcmp(o1z, 'bare')
        % Documented escape hatch (longitudinal tilt, or the explicit 'bare' opt-out): the
        % retired Stage-1 bare-MF diagnostic. It is NOT a result of the requested strict
        % scheme, so the provenance says so and the strict Bc_1z reduction, which anchors on
        % 'ordered' alone, never treats it as an ordered anchor.
        if phase == 1
            phase_1z = 1;  phase_1z_reason = 'bare_escape_hatch';
            static_medium_used = 'n/a_bare_escape';  Sigma0 = pt.Sigma0;
            [ob, okb, eb] = invz_try_solver_call(@() invz_chi_realaxis(ion, T, B, pt, w, copts));
            if okb
                chiz = imag(ob.chi_cc_q(1, :)).';
            else                                     % review F3, as in the 'pm' arm above
                phase_1z = 0;  phase_1z_reason = 'response_failed';
                Sigma0 = NaN;  response_error_id = eb;
            end
        elseif longit
            phase_1z_reason = 'not_attempted_longitudinal';
        else
            phase_1z_reason = 'bare_not_ordered';   % 'bare' needs a moment-form auto state
        end
    else
        % Jensen leg. Under 'ordered_eligible' it may LABEL the column, and only when BOTH
        % acceptance tiers pass; under 'unknown' it is DIAGNOSTIC ONLY, because a failed PM
        % probe is not evidence for order (invz_pm_verdict's whole point).
        so2 = sopts;  so2.ordered_mode = 'jensen';
        [ptj, j_completed, j_error_id] = invz_try_solver_call( ...
            @() invz_solve_point_ordered(ion, T, B, Jnu, so2));
        hmf_seed_out = local_hmf_seed(ptj);
        consistent = j_completed && ~isempty(ptj) && ptj.is_ordered && ptj.converged && ...
                     isfinite(ptj.Sigma0);
        % stable_1z is jensen-only and ABSENT (not NaN) in bare mode: guard on isfield.
        stability_1z = j_completed && ~isempty(ptj) && isfield(ptj, 'stable_1z') && ptj.stable_1z;
        if consistent
            m_1z = getf(ptj, 'm0', NaN);  D_ord = getf(ptj, 'D_uni', NaN);
        end
        % REVIEW F2. The verdict of this solve is recorded on BOTH arms, not just the arm that
        % may act on it. 'unknown' pays for a full jensen solve for diagnostics; discarding what
        % it found is what made the sweep counters report "0 degenerate-doublet" for a B = 0
        % column that IS a degenerate doublet. 'accepted' rather than a local_ordered_reason
        % call on the consistent path: that classifier reduces REJECTION causes, and would
        % flatten a perfectly good ordered point into 'solver_failed'.
        if consistent
            ordered_diag_reason = 'accepted';
        else
            ordered_diag_reason = local_ordered_reason(ptj, j_completed, j_error_id);
        end
        if strcmp(verdict, 'ordered_eligible') && consistent && stability_1z
            phase_1z = 1;  phase_1z_reason = 'ordered';  Sigma0 = ptj.Sigma0;
            [oj, okj, ej] = invz_try_solver_call(@() invz_chi_realaxis(ion, T, B, ptj, w, copts));
            if okj
                chiz = imag(oj.chi_cc_q(1, :)).';
            else                                     % review F3, as in the 'pm' arm above
                phase_1z = 0;  phase_1z_reason = 'response_failed';
                Sigma0 = NaN;  response_error_id = ej;
            end
        elseif strcmp(verdict, 'ordered_eligible')
            if consistent            % closure-consistent but the endpoint is not stable
                phase_1z_reason = 'unstable_endpoint';
            else
                phase_1z_reason = ordered_diag_reason;
            end
        elseif strcmp(pm_probe_status, 'boundary_band')
            phase_1z_reason = 'boundary_indeterminate';   % finite PM point inside the band
        else
            % PM-probe dominance is DELIBERATE and unchanged: a failed PM probe must never vote
            % 'ordered', so the PHASE reason stays a PM-probe reason even when the diagnostic
            % solve above knows more. The knowable cause is not lost -- it is in
            % ordered_diag_reason, which is also what the sweep-summary counters attribute on.
            phase_1z_reason = 'pm_probe_unknown';         % nonconvergence/nonfinite/recoverable
        end
    end

    % (3) auto/RPA overlay: the same branch-specific arithmetic as the historical body, and it
    % NEVER votes on phase_1z -- the two legs are independently phased.
    c0s = copts;  c0s.npass = 1;
    if phase == 1
        [o0s, ok0s, e0s] = invz_try_solver_call( ...
            @() invz_chi_realaxis(ion, T, B, invz_zero_sigma_overlay(pt), w, c0s));
    else
        c0s.si = pt.si;
        pt0s = struct('alpha', 0, 'lambda', [0; 0], 'tl', pt.tl, 'K', []);
        [o0s, ok0s, e0s] = invz_try_solver_call(@() invz_chi_realaxis(ion, T, B, pt0s, w, c0s));
    end
    if ok0s
        chirpa = imag(o0s.chi_cc_q(1, :)).';
    elseif isempty(response_error_id)
        % Review F3, overlay half: RECORD ONLY. A failed overlay evaluation leaves chirpa NaN on
        % a column the driver does not mask (invalid_auto keys on S.phase), so without this the
        % NaN would be causeless -- but S.phase is the AUTO SOLVE's label and must not be flipped
        % by an overlay-evaluation failure, which would silently move Bc_auto. First identifier
        % wins, so a chiz failure above (which does change the column's state) is never
        % overwritten by an overlay one.
        response_error_id = e0s;
    end
    return;
end
% ============================ HISTORICAL 'resummed' BODY ====================================
% Unchanged in order and arithmetic; the only additions are provenance assignments.

[pt, phase] = invz_solve_auto(ion, T, B, Jnu, sopts);

if phase == 1                                     % --- moment-form branch (FM or induced) ---
    % 1/z leg: strict-PM stability probe, TRANSVERSE ONLY (the strict-PM solver's m = 0
    % gate would raise invz:orderedPhase under a longitudinal tilt -- review finding 5).
    ptp = [];
    if B(3) == 0
        [ptp, pm_completed, pm_probe_error_id] = invz_try_solver_call( ...
            @() invz_solve_point(ion, T, B, Jnu, sopts));
        if ~pm_completed, pm_probe_status = 'recoverable_error'; end
    end
    if ~isempty(ptp) && isfield(ptp, 'crit') && isfinite(ptp.crit), crit_pm = ptp.crit; end
    if ~strcmp(pm_probe_status, 'recoverable_error')       % provenance only -- pm_valid below
        pm_probe_status = local_pm_status(ptp, crit_pm, crit_tol);   % is the historical gate
    end
    % crit_pm doubles as the field-guarded crit: pm_valid never touches ptp.crit directly
    % (an early-return ptp may lack the field).
    pm_valid = ~isempty(ptp) && ptp.converged && isfinite(ptp.Sigma0) && ...
               isfinite(crit_pm) && crit_pm > 0;

    if pm_valid                         % --- Bc_1z < |B| < bare boundary: 1/z theory is PM ---
        phase_1z = 2;  Sigma0 = ptp.Sigma0;
        phase_1z_reason = 'pm';  stability_1z = true;
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
            [ptj, j_completed, j_error_id] = invz_try_solver_call( ...
                @() invz_solve_point_ordered(ion, T, B, Jnu, so2));
            hmf_seed_out = local_hmf_seed(ptj);
            if ~isempty(ptj) && ptj.is_ordered && ptj.converged && isfinite(ptj.Sigma0)
                phase_1z = 1;  Sigma0 = ptj.Sigma0;  m_1z = ptj.m0;  D_ord = ptj.D_uni;
                phase_1z_reason = 'ordered';
                ordered_diag_reason = 'accepted';          % provenance only, no decision here
                stability_1z = isfield(ptj, 'stable_1z') && ptj.stable_1z;
                oj = invz_chi_realaxis(ion, T, B, ptj, w, copts);   % jensen si differs from
                chiz = imag(oj.chi_cc_q(1, :)).';                   % the auto pt's -- no sharing
            else                             % phase_1z stays 0 -> chiz column masked, but the
                                             % CAUSE is now recorded instead of inferred
                phase_1z_reason = local_ordered_reason(ptj, j_completed, j_error_id);
                ordered_diag_reason = phase_1z_reason;     % same verdict; recorded on both paths
            end
            pt0 = invz_zero_sigma_overlay(pt);                      % overlay: UNCHANGED auto state
            c0opts = copts;  c0opts.npass = 1;
            o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
            chirpa = imag(o0.chi_cc_q(1, :)).';
        else                                 % 'bare', or a longitudinal tilt under 'jensen':
                                              % Stage-1 diagnostic escape hatch
            phase_1z = 1;  Sigma0 = pt.Sigma0;
            phase_1z_reason = 'bare_escape_hatch';
            static_medium_used = 'n/a_bare_escape';   % not a result of the requested scheme
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
    phase_1z_reason = 'not_attempted_longitudinal';   % neither 1/z leg is defined under a tilt
    if ~isempty(pt) && ~isempty(pt.si) && isfield(pt, 'tl') && ~isempty(pt.tl)
        pt0 = invz_zero_sigma_overlay(pt);
        c0opts = copts;  c0opts.npass = 1;
        [o0, ok0, e0] = invz_try_solver_call(@() invz_chi_realaxis(ion, T, B, pt0, w, c0opts));
        if ok0, chirpa = imag(o0.chi_cc_q(1, :)).';
        else,   response_error_id = e0;    % review F3: provenance only, no decision changes
        end
    end
    if ~isempty(pt) && isfield(pt, 'Sigma0'), Sigma0 = pt.Sigma0; end
    return;
end

% --- transverse paramagnetic side: unchanged historical logic --------------------------
% THE ONLY SITE LEFT IN THIS FILE THAT CAN ABORT THE SWEEP ON A *WHITELISTED RECOVERABLE*
% SIGNAL (review F4, headline qualified by M4; task-15 report SS10.5). It is NOT the only
% aborting site, and the unqualified claim was false: the seven unwrapped invz_chi_realaxis
% calls in this file's historical body (six above, one below) and invz_solve_auto itself abort
% on any NON-whitelisted invz:* identifier -- deliberately, that is the narrowed error policy.
% What is unique HERE is that the invz_twolevel call below sits OUTSIDE every
% invz_try_solver_call boundary and can raise the whitelisted invz:degenerateDoublet
% (invz_twolevel.m: Delta < 1e-4 meV domain floor), which
% aborts the whole parfor sweep rather than masking one column. That is deliberate and must stay:
% wrapping it would WIDEN error absorption on the protected 'resummed' path, the exact opposite
% of this task's direction, and it is the pre-existing HEAD behaviour. It is reached only when
% the auto solve did NOT return a usable phase-2 point at a transverse field -- in practice
% unreachable at B -> 0, where the auto leg orders first -- so the exposure is a transverse
% column whose auto solve fails AND whose doublet is degenerate. Recorded here because the
% report that used to be its only record is gitignored and will not survive a merge.
if phase == 2 && ~isempty(pt), tl0 = pt.tl;  si0 = pt.si;
else, tl0 = invz_twolevel(ion, T, B, struct('Jxx0', Jaa0, 'transverse_mf', tmf));  si0 = []; end
chi0cc = [];
pt0 = struct('alpha', 0, 'lambda', [0; 0], 'tl', tl0, 'K', []);
c0opts = copts;  c0opts.npass = 1;  c0opts.si = si0;
[o0, ok0, e0] = invz_try_solver_call(@() invz_chi_realaxis(ion, T, B, pt0, w, c0opts));
if ok0
    chirpa = imag(o0.chi_cc_q(1, :)).';
    chi0cc = o0.chi0cc_w;                          % share the bare cc with the 1/z call
else
    response_error_id = e0;                        % review F3: provenance only, no decision
end
if ~isempty(pt) && isfield(pt, 'Sigma0'), Sigma0 = pt.Sigma0; end
if ~isempty(pt) && isfield(pt, 'crit') && isfinite(pt.crit), crit_pm = pt.crit; end
% NAME CAVEAT (review F9, provenance only): on a phase == 0 column `pt` is NOT a PM probe -- it
% is whatever non-accepted point invz_solve_auto returned from its own dispatch -- so
% pm_probe_status here reduces that point, not a dedicated strict-PM probe. Only the phase == 1
% branch above and the strict block run a probe worthy of the name. Read this field as "status of
% the PM-side point this column had", and pair it with S.phase before drawing a conclusion.
if ~isempty(pt), pm_probe_status = local_pm_status(pt, crit_pm, crit_tol); end
if phase == 2                                     % --- converged paramagnetic 1/z ---
    if isfinite(crit_pm) && crit_pm > 0
        phase_1z = 2;                             % the auto PM state IS the 1/z leg here
        phase_1z_reason = 'pm';  stability_1z = true;
    else                                          % crit <= 0: spurious below-Bc PM point
        phase_1z_reason = 'unstable_endpoint';    % (documented tensor-side failure mode) --
    end                                           % phase_1z stays 0, column flagged suspect
    copts1 = copts;  copts1.si = pt.si;           % reuse the solve's si (C1 companion: pt is
    if ~isempty(chi0cc), copts1.chi0cc_w = chi0cc; end   % not is_ordered, so without si this
    o = invz_chi_realaxis(ion, T, B, pt, w, copts1);     % would silently rebuild at 'legacy_x'
    chiz = imag(o.chi_cc_q(1, :)).';
end
end

% -------------------------------------------------------------------------------------------
function seed = local_hmf_seed(ptj)
%LOCAL_HMF_SEED Extract only a residual-accepted Jensen predictor continuation state.
seed = [];
if isempty(ptj) || ~isstruct(ptj) || ~isscalar(ptj), return; end
hp = getf(ptj, 'hmf_prof', struct());
candidate = getf(hp, 'hmf_seed_out', []);
if isstruct(candidate) && isscalar(candidate) && ...
        isfield(candidate, 'Sigma') && isfield(candidate, 'K0s') && ...
        isnumeric(candidate.Sigma) && isreal(candidate.Sigma) && ...
        ~isempty(candidate.Sigma) && all(isfinite(candidate.Sigma(:))) && ...
        isnumeric(candidate.K0s) && isreal(candidate.K0s) && ...
        isscalar(candidate.K0s) && isfinite(candidate.K0s)
    seed = struct('Sigma', candidate.Sigma(:), 'K0s', candidate.K0s);
end
end

% -------------------------------------------------------------------------------------------
function label = local_seed_label(supplied, sourceField)
%LOCAL_SEED_LABEL Compact progress-only physical-field seed provenance.
if supplied && isfinite(sourceField)
    label = sprintf('warm<-%0.3fT', sourceField);
else
    label = 'cold';
end
end

% -------------------------------------------------------------------------------------------
function s = local_pm_status(ptp, crit_pm, crit_tol)
%LOCAL_PM_STATUS Narrow PM-probe vocabulary (spec SS4.4), so a masked column's cause is
% recorded rather than inferred. Fatal/unclassified errors never reach here: they rethrow.
if isempty(ptp),                          s = 'not_attempted';
elseif ~isfield(ptp, 'converged') || ~ptp.converged, s = 'nonconverged';
elseif ~isfinite(crit_pm),                s = 'nonfinite';
elseif abs(crit_pm) <= crit_tol,           s = 'boundary_band';
elseif crit_pm > crit_tol,                 s = 'stable';
else,                                      s = 'unstable';
end
end

% -------------------------------------------------------------------------------------------
function r = local_ordered_reason(ptj, completed, error_id)
%LOCAL_ORDERED_REASON Narrow masked-column vocabulary for a REJECTED ordered (jensen) leg.
% The solver already reduces its own node records with a binding precedence
% (invz_hmf_status.m: degenerate > reference-domain > node failure > unresolved), and exports
% the verdict as pt.hmf_status / pt.hmf_prof.status. Reading it is what keeps a domain or
% degeneracy cause from being flattened into a generic 'solver_failed'.
% completed/error_id come from invz_try_solver_call: only whitelisted recoverable identifiers
% can arrive here at all, so anything outside them already rethrew.
% It reduces REJECTION causes only, so an ACCEPTED point must never be passed through it (it
% would come back 'solver_failed'); callers record 'accepted' for that case instead. It also
% feeds S.ordered_diag_reason, which is why it is called on the strict 'unknown' arm too, where
% the phase reason stays PM-probe-dominant but the verdict must still be recorded.
if ~completed
    if strcmp(error_id, 'invz:degenerateDoublet'), r = 'degenerate_doublet';
    else,                                          r = 'solver_failed';
    end
    return;
end
if isempty(ptj), r = 'solver_failed';  return; end
hs = getf(ptj, 'hmf_status', '');
if isempty(hs), hs = getf(getf(ptj, 'hmf_prof', struct()), 'status', ''); end
ms = getf(ptj, 'medium_status', 'not_applicable');
if strcmp(hs, 'degenerate_doublet')
    r = 'degenerate_doublet';
elseif strcmp(hs, 'medium_out_of_domain') || ~any(strcmp(ms, {'ok', 'not_applicable'}))
    r = 'medium_out_of_domain';
else
    r = 'solver_failed';         % node_failed / unresolved / no_bare_order / not accepted
end
end

% Ewald Step-5 Task 7: precomputed-coupling provenance/conflict validation now lives in the
% shared invz_check_coupling_opts.m helper (task-7 review dedup fix), used identically by
% invz_spectra_qpath.m.
