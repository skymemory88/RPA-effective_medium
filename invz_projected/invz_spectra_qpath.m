function S = invz_spectra_qpath(ion, T, B, qpath, w, opts)
%INVZ_SPECTRA_QPATH Ferromagnetic-mode chi''_cc(q, omega) along a q-path at fixed (T, B).
%   The default dispersion follows the UNIFORM ferromagnetic-mode coupling J(q) = v'*Jcc*v
%   (invz_jq_path P.Juni) -- the same effective coupling MF_RPA_Yikai.m uses, whose energy
%   TRENDS reproduce R 2007 Fig 3 (mode softens monotonically along (1,0,0)->(2,0,0)). This
%   is NOT neutron scattering intensity (no eigenvector/form-factor weights); quantitative
%   reproduction is not claimed.
%
%   S = invz_spectra_qpath(ion, T, B, qpath, w) computes chi''_cc along qpath (nq x 3,
%   r.l.u.) at fixed T (K) and field B (scalar or 1x3 vector, T), for the 1/z and bare-RPA theories. The 1/z
%   medium (Sigma, K, lambda) is solved ONCE at (T, B) on the BZ-integration grid
%   (ordered-first via invz_solve_auto); the path susceptibility follows from the single-
%   site response chit(w) via chi(q, w) = chit/(1 - J(q) chit), J(q) from invz_jq_path.
%
%   Demag semantics: the strict-uniform Jshape_cc transform is NOT applied here -- a
%   finite-q probe measures the intrinsic response, and Gamma-equivalent path points use
%   the q -> 0+ limit where the demagnetizing field cancels (R 2007). Spectra can still
%   shift with ion.demag through the transverse applied/internal-field relation (info.Jaa0).
%
%   Returns:
%     S.chiz, S.chirpa     [nw x nq]  1/z and bare-RPA chi''_cc(q, w)
%     S.Epeak, S.Epeak_rpa [1 x nq]   censored, parabolic-refined peak energy per q (NaN at
%                                     a non-positive or boundary maximum)
%     S.Jq      [1 x nq]   selected coupling along the path (meV): uniform FM-mode
%                          projection by default, or a sorted branch if opts.branch is set
%     S.snapped [1 x nq]   true where invz_jq_path replaced the raw truncated sum
%     S.s, S.s_cart        path distance in index (r.l.u.) / Cartesian (Ang^-1) coordinates
%     S.x, S.xlab          plot coordinate + label (varying Miller component for a single-
%                          axis path, e.g. h = 1..2, else falls back to S.s)
%     S.qpath, S.w, S.T, S.Bvec, S.Bmag, S.info, S.phase (1 = moment-form solve, 2 = strict-PM solve)
%     S.transverse_mf      resolved MF mode string (echoes opts.solve_opts.transverse_mf,
%                           default 'legacy_x')
%
%   opts fields (all optional):
%     .grid ([16 16 16]), .dpRng (30), .eta (5e-3)   as in invz_spectra_map
%     .bz_tol (1e-9)     T; longitudinal field threshold for dead band (same as invz_solve_auto)
%     .solve_opts        struct of reserved-field-checked solver overrides (fields J0eff,
%                        Jxx0, hyp are driver-owned and will error if present). transverse_mf
%                        ('legacy_x' | 'none' | 'vector_ab', default 'legacy_x') is a legal
%                        field forwarded to the solvers. Under 'legacy_x' (x-only mean field)
%                        a nonzero b-axis (y) field component is C4-inconsistent and errors
%                        invz:transverseMF; use 'vector_ab' for genuine in-plane rotation.
%     .branch (0)        0 (default) = uniform FM-mode coupling v'*Jcc*v (the physical
%                        single mode); 1..4 = follow that sorted-eigenvalue branch instead
%                        (exploratory; sorted index, NOT a tracked mode identity through
%                        crossings)
%     .snapfac (2.5)     Gamma-limit trust-radius factor (see invz_jq_path)
%     .peak_wmin (0.05)  meV; excludes the low-frequency hyperfine pole from the peak search
%     .Jnu, .info        precomputed BZ-grid branches / info (skips the lattice sum; tests).
%                        Must be supplied TOGETHER (error invz:spectraPrecomputedPartial
%                        otherwise).
%     .dipole, .ewald    opt-in Ewald dipolar backend: forwarded BY PRESENCE into
%                        invz_bz_couplings when computed, and the SAME resolved backend is
%                        forwarded to invz_jq_path (absent => unchanged brute-force default
%                        throughout). .dipole is 'bruteforce' | 'ewald'; .ewald is a complete
%                        {alpha,r_cut,g_cut,boundary} struct, required only with .dipole =
%                        'ewald' (docs/invzp_ewald_integration_map.md Sec.5A/6).
%     .gridConvention, .gridOffset, .gammaPolicy   opt-in BZ quadrature policy, forwarded BY
%                        PRESENCE into invz_bz_couplings when computed (any one present there
%                        switches to the invz_phase1_qgrid route). NEVER forwarded to
%                        invz_jq_path (the q-path has no BZ quadrature of its own).
%     .cache (true)      invz_jq_modes file cache; forwarded BY PRESENCE when computed.
%   The hyperfine manifold is always included; there is deliberately no .hyp option here.
%
%   Precomputed .Jnu/.info + backend/grid strictness: identical rule table to invz_spectra_map
%   (see there for the full contract and error identifiers). ADDITIONALLY here, whenever the
%   (precomputed or computed) BZ info carries complete info.dipole provenance, the q-path's OWN
%   resolved P.dipole (from invz_jq_path) must be isequaln to it -- a mismatch errors
%   invz:spectraPathDipoleMismatch, never a silent fall-through. S.path_dipole = P.dipole is
%   always returned, whether or not that check applied.

if nargin < 6, opts = struct(); end
grid    = getf(opts, 'grid', [16 16 16]);
dpRng   = getf(opts, 'dpRng', 30);
eta     = getf(opts, 'eta', 5e-3);
branch  = getf(opts, 'branch', 0);   % 0 = uniform FM mode (default); 1..4 = sorted branch
snapfac = getf(opts, 'snapfac', 2.5);
wmin    = getf(opts, 'peak_wmin', 0.05);

bztol = getf(opts, 'bz_tol', 1e-9);
sxtra = getf(opts, 'solve_opts', struct());
invz_check_solve_opts(sxtra);
B = invz_field_vec(B);
if abs(B(3)) <= bztol, B(3) = 0; end             % same dead band as invz_solve_auto

tmf = invz_check_transverse_mf(sxtra, B(2));

w = w(:);

% Ewald Step-5 Task 7 (docs/invzp_ewald_integration_map.md Sec.5A): backend/grid options are
% forwarded BY PRESENCE into invz_bz_couplings on the compute branch, and a precomputed
% opts.Jnu/opts.info pair is validated against any EXPLICIT backend/grid-policy request rather
% than trusted blindly -- see invz_check_coupling_opts.m (shared with invz_spectra_map.m,
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
Jaa0 = ion.Jxx0;  if isfield(info, 'Jaa0'), Jaa0 = info.Jaa0; end

% one medium solve at (T, B) -- FM below the (bare-MF) boundary, PM above
sopts = sxtra;
sopts.hyp = true;  sopts.J0eff = Jcc0;  sopts.Jxx0 = Jaa0;  sopts.bz_tol = bztol;
[pt, phase] = invz_solve_auto(ion, T, B, Jnu, sopts);
if phase == 0
    error('invz:noSolution', ...
        ['No converged 1/z solution at T = %.3f K, B = %s T ' ...
         '(near-degenerate doublet, critical band, or non-converged moment branch).'], ...
        T, mat2str(B, 4));
end

% guarded coupling along the path: physical uniform FM mode by default (P.Juni), or an
% exploratory sorted eigenvalue branch when opts.branch is 1..4. Resolve ONE dipolar backend
% for the path (Ewald Step-5 Task 7, docs/invzp_ewald_integration_map.md Sec.5A): an explicit
% spectra-level request is forwarded verbatim (already checked against info above when
% precomputed); otherwise a complete BZ info.dipole naming the Ewald backend is INHERITED so
% the path never silently falls back to brute force under it; an (inherited-or-default)
% bruteforce backend, or a provenance-less legacy info, is left unforwarded here -- bit-
% identical to the pre-Task-7 call. Grid convention/offset/gammaPolicy are BZ-quadrature-only
% and are NEVER forwarded to invz_jq_path.
pathOpts = struct('dpRng', dpRng, 'cache', true, 'snapfac', snapfac);
if hasBackendReq
    pathOpts.dipole = backendReq;
    if strcmp(backendReq, 'ewald'), pathOpts.ewald = eoptsReq; end
elseif chk.has_complete_dipole_provenance(info) && strcmp(info.dipole.backend, 'ewald')
    pathOpts.dipole = 'ewald';
    pathOpts.ewald  = info.dipole.ewald;
end
P  = invz_jq_path(ion, qpath, pathOpts);
if chk.has_complete_dipole_provenance(info) && ~isequaln(P.dipole, info.dipole)
    error('invz:spectraPathDipoleMismatch', ['invz_jq_path resolved a dipole backend/provenance ' ...
        'that does not match the BZ coupling info.dipole provenance; this must never be silently ' ...
        'mixed (compare S.path_dipole against S.info.dipole).']);
end
if branch == 0
    Jq = P.Juni(:).';                 % uniform ferromagnetic-mode coupling v'*Jcc*v
else
    Jq = P.Jnu(:, branch).';          % exploratory: one ascending-sorted branch
end

% path spectra from ONE real-axis evaluation each: Jsel is vectorized over q, and the
% single-site chit(w) is q-independent. Intrinsic: no Jshape here (see header).
% transverse_mf + si: without these, invz_chi_realaxis's strict-PM (phase==2) fallback
% rebuilds the single-ion state at the default 'legacy_x' MF model regardless of what
% the solve above actually used, breaking C4 between a field and its +90 deg rotation
% (review finding C1). si is redundant with the ordered-phase branch (which already
% resolves to pt.si on its own) but harmless there, so pass both unconditionally.
copts = struct('Jsel', Jq, 'eta', eta, 'Jxx0', Jaa0, 'hyp', true, ...
               'transverse_mf', tmf, 'si', pt.si);
o = invz_chi_realaxis(ion, T, B, pt, w, copts);
chiz = imag(o.chi_cc_q).';                        % [nw x nq]

% Under a longitudinal B (|Bz| > bz_tol) invz_solve_auto returns phase 1 or the error
% above already fired -- the strict-PM else-branch below is only reachable transversely.
if phase == 1
    pt0 = invz_zero_sigma_overlay(pt);
else
    tl0 = invz_twolevel(ion, T, B, struct('Jxx0', Jaa0, 'transverse_mf', tmf));
    pt0 = struct('alpha', 0, 'lambda', [0; 0], 'tl', tl0, 'K', []);
end
c0opts = copts;  c0opts.npass = 1;  c0opts.chi0cc_w = o.chi0cc_w;   % share bare cc from the 1/z call
o0 = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
chirpa = imag(o0.chi_cc_q).';

S = struct();
S.qpath = qpath;  S.s = P.s;  S.s_cart = P.s_cart;  S.snapped = P.snapped(:).';
% Plot coordinate S.x: when the path varies along exactly ONE Miller axis (monotonically),
% use that actual component -- e.g. h = 1..2 for the R2007 Fig-3 path -- so the axis shows
% real q positions. (The distance-from-start S.s reads 0..1 for EVERY unit-length window,
% hiding where in the zone the path sits.) Multi-axis or non-monotonic paths fall back to S.s.
span = max(qpath, [], 1) - min(qpath, [], 1);
vary = find(span > 1e-12);
if numel(vary) == 1 && (all(diff(qpath(:, vary)) > 0) || all(diff(qpath(:, vary)) < 0))
    S.x = qpath(:, vary).';
    lbl = {'h', 'k', 'l'};  parts = cell(1, 3);
    for c = 1:3
        if c == vary, parts{c} = lbl{c}; else, parts{c} = sprintf('%g', qpath(1, c)); end
    end
    S.xlab = sprintf('Q = (%s, %s, %s) (r.l.u.)', parts{:});
else
    S.x = P.s;
    S.xlab = sprintf('s along path from Q = [%g %g %g] (index r.l.u.)', qpath(1, :));
end
S.w = w;  S.T = T;  S.Bvec = B;  S.Bmag = norm(B);  S.phase = phase;  S.info = info;  S.Jq = Jq;
S.transverse_mf = tmf;  S.path_dipole = P.dipole;
S.chiz = chiz;  S.chirpa = chirpa;
S.Epeak     = invz_peak_energy(chiz,   w, wmin);
S.Epeak_rpa = invz_peak_energy(chirpa, w, wmin);
end

% Ewald Step-5 Task 7: precomputed-coupling provenance/conflict validation now lives in the
% shared invz_check_coupling_opts.m helper (task-7 review dedup fix), used identically by
% invz_spectra_map.m.
