function S = invz_spectra_qpath(ion, T, B, qpath, w, opts)
%INVZ_SPECTRA_QPATH Ferromagnetic-mode chi''_cc(q, omega) along a q-path at fixed (T, B).
%   The default dispersion follows the UNIFORM ferromagnetic-mode coupling J(q) = v'*Jcc*v
%   (invz_jq_path P.Juni) -- the same effective coupling MF_RPA_Yikai.m uses, whose energy
%   TRENDS reproduce R 2007 Fig 3 (mode softens monotonically along (1,0,0)->(2,0,0)). This
%   is NOT neutron scattering intensity (no eigenvector/form-factor weights); quantitative
%   reproduction is not claimed.
%
%   S = invz_spectra_qpath(ion, T, B, qpath, w) computes chi''_cc along qpath (nq x 3,
%   r.l.u.) at fixed T (K) and scalar transverse field B (T), for the 1/z and bare-RPA theories. The 1/z
%   medium (Sigma, K, lambda) is solved ONCE at (T, B) on the BZ-integration grid
%   (ordered-first via invz_solve_auto); the path susceptibility follows from the single-
%   site response chit(w) via chi(q, w) = chit/(1 - J(q) chit), J(q) from invz_jq_path.
%   The bare-RPA map has its own phase/state selection from the uniform RPA mass and does
%   not reuse the 1/z-selected single-ion state.
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
%     S.qpath, S.w, S.T, S.B, S.info, S.phase (1/z state) and S.phase_rpa
%                                     (independent bare-RPA state; 1=ordered, 2=paramagnet);
%     S.rpa_mass_pm                   bare-PM mass 1-Jcc0*chi0cc0 at this field
%
%   opts fields (all optional):
%     .grid ([16 16 16]), .dpRng (30), .eta (5e-3)   as in invz_spectra_map
%     .solve_opts        solver overrides; J0eff/Jxx0/hyp are driver-owned
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
%                        'ewald'.
%     .gridConvention, .gridOffset, .gammaPolicy   opt-in BZ quadrature policy, forwarded BY
%                        PRESENCE into invz_bz_couplings when computed (any one present there
%                        switches to the invz_phase1_qgrid route). NEVER forwarded to
%                        invz_jq_path (the q-path has no BZ quadrature of its own).
%     .cache (false)     optional invz_jq_modes file cache.
%   The hyperfine manifold is always included; there is deliberately no .hyp option here.
%
%   Precomputed .Jnu/.info + backend/grid strictness: identical rule table to invz_spectra_map
%   (see there for the full contract and error identifiers). ADDITIONALLY here, whenever the
%   (precomputed or computed) BZ info carries complete info.dipole provenance, the q-path's OWN
%   resolved P.dipole (from invz_jq_path) must be isequaln to it -- a mismatch errors
%   invz:spectraPathDipoleMismatch, never a silent fall-through. S.path_dipole = P.dipole is
%   always returned, whether or not that check applied.

if nargin < 6, opts = struct(); end
dpRng   = getf(opts, 'dpRng', 30);
eta     = getf(opts, 'eta', 5e-3);
branch  = getf(opts, 'branch', 0);   % 0 = uniform FM mode (default); 1..4 = sorted branch
snapfac = getf(opts, 'snapfac', 2.5);
wmin    = getf(opts, 'peak_wmin', 0.05);

sxtra = getf(opts, 'solve_opts', struct());
invz_check_solve_opts(sxtra);
if ~(isnumeric(B) && isreal(B) && isscalar(B) && isfinite(B) && B >= 0)
    error('invz:field', 'B must be a finite nonnegative scalar transverse field.');
end

w = w(:);

hasJnuOpt  = isfield(opts, 'Jnu');
hasInfoOpt = isfield(opts, 'info');
if hasJnuOpt ~= hasInfoOpt
    error('invz:spectraPrecomputedPartial', ...
          'opts.Jnu and opts.info must be supplied together (exactly one was present).');
end

if hasJnuOpt && hasInfoOpt
    Jnu = opts.Jnu(:);   info = opts.info;
else
    [Jnu, info] = invz_bz_couplings(ion, opts);
end
Jcc0 = info.Jcc0;
Jaa0 = ion.Jxx0;  if isfield(info, 'Jaa0'), Jaa0 = info.Jaa0; end

% One 1/z medium solve at (T, B); the RPA state is selected independently below.
sopts = sxtra;
sopts.hyp = true;  sopts.J0eff = Jcc0;  sopts.Jxx0 = Jaa0;
[pt, phase] = invz_solve_auto(ion, T, B, Jnu, sopts);
if phase == 0
    error('invz:noSolution', ...
        ['No converged 1/z solution at T = %.3f K, B = %.4g T ' ...
         '(near-degenerate doublet, critical band, or non-converged moment branch).'], ...
        T, B);
end

% guarded coupling along the path: physical uniform FM mode by default (P.Juni), or an
% exploratory sorted eigenvalue branch when opts.branch is 1..4. The path
% inherits an explicit backend, or Ewald provenance from the integration grid.
pathOpts = struct('dpRng', dpRng, 'cache', getf(opts, 'cache', false), 'snapfac', snapfac);
if isfield(opts, 'dipole')
    pathOpts.dipole = opts.dipole;
    if isfield(opts, 'ewald'), pathOpts.ewald = opts.ewald; end
elseif isfield(info, 'dipole') && isstruct(info.dipole) && ...
        isfield(info.dipole, 'backend') && strcmp(info.dipole.backend, 'ewald')
    pathOpts.dipole = 'ewald';
    pathOpts.ewald  = info.dipole.ewald;
end
P  = invz_jq_path(ion, qpath, pathOpts);
if branch == 0
    Jq = P.Juni(:).';                 % uniform ferromagnetic-mode coupling v'*Jcc*v
else
    Jq = P.Jnu(:, branch).';          % exploratory: one ascending-sorted branch
end

% path spectra from ONE real-axis evaluation each: Jsel is vectorized over q, and the
% single-site chit(w) is q-independent. Intrinsic: no Jshape here (see header).
copts = struct('Jsel', Jq, 'eta', eta, 'Jxx0', Jaa0, 'hyp', true, 'si', pt.si);
o = invz_chi_realaxis(ion, T, B, pt, w, copts);
chiz = imag(o.chi_cc_q).';                        % [nw x nq]

% Independent bare-RPA state: its uniform Jcc0 selects ordered versus PM, while
% the vector Jq remains the q-resolved denominator coupling.
c0opts = copts;
c0opts.J0eff = Jcc0;
c0opts.bare_rpa = true;
c0opts.npass = 1;
if isfield(c0opts, 'si'), c0opts = rmfield(c0opts, 'si'); end
o0 = invz_chi_realaxis(ion, T, B, [], w, c0opts);
chirpa = imag(o0.chi_cc_q).';

S = struct();
S.qpath = qpath;  S.s = P.s;  S.s_cart = P.s_cart;  S.snapped = P.snapped(:).';
% Plot coordinate S.x: when the path varies along exactly ONE Miller axis (monotonically),
% use that actual component -- e.g. h = 1..2 for the R2007 Fig-3 path -- so the axis shows
% real q positions. (The distance-from-start S.s reads 0..1 for EVERY unit-length window,
% hiding where in the zone the path sits.) Multi-axis or non-monotonic paths fall back to S.s.
span = max(qpath, [], 1) - min(qpath, [], 1);
vary = find(span > 1e-12);
if isscalar(vary) && (all(diff(qpath(:, vary)) > 0) || all(diff(qpath(:, vary)) < 0))
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
S.w = w;  S.T = T;  S.B = B;  S.phase = phase;
S.phase_rpa = o0.phase_rpa;  S.rpa_mass_pm = o0.rpa_mass_pm;
S.info = info;  S.Jq = Jq;
S.path_dipole = P.dipole;
S.chiz = chiz;  S.chirpa = chirpa;
S.Epeak     = invz_peak_energy(chiz,   w, wmin);
S.Epeak_rpa = invz_peak_energy(chirpa, w, wmin);
end
