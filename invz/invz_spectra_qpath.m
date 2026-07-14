function S = invz_spectra_qpath(ion, T, B, qpath, w, opts)
%INVZ_SPECTRA_QPATH Branch-resolved chi''_cc(q, omega) along a q-path at fixed (T, B).
%   EXPLORATORY: this is a branch susceptibility -- ONE sorted coupling eigenvalue per q
%   fed to the scalar pole formula. It is NOT neutron scattering intensity (no eigenvector/
%   sublattice-interference weights, no polarization factor, no magnetic form factor), and
%   it inherits the known open issues of the medium solve (closed-grid BZ quadrature, the
%   real-axis continuation's possible negative weight, and the bare-MF FM/PM handoff of
%   invz_solve_auto). Energy TRENDS are comparable with R 2007 Fig 3; quantitative
%   reproduction is not claimed.
%
%   S = invz_spectra_qpath(ion, T, B, qpath, w) computes chi''_cc along qpath (nq x 3,
%   r.l.u.) at fixed T (K) and transverse field B (T), for the 1/z and bare-RPA theories.
%   The 1/z medium (Sigma, K, lambda) is solved ONCE at (T, B) on the BZ-integration grid
%   (ordered-first via invz_solve_auto); the path susceptibility then follows from the same
%   single-site response chit(w) via chi(q, w) = chit/(1 - J_nu(q) chit), with J_nu(q) from
%   invz_jq_path (direction-aware Gamma-limit guard; see its header).
%
%   Demag semantics (canonical): the strict-uniform Jshape_cc transform is NOT applied
%   here -- a finite-q probe measures the intrinsic longitudinal response, and at Gamma-
%   equivalent path points the relevant limit is q -> 0+ where the demagnetizing field
%   cancels (R 2007). The spectra CAN still change with ion.demag through the transverse
%   applied/internal-field relation (demag-aware info.Jaa0 -> solver Jxx0).
%
%   Returns:
%     S.chiz, S.chirpa     [nw x nq]  1/z and bare-RPA chi''_cc(q, w)
%     S.Epeak, S.Epeak_rpa [1 x nq]   censored, parabolic-refined peak energy per q (NaN
%                                     when the maximum is non-positive/non-finite or sits
%                                     in the first/last usable bin -- i.e. the true peak
%                                     lies outside the sampled window)
%     S.Jq      [1 x nq]   selected coupling branch along the path (meV)
%     S.snapped [1 x nq]   true where invz_jq_path replaced the raw truncated sum
%     S.s       [1 x nq]   path distance in INDEX (r.l.u.) coordinates
%     S.s_cart  [1 x nq]   path distance in Cartesian reciprocal Ang^-1
%     S.qpath, S.w, S.T, S.B, S.info, S.phase (1 = FM, 2 = PM solve used)
%
%   opts fields (all optional):
%     .grid ([16 16 16]), .dpRng (30), .eta (5e-3)   as in invz_spectra_map
%     .branch (4)        which ascending-sorted coupling branch to follow (sorted index,
%                        NOT a tracked mode identity through crossings)
%     .snapfac (2.5)     Gamma-limit trust-radius factor (see invz_jq_path)
%     .peak_wmin (0.05)  meV; excludes the low-frequency hyperfine pole (R 2007 Fig 2)
%                        from the peak search so Epeak tracks the doublet mode
%     .Jnu, .info        precomputed BZ-grid branches / info (skips the lattice sum; tests)
%   The hyperfine manifold is always included (electronuclear solve + matching real-axis
%   state); there is deliberately no .hyp option on this API.

if nargin < 6, opts = struct(); end
grid    = getf(opts, 'grid', [16 16 16]);
dpRng   = getf(opts, 'dpRng', 30);
eta     = getf(opts, 'eta', 5e-3);
branch  = getf(opts, 'branch', 4);
snapfac = getf(opts, 'snapfac', 2.5);
wmin    = getf(opts, 'peak_wmin', 0.05);

w = w(:);

if isfield(opts, 'Jnu') && isfield(opts, 'info')
    Jnu = opts.Jnu(:);   info = opts.info;
else
    [qc, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', grid, 'range', [-0.5 0.5]);
    qc = qc(any(abs(qc) > 1e-12, 2), :);
    [Jnu, info] = invz_jq_modes(ion, qc, struct('dpRng', dpRng, 'cache', true));
    Jnu = Jnu(:);
end
Jcc0 = info.Jcc0;
Jaa0 = ion.Jxx0;  if isfield(info, 'Jaa0'), Jaa0 = info.Jaa0; end

% one medium solve at (T, B) -- FM below the (bare-MF) boundary, PM above
[pt, phase] = invz_solve_auto(ion, T, B, Jnu, struct('hyp', true, 'J0eff', Jcc0, 'Jxx0', Jaa0));
if phase == 0
    error('invz:noSolution', ...
        ['No converged 1/z solution at T = %.3f K, B = %.3f T ' ...
         '(near-degenerate doublet or critical band).'], T, B);
end

% guarded coupling branches along the path
P  = invz_jq_path(ion, qpath, struct('dpRng', dpRng, 'cache', true, 'snapfac', snapfac));
Jq = P.Jnu(:, branch).';

% path spectra from ONE real-axis evaluation each: Jsel is vectorized over q, and the
% single-site chit(w) is q-independent. Intrinsic: no Jshape here (see header).
copts = struct('Jsel', Jq, 'eta', eta, 'Jxx0', Jaa0, 'hyp', true);
o = invz_chi_realaxis(ion, T, B, pt, w, copts);
chiz = imag(o.chi_cc_q).';                        % [nw x nq]

if phase == 1
    pt0 = struct('alpha', 0, 'alpha_m', 0, 'lambda', [0; 0; 0], 'tl', pt.tl, ...
                 'K', [], 'is_ordered', true, 'si', pt.si);
else
    tl0 = invz_twolevel(ion, T, B, struct('Jxx0', Jaa0));
    pt0 = struct('alpha', 0, 'lambda', [0; 0], 'tl', tl0, 'K', []);
end
c0opts = copts;  c0opts.npass = 1;
o0 = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
chirpa = imag(o0.chi_cc_q).';

S = struct();
S.qpath = qpath;  S.s = P.s;  S.s_cart = P.s_cart;  S.snapped = P.snapped(:).';
S.w = w;  S.T = T;  S.B = B;  S.phase = phase;  S.info = info;  S.Jq = Jq;
S.chiz = chiz;  S.chirpa = chirpa;
S.Epeak     = peak_energy(chiz,   w, wmin);
S.Epeak_rpa = peak_energy(chirpa, w, wmin);
end

% -------------------------------------------------------------------------------------------
function E = peak_energy(chi, w, wmin)
%PEAK_ENERGY per-column peak of chi''(w) at w >= wmin, parabolic-refined; CENSORED (NaN)
% when the maximum is non-positive/non-finite or sits in the first/last usable bin (a
% boundary maximum means the true peak lies outside the sampled window -- reporting the
% grid edge would fabricate a flat dispersion). Assumes uniform w spacing.
E = nan(1, size(chi, 2));
mask = w >= wmin;
wm = w(mask);
if numel(wm) < 3, return; end
dw = wm(2) - wm(1);
for k = 1:size(chi, 2)
    c = chi(mask, k);
    [cmax, i] = max(c);
    if ~isfinite(cmax) || cmax <= 0 || i == 1 || i == numel(c), continue; end
    d = (c(i-1) - c(i+1)) / (2*(c(i-1) - 2*c(i) + c(i+1)));   % parabolic vertex offset
    if ~isfinite(d) || abs(d) > 1, d = 0; end
    E(k) = wm(i) + d*dw;
end
end

% -------------------------------------------------------------------------------------------
function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
