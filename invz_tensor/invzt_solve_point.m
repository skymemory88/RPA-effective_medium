function pt = invzt_solve_point(ion, T, B, lat, opts)
%INVZT_SOLVE_POINT A1 projected-1/z bridge: self-consistent tensor point solve.
%   pt = INVZT_SOLVE_POINT(ion, T, B, lat, opts) solves the paramagnetic 1/z
%   self-energy fixed point at one (T, B) point using the FULL 12x12 tensor RPA
%   (3 Cartesian x 4 sublattice) as the effective-medium lattice engine. This is
%   the A1 dominant-sector bridge (docs/superpowers/plans/2026-07-17-invz-tensor-
%   full.md, Task 6): the tensor RPA carries the transverse (ODD) mediation
%   automatically and retarded, so A1 needs NO invz_chiperp / invz_odd_deltaJ.
%
%   B : scalar (transverse-along-a, historical) or [Bx By Bz] in Tesla.
%   lat : the LOCKED lattice struct from INVZT_JQ_TENSOR (fields Jt [12,12,nq],
%         qvec, w, conv, JtGamma [12,12], info with Jaa0/Jcc0/...). Production
%         builds lat at parts = 'full', so lat.JtGamma is the full-physics Gamma
%         matrix; its ODD (c<->a,b) blocks vanish by C2 symmetry (asserted).
%
%   OPTIONS (getf defaults):
%     mode          'a1'      : A1 dominant-sector bridge (T9/T12 add 'a2'/'a3').
%     Ecut          40        : Matsubara cutoff (meV) for invz_matsubara.
%     hyp           true      : include nuclear I = 7/2 in the single-ion solve.
%     transverse_mf 'legacy_x': transverse mean-field model (forwarded).
%     mix_outer     0.7       : base damped-mixing factor (Anderson depth-1 step).
%     tol_outer     1e-8      : outer-loop convergence tol (max|dSigma|).
%     max_outer     120       : outer-loop iteration cap. Raised from the projected
%                               solver's 60: the tensor A1 map is stiffer near the
%                               ORDERED-side near-singular RPA (a near-Gamma q whose
%                               cc mode eigenvalue approaches the Jcc0 FM resonance).
%                               Well-conditioned (PM-side / legacy-grid) points still
%                               converge in < 20 iterations; only a deep-ordered
%                               fine-grid point (e.g. 0.1 T, no-ODD, 8^3 halfopen)
%                               needs ~72, where plain damped mixing (the projected
%                               reference) does not converge at all.
%     anderson_depth 5        : Anderson history depth (1 == plain damped mixing).
%     Sigma_seed    []        : warm-start seed for Sigma (same length as wn).
%     odd           true      : false zeroes the Cartesian-OFF-diagonal entries of
%                               EVERY sublattice block of a COPY of lat.Jt (the cc
%                               sector then no longer sees the transverse blocks).
%     chi_rest      true      : false drops the crest (non-dominant) part of the
%                               local chi0 -- only the dominant (ground-manifold)
%                               sector is renormalized/mediated.
%     Esplit        0.4653    : dominant/rest split energy (meV), passed to split.
%     chi0_diag     false     : TEST HOOK. Zeroes the cross-Cartesian elements of
%                               the local tensor ctil (and ctil0) before use, so
%                               with odd=false the cc sector EXACTLY decouples
%                               (enables exact-identity gates).
%
%   Jaa0 THREADING (v2/v3, LOCKED): the single-ion transverse mean field and the
%   two-level params receive Jxx0 = lat.info.Jaa0 (NOT the ion.Jxx0 fallback);
%   pt.Jxx0_used records it. This makes projected-interop parity apples-to-apples
%   (pass the same Jxx0 = info.Jaa0 to the projected leg).
%
%   ZERO-FIELD GUARD (v3): mode 'a1' REQUIRES a symmetry-breaking TRANSVERSE
%   field -- invz_twolevel throws invz:degenerateDoublet below 1e-4 meV splitting.
%   The guard is on the TRANSVERSE component, NOT norm(B): a purely longitudinal
%   [0 0 Bz] has nonzero norm but does not split the doublet. For
%   transverse_mf = 'legacy_x' the guard is abs(B(1)) < 1e-6 T; a vector
%   transverse mode guards hypot(B(1), B(2)) < 1e-6 T. Too-small transverse field
%   errors invzt:a1ZeroField (true zero-field physics belongs to A3 / the
%   projected closed form).
%
%   A1 K-BOOKKEEPING (LOCKED, framework 11.8). Per Matsubara frequency n:
%     ctil(:,:,n) = cdom(:,:,n)/(1+Sigma(n)) + crest(:,:,n)   (renormalized chi0;
%                   only the dominant transition group carries the 1+Sigma factor,
%                   Jensen's "dominant renormalized, weak kept at RPA" rule)
%     Gloc(n)  = -[weighted site-diagonal cc average of the 12x12 RPA built from
%                  ctil(:,:,n)]                                 (invzt_gcc_lattice)
%     G0til(n) = -( cdom_cc(n)/(1+Sigma(n)) + crest_cc(n) ) = -real(ctil_cc(n))
%                                                              (site-local, no BZ)
%     K(n)     = 1./Gloc(n) - 1./G0til(n)
%
%   DECOUPLED-LIMIT REDUCTION (verify): treat the whole cc response as dominant
%   (crest_cc = 0) and Cartesian-diagonal (cc decoupled from the transverse
%   sector). Then G0til = -cdom_cc/(1+Sigma) = G0_scalar/(1+Sigma), so
%   1./G0til = (1+Sigma)./G0_scalar, and Gloc reduces to the scalar effective-
%   medium G. Hence
%     K = 1./Gloc - 1./G0til = 1./G - (1+Sigma)./G0
%   which is EXACTLY invz_emt_scalar's relation (G = G0./(1+Sigma+K.*G0) rearranged
%   gives 1./G = (1+Sigma)./G0 + K). The physical A1 no-ODD solve differs from the
%   scalar solve only by named residuals (cross-Cartesian chi0 elements, dominant-
%   mask vs whole-cc division), hence the 2e-3 frozen-baseline gate, not an exact
%   identity.
%
%   CRITICALITY (v3, Hermitian eigendecomposition -- NO sqrtm). With ctil0 the
%   static (n=1, elastic-on) renormalized 3x3 and C12 = kron(eye(4), ctil0) PSD:
%     [U, D] = eig((C12+C12')/2); d = real(diag(D)); d(d < 1e-12) = 0;
%     S = U*diag(sqrt(max(d,0)))*U';  M = eye(12) - S*JtGamma*S;
%     pt.crit = min(real(eig((M+M')/2)))
%   The symmetric active-subspace square root avoids sqrtm's complex/non-Hermitian
%   noise near criticality; M shares the zeros of I - C12*JtGamma on the active
%   subspace. crit > 0 in the PM phase. Clipped mass + active rank are recorded.
%
%   SUM RULE (report-quality): pt.sumrule_rel = |sum(wts.*G)/beta + si.JzJz_fluct|
%   / si.JzJz_fluct, with G = Gloc and the Matsubara cutoff Ecut.
%
%   Inside the ordered phase the paramagnetic fixed point may not exist; a
%   near-singular RPA leaves Gloc/K non-finite and pt.converged false (surfaced,
%   never silently NaN). Always check pt.converged.
%
%   See also INVZT_JQ_TENSOR, INVZT_GCC_LATTICE, INVZT_CHI0_SPLIT, INVZ_MATSUBARA,
%   INVZ_SINGLE_ION, INVZ_TWOLEVEL, INVZ_G, INVZ_LAMBDAS, INVZ_SIGMA,
%   INVZ_SOLVE_POINT (projected reference / interop parity target).
if nargin < 5, opts = struct(); end
mode  = getf(opts, 'mode', 'a1');
Ecut  = getf(opts, 'Ecut', 40);
hyp   = getf(opts, 'hyp', true);
tmf   = getf(opts, 'transverse_mf', 'legacy_x');
mixo  = getf(opts, 'mix_outer', 0.7);
tolo  = getf(opts, 'tol_outer', 1e-8);
maxo  = getf(opts, 'max_outer', 120);
Esplit    = getf(opts, 'Esplit', 0.4653);
chi_rest  = ~isfield(opts, 'chi_rest') || isempty(opts.chi_rest) || ~isequal(opts.chi_rest, false);
chi0_diag = isfield(opts, 'chi0_diag') && ~isempty(opts.chi0_diag) && ~isequal(opts.chi0_diag, false);
odd = ~isfield(opts, 'odd') || isempty(opts.odd) || ~isequal(opts.odd, false);

if ~(ischar(mode) || isstring(mode)) || ~strcmp(char(mode), 'a1')
    error('invzt:mode', ['invzt_solve_point currently implements mode ''a1'' only ' ...
        '(the A1 projected-1/z bridge; T9/T12 add ''a2''/''a3''); got %s.'], ...
        local_str(mode));
end
mode = char(mode);

B = invz_field_vec(B);                     % scalar -> [B 0 0]; 3-vector passes through

% --- zero-field guard on the TRANSVERSE component (v3), BEFORE any single-ion or
%     two-level solve so the error id is clear (not invz:degenerateDoublet) ------
if strcmp(tmf, 'legacy_x')
    trans = abs(B(1));
else
    trans = hypot(B(1), B(2));             % vector transverse mode
end
if trans < 1e-6
    error('invzt:a1ZeroField', ['mode ''a1'' requires a symmetry-breaking ' ...
        'TRANSVERSE field (transverse component %.3e T < 1e-6 T for ' ...
        'transverse_mf = ''%s''): invz_twolevel throws invz:degenerateDoublet ' ...
        'below 1e-4 meV doublet splitting. True zero-field physics belongs to ' ...
        'the projected closed form and to A3.'], trans, tmf);
end

Jxx0 = lat.info.Jaa0;                       % v3: single-ion / two-level get info.Jaa0

% --- (1) frequency grid, single-ion input, chi0 split, two-level response (once)
[wn, wts, beta] = invz_matsubara(T, Ecut);
nwn = numel(wn);
si  = invz_single_ion(ion, T, B, struct('hyp', hyp, 'transverse_mf', tmf, 'Jxx0', Jxx0));
[cdom, crest, mspec] = invzt_chi0_split(si, T, 1i*wn, struct('Esplit', Esplit));
tl  = invz_twolevel(ion, T, B, struct('Jxx0', Jxx0, 'transverse_mf', tmf));
g   = real(invz_g(tl, 1i*wn));

% cc (3,3) branches are real at imaginary Matsubara frequency (conjugate pairs
% cancel; the dominant mask is symmetric, so cdom/crest inherit the reality).
cdom_cc  = real(squeeze(cdom(3,3,:)));                   % [nwn,1]
crest_cc = real(squeeze(crest(3,3,:)));                  % [nwn,1]
chi0cc0  = cdom_cc(1) + crest_cc(1);                     % bare full local cc chi0(0)
if ~chi_rest
    crest    = zeros(size(crest));
    crest_cc = zeros(nwn, 1);
end

% --- odd = false: zero the Cartesian-off-diagonal entries of a COPY of lat.Jt.
%     cart = mod(i-1,3) is the Cartesian index of composite row/col i; keep only
%     same-Cartesian entries (all aa/bb/cc sublattice couplings, drop ca/cb/ab...).
lat_eff = lat;
if ~odd
    cart = mod((0:11).', 3);
    keepmask = (cart == cart.');                         % [12,12] logical
    lat_eff.Jt = lat.Jt .* keepmask;                     % broadcasts over pages
end

% --- C2 assertion: the Gamma-point ODD (c<->a,b) blocks vanish (< 1e-10*Jcc0) ---
JG = lat.JtGamma;
odd_ca = JG(3:3:12, 1:3:12);  odd_cb = JG(3:3:12, 2:3:12);
oddmag = max(abs([odd_ca(:); odd_cb(:)]));
if ~(oddmag < 1e-10*abs(lat.info.Jcc0))
    error('invzt:a1OddGamma', ['lat.JtGamma ODD (c<->a,b) blocks do not vanish ' ...
        '(max %.3e >= 1e-10*|Jcc0| = %.3e): C2 symmetry violated at Gamma.'], ...
        oddmag, 1e-10*abs(lat.info.Jcc0));
end

% --- (2) outer self-consistent Sigma loop: Anderson-accelerated damped mixing ---
%     The base update is the damped-mix step Sigma <- Sigma + mixo*(g(Sigma)-Sigma)
%     (Anderson depth 1); depth mAA > 1 adds history extrapolation. Plain damped
%     Picard is insufficient near the ORDERED-side near-singular RPA (a near-Gamma
%     q whose cc mode eigenvalue approaches Jcc0 -- the FM resonance): there the
%     projected reference solver also fails to converge (grid + ordered-phase
%     property, not a solver defect). Anderson recovers convergence there while
%     leaving well-conditioned (PM-side / legacy-grid) fixed points bit-unchanged.
mAA = getf(opts, 'anderson_depth', 5);
wstate = warning('off', 'MATLAB:rankDeficientMatrix');   % benign in the LS extrapolation
cleaner = onCleanup(@() warning(wstate));
Sigma = zeros(nwn, 1);
if isfield(opts, 'Sigma_seed') && numel(opts.Sigma_seed) == nwn
    Sigma = opts.Sigma_seed(:);
end
denom = @(s) reshape(1 + s, 1, 1, nwn);
converged = false;
Gloc = nan(nwn, 1);  K = nan(nwn, 1);  diag4 = nan(4, nwn);  lam = [NaN; NaN];
sg = struct('alpha', NaN, 'gamma', nan(nwn,1), 'Sigma', nan(nwn,1));
Fhist = cell(1, 0);  Xhist = cell(1, 0);
for outer = 1:maxo
    % --- medium step g(Sigma): renormalized chi0 -> lattice cc -> K -> Sigma' ---
    ctil = cdom ./ denom(Sigma) + crest;                % [3,3,nwn] renormalized chi0
    if chi0_diag
        for n = 1:nwn, ctil(:,:,n) = diag(diag(ctil(:,:,n))); end
    end
    [Gcc, diag4] = invzt_gcc_lattice(ctil, lat_eff);    % weighted site-diagonal cc
    Gloc  = -Gcc(:);                                     % site-local effective medium
    G0til = -(cdom_cc ./ (1 + Sigma) + crest_cc);        % = -real(ctil_cc)
    K = 1 ./ Gloc - 1 ./ G0til;                          % LOCKED K bookkeeping
    lam = invz_lambdas(K, g, wts, beta, [1 2]);
    sg  = invz_sigma(tl, lam, K, g, beta);
    f = sg.Sigma - Sigma;                               % fixed-point residual
    dS = max(abs(f));
    if dS < tolo, converged = true; break; end          % Gloc/K/lam/sg are at the fixed point
    % --- Anderson-accelerated update (depth 1 == plain damped mix) -------------
    Fhist{end+1} = f;  Xhist{end+1} = sg.Sigma;
    if numel(Fhist) > mAA, Fhist(1) = [];  Xhist(1) = []; end
    mk = numel(Fhist);
    if mk == 1
        Sigma = Sigma + mixo*f;
    else
        DF = zeros(nwn, mk-1);  DX = zeros(nwn, mk-1);
        for j = 1:mk-1, DF(:,j) = Fhist{j+1} - Fhist{j};  DX(:,j) = Xhist{j+1} - Xhist{j}; end
        gcoef = DF \ f;                                 % least-squares (min residual)
        Snew  = sg.Sigma - DX*gcoef;
        if any(~isfinite(Snew)), Snew = Sigma + mixo*f; end   % safeguard -> damped mix
        Sigma = Snew;
    end
end

% --- criticality (Hermitian eigendecomposition, rank-clipped PSD square root) ---
ctil0 = cdom(:,:,1) / (1 + Sigma(1)) + crest(:,:,1);     % static renormalized 3x3
if chi0_diag, ctil0 = diag(diag(ctil0)); end
C12  = kron(eye(4), ctil0);
[U, D] = eig((C12 + C12')/2);
d = real(diag(D));
rank_tol = 1e-12;
clip = d < rank_tol;
crit_clipped_mass = sum(abs(d(clip)));
crit_active_rank  = sum(~clip);
d(clip) = 0;
S = U * diag(sqrt(max(d, 0))) * U';
M = eye(size(S,1)) - S*lat.JtGamma*S;
crit = min(real(eig((M + M')/2)));

% --- assemble pt ----------------------------------------------------------------
pt.Sigma0 = Sigma(1);
pt.Sigma  = Sigma;
pt.alpha  = sg.alpha;
pt.lambda = lam;
pt.K = K;
pt.G = Gloc;                                             % G = Gloc (effective medium)
pt.tl = tl;
pt.si = si;
pt.chi0cc0 = chi0cc0;
pt.crit = crit;
pt.crit_clipped_mass = crit_clipped_mass;
pt.crit_active_rank  = crit_active_rank;
pt.sumrule_rel = abs(sum(wts.*Gloc)/beta + si.JzJz_fluct) / si.JzJz_fluct;
pt.converged = converged && all(isfinite(Sigma)) && all(isfinite(K)) && all(isfinite(Gloc));
pt.outer_iters = outer;
% REPORT: sublattice spread of the site-diagonal cc average (S4 symmetry breaking;
% only meaningful on a full symmetric BZ grid, ~0 there; noisy for explicit-q lat).
pt.diag4_spread = max(max(diag4, [], 1) - min(diag4, [], 1));
pt.mode = mode;
pt.odd  = odd;
pt.chi_rest = chi_rest;
pt.mspec = mspec;
pt.Jxx0_used = Jxx0;
% lat provenance (drop the bulk Jt pages; keep hash/conv/JtGamma/info/Jxx0_used).
pt.lat = struct('qvec_hash', hash_vec(lat.qvec(:)), 'conv', lat.conv, ...
    'JtGamma', lat.JtGamma, 'info', lat.info, 'Jxx0_used', Jxx0);
end

% ------------------------------- local helpers ----------------------------------

function h = hash_vec(v)
% Weak grid-identity hash, same formula as invz_cache_key / invzt_qgrid.
v = v(:);
h = sprintf('%dv_%08x', numel(v), typecast(single(sum(v.*(1:numel(v))')), 'uint32'));
end

function s = local_str(x)
if ischar(x) || (isstring(x) && isscalar(x))
    s = char(x);
else
    s = mat2str(x);
end
end
