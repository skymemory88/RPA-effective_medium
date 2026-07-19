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

if ~(ischar(mode) || isstring(mode)) || ~ismember(char(mode), {'a1','a2','a3'})
    error('invzt:mode', ['invzt_solve_point implements modes ''a1'' (A1 projected-1/z ' ...
        'bridge), ''a2'' (A2 direct matrix effective medium) and ''a3'' (A3 genuine ' ...
        'tensor 1/z self-energy from the exact four-point vertex); got %s.'], invzt_str(mode));
end
mode = char(mode);
emt_rank_tol = getf(opts, 'rank_tol', 1e-12);

% nlevels routes the single-ion construction:
%   'three'          -> the explicit three-state toy (invzt_threestate) for ALL modes
%                       (a1/a2/a3), hyperfine EXCLUDED at this rung;
%   'eN' / 'eNxI8'   -> A4 state-space ladder rung (Task 13): the field/MF single-ion
%                       Hamiltonian is built and diagonalized IN the multiplet-complete
%                       reduced CF basis (invzt_rung_basis); 'xI8' tensors with the full
%                       I=7/2 nuclear space. The two-level driver tl comes from the
%                       ELECTRONIC ground doublet of the SAME reduced CF basis (hyp=false),
%                       mirroring invz_twolevel;
%   'std'  (default) -> the full electronuclear invz_single_ion.
nlevels    = getf(opts, 'nlevels', 'std');
nlch       = char(nlevels);
threestate = strcmp(nlch, 'three');
isrung     = ~threestate && ~strcmp(nlch, 'std') && ~isempty(regexp(nlch, '^e\d+(xI8)?$', 'once'));

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
if threestate
    % Explicit three-state toy (hyperfine EXCLUDED). The toy targets invz_twolevel at
    % the DEFAULT Jxx0 (a standalone electronic model), NOT lat.info.Jaa0 -- the doublet
    % match absorbs the real transverse mean field into Delta1. tl derives from the toy.
    ts_opts = struct('far_excited',    getf(opts, 'far_excited', false), ...
                     'chiperp_scale',  getf(opts, 'chiperp_scale', 1), ...
                     'chiperp_target', getf(opts, 'chiperp_target', 11.05), ...
                     'transverse_mf',  getf(opts, 'transverse_mf', 'none'), 'hz', 0);
    si = invzt_threestate(ion, T, B, ts_opts);
    tl = local_tl_from_si(si, T);
elseif isrung
    % A4 ladder rung: build the field/MF single-ion Hamiltonian and diagonalize it IN the
    % multiplet-complete reduced CF basis (Rayleigh-Ritz truncation). si is the reduced
    % electronuclear (xI8) or electronic response the A1/A2/A3 chain dresses. tl (the
    % two-level self-energy driver) is read off the ELECTRONIC ground doublet of the same
    % reduced CF basis (hyp=false), the rung analogue of invz_twolevel -- so as the basis
    % grows (e3 -> e6 -> e17 -> e17xI8 = 136) tl converges to the full electronic doublet.
    rb = invzt_rung_basis(ion, nlch);
    si = local_rung_si(ion, T, B, rb, rb.hyp, tmf, Jxx0);
    if rb.hyp
        rb_el = invzt_rung_basis(ion, rb.base_label);
        si_el = local_rung_si(ion, T, B, rb_el, false, tmf, Jxx0);
        tl = local_tl_from_si(si_el, T);
    else
        tl = local_tl_from_si(si, T);
    end
else
    si = invz_single_ion(ion, T, B, struct('hyp', hyp, 'transverse_mf', tmf, 'Jxx0', Jxx0));
    tl = invz_twolevel(ion, T, B, struct('Jxx0', Jxx0, 'transverse_mf', tmf));
end
[cdom, crest, mspec] = invzt_chi0_split(si, T, 1i*wn, struct('Esplit', Esplit));
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

% --- odd = false: zero the Cartesian-off-diagonal entries of a COPY of lat.Jt
%     (all aa/bb/cc sublattice couplings kept, ca/cb/ab... dropped). Shared
%     rule -- see INVZT_ODD_MASK for the exact semantics.
lat_eff = lat;
if ~odd
    lat_eff.Jt = invzt_odd_mask(lat.Jt);
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

% --- (2) outer self-consistent loop -----------------------------------------------
% A3 self-energy fields (populated only in mode 'a3'; NaN/[] otherwise).
Vmat = [];  chi_til = [];  Sigma_cc_equiv = nan(nwn, 1);
resum_spread_crit = NaN;  eps_el = NaN;  eps_cross = NaN;  a3_active_rank = NaN;
ctil0_add = [];  a3_dress = '';

if strcmp(mode, 'a3')
    % ===== A3: genuine tensor 1/z self-energy from the exact four-point vertex =====
    % invzt_sigma_tensor runs its OWN Vmat fixed-point (Gamma4 precomputed once); it
    % consumes pt.Kmat's full non-Hermitian K AS-IS (constraint 9, never symmetrized).
    st_opts = struct('rank_tol', emt_rank_tol);
    if isfield(opts, 'mix_outer'), st_opts.mix_outer = opts.mix_outer; end
    if isfield(opts, 'tol_outer'), st_opts.tol_outer = opts.tol_outer; end
    if isfield(opts, 'max_outer'), st_opts.max_outer = opts.max_outer; end
    if isfield(opts, 'Vmat_seed'), st_opts.Vmat_seed = opts.Vmat_seed; end
    if isfield(opts, 'dress'),     st_opts.dress     = opts.dress;     end  % 'full'|'dominant' (E1 match)
    st = invzt_sigma_tensor(si, T, lat_eff, wn, beta, st_opts);
    Vmat = st.Vmat;  chi_til = st.chi_til;  Kmat = st.Kmat;  emtinfo = st.emtinfo;
    Gloc = st.Gloc;  converged = st.converged;  outer = st.outer_iters;
    a3_active_rank = st.active_rank;  a3_dress = st.dress;
    diag4 = emtinfo.diag4_cc;
    K = real(squeeze(Kmat(3,3,:)));                     % cc medium kernel (report)
    G0cc = squeeze(st.G0(3,3,:));
    ok = abs(G0cc) > emt_rank_tol;
    Sigma_cc_equiv(ok) = squeeze(Vmat(3,3,ok)) ./ G0cc(ok);   % +V/G0 (v3 POSITIVE; DIAGNOSTIC)
    Sigma = real(Sigma_cc_equiv);
    lam = [NaN; NaN];  sg = struct('alpha', NaN, 'gamma', nan(nwn,1), 'Sigma', Sigma);
    if strcmp(a3_dress, 'dominant')
        % MATCHED truncation to E1: the crit's static renormalized chi uses the SAME
        % dominant-renormalized/rest-bare split as A1/A2 (cdom/(1+Sigma) + crest), driven
        % by A3's dominant self-energy -- NOT the whole-cc Dyson chi_til (which would
        % renormalize the non-dominant crest_cc too, a beyond-E1 effect). This makes
        % A3-dominant a faithful E1 truncation in BOTH the self-energy and the criticality.
        ctil0 = herm_real(cdom(:,:,1) / (1 + real(Sigma_cc_equiv(1))) + crest(:,:,1));
    else
        ctil0 = herm_real(chi_til(:,:,1));             % full A3: static Dyson renormalized 3x3
    end
    ctil0_add = herm_real(st.chi_til_add(:,:,1));       % static additive renormalized 3x3
    % monitors: elastic-sector control (constraint 7) and cross-Cartesian leakage
    eps_el    = beta * abs(K(1)) * si.JzJz_fluct;
    eps_cross = a3_eps_cross(Kmat, emt_rank_tol);
else
% --- a1/a2 outer self-consistent Sigma loop: Anderson-accelerated damped mixing ---
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
Kmat = [];  emtinfo = struct();                          % A2: matrix medium (mode 'a2' only)
sg = struct('alpha', NaN, 'gamma', nan(nwn,1), 'Sigma', nan(nwn,1));
Fhist = cell(1, 0);  Xhist = cell(1, 0);
for outer = 1:maxo
    % --- medium step g(Sigma): renormalized chi0 -> lattice cc -> K -> Sigma' ---
    ctil = cdom ./ denom(Sigma) + crest;                % [3,3,nwn] renormalized chi0
    if chi0_diag
        for n = 1:nwn, ctil(:,:,n) = diag(diag(ctil(:,:,n))); end
    end
    if strcmp(mode, 'a1')
        [Gcc, diag4] = invzt_gcc_lattice(ctil, lat_eff);    % weighted site-diagonal cc
        Gloc  = -Gcc(:);                                     % site-local effective medium
        G0til = -(cdom_cc ./ (1 + Sigma) + crest_cc);        % = -real(ctil_cc)
        K = 1 ./ Gloc - 1 ./ G0til;                          % LOCKED K bookkeeping
    else   % mode 'a2': A2 DIRECT matrix effective medium (K = ctil^-1 - chibar^-1)
        [Kmat, chibar, emtinfo] = invzt_emt_matrix(ctil, lat_eff, ...
            struct('rank_tol', emt_rank_tol));
        K = real(squeeze(Kmat(3,3,:)));                      % +K_cc (LOCKED POSITIVE, v2 sign)
        Gloc = -real(squeeze(chibar(3,3,:)));                % = -chi_bar_cc (= -Gcc)
        diag4 = emtinfo.diag4_cc;                            % per-sublattice cc medium (S4 report)
    end
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

    ctil0 = cdom(:,:,1) / (1 + Sigma(1)) + crest(:,:,1);     % static renormalized 3x3
    if chi0_diag, ctil0 = diag(diag(ctil0)); end
end                                                          % end mode branch

% --- criticality (Hermitian eigendecomposition, rank-clipped PSD square root) ---
% The static renormalized 3x3 ctil0 (A1/A2: cdom/(1+Sigma)+crest; A3: chi_til(:,:,1))
% shares the zeros of I - C12*JtGamma on its active subspace (C12 = kron(eye(4),ctil0)
% PSD). crit > 0 in the PM phase. A3 also reports the Cartesian-Dyson-vs-additive
% resummation SPREAD in crit (constraint 8: the O(1/z^2) method error bar).
rank_tol = 1e-12;
[crit, crit_clipped_mass, crit_active_rank] = invzt_crit_static(ctil0, lat.JtGamma, rank_tol);
if strcmp(mode, 'a3')
    crit_add = invzt_crit_static(ctil0_add, lat.JtGamma, rank_tol);
    resum_spread_crit = crit - crit_add;
end

% --- assemble pt ----------------------------------------------------------------
pt.Sigma0 = Sigma(1);
pt.Sigma  = Sigma;
pt.alpha  = sg.alpha;
pt.lambda = lam;
pt.K = K;
pt.Kmat = Kmat;                                          % A2: [3,3,nwn] matrix medium K (mode 'a2'; [] for 'a1'), for A3
pt.emt = emtinfo;                                         % A2: projector/persub_spread/diag4_cc/Herm diagnostics
pt.G = Gloc;                                             % G = Gloc (effective medium)
pt.tl = tl;
pt.si = si;
pt.chi0cc0 = chi0cc0;
pt.crit = crit;
pt.crit_clipped_mass = crit_clipped_mass;
pt.crit_active_rank  = crit_active_rank;
pt.sumrule_rel = abs(sum(wts.*Gloc)/beta + si.JzJz_fluct) / si.JzJz_fluct;
if strcmp(mode, 'a3')
    % Sigma_cc_equiv is a DIAGNOSTIC that is legitimately NaN where |G0_cc| <= rank_tol
    % (high-frequency slots where the bare cc propagator decays to ~0), so it must NOT
    % gate convergence -- the physical A3 objects (Vmat, chi_til, K, Gloc) do.
    phys_finite = all(isfinite(Vmat(:))) && all(isfinite(chi_til(:))) ...
        && all(isfinite(K)) && all(isfinite(Gloc));
else
    phys_finite = all(isfinite(Sigma)) && all(isfinite(K)) && all(isfinite(Gloc));
end
pt.converged = converged && phys_finite;
pt.outer_iters = outer;
% REPORT: sublattice spread of the site-diagonal cc average (S4 symmetry breaking;
% only meaningful on a full symmetric BZ grid, ~0 there; noisy for explicit-q lat).
pt.diag4_spread = max(max(diag4, [], 1) - min(diag4, [], 1));
pt.mode = mode;
pt.odd  = odd;
pt.chi_rest = chi_rest;
pt.mspec = mspec;
pt.Jxx0_used = Jxx0;
pt.nlevels = char(nlevels);
% --- A3-specific outputs (constraints 2/7/8) ------------------------------------
pt.Sigma_cc_equiv    = Sigma_cc_equiv;      % +V_cc/G0_cc (v3 POSITIVE sign; DIAGNOSTIC ONLY)
pt.Vmat              = Vmat;                 % [3,3,nwn] tensor self-energy correction (V=G0.Sigma)
pt.chi_til           = chi_til;             % [3,3,nwn] Dyson renormalized local chi
pt.resum_spread_crit = resum_spread_crit;   % crit(dyson) - crit(additive) -- O(1/z^2) error bar
pt.eps_el            = eps_el;              % beta*|K_cc(0)|*<dJz^2> elastic-sector control (constraint 7)
pt.eps_cross         = eps_cross;          % cross-Cartesian (c<->a,b) leakage in K, relative
pt.a3_active_rank    = a3_active_rank;      % resummation active-subspace rank
pt.a3_dress          = a3_dress;           % 'full' | 'dominant' (E1-matched truncation)
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


function A = herm_real(A)
% Real symmetric (Hermitian) part -- the static renormalized chi is Hermitian at
% wn = 0 (the gyrotropic anti-Hermitian part vanishes there); used for crit.
A = real((A + A')/2);
end

function e = a3_eps_cross(Kmat, tol)
% Cross-Cartesian (c<->a,b) leakage monitor: the largest c<->(a,b) medium-kernel
% magnitude relative to the largest |K_cc|, over all frequencies. ~0 with odd=false
% (block-diagonal K); O(ODD strength) with odd=true.
Kcc = max(abs(real(squeeze(Kmat(3,3,:)))));
cross = 0;
for pr = {[3 1],[3 2],[1 3],[2 3]}
    ij = pr{1};
    cross = max(cross, max(abs(squeeze(Kmat(ij(1), ij(2), :)))));
end
e = cross / max(Kcc, tol);
end

function tl = local_tl_from_si(si, T)
% Two-level params derived from a three-state toy's ground doublet (mirrors
% invz_twolevel's fields, but read off the toy's E/Mz rather than a fresh single-ion
% solve). Used by A1/A2 when nlevels = 'three'.
C = invz_const();
tl.Delta = si.E(2) - si.E(1);
if tl.Delta < 1e-4
    error('invz:degenerateDoublet', ...
        'Toy doublet splitting %.2e meV too small for the two-level self-energy.', tl.Delta);
end
tl.M2  = abs(si.Mz(1,2))^2;
tl.m   = real(si.Mz(1,1));
if abs(tl.m) > 1e-3
    error('invz:orderedPhase', 'Nonzero toy diagonal moment m=%.3g: outside paramagnetic scope.', tl.m);
end
tl.n01 = tanh(tl.Delta/(2*C.kB*T));
tl.g0  = 2*tl.n01/tl.Delta;
tl.transverse_mf = 'toy';
end

function si = local_rung_si(ion, T, B, rb, hyp, tmf, Jxx0)
% Build the field/MF single-ion Hamiltonian and diagonalize it IN the reduced rung basis
% rb.projector ("build the Hamiltonian in the reduced basis", Task 13). Mirrors
% invz_single_ion's transverse mean-field fixed point (legacy_x / vector_ab / none) with
% every operator projected P'*op*P first (Rayleigh-Ritz in the truncated CF [x nuclear]
% basis). Paramagnetic rung: hz = 0. Returns a COMPLETE si struct (same field surface as
% invz_single_ion), so it drops straight into invz_chi0z / invzt_chi0_split / the vertex.
%
% JzJz_fluct uses the REDUCED-basis intermediate sum sum_a p_a (Mz*Mz)_aa (Mz the reduced
% matrix elements), NOT the full P'*Jz^2*P: the sum rule is a consistency test of the
% reduced response, so the fluctuation it closes against must run over the SAME reduced
% intermediate states the reduced chi0 sums (they coincide at the full basis, where Mz is
% square-unitary). The missing full-vs-reduced fluctuation is the virtual content the
% driver's chi0_virtual_deficit diagnoses.
C = invz_const();
B = invz_field_vec(B);
oJ = stevens_ops(ion.J);
B44s = 0;  if isfield(ion, 'B44s'), B44s = ion.B44s; end
Hcf = ion.B20*oJ.O20 + ion.B40*oJ.O40 + ion.B44*oJ.O44 + B44s*oJ.O44s ...
    + ion.B60*oJ.O60 + ion.B64c*oJ.O64c + ion.B64s*oJ.O64s;
if hyp
    oI = stevens_ops(ion.I);  nI = size(oI.Jz, 1);
    kJ = @(M) kron(M, eye(nI));
    Hhf = ion.A*(kron(oJ.Jx, oI.Jx) + kron(oJ.Jy, oI.Jy) + kron(oJ.Jz, oI.Jz));
else
    kJ = @(M) M;  Hhf = 0;
end
Jx = kJ(oJ.Jx);  Jy = kJ(oJ.Jy);  Jz = kJ(oJ.Jz);
H0full = kJ(Hcf) + Hhf - ion.gL*C.muB*(B(1)*Jx + B(2)*Jy + B(3)*Jz);
% --- project every operator into the reduced rung basis -------------------------------
P = rb.projector;
H0  = P'*H0full*P;   H0  = (H0 + H0')/2;
Jxr = P'*Jx*P;       Jyr = P'*Jy*P;      Jzr = P'*Jz*P;
beta = 1/(C.kB*T);
vecmf  = strcmp(tmf, 'vector_ab');
nonemf = strcmp(tmf, 'none');
hx = 0;  hy = 0;  hz = 0;                        % paramagnetic rung: no longitudinal MF
converged = false;  it = 0;
for it = 1:200                                   % transverse mean-field fixed point (hx[,hy])
    H = H0 - hx*Jxr - hy*Jyr - hz*Jzr;  H = (H + H')/2;
    [Vr, Dr] = eig(H, 'vector');  [E, ixe] = sort(real(Dr));  Vr = Vr(:, ixe);
    p = exp(-beta*(E - E(1)));  p = p/sum(p);
    if nonemf
        hxn = 0;  hyn = 0;
    else
        hxn = Jxx0*(real(diag(Vr'*Jxr*Vr)).'*p);
        hyn = 0;
        if vecmf, hyn = Jxx0*(real(diag(Vr'*Jyr*Vr)).'*p); end
    end
    dmf = max([abs(hxn - hx), abs(hyn - hy)]);
    if dmf < 1e-12, hx = hxn;  hy = hyn;  converged = true;  break; end
    hx = hxn;  hy = hyn;
end
% recompute all fields ONCE from the converged (hx, hy)
H = H0 - hx*Jxr - hy*Jyr - hz*Jzr;  H = (H + H')/2;
[Vr, Dr] = eig(H, 'vector');  [E, ixe] = sort(real(Dr));  Vr = Vr(:, ixe);
p = exp(-beta*(E - E(1)));  p = p/sum(p);
si.E  = E - E(1);
si.V  = Vr;
si.P  = p;
si.Mx = Vr'*Jxr*Vr;  si.My = Vr'*Jyr*Vr;  si.Mz = Vr'*Jzr*Vr;
si.Jexp = [real(diag(si.Mx)).'*p; real(diag(si.My)).'*p; real(diag(si.Mz)).'*p];
si.hx = hx;  si.hy = hy;  si.hz = hz;
jz2 = real(diag(si.Mz*si.Mz)).'*p;               % reduced-basis intermediate sum (see header)
si.JzJz_fluct = jz2 - si.Jexp(3)^2;
si.mf_converged = converged;
si.mf_iters = it;
si.E0 = E(1);
si.transverse_mf = tmf;
si.rung = rb.rung;
si.dim_actual = rb.dim_actual;
end
