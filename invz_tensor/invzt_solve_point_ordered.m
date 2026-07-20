function pt = invzt_solve_point_ordered(ion, T, B, lat, opts)
%INVZT_SOLVE_POINT_ORDERED Tensor a1 solve at one FERROMAGNETIC (T, B) point.
%   pt = INVZT_SOLVE_POINT_ORDERED(ion, T, B, lat, opts) is the ordered-phase
%   counterpart of INVZT_SOLVE_POINT mode 'a1': the single-ion problem is solved
%   with the longitudinal ORDERING mean field (spontaneous m0 = <Jz>, hz =
%   J0z*<Jz> with J0z = lat.info.Jcc0 -- tensor-native threading, NOT ion.J0eff),
%   the self-energy uses the moment form (INVZ_SIGMA_ORDERED, HTML eqs 37-38,
%   lambda_{1,2,3}), and the medium step is the SAME tensor lattice engine as
%   the PM solver (INVZT_GCC_LATTICE 12x12 RPA + LOCKED K bookkeeping).
%
%   MEDIUM (2026-07-20 amendment, WHOLE-CC, no dominant/rest split): the local
%   chi0 fed to the lattice engine is the FULL electronuclear c0 = INVZ_CHI0Z(si,
%   T, i*wn, 'elastic', true), renormalized in one piece (ctil = c0./(1+Sigma)) --
%   mirroring the projected ordered reference (INVZ_SOLVE_POINT_ORDERED), which
%   also renormalizes the full G0, never a dominant-sector subset. The PM
%   solver's E < Esplit content cut (INVZT_CHI0_SPLIT) is NOT used here: measured
%   at T = 0.1 K, Bx = 3.0 T it selects only ndom = 13 states and carries a
%   DRIFTING ~48% of the local weight in "rest" (fdom_cc0 = 0.4762, review P0-3),
%   which made the no-ODD interop parity fail (dSigma0 = 8.6e-3, dalpha_m =
%   5.8e-2 at the T3 dispatch, both > 5e-3) when that rest weight was dropped by
%   chi_rest = false. A fixed content-defined dominant-vertex basis is a3d's
%   concern (Task 7D), not this a1 solver's -- so 'Esplit'/'chi_rest' are REJECTED
%   options here (invzt:orderedSplitKnobs), not silently accepted knobs.
%
%   SCOPE: TRANSVERSE fields only (spontaneous route; invzt:orderedLongitudinal
%   otherwise -- no forced_moment port, 2026-07-19 decision). Modes: 'a1' (moment
%   form, Task 3) and 'a3d' (full-response fixed-rank dominant-VERTEX hybrid,
%   Matsubara-only, Task 7D -- see the a3d branch below); full-dress 'a3' is
%   PERMANENTLY rejected (invzt:orderedMode, 136-state vertex budget-refused).
%   nlevels 'std' only (invzt:orderedNlevels). Option-A parity with the projected
%   ordered solver:
%   m0 is the BARE mean-field order parameter -- it onsets at the MF boundary
%   (~5.0 T at 0.1 K, MEASURED 2026-07-19), well above the tensor crit = 0 QPT
%   (4.65-4.70 T). A converged ordered point is therefore NOT evidence of being
%   in the FM phase; phase assignment belongs to INVZT_SOLVE_AUTO's
%   stability-based rule, never to this solver.
%
%   Mixing: plain damped mixing, CHECK-BEFORE-MIX (same loop ordering as the PM
%   solver, review P2-2): on every exit -- converged or not -- the returned
%   (Sigma, K, G, lambda, alpha, alpha_m) describe the SAME evaluated state.
%
%   OPTIONS (getf defaults): Ecut 40, hyp true, transverse_mf 'legacy_x',
%   mix_outer 0.7, tol_outer 1e-8, max_outer 80, m_tol 1e-2, bz_tol 1e-9,
%   odd true, rank_tol 1e-12, mz_seed/mf_maxit/mf_mix forwarded to
%   invz_single_ion. 'Esplit'/'chi_rest' (the PM dominant/rest split knobs) are
%   REJECTED here (invzt:orderedSplitKnobs) -- see MEDIUM above.
%
%   Returns the INVZT_SOLVE_POINT pt surface plus m0, alpha_m, is_ordered,
%   moment_branch ('spontaneous' | 'none'), J0z_used. Early returns (PM
%   relaxation |m0| <= m_tol; MF non-convergence) carry the projected-parity
%   fixed set: m0, is_ordered=false, converged=false, Sigma0=NaN, crit=NaN,
%   si, tl=[], moment_branch. Always check pt.is_ordered AND pt.converged.
%
%   MODE 'a3d' surface (round-2 P0-1, HONEST). The dense vertex does NOT produce the
%   Jensen moment-form objects, so alpha/alpha_m/lambda/Sigma/Sigma0 = NaN -- consumers
%   read the DIAGNOSTIC pt.Sigma_cc_equiv(1) EXPLICITLY. Extra fields: Vcc [nwn,1]
%   (cc vertex self-energy correction), chi_til [3,3,nwn] (the ONE hybrid response),
%   Kmat, Sigma_cc_equiv, eps_el/c_d (7A elastic-sector control), vb (16-state basis
%   diagnostics), a3d (vertex provenance + reeval seed). opts.reeval = <a3d pt> skips
%   the fixed point and evaluates the map ONCE at the returned state. a3d ignores the
%   moment-form knobs; it forwards mix_outer/tol_outer/max_outer/Vmat_seed/Lmax/tile_nl
%   and the invzt:orderedA3Budget guards to invzt_sigma_tensor.
%
%   See also INVZT_SOLVE_POINT, INVZ_SOLVE_POINT_ORDERED (projected reference),
%   INVZ_SIGMA_ORDERED, INVZ_TWOLEVEL_ORDERED, INVZT_GCC_LATTICE,
%   INVZT_CRIT_STATIC, INVZT_SOLVE_AUTO.
if nargin < 5, opts = struct(); end
mode   = getf(opts, 'mode', 'a1');
% Modes: 'a1' (moment form, Task 3) and 'a3d' (full-response dominant-vertex,
% Matsubara-only, Task 7D). 'a3' (full 136-state dress) stays rejected PERMANENTLY:
% the 136-state vertex is budget-refused (A4 ladder gate).
if ~ismember(char(mode), {'a1', 'a3d'})
    error('invzt:orderedMode', ['invzt_solve_point_ordered implements modes ''a1'' ' ...
        'and ''a3d'' -- got ''%s''. Full-dress ordered ''a3'' is permanently out of ' ...
        'scope: the 136-state vertex is budget-refused.'], char(mode));
end
Ecut   = getf(opts, 'Ecut', 40);
hyp    = getf(opts, 'hyp', true);
tmf    = getf(opts, 'transverse_mf', 'legacy_x');
mixo   = getf(opts, 'mix_outer', 0.7);
tolo   = getf(opts, 'tol_outer', 1e-8);
maxo   = getf(opts, 'max_outer', 80);
mtol   = getf(opts, 'm_tol', 1e-2);
bztol  = getf(opts, 'bz_tol', 1e-9);
rank_tol = getf(opts, 'rank_tol', 1e-12);
odd = ~isfield(opts, 'odd') || isempty(opts.odd) || ~isequal(opts.odd, false);
% 2026-07-20 amendment (P0-3 fallout): the ordered medium is WHOLE-CC (no
% dominant/rest split) -- 'Esplit'/'chi_rest' are the PM solver's split knobs
% and have no meaning here. Fail loud rather than silently ignore them.
if isfield(opts, 'Esplit') || isfield(opts, 'chi_rest')
    error('invzt:orderedSplitKnobs', ['invzt_solve_point_ordered renormalizes the ' ...
        'WHOLE local cc susceptibility (no dominant/rest split, 2026-07-20 amendment): ' ...
        '''Esplit''/''chi_rest'' are PM-solver-only knobs and are not accepted here. The ' ...
        'E < Esplit content cut drifts to a ~48%% ''rest'' weight at ordered 0.1 K points ' ...
        '(review P0-3), which is ill-defined this deep in the ordered phase; a fixed ' ...
        'content-defined dominant-vertex basis is a3d''s concern (Task 7D), not a1''s.']);
end
nlevels = getf(opts, 'nlevels', 'std');
if ~strcmp(char(nlevels), 'std')
    error('invzt:orderedNlevels', ['invzt_solve_point_ordered supports nlevels = ' ...
        '''std'' only (full electronuclear single ion); got ''%s''.'], char(nlevels));
end

B = invz_field_vec(B);
if abs(B(3)) > bztol
    error('invzt:orderedLongitudinal', ['invzt_solve_point_ordered is the ' ...
        'TRANSVERSE (spontaneous-moment) route only: got Bz = %.3g T > bz_tol = ' ...
        '%.3g T. The forced-moment longitudinal route is not ported to the ' ...
        'tensor branch (2026-07-19 scope decision).'], B(3), bztol);
end
B(3) = 0;                                            % dead band: exactly transverse

J0z  = lat.info.Jcc0;                                % tensor-native uniform cc coupling
Jxx0 = lat.info.Jaa0;                                % v3 threading (parity with PM solver)

[wn, wts, beta] = invz_matsubara(T, Ecut);
nwn = numel(wn);

% --- ordered single-ion mean-field fixed point (full electronuclear space) -------
siopts = struct('hyp', hyp, 'order', true, 'J0z', J0z, 'Jxx0', Jxx0, 'transverse_mf', tmf);
for f = {'mz_seed', 'mf_maxit', 'mf_mix'}
    if isfield(opts, f{1}), siopts.(f{1}) = opts.(f{1}); end
end
si = invz_single_ion(ion, T, B, siopts);
m0 = si.Jexp(3);
if ~si.mf_converged
    pt = early_return(m0, si, 'spontaneous');
    return;
end
if abs(m0) <= mtol
    pt = early_return(m0, si, 'none');               % paramagnetic point: use invzt_solve_point
    return;
end

% --- two-level driver at the converged ordering field hz; WHOLE-CC chi0 (2026-07-20:
%     no dominant/rest split -- mirrors the projected ordered reference's full G0) ---
tl = invz_twolevel_ordered(ion, T, B, si.hz, struct('Jxx0', Jxx0, 'transverse_mf', tmf));
g  = real(invz_g(tl, 1i*wn));
c0 = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
c0_cc = real(squeeze(c0(3,3,:)));
chi0cc0 = c0_cc(1);

% --- odd = false mask on a COPY of lat.Jt (same rule as the PM solver) -----------
lat_eff = lat;
if ~odd
    lat_eff.Jt = invzt_odd_mask(lat.Jt);
end

% --- C2 assertion at Gamma (same as PM solver) -----------------------------------
JG = lat.JtGamma;
odd_ca = JG(3:3:12, 1:3:12);  odd_cb = JG(3:3:12, 2:3:12);
oddmag = max(abs([odd_ca(:); odd_cb(:)]));
if ~(oddmag < 1e-10*abs(lat.info.Jcc0))
    error('invzt:a1OddGamma', ['lat.JtGamma ODD (c<->a,b) blocks do not vanish ' ...
        '(max %.3e >= 1e-10*|Jcc0| = %.3e): C2 symmetry violated at Gamma.'], ...
        oddmag, 1e-10*abs(lat.info.Jcc0));
end

if strcmp(char(mode), 'a3d')
    % ============================ a3d branch (Task 7D) ==============================
    %  Full-response, fixed-rank, FIELD-ADAPTED dominant-vertex hybrid, Matsubara-only.
    %  The V fixed point iterates on the COMPLETE HYBRID MAP (7D rework, controller
    %  decision 2026-07-20 -- the isolated-16-state map was REJECTED):
    %
    %      chi_dom_til(V) = dress[chi_dom, V]                   % dominant Dyson dressing (cc)
    %      chi_H(V)       = chi_full + (chi_dom_til(V) - chi_dom)  % dominant dressed, rest bare
    %      K_H(V)         = invzt_emt_matrix(chi_H(V), lat_eff)
    %      V_new(n)       = (1/2beta) sum_l K_H,cc(l) Gamma16_cc;cc(n,l)
    %
    %  The DEFINING INVARIANT: the returned Vcc is GENERATED BY the returned Kmat. The
    %  omitted (n_full - n_vertex) spectator response (incl. transverse/ODD channels)
    %  modifies K_H at the SAME vertex order -- var_share ~ 0.665-0.84 says that feedback
    %  is NOT an O(1/z^2) nicety -- so it must be inside the fixed point, not a post-hoc
    %  EMT of the converged truncated Vcc. Measured effect at the anchor: dVcc(0) ~ 28.5%,
    %  dcrit ~ -0.014 vs the rejected isolated map (test_a3d_map_choice_report).
    %
    %  MECHANIsm (minimal surgery). invzt_sigma_tensor is run on the reduced 16-state
    %  si_vb (so its G0/dressing/vertex are all the 16-state basis) with the EXPLICIT
    %  offset opts.chi_base = chi_full - chi_dom: its own V-loop then adds chi_base to the
    %  dressed 16-state response BEFORE its EMT step, i.e. it IS the complete map above.
    %  st.chi_til == chi_H, st.Kmat == K_H, st.Vmat == the self-generating V -- all from
    %  the SAME converged fixed point. chi_dom = invz_chi0z(si_vb) = -st.G0 exactly.
    %  (Passing FULL si instead would make st.chi_til whole-cc Dyson, NOT dominant-only --
    %  invzt_solve_point documents that at its dominant branch and avoids it.)
    %
    %  tl (built above) is NOT consumed by the a3d self-energy -- it is retained on pt
    %  only for interface parity with the a1 surface (and invzt_chi_realaxis rejection).

    % --- fixed-rank field-adapted 16-state vertex basis (7B) -----------------------
    vb    = invzt_ordered_vertex_basis(ion, T, si, struct());
    si_vb = vb.si_vertex;                          % reduced si-like struct (7B built it)

    % --- the hybrid ingredients: chi_full (all n_full), chi_dom (16-state), offset ---
    chi_full = c0;                                 % bare FULL response (elastic ON, above)
    chi_dom  = invz_chi0z(si_vb, T, 1i*wn, struct('elastic', true));   % bare dominant (= -st.G0)
    chi_base = chi_full - chi_dom;                 % [3,3,nwn] complete-map hybrid offset

    % --- vertex solve on the dominant basis, iterating the COMPLETE map via chi_base --
    st_opts = struct('dress', 'dominant', 'dom_basis', vb, 'rank_tol', rank_tol, ...
                     'chi_base', chi_base);
    for f = {'mix_outer', 'tol_outer', 'max_outer', 'Vmat_seed', 'tile_nl', ...
             'Lmax', 'max_vertex_states', 'max_vertex_work', 'max_vertex_bytes'}
        if isfield(opts, f{1}), st_opts.(f{1}) = opts.(f{1}); end
    end
    % reeval hook (gate 2): seed the converged Vmat and take ONE outer pass of the
    % COMPLETE map, so the returned Vcc/Kmat/chi_til are the map evaluated once at the
    % reeval state (dress + chi_base + EMT + contraction).
    reev = isfield(opts, 'reeval') && ~isempty(opts.reeval) && isstruct(opts.reeval);
    if reev
        rp = opts.reeval;
        if isfield(rp, 'a3d') && isstruct(rp.a3d) && isfield(rp.a3d, 'Vmat')
            st_opts.Vmat_seed = rp.a3d.Vmat;
        end
        st_opts.max_outer = 1;
    end
    st = invzt_sigma_tensor(si_vb, T, lat_eff, wn, beta, st_opts);

    % --- ALL theory objects come from the ONE converged complete-map fixed point -----
    chi_til = st.chi_til;                          % = chi_H (hybrid) [3,3,nwn]
    Kmat    = st.Kmat;                             % = K_H, the medium that GENERATED Vcc
    emtinfo = st.emtinfo;
    Gloc    = st.Gloc(:);                          % -real(chi_bar_cc) of the hybrid medium
    Kcc     = real(squeeze(Kmat(3, 3, :)));

    % --- crit from the RETURNED hybrid chi_til (herm-real static block) -------------
    [crit, crit_clipped_mass, crit_active_rank] = ...
        invzt_crit_static(herm_real(chi_til(:, :, 1)), lat.JtGamma, rank_tol);

    % --- Sigma_cc_equiv (DIAGNOSTIC): vertex self-energy vs the MATCHING dominant G0 -
    G0cc = squeeze(st.G0(3, 3, :));
    ok   = abs(G0cc) > rank_tol;
    Sigma_cc_equiv = nan(nwn, 1);
    Sigma_cc_equiv(ok) = squeeze(st.Vmat(3, 3, ok)) ./ G0cc(ok);   % +V/G0 (v3 POSITIVE)
    Vcc  = squeeze(st.Vmat(3, 3, :));  Vcc = Vcc(:);

    % --- eps_el / c_d: 7A elastic-sector control on the FULL ordered si -------------
    c0_el   = invz_chi0z(si, T, 1i*wn(1), struct('elastic', true));
    c0_inel = invz_chi0z(si, T, 1i*wn(1), struct('elastic', false));
    c_d     = real(c0_el(3, 3, 1) - c0_inel(3, 3, 1)) / beta;
    eps_el  = beta * abs(Kcc(1)) * c_d;

    phys_finite = all(isfinite(chi_til(:))) && all(isfinite(Vcc)) && all(isfinite(Gloc));
    sumrule_rel = abs(sum(wts.*Gloc)/beta + si.JzJz_fluct) / max(abs(si.JzJz_fluct), 1e-12);

    % --- assemble the HONEST a3d surface -------------------------------------------
    %  alpha/alpha_m/lambda/Sigma/Sigma0 are Jensen moment-form objects the dense
    %  vertex does NOT produce -- fabricating them is forbidden (round-2 P0-1). NaN.
    %  Consumers use pt.Sigma_cc_equiv(1) EXPLICITLY (it is a diagnostic, not Sigma0).
    pt = struct();
    pt.m0     = m0;
    pt.is_ordered = true;
    pt.moment_branch = 'spontaneous';
    pt.Sigma0 = NaN;                               % honest: a1 moment-form object only
    pt.Sigma  = nan(nwn, 1);
    pt.alpha  = NaN;
    pt.alpha_m = NaN;
    pt.lambda = nan(3, 1);
    pt.K = Kcc;                                     % hybrid-medium cc kernel (meaningful)
    pt.G = Gloc;
    pt.tl = tl;                                     % parity only (not used by a3d self-energy)
    pt.si = si;
    pt.chi0cc0 = real(chi_full(3, 3, 1));
    pt.crit = crit;
    pt.crit_clipped_mass = crit_clipped_mass;
    pt.crit_active_rank  = crit_active_rank;
    pt.sumrule_rel = sumrule_rel;
    pt.converged = st.converged && phys_finite;
    pt.outer_iters = st.outer_iters;
    pt.diag4_spread = a3d_diag4_spread(emtinfo);
    pt.mode = 'a3d';
    pt.odd  = odd;
    pt.chi_rest = true;                            % WHOLE response carried (parity)
    pt.mspec = [];
    pt.Jxx0_used = Jxx0;
    pt.J0z_used  = J0z;
    pt.nlevels = 'std';
    pt.lat = struct('qvec_hash', hash_vec(lat.qvec(:)), 'conv', lat.conv, ...
        'JtGamma', lat.JtGamma, 'info', lat.info, 'Jxx0_used', Jxx0);
    % --- a3d-specific surface (round-2 P0-1) ---------------------------------------
    pt.Vcc = Vcc;                                  % [nwn,1] cc vertex self-energy correction
    pt.chi_til = chi_til;                          % [3,3,nwn] the ONE hybrid response
    pt.Kmat = Kmat;                                % full non-Hermitian medium kernel
    pt.Sigma_cc_equiv = Sigma_cc_equiv;            % [nwn,1] DIAGNOSTIC (vs dominant G0)
    pt.eps_el = eps_el;
    pt.c_d    = c_d;
    pt.vb = vb;                                    % full basis-diagnostics struct
    pt.a3d = struct('Vmat', st.Vmat, 'G4cc', st.G4cc, 'Lmax', st.Lmax, ...
        'st_converged', st.converged, ...
        'st_outer_iters', st.outer_iters, 'st_active_rank', st.active_rank, ...
        'st_dress', st.dress, 'reeval', reev, ...
        'n_full', vb.n_full, 'n_vertex', vb.n_vertex, 'chi_share', vb.chi_share, ...
        'var_share', vb.var_share, 'p_mass', vb.p_mass, 'gap_kBT', vb.gap_kBT);
    return;
end

% --- outer moment-form loop: damped mixing, CHECK-BEFORE-MIX (P2-2) --------------
Sigma = zeros(nwn, 1);
denom = @(s) reshape(1 + s, 1, 1, nwn);
converged = false;
Gloc = nan(nwn, 1);  K = nan(nwn, 1);  diag4 = nan(4, nwn);
lam = nan(3, 1);  sg = struct('alpha', NaN, 'alpha_m', NaN, 'gamma', nan(nwn,1), 'Sigma', nan(nwn,1));
for outer = 1:maxo
    ctil = c0 ./ denom(Sigma);                       % [3,3,nwn] WHOLE-CC renormalized chi0
    [Gcc, diag4] = invzt_gcc_lattice(ctil, lat_eff);
    Gloc  = -Gcc(:);
    G0til = -(c0_cc ./ (1 + Sigma));
    K = 1 ./ Gloc - 1 ./ G0til;                      % LOCKED K bookkeeping (same as PM)
    lam = invz_lambdas(K, g, wts, beta, [1 2 3]);    % moment form needs lambda_3
    sg  = invz_sigma_ordered(tl, lam, K, g, beta);   % HTML eqs 37-38
    dS  = max(abs(sg.Sigma - Sigma));
    if dS < tolo, converged = true; break; end       % exit BEFORE mixing: state consistent
    Sigma = Sigma + mixo*(sg.Sigma - Sigma);
end

ctil0 = c0(:,:,1) / (1 + Sigma(1));
[crit, crit_clipped_mass, crit_active_rank] = invzt_crit_static(ctil0, lat.JtGamma, rank_tol);

% --- assemble pt (invzt_solve_point surface + ordered extras) --------------------
pt.m0     = m0;
pt.is_ordered = true;
pt.moment_branch = 'spontaneous';
pt.Sigma0 = Sigma(1);
pt.Sigma  = Sigma;
pt.alpha  = sg.alpha;
pt.alpha_m = sg.alpha_m;
pt.lambda = lam;
pt.K = K;
pt.G = Gloc;
pt.tl = tl;
pt.si = si;
pt.chi0cc0 = chi0cc0;
pt.crit = crit;
pt.crit_clipped_mass = crit_clipped_mass;
pt.crit_active_rank  = crit_active_rank;
pt.sumrule_rel = abs(sum(wts.*Gloc)/beta + si.JzJz_fluct) / max(abs(si.JzJz_fluct), 1e-12);
phys_finite = all(isfinite(Sigma)) && all(isfinite(K)) && all(isfinite(Gloc));
pt.converged = converged && phys_finite;
pt.outer_iters = outer;
pt.diag4_spread = max(max(diag4, [], 1) - min(diag4, [], 1));
pt.mode = 'a1';
pt.odd  = odd;
pt.chi_rest = true;   % interface parity with invzt_solve_point: WHOLE response always
                       % renormalized here (2026-07-20 amendment, no dominant/rest split)
pt.mspec = [];         % no chi0 split performed (mspec is a PM-solver-only diagnostic)
pt.Jxx0_used = Jxx0;
pt.J0z_used  = J0z;
pt.nlevels = 'std';
pt.lat = struct('qvec_hash', hash_vec(lat.qvec(:)), 'conv', lat.conv, ...
    'JtGamma', lat.JtGamma, 'info', lat.info, 'Jxx0_used', Jxx0);
end

% ------------------------------- local helpers ----------------------------------
function pt = early_return(m0, si, branch)
%EARLY_RETURN Complete field set for every non-accepted exit (projected parity).
pt = struct('m0', m0, 'is_ordered', false, 'converged', false, 'Sigma0', NaN, ...
            'crit', NaN, 'si', si, 'tl', [], 'moment_branch', branch);
end

function h = hash_vec(v)
% Weak grid-identity hash, same formula as invz_cache_key / invzt_solve_point.
v = v(:);
h = sprintf('%dv_%08x', numel(v), typecast(single(sum(v.*(1:numel(v))')), 'uint32'));
end

function A = herm_real(A)
% Real symmetric (Hermitian) part -- the static renormalized chi is Hermitian at
% wn = 0 (the gyrotropic anti-Hermitian part vanishes there); used for crit. Same
% definition as invzt_solve_point's local herm_real (the crit gate must match).
A = real((A + A')/2);
end

function s = a3d_diag4_spread(emtinfo)
% S4 sublattice spread of the cc effective medium (a3d analog of the a1
% diag4_spread from invzt_gcc_lattice). emtinfo.diag4_cc is the sublattice-resolved
% cc medium [4,nwn]; NaN if the EMT did not report it.
if isstruct(emtinfo) && isfield(emtinfo, 'diag4_cc') && ~isempty(emtinfo.diag4_cc)
    d = real(emtinfo.diag4_cc);
    s = max(max(d, [], 1) - min(d, [], 1));
else
    s = NaN;
end
end
