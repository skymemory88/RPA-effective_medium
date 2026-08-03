function pt = invzt_solve_point_ordered(ion, T, B, lat, opts)
%INVZT_SOLVE_POINT_ORDERED Tensor a1 solve at one FERROMAGNETIC (T, B) point.
%   pt = INVZT_SOLVE_POINT_ORDERED(ion, T, B, lat, opts) is the ordered-phase
%   counterpart of INVZT_SOLVE_POINT mode 'a1': the single-ion problem is solved
%   with the longitudinal ORDERING mean field (spontaneous m0 = <Jz>). Direct
%   calls retain the bare relation hz = J0z*<Jz>, J0z = lat.info.Jcc0. The
%   phase dispatcher additionally passes opts.hmf_J0z from its rejected PM
%   leg, activating the boundary-exact linearized Jensen modified field
%       hz = J0z*chi_tilde_PM(0)/chi_full_PM(0) * <Jz>.
%   That makes the ordered moment vanish at the renormalized PM instability,
%   rather than at the higher bare-MF boundary, and prevents a spectral jump at
%   the phase handoff. The nonzero scalar root is bracketed directly (fixed-hz
%   single-ion solves), avoiding critical slowing of Picard MF iteration.
%   the self-energy uses the moment form (INVZ_SIGMA_ORDERED, HTML eqs 37-38,
%   lambda_{1,2,3}), and the medium step is the SAME tensor lattice engine as
%   the PM solver (INVZT_GCC_LATTICE 12x12 RPA + LOCKED K bookkeeping).
%
%   MEDIUM (2026-07-20 amendment, WHOLE-CC, no dominant/rest split): the local
%   chi0 fed to the lattice engine is the FULL electronuclear c0 = INVZ_CHI0Z(si,
%   T, i*wn, 'elastic', true), renormalized in one piece (ctil = c0./(1+Sigma)) --
%   mirroring the projected ordered reference (INVZ_SOLVE_POINT_ORDERED), which
%   also renormalizes the full G0, never a dominant-sector subset. The legacy
%   diagnostic E < Esplit content cut (INVZT_CHI0_SPLIT) is NOT used here:
%   measured at T = 0.1 K, Bx = 3.0 T it selects only ndom = 13 states, with only
%   ~48% of the local weight in "dominant" and ~52% in "rest"
%   (fdom_cc0 = 0.4762, review P0-3),
%   which made the no-ODD interop parity fail (dSigma0 = 8.6e-3, dalpha_m =
%   5.8e-2 at the T3 dispatch, both > 5e-3) when that rest weight was dropped by
%   chi_rest = false. A fixed content-defined dominant-vertex basis is a3d's
%   concern (Task 7D), not this a1 solver's -- so 'Esplit'/'chi_rest' are REJECTED
%   options here (invzt:orderedSplitKnobs), not silently accepted knobs.
%   WP0.5 later measured the CURRENT fixed-rank selector at fdom_cc0=0.985--0.997
%   with <0.9% frozen split/whole static shifts. WP1 nevertheless identifies the
%   fixed-rank split as the framework-consistent scalar-A1 candidate; that proof
%   work is side-effect-free and is NOT wired into this production solver yet.
%
%   SCOPE: TRANSVERSE fields only (spontaneous route; invzt:orderedLongitudinal
%   otherwise -- no forced_moment port, 2026-07-19 decision). Modes: 'a1' (moment
%   form, Task 3) and 'a3d' (full-response fixed-rank dominant-VERTEX hybrid,
%   Matsubara-only, Task 7D -- see the a3d branch below); full-dress 'a3' is
%   PERMANENTLY rejected (invzt:orderedMode, 136-state vertex budget-refused).
%   nlevels 'std' only (invzt:orderedNlevels). Direct calls without
%   opts.hmf_J0z retain the historical bare-MF state for a3d/parity work.
%   a3d additionally accepts n_vertex (default 16), with max_vertex_states as
%   an independent workload guard, and evaluation_kind='one_step_reevaluation'
%   for an explicitly seeded diagnostic map evaluation.
%   INVZT_SOLVE_AUTO supplies hmf_J0z for production a1 phase-crossing runs.
%
%   Mixing: plain damped mixing, CHECK-BEFORE-MIX (same loop ordering as the PM
%   solver, review P2-2): on every exit -- converged or not -- the returned
%   (Sigma, K, G, lambda, alpha, alpha_m) describe the SAME evaluated state.
%
%   OPTIONS (getf defaults): Ecut 40, hyp true, transverse_mf 'legacy_x',
%   mix_outer 0.7, tol_outer 1e-8, max_outer 80, m_tol 1e-2, bz_tol 1e-9,
%   odd true, rank_tol 1e-12, mz_seed/mf_maxit/mf_mix forwarded to
%   invz_single_ion. hmf_J0z/hmf_sigma0 (normally internal, supplied by
%   INVZT_SOLVE_AUTO) activates the linearized modified-field relation above;
%   Sigma_seed may seed the moment-form fixed point. 'Esplit'/'chi_rest' (the
%   PM dominant/rest split knobs) are REJECTED here (invzt:orderedSplitKnobs).
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

% --- ordered single-ion mean-field state (full electronuclear space) -------------
% Auto-dispatched a1 points use the boundary-exact linearized modified field.
% Direct calls (including the expensive a3d validation path) retain the historical
% bare-MF fixed point unless hmf_sigma0 is explicitly supplied.
J0z_mf = J0z;
hmf_sigma0 = NaN;
hmf_residual = NaN;
if isfield(opts, 'hmf_J0z') && ~isempty(opts.hmf_J0z)
    if ~strcmp(char(mode), 'a1')
        error('invzt:hmfMode', 'opts.hmf_sigma0 is supported for ordered mode ''a1'' only.');
    end
    hmf_sigma0 = getf(opts, 'hmf_sigma0', NaN);
    if ~(isscalar(hmf_sigma0) && isreal(hmf_sigma0) ...
            && (isnan(hmf_sigma0) || (isfinite(hmf_sigma0) && (1 + hmf_sigma0) > 0)))
        error('invzt:hmfSigma0', ...
            'opts.hmf_sigma0 must be a finite real scalar with 1+Sigma0 > 0.');
    end
    J0z_mf = opts.hmf_J0z;
    if ~(isscalar(J0z_mf) && isreal(J0z_mf) && isfinite(J0z_mf) && J0z_mf > 0)
        error('invzt:hmfJ0z', 'opts.hmf_J0z must be a finite positive real scalar.');
    end
    [si, has_order, hmf_residual] = modified_mf_state( ...
        ion, T, B, hyp, Jxx0, tmf, J0z_mf, opts);
    if ~has_order
        pt = early_return(si.Jexp(3), si, 'none');
        pt.J0z_mf = J0z_mf;
        pt.hmf_sigma0 = hmf_sigma0;
        pt.hmf_residual = hmf_residual;
        return;
    end
else
    siopts = struct('hyp', hyp, 'order', true, 'J0z', J0z, ...
                    'Jxx0', Jxx0, 'transverse_mf', tmf);
    for f = {'mz_seed', 'mf_maxit', 'mf_mix'}
        if isfield(opts, f{1}), siopts.(f{1}) = opts.(f{1}); end
    end
    si = invz_single_ion(ion, T, B, siopts);
end
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
    %  The DEFINING INVARIANT: the returned Vnext is GENERATED BY the returned Kmat. The
    %  omitted (n_full - n_vertex) spectator response (incl. transverse/ODD channels)
    %  modifies K_H at the SAME vertex order, so it must be inside the fixed point,
    %  not a post-hoc EMT of the converged truncated Vcc.
    %
    %  MECHANIsm (minimal surgery). invzt_sigma_tensor is run on the reduced 16-state
    %  si_vb (so its G0/dressing/vertex are all the 16-state basis) with the EXPLICIT
    %  offset opts.chi_base = chi_full - chi_dom: its own V-loop then adds chi_base to the
    %  dressed 16-state response BEFORE its EMT step, i.e. it IS the complete map above.
    %  st.chi_til == chi_H and st.Kmat == K_H at st.Vmat; st.Vnext is generated from
    %  that SAME medium (and differs from st.Vmat only by the exposed residual at a
    %  converged fixed point). chi_dom = invz_chi0z(si_vb) = -st.G0 exactly.
    %  (Passing FULL si instead would make st.chi_til whole-cc Dyson, NOT dominant-only --
    %  invzt_solve_point documents that at its dominant branch and avoids it.)
    %
    %  tl (built above) is NOT consumed by the a3d self-energy -- it is retained on pt
    %  only for interface parity with the a1 surface (and invzt_chi_realaxis rejection).

    % --- fixed-rank field-adapted 16-state vertex basis (7B) -----------------------
    vb_opts = struct();
    if isfield(opts, 'n_vertex'), vb_opts.n_vertex = opts.n_vertex; end
    if isfield(opts, 'vb_prev'),  vb_opts.vb_prev  = opts.vb_prev;  end
    vb    = invzt_ordered_vertex_basis(ion, T, si, vb_opts);
    si_vb = vb.si_vertex;                          % reduced si-like struct (7B built it)

    % --- the hybrid ingredients: chi_full (all n_full), chi_dom (16-state), offset ---
    chi_full = c0;                                 % bare FULL response (elastic ON, above)
    chi_dom  = invz_chi0z(si_vb, T, 1i*wn, struct('elastic', true));   % bare dominant (= -st.G0)
    chi_base = chi_full - chi_dom;                 % [3,3,nwn] complete-map hybrid offset

    % --- vertex solve on the dominant basis, iterating the COMPLETE map via chi_base --
    st_opts = struct('dress', 'dominant', 'dom_basis', vb, 'rank_tol', rank_tol, ...
                     'chi_base', chi_base, 'requested_vertex_rank', vb.n_vertex);
    for f = {'mix_outer', 'tol_outer', 'max_outer', 'Vmat_seed', 'tile_nl', ...
             'Lmax', 'max_vertex_states', 'max_vertex_work', 'max_vertex_bytes', ...
             'anderson_depth', 'evaluation_kind'}
        if isfield(opts, f{1}), st_opts.(f{1}) = opts.(f{1}); end
    end
    % reeval hook (gate 2): seed the converged Vmat and take ONE outer pass of the
    % COMPLETE map. Returned Vmat/Kmat/chi_til share the reeval state and Vnext is the
    % generated output (dress + chi_base + EMT + contraction).
    reev = isfield(opts, 'reeval') && ~isempty(opts.reeval) && isstruct(opts.reeval);
    if reev
        rp = opts.reeval;
        if ~(isfield(rp, 'a3d') && isstruct(rp.a3d) && isfield(rp.a3d, 'Vmat'))
            error('invzt:a3dReevalSeed', ...
                'opts.reeval must contain a3d.Vmat from an A3d point.');
        end
        st_opts.Vmat_seed = rp.a3d.Vmat;
        st_opts.max_outer = 1;
        st_opts.evaluation_kind = 'one_step_reevaluation';
    end
    st = invzt_sigma_tensor(si_vb, T, lat_eff, wn, beta, st_opts);

    % --- ALL theory objects come from the ONE converged complete-map fixed point -----
    chi_til = st.chi_til;                          % = chi_H (hybrid) [3,3,nwn]
    Kmat    = st.Kmat;                             % = K_H, the medium that generated st.Vnext
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
    pt.map_converged = st.converged && phys_finite;
    pt.converged = pt.map_converged ...
        && ~strcmp(st.evaluation_kind, 'one_step_reevaluation');
    pt.outer_iters = st.outer_iters;
    pt.outer_residual = st.residual_inf;
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
    pt.Vmat = st.Vmat;                             % [3,3,nwn] evaluated A3d state
    pt.Vcc = Vcc;                                  % [nwn,1] evaluated cc component
    pt.Vnext = st.Vnext;                           % [3,3,nwn] generated from pt.Kmat
    pt.a3_residual = st.residual;
    pt.a3_residual_scaled = st.residual_scaled;
    pt.a3_residual_component_max = st.residual_component_max;
    pt.a3_residual_frequency_max = st.residual_frequency_max;
    pt.evaluation_kind = st.evaluation_kind;
    pt.chi_til = chi_til;                          % [3,3,nwn] the ONE hybrid response
    pt.Kmat = Kmat;                                % full non-Hermitian medium kernel
    pt.Sigma_cc_equiv = Sigma_cc_equiv;            % [nwn,1] DIAGNOSTIC (vs dominant G0)
    pt.eps_el = eps_el;
    pt.c_d    = c_d;
    pt.vb = vb;                                    % full basis-diagnostics struct
    pt.a3d = struct('Vmat', st.Vmat, 'Vnext', st.Vnext, ...
        'residual', st.residual, 'residual_inf', st.residual_inf, ...
        'residual_scaled', st.residual_scaled, ...
        'residual_component_max', st.residual_component_max, ...
        'residual_frequency_max', st.residual_frequency_max, ...
        'evaluation_kind', st.evaluation_kind, ...
        'terminal_co_state_consistent', st.terminal_co_state_consistent, ...
        'vertex', st.vertex, 'G4cc', st.G4cc, 'Lmax', st.Lmax, ...
        'st_converged', st.converged, ...
        'st_outer_iters', st.outer_iters, 'st_active_rank', st.active_rank, ...
        'st_dress', st.dress, 'reeval', reev, ...
        'n_full', vb.n_full, 'n_vertex', vb.n_vertex, 'chi_share', vb.chi_share, ...
        'var_share', vb.var_share, 'p_mass', vb.p_mass, 'gap_kBT', vb.gap_kBT);
    return;
end

% --- outer moment-form loop: Anderson-accelerated mixing, CHECK-BEFORE-MIX -------
Sigma = zeros(nwn, 1);
if isfield(opts, 'Sigma_seed') && ~isempty(opts.Sigma_seed)
    if numel(opts.Sigma_seed) ~= nwn || any(~isfinite(opts.Sigma_seed(:)))
        error('invzt:orderedSigmaSeed', ...
            'opts.Sigma_seed must contain %d finite entries.', nwn);
    end
    Sigma = real(opts.Sigma_seed(:));
end
denom = @(s) reshape(1 + s, 1, 1, nwn);
converged = false;
Gloc = nan(nwn, 1);  K = nan(nwn, 1);  diag4 = nan(4, nwn);
lam = nan(3, 1);  sg = struct('alpha', NaN, 'alpha_m', NaN, 'gamma', nan(nwn,1), 'Sigma', nan(nwn,1));
mAA = getf(opts, 'anderson_depth', 5);
Fhist = cell(1, 0);  Xhist = cell(1, 0);
wstate = warning('off', 'MATLAB:rankDeficientMatrix');
cleaner = onCleanup(@() warning(wstate));
dS = Inf;
for outer = 1:maxo
    ctil = c0 ./ denom(Sigma);                       % [3,3,nwn] WHOLE-CC renormalized chi0
    [Gcc, diag4] = invzt_gcc_lattice(ctil, lat_eff);
    Gloc  = -Gcc(:);
    G0til = -(c0_cc ./ (1 + Sigma));
    K = 1 ./ Gloc - 1 ./ G0til;                      % LOCKED K bookkeeping (same as PM)
    lam = invz_lambdas(K, g, wts, beta, [1 2 3]);    % moment form needs lambda_3
    sg  = invz_sigma_ordered(tl, lam, K, g, beta);   % HTML eqs 37-38
    f = sg.Sigma - Sigma;
    dS  = max(abs(f));
    if dS < tolo, converged = true; break; end       % exit BEFORE mixing: state consistent
    if outer == maxo, break; end                     % rejected exits retain one evaluated state
    Fhist{end+1} = f;  Xhist{end+1} = sg.Sigma; %#ok<AGROW>
    if numel(Fhist) > mAA, Fhist(1) = []; Xhist(1) = []; end
    mk = numel(Fhist);
    if mk == 1
        Sigma = Sigma + mixo*f;
    else
        DF = zeros(nwn, mk-1);  DX = zeros(nwn, mk-1);
        for j = 1:mk-1
            DF(:,j) = Fhist{j+1} - Fhist{j};
            DX(:,j) = Xhist{j+1} - Xhist{j};
        end
        gcoef = DF \ f;
        Snew = sg.Sigma - DX*gcoef;
        if any(~isfinite(Snew)), Snew = Sigma + mixo*f; end
        Sigma = Snew;
    end
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
pt.outer_residual = dS;
pt.diag4_spread = max(max(diag4, [], 1) - min(diag4, [], 1));
pt.mode = 'a1';
pt.odd  = odd;
pt.chi_rest = true;   % interface parity with invzt_solve_point: WHOLE response always
                       % renormalized here (2026-07-20 amendment, no dominant/rest split)
pt.mspec = [];         % no chi0 split performed (mspec is a PM-solver-only diagnostic)
pt.Jxx0_used = Jxx0;
pt.J0z_used  = J0z;
pt.J0z_mf = J0z_mf;
pt.hmf_sigma0 = hmf_sigma0;
pt.hmf_residual = hmf_residual;
pt.nlevels = 'std';
pt.lat = struct('qvec_hash', hash_vec(lat.qvec(:)), 'conv', lat.conv, ...
    'JtGamma', lat.JtGamma, 'info', lat.info, 'Jxx0_used', Jxx0);
end

% ------------------------------- local helpers ----------------------------------
function pt = early_return(m0, si, branch)
%EARLY_RETURN Complete field set for every non-accepted exit (projected parity).
pt = struct('m0', m0, 'is_ordered', false, 'converged', false, 'Sigma0', NaN, ...
            'crit', NaN, 'si', si, 'tl', [], 'moment_branch', branch, ...
            'outer_iters', 0, 'outer_residual', NaN);
end

function [si, has_order, resid] = modified_mf_state(ion, T, B, hyp, Jxx0, tmf, J0z_mf, opts)
%MODIFIED_MF_STATE Nonzero solution of hz = J0z_mf*<Jz>(hz).
% Fixed-hz single-ion solves remove the critical slowing of the ordinary
% hz-Picard loop. The slope test at hz=0 rejects the PM solution without ever
% allowing fzero to fall onto the symmetry-enforced trivial root.
siopts = struct('hyp', hyp, 'hz_fixed', 0, 'Jxx0', Jxx0, 'transverse_mf', tmf);
for f = {'mf_maxit', 'mf_mix'}
    if isfield(opts, f{1}), siopts.(f{1}) = opts.(f{1}); end
end
si0 = invz_single_ion(ion, T, B, siopts);
c00 = invz_chi0z(si0, T, 0, struct('elastic', true));
slope0 = 1 - J0z_mf*real(c00(3,3,1));
if ~si0.mf_converged || ~(slope0 < 0)
    si = si0;
    has_order = false;
    resid = 0;
    return;
end

% A positive nonzero root is bracketed entirely inside the physical saturated-
% moment bound.  Halving from the upper endpoint resolves an arbitrarily narrow
% negative interval near the continuous boundary without leaving that domain.
hhi = J0z_mf*max(abs(ion.J), 1);
fhi = mf_residual(hhi, ion, T, B, siopts, J0z_mf);
if ~(isfinite(fhi) && fhi > 0)
    error('invzt:hmfBracket', ...
        'Could not establish the physical upper modified-MF bracket (fhi=%.3g).', fhi);
end
b = hhi;
fb = fhi;
bracketed = false;
for k = 1:80
    a = 0.5*b;
    fa = mf_residual(a, ion, T, B, siopts, J0z_mf);
    if isfinite(fa) && fa < 0
        bracketed = true;
        break;
    end
    b = a;
    fb = fa;
end
if ~bracketed || ~(isfinite(fb) && fb >= 0)
    error('invzt:hmfBracket', ...
        'Could not bracket the nonzero modified-MF root inside the physical field bound.');
end

hz = fzero(@(h) mf_residual(h, ion, T, B, siopts, J0z_mf), [a b], ...
           optimset('TolX', 1e-12, 'Display', 'off'));
siopts.hz_fixed = hz;
si = invz_single_ion(ion, T, B, siopts);
resid = hz - J0z_mf*si.Jexp(3);
if ~(si.mf_converged && isfinite(resid) && abs(resid) < 1e-10)
    error('invzt:hmfNotConverged', ...
        'Modified-MF root residual %.6g meV failed the 1e-10 meV gate.', resid);
end
has_order = abs(si.Jexp(3)) > 0;
end

function r = mf_residual(hz, ion, T, B, siopts, J0z_mf)
siopts.hz_fixed = hz;
s = invz_single_ion(ion, T, B, siopts);
if ~s.mf_converged
    error('invzt:hmfMFNotConverged', ...
        'Fixed-hz transverse mean field failed at hz = %.6g meV.', hz);
end
r = hz - J0z_mf*s.Jexp(3);
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
