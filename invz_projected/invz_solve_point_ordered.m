function pt = invz_solve_point_ordered(ion, T, Bx, Jnu_flat, opts)
%INVZ_SOLVE_POINT_ORDERED Self-consistent 1/z solution at one FERROMAGNETIC (T, Bx) point.
% Bx: scalar (transverse, historical) or [Bx By Bz] vector (T).
% Ordered-phase counterpart of invz_solve_point: the single-ion problem is solved with the
% longitudinal ORDERING mean field (spontaneous moment m0 = <Jz>), the self-energy uses the
% elastic-sector form (invz_sigma_ordered, framework SS9.2, J 2.26-2.27), and the whole thing is iterated
% against the effective medium exactly as in the paramagnet.
%
% opts.forced_moment (logical, default false): when true (set by invz_solve_auto's
% longitudinal route) the moment is treated as FIELD-INDUCED rather than spontaneous -- the
% |m0| > m_tol gate is bypassed, a sign-aware seed/one mirrored retry enforces alignment with
% the applied Bz, and a non-converged mean-field loop is itself an early-return condition.
% forced_moment with Bx(3) = 0 skips the alignment check (no field sign to align to).
% opts.mz_seed / opts.mf_maxit / opts.mf_mix forward to invz_single_ion (diagnostics/tests).
%
% Returns pt.is_ordered: strictly "this point uses the moment-form self-energy", true for
% EITHER a spontaneous FM moment (|m0| > m_tol) or a forced_moment field-induced solve; FALSE
% on the spontaneous-mode paramagnetic early return AND on the forced-moment failure paths (MF non-convergence; persistent branch misalignment) -- acceptance always also requires pt.converged.
% When ordered, pt carries the same fields as
% invz_solve_point plus pt.m0 (order parameter), pt.alpha_m, pt.si (the ordered single-ion
% struct), pt.is_ordered = true, and pt.moment_branch ∈ {'spontaneous','field_induced','none'}
% ('none' only on the spontaneous-mode paramagnetic early return). EVERY return path (early
% or converged) carries the full field set: m0, is_ordered, converged, Sigma0, crit, si, tl,
% moment_branch -- tl = [] on early returns (no two-level params were built).
%
% SCOPE NOTE (Option A): for the spontaneous route, m0 is the bare mean-field order parameter
% (the applied-field/H_MF self-consistency (framework SS9.3, J 2.31-2.33) is deferred), so it onsets at the
% mean-field boundary, slightly above the true 1/z boundary; the gap matters only near B_c,
% not deep in the ordered phase.
%
% ODD extension (T1.4, opt-in): same contract as invz_solve_point — opts.odd = true
% requires opts.odd_blocks (Vca/Vcb/Vcc/Jcc0 UNSHIFTED, from invz_odd_blocks) and
% Jnu_flat = [] ('invz:oddArgs' otherwise). The cc modes are rebuilt from
% Vcc + deltaJ(T,Bx) and the explicit -d (E5) is applied ONCE to J0eff below, which
% here feeds BOTH the single-ion ordering mean field (siopts.J0z) and pt.crit: the
% shifted value is the physical uniform coupling, so both uses receive it. chi_perp
% is evaluated at the PARAMAGNETIC-MF single-ion state by design (T1.2: Van Vleck
% dominated, insensitive to the cc order parameter — never solved self-consistently
% with the moment). pt gains pt.odd = struct('d','Xp') on the full solve path only
% (early returns keep their fixed field set). Flag off: byte-identical pre-ODD path.
%
% Jensen ordered mode (Stage-2 task 4, opt-in): opts.ordered_mode = 'bare' (default) runs the
% SCOPE-NOTE bare mean-field order parameter above, unchanged. opts.ordered_mode = 'jensen'
% instead imposes the H_MF applied-field/self-consistency root (invz_hmf_ordered, framework
% SS9.3, J 2.31-2.33) as a FIXED longitudinal field (siopts.hz_fixed, order = false -- P0-1:
% the ordering update is NOT re-applied on top of the imposed root) and closes the final
% Sigma<->EMT loop with the SAME static-sector elastic closure used by that machinery's nodes
% (invz_emt_static_ordered), so pt.G(1) is the CLOSED ELASTIC static function -- NOT the
% ordinary-Dyson value invz_emt_scalar alone would produce at omega = 0. In jensen mode pt
% additionally carries pt.ordered_mode, pt.hmf (the refined root, meV), pt.hmf_prof (the
% invz_hmf_ordered profile), pt.D_uni = 1 + (J0eff - K(1))*pt.G(1) (the ordered static
% inverse response AT THE FINAL STATE -- the pole observable), and pt.final_resid (the
% closure-residual revalidation of the exported (Sigma, K, lambda) tuple, required <
% opts.tol_outer for pt.converged -- the tuple is thereby closure-consistent to the STATED
% OUTER TOLERANCE, not exactly self-consistent). pt.is_ordered is gated on H_MF ROOT
% EXISTENCE (isfinite(hmf_star)), NOT on m_tol, with a paramagnetic early return (existing
% struct shape, pt.hmf_status carrying the invz_hmf_ordered status) when no root exists.
% jensen mode is TRANSVERSE/SPONTANEOUS only: opts.forced_moment true, or an explicit
% longitudinal field |Bx(3)| > opts.bz_tol (default 1e-9), raises 'invz:orderedMode'.
% opts.cold_acceleration = 'signed_aitken1' enables the safeguarded, default-off resummed
% cold-start experiment in the shared node solver; absent/'none' preserves the ordinary
% iteration exactly. It changes no equation or final A--D/all-node acceptance gate.
% CONTRACT (P2-G): pt.crit KEEPS its historical ordinary-Dyson definition
% (1 + Sigma0 - J0eff*chi0cc0) as a legacy diagnostic in EITHER mode -- it is NOT the ordered
% pole mass below the boundary; pt.D_uni is. Callers must not conflate the two.
%
% Static medium (spec SS4.2, opt-in): opts.static_medium = 'resummed' (DEFAULT, legacy and
% numerically bit-identical) | 'strict_1z_dyson_ref' | 'strict_1z_bare_ref', with
% opts.ref_margin (1e-6) the reference-denominator floor. The scheme is resolved ONCE at this
% entry by invz_check_static_medium and stamped into BOTH legs (opts.emt for the dynamic slot,
% opts.emt_static for the ordered static sector) plus the invz_hmf_ordered context, so the two
% sectors can never run different truncation orders; setting it per leg is a CONFLICT
% ('invz:staticMedium'), never an override. Under a strict scheme a non-vector (retarded
% [nJ,nw]) Jnu_flat and opts.odd_retarded(_exact) are both rejected with 'invz:staticMedium':
% there is no retarded ordered path here, and silently solving the static multiset instead
% would be indistinguishable from a retarded solve in the output. pt ALWAYS carries
% pt.static_medium, pt.Jmom (invz_coupling_moments of this point's final spectrum),
% pt.medium_status, pt.medium_denom and pt.medium_margin (= denom - floor, the DISTANCE TO THE
% FLOOR, not the denominator) on EVERY return path, early ones included.
%
% TWO-TIER ACCEPTANCE in jensen mode (spec SS1): pt.converged KEEPS its existing meaning --
% info.accepted, the closure-consistency tier -- so no existing consumer shifts meaning under
% it. ENDPOINT STABILITY is the separate pt.stable_1z: the checker's own res.stability verdict
% (already carrying the frozen crit_tol/D_tol/Dq_tol) for the accepted root, AND agreement
% between that node's crit and the profile's crit_star. The two must not be collapsed --
% intermediate path nodes are the unstable Landau interval by construction, so folding
% stability into acceptance would re-mask the ordered phase. jensen mode additionally exports
% pt.crit_1z (hmf profile's crit_star), pt.Dq_min, the endpoint's omitted-term ratios
% pt.omit_mu3/.omit_cubic/.omit_max, and pt.path_omit_max -- the FAIL-CLOSED maximum of the
% solved path's omit_max (empty or any NaN node => NaN; Inf dominates), the quantity the frozen
% prereg's omit_promote gate reads.
% FIELD-PRESENCE CONTRACT ACROSS MODES: those jensen-only members -- stable_1z, crit_1z,
% Dq_min, D_uni, omit_mu3, omit_cubic, omit_max, path_omit_max -- are ABSENT ENTIRELY in bare
% mode (isfield false, not NaN-valued), so consumers must test isfield before reading them and
% must not infer bare-mode results from a NaN. By contrast static_medium, medium_status,
% medium_denom, medium_margin and Jmom are present in BOTH modes on EVERY return path.
if nargin < 5, opts = struct(); end
Ecut  = getf(opts, 'Ecut', 40);
hyp   = getf(opts, 'hyp', true);
J0eff = getf(opts, 'J0eff', ion.J0eff);
Jxx0  = getf(opts, 'Jxx0', ion.Jxx0);
tmf   = getf(opts, 'transverse_mf', 'legacy_x');
mixo  = getf(opts, 'mix_outer', 0.7);
tolo  = getf(opts, 'tol_outer', 1e-8);
maxo  = getf(opts, 'max_outer', 200);
cold_accel = getf(opts, 'cold_acceleration', 'none');
mtol  = getf(opts, 'm_tol', 1e-2);
eopts = getf(opts, 'emt', struct());
eso_pub = getf(opts, 'emt_static', struct());
% Static-medium scheme: resolved ONCE at this public entry by the sole authority and stamped
% into BOTH leg option structs (spec SS4.2), so the dynamic PM slot and the ordered static
% sector can never run different truncation orders. Absent => 'resummed' => numerically
% identical to the pre-strict path. Validation is idempotent by design, so the stamped structs
% are forwarded verbatim into invz_hmf_ordered, which resolves the scheme again.
[sm, eopts, eso_pub] = invz_check_static_medium(opts, eopts, eso_pub);
Bx = invz_field_vec(Bx);                       % scalar -> [Bx 0 0]; 3-vector passes through
fmom = getf(opts, 'forced_moment', false);

% --- ODD diversion (T1.4): strictly additive and opt-in; everything else in
% this function is the pre-ODD code path, byte-untouched when the flag is off.
oddOn = isfield(opts, 'odd') && ~isempty(opts.odd) && ~isequal(opts.odd, false);
if oddOn
    ob = getf(opts, 'odd_blocks', []);
    if ~(isstruct(ob) && isscalar(ob) && all(isfield(ob, {'Vca', 'Vcb', 'Vcc', 'Jcc0'})))   % Jcc0: contract/audit field, required for caller-side consistency; unused here (J0eff comes via opts)
        error('invz:oddArgs', ['opts.odd = true requires opts.odd_blocks = struct(' ...
            '''Vca'',''Vcb'',''Vcc'',''Jcc0'') precomputed once by the caller from ' ...
            'invz_odd_blocks (P0.4: no disk/cache reads inside solver loops; Jcc0 UNSHIFTED).']);
    end
    if ~isempty(Jnu_flat)
        error('invz:oddArgs', ['opts.odd = true requires Jnu_flat = []: the cc modes are ' ...
            'rebuilt here from odd_blocks + deltaJ, and a caller-supplied baseline Jnu ' ...
            'would silently override the rebuild.']);
    end
    % chi_perp at the shared single-ion option set (T1.2); its si is the
    % PARAMAGNETIC-MF state and is NOT reused here (the ordered solve below
    % needs the ordering-MF state, a different siopts).
    Xp = invz_chiperp(ion, T, Bx, struct('hyp', hyp, 'Jxx0', Jxx0, 'transverse_mf', tmf));
    [dJ, d] = invz_odd_deltaJ(ob.Vca, ob.Vcb, Xp);
    Jnu_odd = invz_odd_modes(ob.Vcc, dJ);      % values-only kernel (shared)
    Jnu_flat = Jnu_odd(:);
    % E5 uniform shift, applied HERE exactly once (T1.3 bookkeeping rule: the grid
    % matrices' diagonal already carries -d via E4; J0eff carries the explicit -d;
    % NO other q = 0 handling). Callers pass the UNSHIFTED info.Jcc0 as opts.J0eff.
    % That unshifted value comes from invz_odd_blocks' infoB.Jcc0 (or the flag-off
    % invz_jq_modes info.Jcc0) -- NOT from invz_jq_modes' opts.odd path, whose
    % exported info.Jcc0 is ALREADY shifted by -d.
    % The shifted J0eff is the physical uniform coupling: it seeds siopts.J0z AND
    % enters pt.crit below (both uses receive the same single shift).
    J0eff = J0eff - d;
end

% Moments AFTER the point's coupling spectrum is resolved: the ODD branch above rebuilds
% Jnu_flat from odd_blocks + deltaJ, so taking them earlier would describe the wrong multiset.
Jmom = invz_coupling_moments(Jnu_flat);
eopts.Jmom = Jmom;  eso_pub.Jmom = Jmom;
if sm.is_strict && ~isvector(Jnu_flat)
    error('invz:staticMedium', ['strict ordered/Jensen mode does not support the [nJ,nw] ' ...
        'retarded coupling matrix in this phase; PM strict mode remains supported.']);
end
% opts.odd_retarded / opts.odd_retarded_exact are SILENTLY IGNORED by this solver (they appear
% nowhere else in it -- there is no retarded ordered path, and no invented 'ordered_retarded'
% option). Under a strict scheme that silence would become load-bearing: the caller would get a
% strict ordered solve built on the static multiset while believing it retarded. Rejected here,
% after the options are actually resolved and before any HMF work starts.
if sm.is_strict
    for f = {'odd_retarded', 'odd_retarded_exact'}
        if isfield(opts, f{1}) && ~isempty(opts.(f{1})) && ~isequal(opts.(f{1}), false)
            error('invz:staticMedium', ['opts.%s is not supported under static_medium ''%s'': ' ...
                'this solver has no retarded ordered path (the flag is silently ignored on the ' ...
                '''resummed'' path), so a strict ordered solve would silently describe the ' ...
                'STATIC coupling multiset. Use the PM leg for retarded strict solves.'], ...
                f{1}, sm.scheme);
        end
    end
end

[wn, wts, beta] = invz_matsubara(T, Ecut);

% Ordered mean-field solve (full electronuclear space): spontaneous moment m0 and field hz.
% J0z is the SAME cc coupling J(0) used by the criticality and the RPA/1z denominator.
% forced_moment (spec 2026-07-16): with an explicit longitudinal Bx(3) the moment is
% field-induced -- the spontaneous |m0| > mtol gate is bypassed and branch alignment
% with the applied Bz is enforced (sign-aware seed + one mirrored retry).
siopts = struct('hyp', hyp, 'order', true, 'J0z', J0eff, 'Jxx0', Jxx0, 'transverse_mf', tmf);
for f = {'mz_seed', 'mf_maxit', 'mf_mix'}                  % diagnostic pass-throughs (tests)
    if isfield(opts, f{1}), siopts.(f{1}) = opts.(f{1}); end
end

omode = getf(opts, 'ordered_mode', 'bare');
if ~any(strcmp(omode, {'bare', 'jensen'}))
    error('invz:orderedMode', 'ordered_mode must be ''bare'' or ''jensen''.');
end
if strcmp(omode, 'jensen')
    if fmom || abs(Bx(3)) > getf(opts, 'bz_tol', 1e-9)
        error('invz:orderedMode', 'ordered_mode ''jensen'' is transverse/spontaneous only.');
    end
    hopts = opts;                                    % FULL numerical context (P1-6) ...
    hopts.J0eff = J0eff;                             % ... with the ODD-shifted coupling
    hopts.Jmom = Jmom;                               % ... the resolved moments (no re-derive)
    hopts.static_medium = sm.scheme;                 % ... and the resolved scheme
    hopts.emt = eopts;  hopts.emt_static = eso_pub;  % ... stamped (validation is idempotent)
    for f = {'ordered_mode', 'forced_moment'}        % ... and mode fields stripped
        if isfield(hopts, f{1}), hopts = rmfield(hopts, f{1}); end
    end
    [hstar, hprof] = invz_hmf_ordered(ion, T, Bx, Jnu_flat, hopts);
    if ~isfinite(hstar)
        si = invz_single_ion(ion, T, Bx, struct('hyp', hyp, 'hz_fixed', 0, ...
                                                'Jxx0', Jxx0, 'transverse_mf', tmf));
        pt = early_return(0, si, 'none', sm, Jmom);  % paramagnetic: PM leg owns this field
        pt.ordered_mode = omode;  pt.hmf_status = hprof.status;  pt.hmf_prof = hprof;
        % Preserve the reducer's deterministic binding cause at point level. The shared early
        % builder supplied safe defaults; override them only when a failed/domain node exists.
        if isfield(hprof, 'status_detail') && isstruct(hprof.status_detail) && ...
                isscalar(hprof.status_detail) && ...
                isfield(hprof.status_detail, 'binding_node') && ...
                isstruct(hprof.status_detail.binding_node) && ...
                isscalar(hprof.status_detail.binding_node)
            binding = hprof.status_detail.binding_node;
            pt.medium_status = binding.medium_status;
            pt.medium_denom = binding.ref_denom;
            pt.medium_margin = binding.ref_margin;
        end
        % Jensen-only members, blanked so the jensen exit schema is uniform. This is the ONLY
        % early return reachable in jensen mode (the other three sit behind forced_moment,
        % which jensen rejects, and behind ~is_ordered, which jensen overrides to true), so
        % there is exactly one blank site and one populated site -- no drift surface. No root
        % was accepted, so there is no endpoint to classify and no solved path to reduce; NaN
        % is also the fail-closed value of path_omit_max. A domain/degenerate reason is NOT
        % lost here -- pt.hmf_status carries it.
        pt.D_uni = NaN;  pt.crit_1z = NaN;  pt.Dq_min = NaN;  pt.stable_1z = false;
        pt.omit_mu3 = NaN;  pt.omit_cubic = NaN;  pt.omit_max = NaN;  pt.path_omit_max = NaN;
        if strcmp(hprof.status, 'unresolved')
            pt.converged = false;                    % round-3 P0-3: ordering was PREDICTED
        end                                          % but unbracketed -- NOT a PM verdict;
        return;                                      % the map masks this column
    end
    siopts.hz_fixed = hstar;                         % impose the jensen molecular field ...
    siopts.order = false;                            % ... WITHOUT the ordering update (P0-1)
end
si = invz_single_ion(ion, T, Bx, siopts);
branch = 'spontaneous';  if fmom, branch = 'field_induced'; end
if fmom && Bx(3) ~= 0 && si.mf_converged && abs(si.Jexp(3)) > 1e-10 && sign(si.Jexp(3)) ~= sign(Bx(3))
    % converged onto the metastable anti-aligned branch: one mirrored retry
    siopts.mz_seed = -sign(si.Jexp(3));
    si2 = invz_single_ion(ion, T, Bx, siopts);
    if si2.mf_converged && sign(si2.Jexp(3)) == sign(Bx(3))
        si = si2;
    else
        warning('invz:branchMismatch', ...
            'Anti-aligned moment persists at Bz = %.3g T after mirrored retry.', Bx(3));
        pt = early_return(si.Jexp(3), si, branch, sm, Jmom);
        return;
    end
end
if fmom && ~si.mf_converged
    pt = early_return(si.Jexp(3), si, branch, sm, Jmom);   % MF gate (second review finding 6)
    return;
end
m0 = si.Jexp(3);
pt.m0 = m0;
pt.is_ordered = fmom || abs(m0) > mtol;
if strcmp(omode, 'jensen'), pt.is_ordered = true; end      % root existence gates jensen (P1-4)
if ~pt.is_ordered
    pt = early_return(m0, si, 'none', sm, Jmom);           % paramagnetic point: use invz_solve_point
    return;
end

c0  = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
G0  = -real(squeeze(c0(3,3,:)));               % full electronuclear cc, ordered moment included
tl  = invz_twolevel_ordered(ion, T, Bx, si.hz, struct('Jxx0', Jxx0, 'transverse_mf', tmf));
g   = real(invz_g(tl, 1i*wn));

Sigma = zeros(size(wn));  K = zeros(size(wn));
converged = false;
if strcmp(omode, 'jensen')
    % Static-sector closure insertion (Task 3's eval_node statement order, P1-F/P1-A):
    % full-electronuclear split weights from the FINAL si, mode-switched chain rule.
    eso = eso_pub;      % the scheme/ref_margin/Jmom-stamped struct resolved at the entry
    eso.warn = false;   % solver gates on so.converged; suppress the per-node console flood
    c0i = invz_chi0z(si, T, 1i*wn(1), struct('elastic', false));   % static inelastic only
    G0inel0 = -real(c0i(3,3,1));                                   % fixed-Hamiltonian slot
    X = real(c0(:, :, 1));                                         % static chi tensor (chi=-G)
    switch tmf
        case 'none'
            fb = 0;
        case 'legacy_x'
            fb = X(3, 1) * (Jxx0 / (1 - Jxx0*X(1, 1))) * X(1, 3);
        case 'vector_ab'
            t = [1 2];
            fb = X(3, t) * (Jxx0 * ((eye(2) - Jxx0*X(t, t)) \ X(t, 3)));
        otherwise
            error('invz:transverseMF', 'unknown transverse_mf ''%s''', tmf);
    end
    G0bare0 = -(X(3, 3) + fb);
    G0el0   = G0bare0 - G0inel0;                                   % elastic + feedback (SS4a)
    % --- stage-2c task 1b-ii-B: ONE call to the shared, checker-gated node solver --------
    % (replaces the old inline cold-init + 7-step loop + post-loop refresh + ad hoc ctol
    % gate: invz_ordered_node_solve runs that SAME map verbatim and gates acceptance on the
    % complete residual checker, invz_ordered_residual -- see both files' headers.)
    node = struct('tl', tl, 'G0', G0, 'g', g, 'wts', wts, 'wn', wn, 'beta', beta, ...
        'J0eff', J0eff, 'G0inel0', G0inel0, 'G0el0', G0el0, 'G0bare0', G0bare0, ...
        'eso', eso, 'eopts', eopts, 'Jnu_flat', Jnu_flat, 'Jmom', Jmom);
    % Jmom (resolved ONCE above) is REQUIRED by invz_ordered_node_solve whenever node.eso
    % selects a strict scheme ('invz:nodeSolveNode' otherwise); it is threaded into both EMT
    % leaves there so neither silently re-derives it. Harmless-absent under 'resummed'.
    sopts = struct('mix_outer', mixo, 'max_outer', maxo, 'tol_outer', tolo, ...
        'cold_retry', true, 'trace', false, 'cold_acceleration', cold_accel);
    [state, info] = invz_ordered_node_solve(node, [], sopts);   % COLD: this leg always
                                                                 % starts from Sigma=0/K0s=0
    Sigma = state.Sigma;  K = state.K;  lam = state.lam;
    med = info.med;  med.G(1) = info.so.Gstat;         % elastic static function (P1-D),
                                                        % not the ordinary Dyson value
    sg  = invz_sigma_ordered(tl, lam, K, g, beta);     % recompute for pt.alpha/pt.alpha_m --
                                                        % info does not expose sg
    converged = info.loop_converged;   % in-loop verdict -- diagnostic only; pt.converged
    outer     = info.outer_iters;      % below is checker-gated on info.accepted instead
else
    for outer = 1:maxo
        eopts.K0 = K;
        med = invz_emt_scalar(G0, Sigma, Jnu_flat, eopts);
        K   = med.K;
        if ~any(strcmp(med.medium_status, {'ok', 'not_applicable'}))
            % Mirrors invz_solve_point's PM domain halt verbatim (same two whitelisted strings,
            % same lam/sg reset, same break into the COMMON export block below). 'bare' is the
            % DEFAULT ordered_mode and is reachable from invz_solve_auto, so a strict-scheme
            % reference/closure domain event must halt HERE, before invz_lambdas and
            % invz_sigma_ordered consume the invalid medium: without this the NaN Sigma feeds
            % the next iteration's reference, the loop burns max_outer iterations on NaNs, and
            % the exported medium_status degrades from the TRUE first cause (e.g.
            % ref_denom_small) to 'nonfinite' -- provenance naming the wrong condition.
            % INERT under the default 'resummed' scheme, where invz_emt_scalar's whole strict
            % block is gated off (:73) and medium_status is always 'not_applicable'.
            % lam/sg are reset (not left at the previous iterate) so the exported point carries
            % no stale lambda/alpha for a medium that has no solution; the export block then
            % produces the COMPLETE field set with pt.converged false (med.converged is false
            % here by construction, K(1)/G(1) = NaN) and pt.medium_status the exact domain
            % string. lam is 3-long and sg carries alpha_m here -- the ordered leg's shapes.
            lam = nan(3,1);
            sg  = struct('Sigma', nan(size(Sigma)), 'alpha', NaN, 'alpha_m', NaN);
            break;
        end
        lam = invz_lambdas(K, g, wts, beta, [1 2 3]);
        sg  = invz_sigma_ordered(tl, lam, K, g, beta);
        dS  = invz_finite_max_abs(sg.Sigma, Sigma);
        Sigma = Sigma + mixo*(sg.Sigma - Sigma);
        if dS < tolo, converged = true; break; end
    end
end
pt.Sigma0 = Sigma(1);  pt.alpha = sg.alpha;  pt.alpha_m = sg.alpha_m;  pt.lambda = lam;
pt.K = K;  pt.G = med.G;  pt.Sigma = Sigma;  pt.tl = tl;  pt.si = si;
pt.chi0cc0 = -G0(1);
pt.crit = 1 + pt.Sigma0 - J0eff*pt.chi0cc0;
pt.sumrule_rel = abs(sum(wts.*med.G)/beta + si.JzJz_fluct) / max(abs(si.JzJz_fluct), 1e-12);
pt.converged = converged && med.converged;
pt.outer_iters = outer;
pt.moment_branch = branch;
% --- static-medium provenance (spec SS4.2): stamped in the COMMON export, so BOTH ordered
% modes carry the same members. Here the three medium_* values describe the DYNAMIC medium
% invz_emt_scalar returned; in jensen mode the block below overwrites them with the ordered
% STATIC sector's own record, which is the sector that actually owns slot 1 there. `medium.ref`
% is [] under 'resummed' and getf is safe on that, so both numbers read NaN.
pt.static_medium = sm.scheme;
pt.Jmom = Jmom;
ref = getf(getf(med, 'medium', struct()), 'ref', struct());
pt.medium_status = getf(med, 'medium_status', 'not_applicable');
pt.medium_denom  = getf(ref, 'denom', NaN);
pt.medium_margin = getf(ref, 'margin', NaN);   % distance to floor, NOT the denominator
if oddOn
    pt.odd = struct('d', d, 'Xp', Xp);         % T1.4 diagnostics (absent when flag off)
end
if strcmp(omode, 'jensen')
    % stage-2c task 1b-ii-B: acceptance = the complete four-block checker verdict, NOT the
    % old in-loop dS/sout.converged + post-loop staticok/ctol re-gate (both folded away --
    % info.accepted already IS invz_ordered_residual's res.accepted for this exported state).
    % final_resid is diagnostic-only from here on: exposed as block C (the derived lam/Sigma
    % chain from the exported K) -- the EXACT quantity the old final_resid computed
    % (production's own invz_ordered_residual.m docstring: "this is production's existing
    % final_resid, named ... here" -- the two must stay byte-equal).
    pt.final_resid = info.res.blockC.resid;
    pt.converged = info.accepted;
    pt.ordered_mode = omode;  pt.hmf = hstar;  pt.hmf_prof = hprof;
    pt.D_uni = info.so.D_uni;                        % pole observable AT THE FINAL STATE
                                                      % (checker's own value; algebraically
                                                      % identical to 1+(J0eff-K(1))*med.G(1))
    % TWO-TIER EXPORT (spec SS1): pt.converged keeps its existing meaning (info.accepted, the
    % consistency tier) so no existing consumer shifts meaning under it. Endpoint stability is
    % a SEPARATE field -- collapsing them would re-mask the ordered phase, because intermediate
    % path nodes are the unstable Landau interval by construction. The classification itself is
    % the checker's (res.stability already applied the frozen crit_tol/D_tol/Dq_tol); it is
    % never rebuilt here from raw > 0 comparisons. The crit agreement term only confirms that
    % the exported endpoint IS the root the profile accepted.
    pt.crit_1z = hprof.crit_star;
    pt.Dq_min = info.res.stability.Dq_min;
    pt.stable_1z = pt.converged && info.res.stability.pass && ...
                   isfinite(hprof.crit_star) && ...
                   abs(hprof.crit_star-info.res.stability.crit) <= ...
                       getf(opts,'crit_tol',1e-6);
    pt.omit_mu3 = getf(info.so, 'omit_mu3', NaN);
    pt.omit_cubic = getf(info.so, 'omit_cubic', NaN);
    pt.omit_max = getf(info.so, 'omit_max', NaN);
    % path_omit_max FAILS CLOSED (plan-owner ruling 2026-07-26, mirroring the task-3 omit_max
    % ruling and DELIBERATELY replacing this brief's isfinite-filtered version): it is a
    % load-bearing promotion gate (the frozen predicate compares max(omit_max) over the
    % solved path against omit_promote = 0.10 -- quoted in
    % docs/invzp_strict_medium_gate0_report.md SS1), so a
    % corrupted node must never be quietly dropped from the maximum. Inf DOMINATES, NaN POISONS.
    if isempty(hprof.omit_max)
        pt.path_omit_max = NaN;
    elseif any(isnan(hprof.omit_max))
        pt.path_omit_max = NaN;      % fail closed: a corrupted node
                                     % must not be silently dropped
    else
        pt.path_omit_max = max(hprof.omit_max);   % Inf dominates,
    end                                           % per the frozen
                                                  % zero-denominator
                                                  % convention
    % The ordered STATIC sector owns slot 1 here, so its own medium record supersedes the
    % dynamic one stamped in the common block above.
    ref = getf(getf(info, 'medium', struct()), 'ref', struct());
    pt.medium_status = getf(info, 'medium_status', 'not_applicable');
    pt.medium_denom = getf(ref, 'denom', NaN);
    pt.medium_margin = getf(ref, 'margin', NaN);
    % CONTRACT (P2-G): pt.crit keeps its historical ordinary-Dyson definition and is NOT
    % the ordered pole mass below the boundary -- pt.D_uni is (see docstring).
    if abs(pt.si.hz - hstar) > 1e-12
        error('invz:hzFixed', 'jensen final state did not hold hmf: %.6g vs %.6g', pt.si.hz, hstar);
    end
end
end

% -------------------------------------------------------------------------------------------
function pt = early_return(m0, si, branch, sm, Jmom)
%EARLY_RETURN Complete field set for every non-accepted exit (spec: callers never
% probe a missing member; tl = [] flags "no two-level params were built").
% The static-medium provenance (spec SS4.2) is stamped HERE, inside the shared builder,
% rather than repeated as five post-call assignments at each of the four call sites --
% repeating them is exactly how a member goes missing on one of them. No medium was ever
% solved on any of these paths, so medium_status is 'not_applicable' and both reference
% numbers are NaN; a jensen exit that failed for a domain reason still reports it, through
% pt.hmf_status, which that call site sets.
pt = struct('m0', m0, 'is_ordered', false, 'converged', false, 'Sigma0', NaN, ...
            'crit', NaN, 'si', si, 'tl', [], 'moment_branch', branch, ...
            'static_medium', sm.scheme, 'Jmom', Jmom, ...
            'medium_status', 'not_applicable', 'medium_denom', NaN, 'medium_margin', NaN);
end
