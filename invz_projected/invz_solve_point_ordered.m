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
% CONTRACT (P2-G): pt.crit KEEPS its historical ordinary-Dyson definition
% (1 + Sigma0 - J0eff*chi0cc0) as a legacy diagnostic in EITHER mode -- it is NOT the ordered
% pole mass below the boundary; pt.D_uni is. Callers must not conflate the two.
if nargin < 5, opts = struct(); end
Ecut  = getf(opts, 'Ecut', 40);
hyp   = getf(opts, 'hyp', true);
J0eff = getf(opts, 'J0eff', ion.J0eff);
Jxx0  = getf(opts, 'Jxx0', ion.Jxx0);
tmf   = getf(opts, 'transverse_mf', 'legacy_x');
mixo  = getf(opts, 'mix_outer', 0.7);
tolo  = getf(opts, 'tol_outer', 1e-8);
maxo  = getf(opts, 'max_outer', 200);
mtol  = getf(opts, 'm_tol', 1e-2);
eopts = getf(opts, 'emt', struct());
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
    for f = {'ordered_mode', 'forced_moment'}        % ... and mode fields stripped
        if isfield(hopts, f{1}), hopts = rmfield(hopts, f{1}); end
    end
    [hstar, hprof] = invz_hmf_ordered(ion, T, Bx, Jnu_flat, hopts);
    if ~isfinite(hstar)
        si = invz_single_ion(ion, T, Bx, struct('hyp', hyp, 'hz_fixed', 0, ...
                                                'Jxx0', Jxx0, 'transverse_mf', tmf));
        pt = early_return(0, si, 'none');            % paramagnetic: PM leg owns this field
        pt.ordered_mode = omode;  pt.hmf_status = hprof.status;
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
        pt = early_return(si.Jexp(3), si, branch);
        return;
    end
end
if fmom && ~si.mf_converged
    pt = early_return(si.Jexp(3), si, branch);             % MF gate (second review finding 6)
    return;
end
m0 = si.Jexp(3);
pt.m0 = m0;
pt.is_ordered = fmom || abs(m0) > mtol;
if strcmp(omode, 'jensen'), pt.is_ordered = true; end      % root existence gates jensen (P1-4)
if ~pt.is_ordered
    pt = early_return(m0, si, 'none');                     % paramagnetic point: use invz_solve_point
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
    eso = getf(opts, 'emt_static', struct());
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
    K0s = 0;  lam = [0; 0; 0];
    for outer = 1:maxo
        % (1) dynamic sector -- MIRRORS the bare loop's emt call verbatim
        eopts.K0 = K;
        med = invz_emt_scalar(G0, Sigma, Jnu_flat, eopts);
        K   = med.K;
        % (2) static sector (P0-2/P0-A), threaded opts (P1-F):
        [K0s, ~, sout] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), Jnu_flat, K0s, ...
                                                 beta, J0eff, G0inel0, G0el0, eso);
        K(1) = K0s;
        % (3)-(5) lambdas, ordered Sigma, damped mix -- MIRRORS the bare loop's statements
        lam = invz_lambdas(K, g, wts, beta, [1 2 3]);
        sg  = invz_sigma_ordered(tl, lam, K, g, beta);
        dS  = max(abs(sg.Sigma - Sigma));
        Sigma = Sigma + mixo*(sg.Sigma - Sigma);
        if dS < tolo && sout.converged, converged = true; break; end
    end
    % Final post-loop static refresh (round-4 P1-B / round-5 P1-B): KEEPS its newly
    % closed K0 (computed on the just-mixed Sigma(1)), written to K(1) before packing;
    % gated on its OWN convergence/residual, folded into pt.converged below.
    [K0s, Gstat, so] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), Jnu_flat, K0s, ...
                                               beta, J0eff, G0inel0, G0el0, eso);
    K(1) = K0s;
    ctol = getf(eso, 'resid_tol', 1e-10);
    staticok = so.converged && isfinite(so.resid) && so.resid < ctol;
    % Residual-only lambda/Sigma revalidation (round-5 P2): the exported tuple is
    % closure-consistent to the STATED OUTER TOLERANCE, not exactly self-consistent.
    lam_check   = invz_lambdas(K, g, wts, beta, [1 2 3]);
    Sigma_check = invz_sigma_ordered(tl, lam_check, K, g, beta);
    final_resid = max(abs(Sigma_check.Sigma - Sigma));
    med.G(1) = Gstat;                                  % elastic static function (P1-D),
                                                        % not the ordinary Dyson value
else
    for outer = 1:maxo
        eopts.K0 = K;
        med = invz_emt_scalar(G0, Sigma, Jnu_flat, eopts);
        K   = med.K;
        lam = invz_lambdas(K, g, wts, beta, [1 2 3]);
        sg  = invz_sigma_ordered(tl, lam, K, g, beta);
        dS  = max(abs(sg.Sigma - Sigma));
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
if oddOn
    pt.odd = struct('d', d, 'Xp', Xp);         % T1.4 diagnostics (absent when flag off)
end
if strcmp(omode, 'jensen')
    pt.final_resid = final_resid;
    pt.converged = pt.converged && staticok && (pt.final_resid < tolo);
    pt.ordered_mode = omode;  pt.hmf = hstar;  pt.hmf_prof = hprof;
    pt.D_uni = 1 + (J0eff - K(1))*med.G(1);          % pole observable AT THE FINAL STATE
    % CONTRACT (P2-G): pt.crit keeps its historical ordinary-Dyson definition and is NOT
    % the ordered pole mass below the boundary -- pt.D_uni is (see docstring).
    if abs(pt.si.hz - hstar) > 1e-12
        error('invz:hzFixed', 'jensen final state did not hold hmf: %.6g vs %.6g', pt.si.hz, hstar);
    end
end
end

% -------------------------------------------------------------------------------------------
function pt = early_return(m0, si, branch)
%EARLY_RETURN Complete field set for every non-accepted exit (spec: callers never
% probe a missing member; tl = [] flags "no two-level params were built").
pt = struct('m0', m0, 'is_ordered', false, 'converged', false, 'Sigma0', NaN, ...
            'crit', NaN, 'si', si, 'tl', [], 'moment_branch', branch);
end
