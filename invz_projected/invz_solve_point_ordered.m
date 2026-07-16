function pt = invz_solve_point_ordered(ion, T, Bx, Jnu_flat, opts)
%INVZ_SOLVE_POINT_ORDERED Self-consistent 1/z solution at one FERROMAGNETIC (T, Bx) point.
% Bx: scalar (transverse, historical) or [Bx By Bz] vector (T).
% Ordered-phase counterpart of invz_solve_point: the single-ion problem is solved with the
% longitudinal ORDERING mean field (spontaneous moment m0 = <Jz>), the self-energy uses the
% elastic-sector form (invz_sigma_ordered, HTML eqs 37-38), and the whole thing is iterated
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
% (the applied-field/H_MF self-consistency of HTML eqs 41-43 is deferred), so it onsets at the
% mean-field boundary, slightly above the true 1/z boundary; the gap matters only near B_c,
% not deep in the ordered phase.
if nargin < 5, opts = struct(); end
Ecut  = 40;   if isfield(opts,'Ecut'),      Ecut  = opts.Ecut;      end
hyp   = true; if isfield(opts,'hyp'),       hyp   = opts.hyp;       end
J0eff = ion.J0eff; if isfield(opts,'J0eff'), J0eff = opts.J0eff;    end
Jxx0  = ion.Jxx0;  if isfield(opts,'Jxx0'), Jxx0 = opts.Jxx0; end
mixo  = 0.7;  if isfield(opts,'mix_outer'), mixo  = opts.mix_outer; end
tolo  = 1e-8; if isfield(opts,'tol_outer'), tolo  = opts.tol_outer; end
maxo  = 80;   if isfield(opts,'max_outer'), maxo  = opts.max_outer; end
mtol  = 1e-2; if isfield(opts,'m_tol'),     mtol  = opts.m_tol;     end
eopts = struct(); if isfield(opts,'emt'), eopts = opts.emt; end
Bx = invz_field_vec(Bx);                       % scalar -> [Bx 0 0]; 3-vector passes through
fmom = isfield(opts,'forced_moment') && opts.forced_moment;

[wn, wts, beta] = invz_matsubara(T, Ecut);

% Ordered mean-field solve (full electronuclear space): spontaneous moment m0 and field hz.
% J0z is the SAME cc coupling J(0) used by the criticality and the RPA/1z denominator.
% forced_moment (spec 2026-07-16): with an explicit longitudinal Bx(3) the moment is
% field-induced -- the spontaneous |m0| > mtol gate is bypassed and branch alignment
% with the applied Bz is enforced (sign-aware seed + one mirrored retry).
siopts = struct('hyp', hyp, 'order', true, 'J0z', J0eff, 'Jxx0', Jxx0);
for f = {'mz_seed', 'mf_maxit', 'mf_mix'}                  % diagnostic pass-throughs (tests)
    if isfield(opts, f{1}), siopts.(f{1}) = opts.(f{1}); end
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
if ~pt.is_ordered
    pt = early_return(m0, si, 'none');                     % paramagnetic point: use invz_solve_point
    return;
end

c0  = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
G0  = -real(squeeze(c0(3,3,:)));               % full electronuclear cc, ordered moment included
tl  = invz_twolevel_ordered(ion, T, Bx, si.hz, struct('Jxx0', Jxx0));
g   = real(invz_g(tl, 1i*wn));

Sigma = zeros(size(wn));  K = zeros(size(wn));
converged = false;
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
pt.Sigma0 = Sigma(1);  pt.alpha = sg.alpha;  pt.alpha_m = sg.alpha_m;  pt.lambda = lam;
pt.K = K;  pt.G = med.G;  pt.Sigma = Sigma;  pt.tl = tl;  pt.si = si;
pt.chi0cc0 = -G0(1);
pt.crit = 1 + pt.Sigma0 - J0eff*pt.chi0cc0;
pt.sumrule_rel = abs(sum(wts.*med.G)/beta + si.JzJz_fluct) / max(abs(si.JzJz_fluct), 1e-12);
pt.converged = converged && med.converged;
pt.outer_iters = outer;
pt.moment_branch = branch;
end

% -------------------------------------------------------------------------------------------
function pt = early_return(m0, si, branch)
%EARLY_RETURN Complete field set for every non-accepted exit (spec: callers never
% probe a missing member; tl = [] flags "no two-level params were built").
pt = struct('m0', m0, 'is_ordered', false, 'converged', false, 'Sigma0', NaN, ...
            'crit', NaN, 'si', si, 'tl', [], 'moment_branch', branch);
end
