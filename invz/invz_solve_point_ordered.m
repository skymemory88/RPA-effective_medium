function pt = invz_solve_point_ordered(ion, T, Bx, Jnu_flat, opts)
%INVZ_SOLVE_POINT_ORDERED Self-consistent 1/z solution at one FERROMAGNETIC (T, Bx) point.
% Ordered-phase counterpart of invz_solve_point: the single-ion problem is solved with the
% longitudinal ORDERING mean field (spontaneous moment m0 = <Jz>), the self-energy uses the
% elastic-sector form (invz_sigma_ordered, HTML eqs 37-38), and the whole thing is iterated
% against the effective medium exactly as in the paramagnet.
%
% Returns pt.is_ordered: FALSE if the (T,Bx) point has no spontaneous moment (|m0| < m_tol) --
% the caller should then fall back to invz_solve_point (paramagnetic). When ordered, pt carries
% the same fields as invz_solve_point plus pt.m0 (order parameter), pt.alpha_m, pt.si (the
% ordered single-ion struct), and pt.is_ordered = true.
%
% SCOPE NOTE (Option A): m0 is the bare mean-field order parameter (the applied-field/H_MF
% self-consistency of HTML eqs 41-43 is deferred), so it onsets at the mean-field boundary,
% slightly above the true 1/z boundary; the gap matters only near B_c, not deep in the
% ordered phase.
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

[wn, wts, beta] = invz_matsubara(T, Ecut);

% Ordered mean-field solve (full electronuclear space): spontaneous moment m0 and field hz.
% J0z is the SAME cc coupling J(0) used by the criticality (line below) and the RPA/1z
% denominator (invz_chi_realaxis Jsel), so the mean-field ordering, the 1/z instability, and
% the RPA soft mode all read one J(0) instead of the hardcoded invz_single_ion default.
si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', hyp, 'order', true, 'J0z', J0eff, 'Jxx0', Jxx0));
m0 = si.Jexp(3);
pt.m0 = m0;
pt.is_ordered = abs(m0) > mtol;
if ~pt.is_ordered
    pt.converged = false;  pt.Sigma0 = NaN;  pt.crit = NaN;  pt.si = si;
    return;                                    % paramagnetic point: use invz_solve_point
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
end
