function pt = invz_solve_point(ion, T, Bx, Jnu_flat, opts)
%INVZ_SOLVE_POINT Self-consistent 1/z solution at one paramagnetic (T, Bx) point.
% Bx: scalar (transverse, historical) or [Bx By Bz] vector (T).
% Outer loop: EMT K (Task 7) <-> lambda_p (Task 8) <-> Sigma, at fixed single-ion input.
% Inside the ordered phase the paramagnetic EMT fixed point does not exist; outputs may be non-finite and pt.converged false — invz_critical relies on this as the ordered-phase signal. Always check pt.converged.
if nargin < 5, opts = struct(); end
Ecut  = 40;   if isfield(opts,'Ecut'),      Ecut  = opts.Ecut;      end
hyp   = true; if isfield(opts,'hyp'),       hyp   = opts.hyp;       end
J0eff = ion.J0eff; if isfield(opts,'J0eff'), J0eff = opts.J0eff;    end
Jxx0  = ion.Jxx0;  if isfield(opts,'Jxx0'), Jxx0 = opts.Jxx0; end
mixo  = 0.7;  if isfield(opts,'mix_outer'), mixo  = opts.mix_outer; end
tolo  = 1e-8; if isfield(opts,'tol_outer'), tolo  = opts.tol_outer; end
maxo  = 60;   if isfield(opts,'max_outer'), maxo  = opts.max_outer; end
eopts = struct(); if isfield(opts,'emt'), eopts = opts.emt; end
Bx = invz_field_vec(Bx);                       % scalar -> [Bx 0 0]; 3-vector passes through

[wn, wts, beta] = invz_matsubara(T, Ecut);
si  = invz_single_ion(ion, T, Bx, struct('hyp', hyp, 'Jxx0', Jxx0));
c0  = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
G0  = -real(squeeze(c0(3,3,:)));                 % full (electro)nuclear cc Green function
tl  = invz_twolevel(ion, T, Bx, struct('Jxx0', Jxx0));   % electronic two-level params for Sigma
g   = real(invz_g(tl, 1i*wn));

% Warm-start seed (finding #6): a converged neighbouring point's Sigma speeds up
% the outer fixed point. It must be the SAME length as wn (same T), and only
% changes the iteration path, not the converged fixed point; ignored otherwise.
% (Only Sigma is seeded: the closed-form EMT solves K directly and ignores any K
% seed, so seeding K had no effect.)
Sigma = zeros(size(wn));
if isfield(opts,'Sigma_seed') && numel(opts.Sigma_seed) == numel(wn), Sigma = opts.Sigma_seed(:); end
converged = false;
for outer = 1:maxo
    med = invz_emt_scalar(G0, Sigma, Jnu_flat, eopts);
    K   = med.K;
    lam = invz_lambdas(K, g, wts, beta, [1 2]);
    sg  = invz_sigma(tl, lam, K, g, beta);
    dS  = max(abs(sg.Sigma - Sigma));
    Sigma = Sigma + mixo*(sg.Sigma - Sigma);
    if dS < tolo, converged = true; break; end
end
pt.Sigma0 = Sigma(1);  pt.alpha = sg.alpha;  pt.lambda = lam;
pt.K = K;  pt.G = med.G;  pt.Sigma = Sigma;  pt.tl = tl;  pt.si = si;
pt.chi0cc0 = -G0(1);
pt.crit = 1 + pt.Sigma0 - J0eff*pt.chi0cc0;
pt.sumrule_rel = abs(sum(wts.*med.G)/beta + si.JzJz_fluct) / si.JzJz_fluct;
pt.converged = converged && med.converged;
pt.outer_iters = outer;
end
