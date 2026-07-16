function tc = invz_critical_T(ion, Bx, Jnu_flat, opts)
%INVZ_CRITICAL_T Critical temperature at fixed transverse field Bx (paramagnetic side).
% Root of pt.crit = 1 + Sigma(0) - J(0)*chi0_cc(0) over temperature at fixed
% field. crit is positive in the paramagnet (high T) and crosses zero at the
% boundary; inside the ordered phase (low T) the paramagnetic EMT fixed point
% does not exist and invz_solve_point returns non-finite / non-converged crit.
% Use where the boundary is nearly parallel to the field axis (T near the
% zero-field Tc0, ~1.74-1.78 K): a fixed-B cut crosses it transversally, where
% the fixed-T cut of invz_critical is ill-conditioned.
%
% Classify from CONVERGED points only: near the boundary the solver suffers
% critical slowing down and returns non-finite/non-monotone crit that a naive
% classifier would misread as ordered. Samples crit on a grid across the
% window, keeps only converged/finite points, and interpolates the highest-T
% sign change among them.
%
% opts.window = [Tlo Thi] (K): explicit search window.  If omitted, an
%   adaptive window is used: top anchored at Tc0+0.05 (Tc0 from opts.Tc0 or
%   computed via invz_critical_T0field), spanning 0.5 K down, and slid up/down
%   until it brackets a converged crossing.
% opts.Tc0    zero-field Tc anchor (K), to avoid recomputing it per field.
% opts.width  adaptive-window width (K, default 0.5).
% opts.gridstep coarse-grid step (K, default 1/30 K ~ 0.033).
% opts.tol    crossing refinement tolerance (K, default 0.005).
% Remaining opts pass through to invz_solve_point.
%
% Re-entrance: if more than one converged sign change is found, invz:multipleCrossings
% is raised and the highest-T crossing is returned.
%
% Small-B caveat: below ~0.5 T the doublet is near-degenerate and few points
% converge; if the whole window fails to bracket a crossing the function
% errors rather than guessing.
if nargin < 4, opts = struct(); end
width = 0.5;    if isfield(opts,'width'),    width = opts.width;    end
gstep = 1/30;   if isfield(opts,'gridstep'), gstep = opts.gridstep; end
tol   = 0.005;  if isfield(opts,'tol'),      tol   = opts.tol;      end
Tmin  = 0.02;                                   % single-ion solve floor

if isfield(opts,'window')
    Tlo = opts.window(1);  Thi = opts.window(2);
else
    Tc0 = adaptive_anchor(ion, Jnu_flat, opts);
    Thi = Tc0 + 0.05;  Tlo = Thi - width;
end

for slide = 0:8
    Tlo = max(Tlo, Tmin);
    ng  = max(5, round((Thi - Tlo)/gstep) + 1);
    Tg  = linspace(Tlo, Thi, ng);
    c   = nan(1, ng);  ok = false(1, ng);
    for i = 1:ng
        [c(i), ok(i)] = crit_at(ion, Tg(i), Bx, Jnu_flat, opts);
    end
    Tv = Tg(ok);  cv = c(ok);                   % converged, finite: the voters
    if numel(cv) >= 2
        sc = find(diff(sign(cv)) ~= 0);         % ordered(-) <-> para(+) transitions
        if numel(sc) > 1
            warning('invz:multipleCrossings', ...
                ['Bx = %.3f T: %d converged sign changes in [%.3f, %.3f] K ' ...
                 '(possible re-entrance); returning the highest-T crossing.'], ...
                Bx, numel(sc), Tlo, Thi);
        end
        up = sc(sign(cv(sc)) < 0 & sign(cv(sc+1)) > 0);   % low-T ordered -> high-T para
        if ~isempty(up)
            k  = up(end);                       % highest-T ordered->para crossing
            tc = refine_crossing(ion, Bx, Jnu_flat, opts, ...
                                 Tv(k), cv(k), Tv(k+1), cv(k+1), tol);
            return;
        end
    end
    % No converged crossing in this window: slide toward where it must be.
    % (Check isempty first: all([]>0) is true in MATLAB.)
    if isempty(cv)                              % nothing converged: keep top, grow down
        Tlo = Tlo - width;
    elseif all(cv > 0)                          % window all paramagnet: Tc below
        Thi = Tlo;  Tlo = Tlo - width;
    elseif all(cv < 0)                          % window all ordered: Tc above
        Tlo = Thi;  Thi = Thi + width;
    else
        break;                                  % mixed signs but no ord->para: give up
    end
end
error('invz:bracket', ...
    'Bx = %.3f T: no converged paramagnet/ordered crossing found (last window [%.3f, %.3f] K).', ...
    Bx, Tlo, Thi);
end

% ------------------------------------------------------------------------
function [c, ok] = crit_at(ion, T, Bx, Jf, opts)
% crit at (T,Bx) and whether it is a trustworthy (converged, finite) verdict.
try
    pt = invz_solve_point(ion, T, Bx, Jf, opts);
    c  = pt.crit;
    ok = pt.converged && isfinite(c);
catch
    c = NaN;  ok = false;                       % e.g. invz:degenerateDoublet
end
end

% ------------------------------------------------------------------------
function tc = refine_crossing(ion, Bx, Jf, opts, Ta, ca, Tb, cb, tol)
% Regula-falsi on CONVERGED crit between bracket ends (ca<0 at Ta, cb>0 at Tb).
% A non-converged interior sample is skipped (it does not corrupt the bracket);
% if the interior will not converge, fall back to linear interpolation.
for r = 1:20
    if Tb - Ta < tol, break; end
    Tm = Ta - ca*(Tb - Ta)/(cb - ca);          % false-position estimate
    Tm = min(max(Tm, Ta + 0.05*(Tb-Ta)), Tb - 0.05*(Tb-Ta));
    [cm, okm] = crit_at(ion, Tm, Bx, Jf, opts);
    if ~okm
        Tm = 0.5*(Ta + Tb);                     % try the midpoint instead
        [cm, okm] = crit_at(ion, Tm, Bx, Jf, opts);
        if ~okm, break; end                     % interior won't converge: interpolate
    end
    if cm <= 0, Ta = Tm; ca = cm; else, Tb = Tm; cb = cm; end
end
tc = Ta - ca*(Tb - Ta)/(cb - ca);              % final linear interpolation
end

% ------------------------------------------------------------------------
function Tc0 = adaptive_anchor(ion, Jf, opts)
% Zero-field Tc used to anchor the adaptive window top.
if isfield(opts,'Tc0')
    Tc0 = opts.Tc0;  return;
end
J0 = ion.J0eff;  if isfield(opts,'J0eff'), J0 = opts.J0eff; end
Sc = invz_sigma_crit(J0, Jf);
Tc0 = invz_critical_T0field(ion, Sc, J0);
end
