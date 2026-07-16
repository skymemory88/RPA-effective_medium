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
width = getf(opts, 'width',    0.5);
gstep = getf(opts, 'gridstep', 1/30);
tol   = getf(opts, 'tol',      0.005);
Tmin  = 0.02;                                   % single-ion solve floor
f = @(T) invz_crit_at(ion, T, Bx, Jnu_flat, opts);

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
        [c(i), ok(i)] = f(Tg(i));
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
            tc = invz_refine_crossing(f, Tv(k), cv(k), Tv(k+1), cv(k+1), tol);
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
function Tc0 = adaptive_anchor(ion, Jf, opts)
% Zero-field Tc used to anchor the adaptive window top.
if isfield(opts,'Tc0')
    Tc0 = opts.Tc0;  return;
end
J0 = getf(opts, 'J0eff', ion.J0eff);
Sc = invz_sigma_crit(J0, Jf);
Tc0 = invz_critical_T0field(ion, Sc, J0);
end
