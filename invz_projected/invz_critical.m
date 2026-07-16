function bx = invz_critical(ion, T, Jnu_flat, opts)
%INVZ_CRITICAL Critical transverse field at temperature T (paramagnetic side).
% Root of pt.crit = 1 + Sigma0 - J(0)*chi0_cc(0) over field at fixed T: crit is
% positive in the paramagnet (high B) and crosses zero at the boundary; the
% ordered side (low B) is where the paramagnetic EMT fixed point ceases to
% exist and invz_solve_point returns non-finite / non-converged crit.
%
% Classify from CONVERGED points only: near the boundary the paramagnetic
% self-consistency suffers critical slowing down and can return non-finite
% crit that would otherwise look like a false ordered/paramagnet crossing.
% Scans the field DOWN from the window top, taking the first ordered/paramagnet
% crossing between converged points; non-converged points get no vote.
%
% opts.window = [Blo Bhi] (T, default [2 7]): Bhi must be a converged paramagnet
%   (crit>0) at this T; Blo the low-field floor of the search.
% opts.fieldstep coarse scan step (T, default 0.20).
% opts.tol       crossing tolerance (T, default 0.02).
% Remaining opts pass through to invz_solve_point.
if nargin < 4, opts = struct(); end
win  = [2 7];  if isfield(opts,'window'),    win  = opts.window;    end
tol  = 0.02;   if isfield(opts,'tol'),       tol  = opts.tol;       end
step = 0.20;   if isfield(opts,'fieldstep'), step = opts.fieldstep; end
gap_steps = 6;                       % how far below a non-converged patch to hunt

f = @(B) crit_at(ion, T, B, Jnu_flat, opts);

[cTop, okTop] = f(win(2));
assert(okTop && cTop > 0, 'invz:bracket', ...
    ['Top of field window B=%.2f T is not a converged paramagnet at T=%.3f K ' ...
     '(crit=%.3g, converged=%d): raise the window top.'], win(2), T, cTop, okTop);

Bpara = win(2);  cpara = cTop;        % last confirmed converged paramagnet
B = win(2) - step;
while B >= win(1) - 1e-9
    [c, ok] = f(B);
    if ok && c > 0
        Bpara = B;  cpara = c;  B = B - step;  continue;   % descend through paramagnet
    end
    if ok && c <= 0
        bx = refine_crossing(f, B, c, Bpara, cpara, tol);   % clean converged bracket
        return;
    end
    % Non-converged at B: hunt a few steps down for a converged verdict so the
    % crossing is interpolated between converged points across the failure band.
    found = false;
    Bh = B - step;
    hardStop = max(win(1), B - gap_steps*step);
    while Bh >= hardStop - 1e-9
        [ch, okh] = f(Bh);
        if okh && ch <= 0
            bx = refine_crossing(f, Bh, ch, Bpara, cpara, tol);   % bracket spans the gap
            return;
        elseif okh && ch > 0
            Bpara = Bh;  cpara = ch;  B = Bh - step;  found = true;  break;   % non-conv island
        end
        Bh = Bh - step;
    end
    if ~found
        warning('invz:orderedSideNoConverge', ...
            ['T = %.3f K: the ordered side near the boundary did not converge; ' ...
             'returning the paramagnet-edge field estimate.'], T);
        bx = para_edge(f, B, Bpara, tol);
        return;
    end
end
error('invz:bracket', ...
    'T = %.3f K: no ordered/paramagnet crossing in field window [%.2f, %.2f] T.', ...
    T, win(1), win(2));
end

% ------------------------------------------------------------------------
function [c, ok] = crit_at(ion, T, B, Jf, opts)
% crit at (T,B) and whether it is a trustworthy (converged, finite) verdict.
try
    pt = invz_solve_point(ion, T, B, Jf, opts);
    c  = pt.crit;
    ok = pt.converged && isfinite(c);
catch
    c = NaN;  ok = false;                     % e.g. invz:degenerateDoublet
end
end

% ------------------------------------------------------------------------
function bx = refine_crossing(f, Ba, ca, Bb, cb, tol)
% Regula-falsi between CONVERGED bracket ends (ca<=0 at Ba, cb>0 at Bb, Ba<Bb).
% A non-converged interior sample is skipped; if the interior will not converge,
% fall back to linear interpolation across the bracket.
for r = 1:20
    if Bb - Ba < tol, break; end
    Bm = Ba - ca*(Bb - Ba)/(cb - ca);
    Bm = min(max(Bm, Ba + 0.05*(Bb-Ba)), Bb - 0.05*(Bb-Ba));
    [cm, okm] = f(Bm);
    if ~okm
        Bm = 0.5*(Ba + Bb);
        [cm, okm] = f(Bm);
        if ~okm, break; end
    end
    if cm <= 0, Ba = Bm; ca = cm; else, Bb = Bm; cb = cm; end
end
bx = Ba - ca*(Bb - Ba)/(cb - ca);
end

% ------------------------------------------------------------------------
function bx = para_edge(f, Blo, Bhi, tol)
% Fallback when the ordered side will not converge: bisect on "converged
% paramagnet" to locate the lower edge of the paramagnet (Blo not-para, Bhi para).
for r = 1:40
    if Bhi - Blo < tol, break; end
    Bm = 0.5*(Blo + Bhi);
    [cm, okm] = f(Bm);
    if okm && cm > 0, Bhi = Bm; else, Blo = Bm; end
end
bx = 0.5*(Blo + Bhi);
end
