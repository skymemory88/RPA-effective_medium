function bx = invz_refine_crossing(f, Ba, ca, Bb, cb, tol)
%INVZ_REFINE_CROSSING Regula-falsi between CONVERGED bracket ends (ca<=0 at Ba, cb>0 at Bb, Ba<Bb).
% A non-converged interior sample is skipped by retrying at the midpoint; if the
% interior will not converge either, falls back to linear interpolation across
% the bracket. f is a closure over whatever scalar is being bracketed (e.g.
% @(B) invz_crit_at(ion, T, B, Jnu_flat, opts)), returning [value, converged].
% Shared by invz_critical (field root) and invz_critical_T (temperature root).
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
