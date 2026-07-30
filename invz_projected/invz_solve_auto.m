function [pt, phase, diag] = invz_solve_auto(ion, T, Bx, Jnu_flat, opts)
%INVZ_SOLVE_AUTO Jensen-ordered solve with a paramagnetic fallback.
% phase: 1 = accepted ordered state, 2 = accepted stable paramagnet,
%        0 = no accepted state.
if nargin < 5, opts = struct(); end
pt = [];
phase = 0;
diag = struct('ordered_err', '', 'para_err', '');

try
    ordered = invz_solve_point_ordered(ion, T, Bx, Jnu_flat, opts);
    if ordered.is_ordered && ordered.converged && isfinite(ordered.Sigma0)
        pt = ordered;
        phase = 1;
        return;
    end
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    diag.ordered_err = err.identifier;
end

try
    para = invz_solve_point(ion, T, Bx, Jnu_flat, opts);
    if para.converged && isfinite(para.Sigma0) && isfinite(para.crit) && para.crit > 0
        pt = para;
        phase = 2;
    end
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    diag.para_err = err.identifier;
end
end
