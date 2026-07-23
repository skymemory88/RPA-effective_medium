function [pass, diff_val, tol] = invz_phase1_gate(kind, v1, v2, J_ref)
%INVZ_PHASE1_GATE Frozen Phase-1 pairwise agreement-gate formulas (docs/invzp_phase1_quadrature_
% prereg.md item 5, "Normalized-shape gate" / "Energy-scalar gate"), shared verbatim by item 5's
% refinement gate and item 6's offset-sensitivity gate ("apply the pass/fail agreement gate
% (item-5 tolerances)"). Pure arithmetic, no I/O, no state.
%
% kind == 'shape'  (mean/variance/min/max/quantiles, ALREADY NORMALIZED: s = stat/J_ref, or
%                   stat/J_ref^2 for variance -- normalization happens in the CALLER, this
%                   function only applies the tolerance to the two already-normalized numbers):
%     tol  = 1e-6 + 1e-3*max(|v1|,|v2|)
%     pass = |v2 - v1| <= tol
%
% kind == 'energy' (J0eff, Jcc0, max(Jnu) -- RAW meV, NOT pre-normalized; J_ref enters the
%                   tolerance formula itself, not the compared values):
%     tol  = 1e-6*J_ref + 1e-4*max(|v1|,|v2|)
%     pass = |v2 - v1| <= tol
%
% OUTPUTS: pass (logical), diff_val = |v2-v1|, tol (the threshold actually applied) -- report-only
% triple so a caller can show the numeric margin, not just the boolean.
if ~(isnumeric(v1) && isscalar(v1) && isfinite(v1))
    error('invz:phase1Config', 'invz_phase1_gate: v1 must be a finite scalar; got %s.', mat2str(v1));
end
if ~(isnumeric(v2) && isscalar(v2) && isfinite(v2))
    error('invz:phase1Config', 'invz_phase1_gate: v2 must be a finite scalar; got %s.', mat2str(v2));
end
diff_val = abs(v2 - v1);
switch kind
    case 'shape'
        tol = 1e-6 + 1e-3 * max(abs(v1), abs(v2));
    case 'energy'
        if ~(isnumeric(J_ref) && isscalar(J_ref) && isfinite(J_ref) && J_ref > 0)
            error('invz:phase1Config', 'invz_phase1_gate: J_ref must be a finite positive scalar for kind=''energy''; got %s.', mat2str(J_ref));
        end
        tol = 1e-6 * J_ref + 1e-4 * max(abs(v1), abs(v2));
    otherwise
        error('invz:phase1Config', 'invz_phase1_gate: kind must be ''shape'' or ''energy''; got %s.', invz_phase1_safestr(kind));
end
pass = diff_val <= tol;
end

function s = invz_phase1_safestr(x)
if ischar(x), s = x; elseif isstring(x), s = char(x); else, s = mat2str(x); end
end
