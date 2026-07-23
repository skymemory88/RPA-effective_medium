function res = invz_phase1_refinement_gate(kind, v1, v2, v3, J_ref)
%INVZ_PHASE1_REFINEMENT_GATE Item 5's three-rung refinement determination (docs/invzp_phase1_
% quadrature_prereg.md item 5, "Refinement"): given a statistic's value at three successive
% resolution rungs (v1=coarsest, v2=middle, v3=finest -- e.g. grid N=12/16/20 at fixed dpRng, or
% cutoff dpRng=30/40/50 at fixed N), gate CONVERGENCE on the FINEST comparison only (v2 vs v3, via
% invz_phase1_gate) AND additionally require the spread not to grow relative to the preceding
% step: |v3-v2| <= |v2-v1| ("final-step spread <= preceding-step spread" -- "Requiring every step
% to individually meet the final tolerance would wrongly demand the coarsest rung already be
% converged.").
%
% kind      'shape' (v's already normalized: stat/J_ref or stat/J_ref^2) | 'energy' (v's raw meV).
% v1,v2,v3  the statistic's value at the coarse/mid/fine rungs (finite scalars).
% J_ref     required for kind=='energy' (passed through to invz_phase1_gate); pass [] for 'shape'.
%
% OUTPUT res (struct):
%   .step_coarse_mid       |v2-v1|, the "preceding step" spread.
%   .step_mid_fine         |v3-v2|, the "final step" spread (== the finest-comparison difference).
%   .finest_pass           logical, from invz_phase1_gate(kind,v2,v3,J_ref).
%   .finest_diff, .finest_tol   the invz_phase1_gate diff/tol pair for the finest comparison.
%   .spread_nonincreasing  logical, step_mid_fine <= step_coarse_mid (a tiny relative epsilon,
%                          1e-9*max(1,step_coarse_mid), absorbs floating-point roundoff sitting at
%                          near-equality -- this loosens only the "<=" comparison's own numerical
%                          noise floor, not either compared value or any frozen tolerance).
%   .pass                  finest_pass && spread_nonincreasing -- item 5's full refinement PASS
%                          for this statistic/rung-triple.
if nargin < 5, J_ref = []; end
if ~(isnumeric(v1) && isscalar(v1) && isfinite(v1) && isnumeric(v2) && isscalar(v2) && isfinite(v2) ...
        && isnumeric(v3) && isscalar(v3) && isfinite(v3))
    error('invz:phase1Config', 'invz_phase1_refinement_gate: v1,v2,v3 must be finite scalars.');
end
[finest_pass, finest_diff, finest_tol] = invz_phase1_gate(kind, v2, v3, J_ref);
step_coarse_mid = abs(v2 - v1);
step_mid_fine   = abs(v3 - v2);
EPS_REL = 1e-9 * max(1, step_coarse_mid);

res.step_coarse_mid = step_coarse_mid;
res.step_mid_fine   = step_mid_fine;
res.finest_pass = finest_pass;
res.finest_diff = finest_diff;
res.finest_tol  = finest_tol;
res.spread_nonincreasing = step_mid_fine <= step_coarse_mid + EPS_REL;
res.pass = res.finest_pass && res.spread_nonincreasing;
end
