function [resolved, detail] = invz_task2_ladder_ok(results_over_ladder)
%INVZ_TASK2_LADDER_OK Mesh size / offset / dpRng convergence gate (stage-2c task 2a; prereg
% docs/invzp_task2_prereg.md Sec. D, FROZEN 2026-07-23 -- values below cited VERBATIM).
%
% Requires BOTH, across every rung of the accepted lattice ladder (grid size, grid offset,
% dpRng), at EVERY physical field:
%   (1) IDENTICAL stability classification (the discrete class of Sec. A is the SAME at
%       each field across every rung), AND
%   (2) numeric convergence of s (and its components D_uni, min(Dq)) under the combined
%       test |Delta| <= AbsTol + RelTol*max(|.|,|.|), AbsTol = 1e-6, RelTol = 1e-4
%       (FROZEN, dimensionless -- s/D_uni/Dq_min are all O(1) static-pole quantities).
% A class disagreement (e.g. a sign flip) OR a numeric non-convergence anywhere returns
% resolved = false -- prereg Sec. F: this is the LATTICE/MESH-UNRESOLVED verdict, a hard
% stop preceding 3A/3B, "NOT evidence for 3A or 3B."
%
% results_over_ladder (struct array, one element per ladder rung -- e.g. {grid, offset,
%   dpRng} combinations -- each covering the SAME set of physical fields, in the SAME
%   order, so cross-rung comparison is field-aligned):
%   .class   (1xN cellstr)  the Sec.-A class string at each of the N physical fields
%   .s       (1xN double)   the stability scalar at each field
%   .D_uni   (1xN double)   the component at each field
%   .Dq_min  (1xN double)   the component at each field
%   .label   (char, OPTIONAL) free-form rung tag, echoed into detail for reporting only.
%
% resolved (logical). detail (struct):
%   .class_disagree    (1xN logical) which fields have a class mismatch across rungs
%   .numeric_disagree  (1xN logical) which fields fail the |Delta| convergence test
%                      (checked independently for s, D_uni, Dq_min; any one failing marks
%                      the field)
%   .worst             struct('s','D_uni','Dq_min') -- the largest (diff - tol) margin seen
%                      across all rung pairs for that quantity (>0 => some pair failed)
%   .n_rungs, .n_fields
AbsTol = 1e-6;  RelTol = 1e-4;   % prereg Sec. D, FROZEN -- do not tune.

nR = numel(results_over_ladder);
if nR < 2
    error('invz:task2LadderOk', 'need >= 2 ladder rungs to assess convergence; got %d.', nR);
end
nF = numel(results_over_ladder(1).class);
for r = 2:nR
    if numel(results_over_ladder(r).class) ~= nF
        error('invz:task2LadderOk', ...
            'rung %d has %d fields; rung 1 has %d -- ladder rungs must be field-aligned.', ...
            r, numel(results_over_ladder(r).class), nF);
    end
end

class_disagree   = false(1, nF);
numeric_disagree = false(1, nF);
worst = struct('s', -Inf, 'D_uni', -Inf, 'Dq_min', -Inf);
for f = 1:nF
    cls1 = results_over_ladder(1).class{f};
    for r = 2:nR
        if ~strcmp(results_over_ladder(r).class{f}, cls1)
            class_disagree(f) = true;
        end
    end
    for qn = {'s', 'D_uni', 'Dq_min'}
        q = qn{1};
        vals = arrayfun(@(rr) rr.(q)(f), results_over_ladder);
        for r1 = 1:nR
            for r2 = (r1+1):nR
                a = vals(r1);  b = vals(r2);
                tol = AbsTol + RelTol*max(abs(a), abs(b));
                if isfinite(a) && isfinite(b)
                    d = abs(a - b);
                    if d > tol, numeric_disagree(f) = true; end
                    worst.(q) = max(worst.(q), d - tol);
                else
                    numeric_disagree(f) = true;    % non-finite anywhere in the ladder: unresolved
                    worst.(q) = Inf;
                end
            end
        end
    end
end

resolved = ~any(class_disagree) && ~any(numeric_disagree);
detail = struct('class_disagree', class_disagree, 'numeric_disagree', numeric_disagree, ...
                'worst', worst, 'n_rungs', nR, 'n_fields', nF);
end
