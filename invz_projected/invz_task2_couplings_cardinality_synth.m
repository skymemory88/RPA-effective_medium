function [Jnu_flat, detail] = invz_task2_couplings_cardinality_synth(base_vals, n_total)
%INVZ_TASK2_COUPLINGS_CARDINALITY_SYNTH Deterministic cardinality-matched duplication of a
% small synthetic coupling fixture up to a target cardinality (stage-2c task 2b-driver;
% prereg docs/invzp_task2_prereg.md Sec. G item 4, brief G4: "duplicate the 24-pt synthetic
% distribution Jnu=linspace(-2e-3,6e-3,24) up to the real cardinality (16384) ... isolates
% cardinality from distribution shape. Lattice-diagnostic; not an existence endpoint.").
%
% Construction (NO RNG): tile base_vals as many WHOLE times as fit within n_total
% (n_rep = floor(n_total/n_base)), then pad with the FIRST n_rem = n_total - n_rep*n_base
% entries of base_vals to reach EXACTLY n_total. Fully deterministic: the same base_vals/
% n_total always produces bit-identical output, and the construction uses only integer
% division/indexing, never a random draw.
%
% base_vals  numeric vector (e.g. the pinned 24-pt fixture linspace(-2e-3,6.0e-3,24)).
% n_total    target cardinality (e.g. 16384, the real 16^3-grid x 4-branch flattened count).
%
% Jnu_flat  [n_total x 1] column vector containing ONLY values from base_vals (no new values
%           introduced), each appearing floor(n_total/n_base) or floor(n_total/n_base)+1 times.
% detail    struct('n_base','n_total','n_rep','n_rem') -- documents the exact tiling/padding
%           split used, for the caller/test to assert against directly.
base_vals = base_vals(:);
n_base = numel(base_vals);
if n_base < 1
    error('invz:task2Config', 'invz_task2_couplings_cardinality_synth: base_vals must be non-empty.');
end
if ~(isscalar(n_total) && isfinite(n_total) && n_total == round(n_total) && n_total >= n_base)
    error('invz:task2Config', ['invz_task2_couplings_cardinality_synth: n_total must be an ' ...
        'integer >= numel(base_vals) (%d); got %s.'], n_base, mat2str(n_total));
end
n_rep = floor(n_total / n_base);
n_rem = n_total - n_rep * n_base;
Jnu_flat = [repmat(base_vals, n_rep, 1); base_vals(1:n_rem)];
detail = struct('n_base', n_base, 'n_total', n_total, 'n_rep', n_rep, 'n_rem', n_rem);
end
