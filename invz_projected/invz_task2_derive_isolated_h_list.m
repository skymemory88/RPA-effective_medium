function h_list = invz_task2_derive_isolated_h_list(swept_out)
%INVZ_TASK2_DERIVE_ISOLATED_H_LIST The h-values an 'isolated' companion cell should re-solve,
% derived from its 'swept' sibling's OWN checkpointed output (stage-2c task 2b-driver; prereg
% docs/invzp_task2_prereg.md Sec. G item 1, brief G1/G2: "isolated vs swept ... node converges
% alone but not in continuation?"; invz_task2_run_config.m's own run_isolated_mode docstring:
% "re-solves the SAME h a swept run visited, isolated/cold, via exactly this path").
%
% Rule: h_list = unique([swept_out.nodes.h], 'stable') -- EVERY h the swept profile visited
% (predictor/sweep/extend/redensify/bisect/root phases alike, per invz_hmf_ordered.m's own
% phase taxonomy), de-duplicated, in first-visited order. This is a judgement call flagged
% for review (the prereg/brief specify the COMPARISON, not a literal h-selection formula), but
% it is the only choice that answers the question "does the SAME node converge alone" in
% EVERY case -- including the (expected, per stage2c-context.md/the invzp-jensen-realcoupling-
% nonconvergence memory) case where the swept profile MASKS (hstar non-finite): using the
% full node list rather than just hstar means this function needs no special-casing for a
% masked vs. converged swept run, and never silently produces an empty/undefined h_list.
%
% swept_out  the `out` struct invz_task2_run_config(cfg) returns for a cfg.mode='swept' cell.
%            Errors ('invz:task2Matrix') if swept_out.nodes is empty (e.g. the swept cell
%            itself was an absorbed-error failed-config record, out.status='failed') --
%            callers (invz_task2_resolve_cell_cfg.m) check for this BEFORE calling here, so
%            reaching this error indicates a driver bug, not an expected physical outcome.
%
% h_list  column vector of h-values, in first-visited order (NOT sorted by value -- the
%         VISIT order is itself part of the isolated-vs-swept provenance a human/2b-report
%         may want to see, e.g. which node index in the swept trace each isolated companion
%         corresponds to).
if isempty(swept_out.nodes)
    error('invz:task2Matrix', ['invz_task2_derive_isolated_h_list: swept_out.nodes is empty -- ' ...
        'cannot derive an isolated h_list from a swept cell with no node records (check ' ...
        'out.status/out.err_id first).']);
end
h_all = [swept_out.nodes.h];
h_list = unique(h_all(:), 'stable');
end
