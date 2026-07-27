function [status, status_detail] = invz_hmf_status(pred, nodes, slope_pred, F, fallback_status)
%INVZ_HMF_STATUS Reduce Jensen node records with binding reason precedence.
% degenerate > reference-domain > node failure > unresolved > ok.
% G = -chi (meV^-1), ferromagnetic positive J.
%
% Pure precedence reducer (task 13): pred (the h = 0 predictor's per-node record) and
% nodes (a struct array of every other per-node record evaluated so far -- profile,
% extension, redensification, bisection iterate, or the final root, whichever the caller
% has in hand at its own early-return site) are both the SAME fixed schema
% (invz_hmf_ordered.m's blank_node_record), so they concatenate into one array and every
% record's .term_reason/.medium_status/.accepted can be scanned uniformly. slope_pred/F
% are the h = 0 mass and the profile's F = h0 - J0eff*m array, used ONLY for the
% floor-hit flavor of 'unresolved' (ordering predicted, no bracket found).
%
% Precedence, binding (highest first):
%   degenerate_doublet    any record's term_reason is 'degenerate_doublet' (the
%                         two-level Delta < 1e-4 meV domain screen, invz_twolevel_ordered
%                         opts.domain_policy = 'return')
%   medium_out_of_domain  any record's medium_status is neither 'ok' nor 'not_applicable'
%                         (a strict-scheme reference/closure domain event)
%   node_failed           any record did not converge/close (.accepted false), and
%                         neither of the two domain reasons above applies
%   unresolved            ordering was PREDICTED (slope_pred < 0) but no bracket was
%                         found above hmin_abs (all(F >= 0)), and nothing above applies
%   ok                    none of the above
%
% The optional second output is deliberately empty on 'ok'. On every non-ok outcome it is a
% compact, fixed-schema failure ledger plus one deterministic binding node. The first output and
% its precedence are unchanged. fallback_status is used only by the caller's root-refinement
% exhaustion path, whose local 'unresolved' verdict is not inferable from slope_pred/F alone.
if nargin < 5, fallback_status = ''; end
allnodes = [pred, nodes];
reasons = {allnodes.term_reason};
medium  = {allnodes.medium_status};
if any(strcmp(reasons, 'degenerate_doublet'))
    status = 'degenerate_doublet';
elseif any(~cellfun(@(s) any(strcmp(s, {'ok','not_applicable'})), medium))
    status = 'medium_out_of_domain';
elseif any(~[allnodes.accepted])
    status = 'node_failed';
elseif slope_pred < 0 && all(F >= 0)
    status = 'unresolved';
else
    status = 'ok';
end
if strcmp(status, 'ok') && ~isempty(fallback_status)
    status = fallback_status;
end

status_detail = [];
if strcmp(status, 'ok'), return; end

accepted = logical([allnodes.accepted]);
is_domain = ~cellfun(@(s) any(strcmp(s, {'ok','not_applicable'})), medium);
is_degenerate = strcmp(reasons, 'degenerate_doublet');
is_failed = ~accepted | is_domain | is_degenerate;

% Bracket membership is diagnostic only. F describes nodes(1:numel(F)); any extra node passed
% by the caller is a bisection/final-root candidate and carries its own .in_bracket flag.
in_bracket = false(1, numel(allnodes));
nf = min(numel(F), numel(nodes));
if nf > 1
    Fv = F(1:nf);
    crossings = find(isfinite(Fv(1:end-1)) & isfinite(Fv(2:end)) & ...
        ((Fv(1:end-1) <= 0 & Fv(2:end) >= 0) | ...
         (Fv(1:end-1) >= 0 & Fv(2:end) <= 0)));
    in_bracket(1 + crossings) = true;
    in_bracket(1 + crossings + 1) = true;
end
for k = 1:numel(allnodes)
    in_bracket(k) = in_bracket(k) || logical(field_or(allnodes(k), 'in_bracket', false));
end

% Binding precedence refines the public status precedence without changing it:
% first degenerate, first medium-domain event, first failed non-finite residual, otherwise
% the largest normalized failed residual (MATLAB max's first-index tie break is deterministic).
binding_index = find(is_degenerate, 1, 'first');
if isempty(binding_index), binding_index = find(is_domain, 1, 'first'); end
if isempty(binding_index)
    resid = nan(numel(allnodes), 5);
    for k = 1:numel(allnodes)
        resid(k,:) = [field_or(allnodes(k), 'resid_A', NaN), ...
                      field_or(allnodes(k), 'resid_B', NaN), ...
                      field_or(allnodes(k), 'resid_C', NaN), ...
                      field_or(allnodes(k), 'resid_D', NaN), ...
                      field_or(allnodes(k), 'resid_static', NaN)];
    end
    nonfinite_failed = is_failed(:) & any(~isfinite(resid), 2);
    binding_index = find(nonfinite_failed, 1, 'first');
end
if isempty(binding_index)
    score = -inf(1, numel(allnodes));
    for k = find(is_failed)
        score(k) = field_or(allnodes(k), 'resid_norm', -Inf);
        if ~isfinite(score(k)), score(k) = -Inf; end
    end
    [best, binding_index] = max(score);
    if ~isfinite(best), binding_index = find(is_failed, 1, 'first'); end
end

failed_indices = find(is_failed);
compact = repmat(blank_detail_node(), 1, numel(failed_indices));
for k = 1:numel(failed_indices)
    idx = failed_indices(k);
    compact(k) = compact_node(allnodes(idx), idx, idx == 1, in_bracket(idx));
end
if isempty(binding_index)
    binding_node = struct([]);
else
    binding_node = compact_node(allnodes(binding_index), binding_index, ...
                                binding_index == 1, in_bracket(binding_index));
end
status_detail = struct('n_nodes', numel(allnodes), 'n_accepted', nnz(accepted), ...
    'predictor_failed', ~accepted(1), 'binding_index', binding_index, ...
    'binding_node', binding_node, 'nodes', compact);
end

function out = compact_node(node, node_id, is_predictor, in_bracket)
out = blank_detail_node();
out.node_id = node_id;
out.h = field_or(node, 'h', NaN);
out.accepted = logical(field_or(node, 'accepted', false));
out.outer_iters = field_or(node, 'outer_iters', 0);
out.resid_A = field_or(node, 'resid_A', NaN);
out.resid_B = field_or(node, 'resid_B', NaN);
out.resid_C = field_or(node, 'resid_C', NaN);
out.resid_D = field_or(node, 'resid_D', NaN);
out.resid_static = field_or(node, 'resid_static', NaN);
out.medium_status = field_or(node, 'medium_status', 'not_applicable');
out.term_reason = field_or(node, 'term_reason', 'not_evaluated');
out.Dq_min = field_or(node, 'Dq_min', NaN);
out.Dq_abs_min = field_or(node, 'Dq_abs_min', NaN);
out.D_uni = field_or(node, 'D_uni', NaN);
out.ref_denom = field_or(node, 'ref_denom', NaN);
out.ref_margin = field_or(node, 'ref_margin', NaN);
out.is_predictor = is_predictor;
out.in_bracket = in_bracket;
end

function node = blank_detail_node()
node = struct('node_id', NaN, 'h', NaN, 'accepted', false, 'outer_iters', 0, ...
    'resid_A', NaN, 'resid_B', NaN, 'resid_C', NaN, 'resid_D', NaN, ...
    'resid_static', NaN, 'medium_status', 'not_applicable', ...
    'term_reason', 'not_evaluated', 'Dq_min', NaN, 'Dq_abs_min', NaN, ...
    'D_uni', NaN, 'ref_denom', NaN, 'ref_margin', NaN, ...
    'is_predictor', false, 'in_bracket', false);
end

function value = field_or(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), value = s.(name); else, value = default; end
end
