function status = invz_hmf_status(pred, nodes, slope_pred, F)
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
end
