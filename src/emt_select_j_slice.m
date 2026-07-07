function J_q = emt_select_j_slice(Jq_all, idx, n_cvar)
% EMT_SELECT_J_SLICE Select J(q) slice for one scan point.

if ndims(Jq_all) == 3
    J_q = Jq_all;
    return;
end

if ndims(Jq_all) ~= 4
    error('emt_select_j_slice:badDims', 'Jq must be 3-D or 4-D.');
end

dims = size(Jq_all);
if dims(3) == n_cvar
    J_q = squeeze(Jq_all(:,:,idx,:));
elseif dims(4) == n_cvar
    J_q = squeeze(Jq_all(:,:,:,idx));
else
    error('emt_select_j_slice:badShape', ...
        'Cannot infer cVar dimension in Jq size [%s].', num2str(dims));
end

if ndims(J_q) == 2
    J_q = reshape(J_q, 3, 3, 1);
end

end
