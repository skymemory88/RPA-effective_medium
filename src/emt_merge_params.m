function params = emt_merge_params(defaults, overrides)
% EMT_MERGE_PARAMS Merge user overrides into default EMT parameters.

params = defaults;
if nargin < 2 || isempty(overrides)
    return;
end

if isfield(overrides, 'mixing_alpha') && ~isfield(overrides, 'mix_alpha')
    overrides.mix_alpha = overrides.mixing_alpha;
end
if isfield(overrides, 'G_damp') && ~isfield(overrides, 'mix_alpha')
    overrides.mix_alpha = max(0, min(1, 1 - overrides.G_damp));
end
if isfield(overrides, 'trackA') && ~isfield(overrides, 'track_a')
    overrides.track_a = overrides.trackA;
end

params = merge_struct_recursive(params, overrides);
params.mix_alpha = max(min(params.mix_alpha, 1), 1e-3);
params.min_mix_alpha = max(min(params.min_mix_alpha, params.mix_alpha), 1e-3);

if ~isfield(params, 'track_a') || isempty(params.track_a)
    params.track_a = defaults.track_a;
else
    params.track_a = merge_struct_recursive(defaults.track_a, params.track_a);
end

end

function out = merge_struct_recursive(base, update)
    out = base;
    fields = fieldnames(update);
    for ifield = 1:numel(fields)
        key = fields{ifield};
        val = update.(key);
        if isstruct(val) && isfield(out, key) && isstruct(out.(key))
            out.(key) = merge_struct_recursive(out.(key), val);
        else
            out.(key) = val;
        end
    end
end
