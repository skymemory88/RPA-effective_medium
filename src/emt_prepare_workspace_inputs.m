function inputs = emt_prepare_workspace_inputs(data)
% EMT_PREPARE_WORKSPACE_INPUTS Normalize workspace data for EMT solver.
%
% Required fields in data:
%   scanMode, freq_total, chi0, qvec, Jq
% Optional fields:
%   chiq, cVar, fields, temp, dscrt_var, use_single_ion_seed

required = {'scanMode', 'freq_total', 'chi0', 'qvec', 'Jq'};
for ii = 1:numel(required)
    key = required{ii};
    if ~isfield(data, key) || isempty(data.(key))
        error('emt_prepare_workspace_inputs:missingField', ...
            'Missing required input field: %s', key);
    end
end

inputs = struct();
inputs.scan_mode = normalize_scan_mode(data.scanMode);
inputs.freq = data.freq_total(:);
inputs.qvec = data.qvec;
inputs.Jq = data.Jq;
inputs.jq_report = validate_jq_input(data.Jq);
if ~inputs.jq_report.valid
    error('emt_prepare_workspace_inputs:badJq', inputs.jq_report.message);
end
if inputs.jq_report.max_asymmetry > 1e-8
    warning('emt_prepare_workspace_inputs:jAsymmetry', ...
        'Jq asymmetry is high (%.3e).', inputs.jq_report.max_asymmetry);
end
inputs.n_q = size(data.qvec, 1);

if isfield(data, 'use_single_ion_seed')
    use_single_ion_seed = logical(data.use_single_ion_seed);
else
    use_single_ion_seed = false;
end
inputs.use_single_ion_seed = use_single_ion_seed;

inputs.cVar = resolve_cvar(data, inputs.scan_mode);
inputs.dscrt_value = resolve_dscrt_value(data, inputs.scan_mode);
inputs.chi_seed = build_seed_tensor(data, use_single_ion_seed, inputs.n_q);
seed_nq = size(inputs.chi_seed, 5);
if seed_nq ~= inputs.n_q
    if seed_nq == 1
        inputs.chi_seed = repmat(inputs.chi_seed, 1, 1, 1, 1, inputs.n_q);
    else
        error('emt_prepare_workspace_inputs:seedQMismatch', ...
            'chi seed q-dimension (%d) does not match qvec length (%d).', ...
            seed_nq, inputs.n_q);
    end
end

n_cvar_seed = size(inputs.chi_seed, 4);
if numel(inputs.cVar) ~= n_cvar_seed
    n_keep = min(numel(inputs.cVar), n_cvar_seed);
    warning('emt_prepare_workspace_inputs:cVarMismatch', ...
        ['cVar length (%d) and seed dimension (%d) differ. ', ...
         'Using first %d points.'], numel(inputs.cVar), n_cvar_seed, n_keep);
    inputs.cVar = inputs.cVar(1:n_keep);
    inputs.chi_seed = inputs.chi_seed(:,:,:,1:n_keep,:);
end

inputs.n_cVar = numel(inputs.cVar);
inputs.n_omega = size(inputs.chi_seed, 3);

end

function mode = normalize_scan_mode(mode_in)
    if isstring(mode_in)
        mode_in = char(mode_in);
    end
    mode = lower(strtrim(mode_in));
    if strcmp(mode, 'temperature')
        mode = 'temp';
    end
    if ~strcmp(mode, 'field') && ~strcmp(mode, 'temp')
        error('emt_prepare_workspace_inputs:badScanMode', ...
            'Unsupported scanMode: %s', mode_in);
    end
end

function cvar = resolve_cvar(data, mode)
    if strcmp(mode, 'field')
        if isfield(data, 'fields') && ~isempty(data.fields)
            cvar = data.fields;
            return;
        end
    else
        if isfield(data, 'temp') && ~isempty(data.temp)
            cvar = data.temp;
            return;
        end
    end

    if isfield(data, 'cVar') && ~isempty(data.cVar)
        cvar = data.cVar;
        return;
    end

    if strcmp(mode, 'field')
        error('emt_prepare_workspace_inputs:missingFieldAxis', ...
            'Provide fields or cVar for field scan mode.');
    end

    error('emt_prepare_workspace_inputs:missingTempAxis', ...
        'Provide temp or cVar for temp scan mode.');
end

function dscrt_value = resolve_dscrt_value(data, mode)
    dscrt_value = NaN;

    if strcmp(mode, 'field')
        if isfield(data, 'temp') && ~isempty(data.temp)
            dscrt_value = data.temp;
            return;
        end
    else
        if isfield(data, 'fields') && ~isempty(data.fields)
            dscrt_value = data.fields;
            return;
        end
    end

    if isfield(data, 'dscrt_var') && ~isempty(data.dscrt_var)
        dscrt_value = data.dscrt_var;
    end
end

function chi_seed = build_seed_tensor(data, use_single_ion_seed, n_q)
    if use_single_ion_seed
        chi0 = data.chi0;
        if ndims(chi0) == 4
            chi_seed = repmat(chi0, 1, 1, 1, 1, n_q);
        elseif ndims(chi0) == 5
            chi_seed = chi0;
        else
            error('emt_prepare_workspace_inputs:badChi0Dims', ...
                'chi0 must be 4-D or 5-D array.');
        end
        return;
    end

    if ~isfield(data, 'chiq') || isempty(data.chiq)
        error('emt_prepare_workspace_inputs:missingChiq', ...
            'chiq is required unless use_single_ion_seed is true.');
    end

    chiq = data.chiq;
    if ndims(chiq) == 4
        chi_seed = reshape(chiq, size(chiq, 1), size(chiq, 2), ...
            size(chiq, 3), size(chiq, 4), 1);
    else
        chi_seed = chiq;
    end
end
