% effective_medium.m
% Entry script for rewritten EMT backbone.
%
% Required workspace variables:
%   scanMode, freq_total, chi0, qvec, Jq
%
% Optional workspace variables:
%   chiq, cVar, fields, temp, dscrt_var
%   USE_SINGLE_ION_SEED (logical)
%   EMT_PARAMS (struct, see src/emt_default_params.m)
%
% Outputs written to workspace:
%   chi_emt, K_emt, G_local_emt, G_q_emt, effective_medium_info
%   effective_medium_validation, effective_medium_inputs
%
% Notes:
% - This implementation uses a modular EMT solver in src/.
% - Optional Track-A renormalization (lambda/X) is available via
%   EMT_PARAMS.track_a.enabled = true.
% - The Track-A X model is configurable (moment-based ansatz by default).

required = {'scanMode', 'freq_total', 'chi0', 'qvec', 'Jq'};
for ii = 1:numel(required)
    if ~exist(required{ii}, 'var')
        error('effective_medium:missingInput', ...
            'Missing required workspace variable: %s', required{ii});
    end
end

repo_src = fullfile(pwd, 'src');
if exist(repo_src, 'dir') == 7
    if isempty(strfind([path pathsep], [repo_src pathsep]))
        addpath(repo_src);
    end
else
    error('effective_medium:missingSrc', 'src/ directory not found.');
end

data = struct();
data.scanMode = scanMode;
data.freq_total = freq_total;
data.chi0 = chi0;
data.qvec = qvec;
data.Jq = Jq;

if exist('chiq', 'var')
    data.chiq = chiq;
end
if exist('cVar', 'var')
    data.cVar = cVar;
end
if exist('fields', 'var')
    data.fields = fields;
end
if exist('temp', 'var')
    data.temp = temp;
end
if exist('dscrt_var', 'var')
    data.dscrt_var = dscrt_var;
end
if exist('USE_SINGLE_ION_SEED', 'var')
    data.use_single_ion_seed = logical(USE_SINGLE_ION_SEED);
else
    data.use_single_ion_seed = false;
end

if exist('EMT_PARAMS', 'var')
    params = EMT_PARAMS;
else
    params = struct();
end

outputs = emt_run_from_workspace(data, params);

chi_emt = outputs.chi_emt;
K_emt = outputs.K_emt;
G_local_emt = outputs.G_local_emt;
G_q_emt = outputs.G_q_emt;

effective_medium_info = outputs.info;
effective_medium_validation = outputs.validation;
effective_medium_inputs = outputs.inputs;

if isfield(params, 'enable_plot')
    enable_plot = logical(params.enable_plot);
else
    enable_plot = true;
end

if enable_plot && exist('plot_comparison', 'file') == 2
    plot_comparison( ...
        outputs.inputs.cVar, outputs.inputs.freq, outputs.inputs.chi_seed, ...
        chi_emt, outputs.inputs.scan_mode, outputs.inputs.dscrt_value);
end
