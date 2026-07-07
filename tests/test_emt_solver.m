function tests = test_emt_solver
% TEST_EMT_SOLVER Function-based tests for rewritten EMT backbone.

tests = functiontests(localfunctions);

end

function setupOnce(~)
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
end

function test_noninteracting_limit(testCase)
data = synthetic_data();
data.Jq = zeros(3, 3, size(data.qvec, 1));
data.use_single_ion_seed = true;

params = emt_default_params();
params.max_iter = 60;
params.tol = 1e-10;
params.mix_alpha = 0.25;
params.use_neighbor_seed = false;
params.enable_plot = false;
params.store_gq = false;

outputs = emt_run_from_workspace(data, params);
chi_ref = data.chi0;

verifyLessThan(testCase, max(abs(outputs.chi_emt(:) - chi_ref(:))), 1e-7);
verifyTrue(testCase, all(outputs.converged));
end

function test_k_symmetry(testCase)
data = synthetic_data();
data.use_single_ion_seed = true;

params = emt_default_params();
params.max_iter = 40;
params.tol = 1e-7;
params.mix_alpha = 0.2;
params.use_neighbor_seed = false;
params.enable_plot = false;

outputs = emt_run_from_workspace(data, params);

max_asym = 0;
for ic = 1:size(outputs.K_emt, 4)
    for iw = 1:size(outputs.K_emt, 3)
        k = outputs.K_emt(:,:,iw,ic);
        max_asym = max(max_asym, norm(k - k', 'fro'));
    end
end

verifyLessThan(testCase, max_asym, 1e-7);
end

function test_scanmode_temperature_alias(testCase)
data = synthetic_data();
data.scanMode = 'temperature';
data.temp = [0.4 0.8];
data.fields = 0.35;

params = emt_default_params();
params.max_iter = 5;
params.enable_plot = false;

outputs = emt_run_from_workspace(data, params);
verifyEqual(testCase, outputs.info.scan_mode, 'temp');
end

function data = synthetic_data()
qvec = [
    0.0 0.0 0.0;
    0.1 0.0 0.0;
    0.0 0.1 0.0;
    0.0 0.0 0.1;
    0.1 0.1 0.0
];

n_w = 12;
n_c = 2;
freq = linspace(0, 4, n_w);
chi0 = zeros(3, 3, n_w, n_c);
for ic = 1:n_c
    for iw = 1:n_w
        base = 0.15 + 0.02 * ic + 0.01 * iw;
        chi0(:,:,iw,ic) = [
            base + 0.01i, 0, 0;
            0, 0.8 * base + 0.02i, 0;
            0, 0, 0.6 * base + 0.03i
        ];
    end
end

Jq = zeros(3, 3, size(qvec, 1));
for iq = 1:size(qvec, 1)
    x = qvec(iq, 1);
    y = qvec(iq, 2);
    z = qvec(iq, 3);
    Jq(:,:,iq) = [
        0.04 + 0.01 * x, 0.004, 0;
        0.004, 0.03 + 0.01 * y, 0.003;
        0, 0.003, 0.02 + 0.01 * z
    ];
    Jq(:,:,iq) = 0.5 * (Jq(:,:,iq) + Jq(:,:,iq).');
end

data = struct();
data.scanMode = 'field';
data.freq_total = freq;
data.chi0 = chi0;
data.chiq = repmat(chi0, 1, 1, 1, 1, size(qvec, 1));
data.qvec = qvec;
data.Jq = Jq;
data.fields = [0.2 0.4];
data.temp = 0.6;
end
