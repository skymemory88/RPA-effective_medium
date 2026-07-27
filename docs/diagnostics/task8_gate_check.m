ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));

beta = 1/(0.0862*0.31);
tl = struct('Delta', 0.02, 'M2', 0.8, 'm', 0.6, 'n01', tanh(0.02*beta/2), 'g0', 1);
tl.g0  = 2*tl.n01/tl.Delta;
lam = [0.01; 0.02];  Sigma0 = 0.05;
Jnu = linspace(-2e-3, 6e-3, 24).';
K0_seed = 0;  J0eff = 6.42444e-3;  G0inel0 = -300;  G0el0 = -20;

% --- strict path: out.Gtil0/out.r must equal a direct stable_form=true call at the same K0 ---
os = struct('warn', false, 'static_medium', 'strict_1z_dyson_ref');
[K0s, ~, outs] = invz_emt_static_ordered(tl, lam, Sigma0, Jnu, K0_seed, beta, J0eff, G0inel0, G0el0, os);
[~, go_true]  = invz_gstat_ordered(tl, lam, K0s, Sigma0, beta, G0inel0, G0el0, struct('stable_form', true));
[~, go_false] = invz_gstat_ordered(tl, lam, K0s, Sigma0, beta, G0inel0, G0el0, struct('stable_form', false));
fprintf('STRICT out.Gtil0 = %.17g\n', outs.Gtil0);
fprintf('STRICT direct stable_form=true  Gtil0 = %.17g  (match AbsTol0 = %d)\n', ...
    go_true.Gtil0, isequal(outs.Gtil0, go_true.Gtil0));
fprintf('STRICT direct stable_form=false Gtil0 = %.17g  (match AbsTol0 = %d)\n', ...
    go_false.Gtil0, isequal(outs.Gtil0, go_false.Gtil0));
fprintf('STRICT out.r = %.17g\n', outs.r);
fprintf('STRICT direct stable_form=true  r = %.17g  (match AbsTol0 = %d)\n', ...
    go_true.r, isequal(outs.r, go_true.r));

% --- resummed/legacy path: out.Gtil0/out.r must equal a direct PLAIN 7-arg call (stable_form
% default false), and must NOT be produced by an opts struct at all (arg count check below) ---
or_ = struct('warn', false);   % no static_medium field: default 'resummed'
[K0r, ~, outr] = invz_emt_static_ordered(tl, lam, Sigma0, Jnu, K0_seed, beta, J0eff, G0inel0, G0el0, or_);
[~, go_plain] = invz_gstat_ordered(tl, lam, K0r, Sigma0, beta, G0inel0, G0el0);   % 7-arg, legacy
fprintf('\nRESUMMED out.Gtil0 = %.17g\n', outr.Gtil0);
fprintf('RESUMMED direct 7-arg (legacy) Gtil0 = %.17g  (match AbsTol0 = %d)\n', ...
    go_plain.Gtil0, isequal(outr.Gtil0, go_plain.Gtil0));
fprintf('RESUMMED out.r = %.17g\n', outr.r);
fprintf('RESUMMED direct 7-arg (legacy) r = %.17g  (match AbsTol0 = %d)\n', ...
    go_plain.r, isequal(outr.r, go_plain.r));

% --- static source-code gate: count nargin()==8 call sites of invz_gstat_ordered inside the
% strict branch only (grep-level confirmation, printed for the record) ---
src = fileread(fullfile(ROOT, 'invz_projected', 'invz_emt_static_ordered.m'));
n8 = numel(regexp(src, 'invz_gstat_ordered\([^)]*stable_form', 'match'));
fprintf('\nSOURCE occurrences of an invz_gstat_ordered(...) call carrying stable_form: %d\n', n8);
