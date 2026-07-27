ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));

% Exact fixture from test_out_of_domain_reference_is_a_status
beta = 1/(0.0862*0.31);
tl = struct('Delta', 0.02, 'M2', 0.8, 'm', 0.6, 'n01', tanh(0.02*beta/2), 'g0', 1);
tl.g0  = 2*tl.n01/tl.Delta;
lam = [0.01; 0.02];  Sigma0 = -1;                     % 1 + Sigma0 = 0 -> ref_denom_nonpositive
Jnu = linspace(-2e-3, 6e-3, 24).';
K0_seed = 0;  J0eff = 6.42444e-3;  G0inel0 = -300;  G0el0 = -20;
o = struct('warn', false, 'static_medium', 'strict_1z_dyson_ref');

[K0, Gstat, out] = invz_emt_static_ordered(tl, lam, Sigma0, Jnu, K0_seed, beta, J0eff, ...
                                           G0inel0, G0el0, o);

fprintf('MEASURED medium_status = %s\n', out.medium_status);
fprintf('MEASURED K0 = %.17g\n', K0);
fprintf('MEASURED Gstat = %.17g\n', Gstat);
fprintf('MEASURED D_uni = %.17g\n', out.D_uni);
fprintf('MEASURED Dq_min = %.17g\n', out.Dq_min);
fprintf('MEASURED Dq_max = %.17g\n', out.Dq_max);
fprintf('MEASURED Dq_neg_count = %d\n', out.Dq_neg_count);
fprintf('MEASURED resid = %.17g\n', out.resid);
fprintf('MEASURED converged = %d\n', out.converged);
fprintf('MEASURED iters = %d\n', out.iters);
fprintf('MEASURED class(Dq_neg_count) = %s\n', class(out.Dq_neg_count));
