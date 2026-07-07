function params = emt_default_params()
% EMT_DEFAULT_PARAMS Return default parameters for EMT solver.

params = struct();
params.frequency_domain = 'real';
params.eta = 1e-3;
params.mix_alpha = 0.15;
params.tol = 1e-6;
params.max_iter = 1000;
params.reg_epsilon = 1e-8;
params.rcond_tol = 1e-10;
params.verbose = false;
params.use_neighbor_seed = true;
params.closure_tol = 1e-3;
params.enable_plot = true;
params.store_gq = false;
params.init_k = 'one_pass';
params.backoff_factor = 0.5;
params.backoff_trigger = 1.05;
params.min_mix_alpha = 0.02;

track_a = struct();
track_a.enabled = false;
track_a.component = 1;
track_a.apply_mode = 'channel';
track_a.model = 'moment_ratio';
track_a.beta = 1.0;
track_a.M2 = 1.0;
track_a.n01 = 1.0;
track_a.n0 = [];
track_a.n1 = [];
track_a.p_list = 1:4;
track_a.frequency_weights = [];
track_a.clamp_scale = 10.0;
track_a.linear_a = [0 1 0 0 0];
track_a.linear_gamma = [0 0 1 0 0];
track_a.custom_x_function = [];
params.track_a = track_a;

end
