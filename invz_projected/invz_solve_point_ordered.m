function pt = invz_solve_point_ordered(ion, T, Bx, Jnu_flat, opts)
%INVZ_SOLVE_POINT_ORDERED Jensen-consistent ordered 1/z solution at one point.
% Bx is a scalar transverse field along the crystallographic a axis (T).
% The applied-field/H_MF relation and ordered static EMT closure are always
% enforced; the retired bare-order-parameter route is not supported.
% The local functions below own both ordered-only stages so the complete
% ordered solve can be inspected and debugged from this one file.
%
% ODD (opt-in): opts.odd=true requires opts.odd_blocks with Vca/Vcb/Vcc/Jcc0
% from invz_odd_blocks and Jnu_flat=[]. The temperature/field-dependent modes
% and uniform -d shift are rebuilt here.
%
% Controlled approximation (opt-in):
%   opts.hmf_integral_mode = 'missing_area_approx'
%   opts.hmf_missing_area_factor = scalar >= 0.5
% uses only the terminal contiguous certified high-h component and sets
% A=factor*h_e*r(h_e).  The supported branch is explicitly the
% Picard-attracting contiguous high-h component.  This is an approximation,
% not a thermodynamic branch selector.  The default remains 'full_profile'.
%
% Diagnostic continuation (opt-in, no default effect):
%   opts.hmf_profile_state_seed = struct('Sigma',...,'lambda',...,'K0',...)
%   opts.hmf_profile_sweep_direction = 'descending'
% initializes the ordered iteration from a previously accepted state and
% evaluates the profile from high to low h. Every acceptance gate remains
% binding; a seed is never accepted as a result.
%   opts.hmf_profile_sweep_direction = 'descending_local'
% instead follows the same high-h component from high to low using only
% target-field states. It needs no cross-field seed and is the consistent
% field-map route when cold ascending profiles change their certified edge.
%
% Hyperfine-off boundary-linearized state (opt-in):
%   opts.hmf_J0z = finite positive scalar
% replaces only the nonlocal Jensen H_MF integral by the bracketed solution
% hmf = hmf_J0z*<Jz>.  invz_solve_auto supplies hmf_J0z from the rejected
% paramagnetic 1/z fixed point.  This approximation is deliberately rejected
% for hyp=true; the established electro-nuclear Jensen route is unchanged.
if nargin < 5, opts = struct(); end
if any(isfield(opts, {'odd_tier2','tier2','tol_tier2','max_tier2'}))
    error('invz:removedOption', 'The incomplete ODD Tier-2 route has been removed.');
end

if ~(isnumeric(Bx) && isreal(Bx) && isscalar(Bx) && isfinite(Bx))
    error('invz:field', 'Bx must be a finite real scalar transverse field.');
end
Bvec = [Bx 0 0];

Ecut  = getf(opts, 'Ecut', 40);
hyp   = getf(opts, 'hyp', true);
J0eff = getf(opts, 'J0eff', ion.J0eff);
Jxx0  = getf(opts, 'Jxx0', ion.Jxx0);
mixo  = getf(opts, 'mix_outer', 0.7);
tolo  = getf(opts, 'tol_outer', 1e-8);
maxo  = getf(opts, 'max_outer', 200);
eopts = getf(opts, 'emt', struct());

oddOn = isfield(opts, 'odd') && ~isempty(opts.odd) && ~isequal(opts.odd, false);
if oddOn
    ob = getf(opts, 'odd_blocks', []);
    if ~(isstruct(ob) && isscalar(ob) && all(isfield(ob, {'Vca', 'Vcb', 'Vcc', 'Jcc0'})))
        error('invz:oddArgs', ['opts.odd=true requires opts.odd_blocks with ' ...
            'Vca/Vcb/Vcc/Jcc0 from invz_odd_blocks.']);
    end
    if ~isempty(Jnu_flat)
        error('invz:oddArgs', 'opts.odd=true requires Jnu_flat=[]; modes are rebuilt internally.');
    end
    Xp = invz_chiperp(ion, T, Bx, struct('hyp', hyp, 'Jxx0', Jxx0));
    [dJ, d] = invz_odd_deltaJ(ob.Vca, ob.Vcb, Xp);
    Jnu_flat = reshape(invz_odd_modes(ob.Vcc, dJ), [], 1);
    J0eff = J0eff - d;
end

[wn, wts, beta] = invz_matsubara(T, Ecut);
[state_seed, state_seeded] = ordered_profile_state_seed(opts,numel(wn));

hopts = opts;
hopts.J0eff = J0eff;
for f = {'ordered_mode', 'forced_moment', 'transverse_mf', 'bz_tol'}
    if isfield(hopts, f{1}), hopts = rmfield(hopts, f{1}); end
end
linearized_hmf = isfield(opts, 'hmf_J0z');
if linearized_hmf
    if hyp
        error('invz:hmfLinearizedHyperfine', ...
            ['hmf_J0z is restricted to the explicit hyperfine-off route; ' ...
             'the electro-nuclear Jensen route must not be changed.']);
    end
    J0z_mf = opts.hmf_J0z;
    if ~(isnumeric(J0z_mf) && isreal(J0z_mf) && isscalar(J0z_mf) && ...
            isfinite(J0z_mf) && J0z_mf > 0)
        error('invz:hmfLinearizedCoupling', ...
            'hmf_J0z must be a finite positive real scalar.');
    end
    hmf_sigma0 = getf(opts, 'hmf_sigma0', NaN);
    if ~(isnumeric(hmf_sigma0) && isreal(hmf_sigma0) && ...
            isscalar(hmf_sigma0) && isfinite(hmf_sigma0) && ...
            (1 + hmf_sigma0) > 0)
        error('invz:hmfLinearizedSigma', ...
            'hmf_sigma0 must be finite with 1+hmf_sigma0 > 0.');
    end
    [si, mfphase, ~, mfout] = invz_bare_mf_state( ...
        ion, T, Bx, J0z_mf, Jxx0, false);
    hprof = linearized_hmf_profile(J0z_mf, hmf_sigma0, mfphase, mfout);
    if mfphase ~= 1
        pt = early_return(si, hprof);
        return;
    end
    hstar = si.hz;
else
    [hstar, hprof] = invz_hmf_ordered(ion, T, Bx, Jnu_flat, hopts);
    if ~isfinite(hstar)
        si = invz_single_ion(ion, T, Bvec, struct('hyp', hyp, 'hz_fixed', 0, 'Jxx0', Jxx0));
        pt = early_return(si, hprof);
        return;
    end

    si = invz_single_ion(ion, T, Bvec, ...
        struct('hyp', hyp, 'hz_fixed', hstar, 'Jxx0', Jxx0));
end
m0 = si.Jexp(3);

c0 = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
G0 = -real(squeeze(c0(3,3,:)));
tl = invz_twolevel_ordered(ion, T, Bx, si.hz, struct('Jxx0', Jxx0));
g  = real(invz_g(tl, 1i*wn));

% Full-electronuclear static split for the Jensen elastic closure.
eso = getf(opts, 'emt_static', struct());
eso.warn = false;
c0i = invz_chi0z(si, T, 1i*wn(1), struct('elastic', false));
G0inel0 = -real(c0i(3,3,1));
X = real(c0(:, :, 1));
feedback = X(3,1) * (Jxx0 / (1 - Jxx0*X(1,1))) * X(1,3);
G0bare0 = -(X(3,3) + feedback);
G0el0 = G0bare0 - G0inel0;

Sigma = zeros(size(wn));
K = zeros(size(wn));
K0s = 0;
lam = [0; 0; 0];
if state_seeded
    Sigma = state_seed.Sigma;
    K0s = state_seed.K0;
    lam = state_seed.lambda;
    if isfield(hprof, 'target_carrier_fallback_used') && ...
            hprof.target_carrier_fallback_used && ...
            isfield(hprof, 'Sigma_star') && ...
            numel(hprof.Sigma_star) == numel(wn) && ...
            all(isfinite(hprof.Sigma_star(:))) && ...
            isfield(hprof, 'K0_star') && isfinite(hprof.K0_star)
        % The legacy source carrier failed and H_MF refinement converged only
        % through the certified target-field fallback. Keep that carrier for
        % the final ordered closure; legacy-success retries remain unchanged.
        Sigma = hprof.Sigma_star(:);
        K0s = hprof.K0_star;
    end
elseif isfield(hprof, 'profile_sweep_direction') && ...
        strcmp(hprof.profile_sweep_direction, 'descending_local') && ...
        isfield(hprof, 'Sigma_star') && ...
        numel(hprof.Sigma_star) == numel(wn) && ...
        all(isfinite(hprof.Sigma_star(:))) && ...
        isfield(hprof, 'K0_star') && isfinite(hprof.K0_star)
    % The local descending profile has already converged the target-field
    % root without any cross-field state. Preserve that certified carrier in
    % the final ordered closure so the profile and exported state use one path.
    Sigma = hprof.Sigma_star(:);
    K0s = hprof.K0_star;
elseif linearized_hmf && isfield(opts, 'Sigma_seed') && ...
        numel(opts.Sigma_seed) == numel(wn) && all(isfinite(opts.Sigma_seed(:)))
    Sigma = opts.Sigma_seed(:);
end
converged = false;
med = struct('G', nan(size(wn)), 'converged', false);
sg = struct('alpha', NaN, 'alpha_m', NaN);
for outer = 1:maxo
    eopts.K0 = K;
    med = invz_emt_scalar(G0, Sigma, Jnu_flat, eopts);
    K = med.K;

    [K0s, ~, sout] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), ...
        Jnu_flat, K0s, beta, J0eff, G0inel0, G0el0, eso);
    if ~sout.converged
        break; % the outer map is undefined without a certified physical static root
    end
    K(1) = K0s;

    lam = invz_lambdas(K, g, wts, beta, [1 2 3]);
    sg = invz_sigma_ordered(tl, lam, K, g, beta);
    dS = max(abs(sg.Sigma - Sigma));
    Sigma = Sigma + mixo*(sg.Sigma - Sigma);
    if dS < tolo && sout.converged
        converged = true;
        break;
    end
end

% Refresh the static closure at the exported self-energy.
[K0s, Gstat, so] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), ...
    Jnu_flat, K0s, beta, J0eff, G0inel0, G0el0, eso);
K(1) = K0s;
ctol = getf(eso, 'resid_tol', 1e-10);
staticok = so.converged && isfinite(so.resid) && so.resid < ctol;

lam_check = invz_lambdas(K, g, wts, beta, [1 2 3]);
Sigma_check = invz_sigma_ordered(tl, lam_check, K, g, beta);
final_resid = max(abs(Sigma_check.Sigma - Sigma));
med.G(1) = Gstat;

pt.m0 = m0;
pt.is_ordered = true;
pt.Sigma0 = Sigma(1);
pt.alpha = sg.alpha;
pt.alpha_m = sg.alpha_m;
pt.lambda = lam;
pt.K = K;
pt.G = med.G;
pt.Sigma = Sigma;
pt.tl = tl;
pt.si = si;
pt.chi0cc0 = -G0(1);
pt.crit = 1 + pt.Sigma0 - J0eff*pt.chi0cc0;
pt.sumrule_rel = abs(sum(wts.*med.G)/beta + si.JzJz_fluct) ...
    / max(abs(si.JzJz_fluct), 1e-12);
pt.final_resid = final_resid;
pt.converged = converged && med.converged && staticok && final_resid < tolo;
pt.outer_iters = outer;
pt.hmf = hstar;
pt.hmf_prof = hprof;
pt.hmf_status = hprof.status;
pt.production_approximation = hprof.production_approximation;
pt.approximation_only = hprof.approximation_only;
pt.D_uni = 1 + (J0eff - K(1))*med.G(1);
pt.static_status = so.status;
if linearized_hmf
    pt.hmf_J0z = J0z_mf;
    pt.hmf_sigma0 = hmf_sigma0;
    pt.hmf_residual = hstar - J0z_mf*m0;
    if ~(isfinite(pt.hmf_residual) && abs(pt.hmf_residual) < 1e-10)
        error('invz:hmfLinearizedResidual', ...
            'Boundary-linearized H_MF residual %.6g meV failed.', pt.hmf_residual);
    end
end
if oddOn
    pt.odd = struct('d', d, 'Xp', Xp);
end
if abs(pt.si.hz - hstar) > 1e-12
    error('invz:hzFixed', 'Jensen final state did not hold hmf: %.6g vs %.6g.', pt.si.hz, hstar);
end
end

function prof = linearized_hmf_profile(J0z_mf, sigma0, mfphase, mfout)
% Minimal profile-shaped diagnostics for the explicitly approximate route.
if mfphase == 1
    status = 'ok_linearized_pm_handoff';
else
    status = 'no_linearized_order';
end
prof = struct('status', status, 'integral_mode', 'linearized_pm_handoff', ...
    'approximation_only', true, 'production_approximation', true, ...
    'missing_area', NaN, 'hmf_J0z', J0z_mf, 'hmf_sigma0', sigma0, ...
    'bare_mf_phase', mfphase, 'bare_mf_mass_pm', mfout.mass_pm, ...
    'bare_mf_method', mfout.method, ...
    'bare_mf_root_residual', mfout.root_residual);
end

function pt = early_return(si, hprof)
pt = struct('m0', 0, 'is_ordered', false, 'converged', false, ...
    'Sigma0', NaN, 'crit', NaN, 'si', si, 'tl', [], ...
    'hmf', NaN, 'hmf_prof', hprof, 'hmf_status', hprof.status, ...
    'production_approximation',hprof.production_approximation, ...
    'approximation_only',hprof.approximation_only, ...
    'D_uni', NaN, 'final_resid', NaN, 'static_status', 'not_evaluated');
end

function [hmf_star, prof] = invz_hmf_ordered(ion, T, Bx, Jnu_flat, opts)
%INVZ_HMF_ORDERED Jensen applied-field/H_MF self-consistency, spontaneous root (SS9.3, J 2.31-2.33).
%   h0(hmf) = int_0^hmf r(h') dh',   r = G0(0;h')/Gtil0(0;h')
% with Gtil0 built on the STATIC-CLOSURE K0 (invz_emt_static_ordered, P0-2), evaluated on
% fixed-field single-ion states (hz_fixed WITHOUT order -- P0-1, invz:hzFixed asserted).
% Spontaneous condition (zero applied longitudinal field): h0(hmf) = J0eff*<Jz>(hmf); the
% nonzero root is bracketed on a GEOMETRIC profile clustered at 0 (P1-4) and refined by
% bisection with direct node evaluations to opts.tol_root. F(h)/h -> crit as h -> 0+
% (SS5), returned as prof.slope0. Returns NaN when no nonzero root exists, or when the
% separate bare (order=true) bracketing solve does not order.
%
% Status contract (round-5 P2, binding): 'unresolved' means ordering was PREDICTED
% (slope0 < 0) but no bracket was found above hmin_abs, OR the bracket did not refine
% to tol_root. 'node_failed' means ANY node evaluated along the way -- the h=0
% predictor, a profile node, a bisection iterate, or the final root evaluation --
% failed to converge/close. Both map to a NaN hmf_star and MUST be read by callers as
% converged = false, never as a PM label.
if nargin < 5, opts = struct(); end
if ~(isnumeric(Bx) && isreal(Bx) && isscalar(Bx) && isfinite(Bx))
    error('invz:field', 'Bx must be a finite real scalar transverse field.');
end
Bvec = [Bx 0 0];
J0eff = opts.J0eff;                                  % required, no default (caller-owned)
Jxx0  = getf(opts, 'Jxx0', ion.Jxx0);
hyp   = getf(opts, 'hyp', true);
nH    = getf(opts, 'nH', 33);
hfac  = getf(opts, 'hmax_fac', 1.25);
hfrac = getf(opts, 'hmin_frac', 1e-3);
trt   = getf(opts, 'tol_root', 1e-3);
mixo  = getf(opts, 'mix_outer', 0.7);
tolo  = getf(opts, 'tol_outer', 1e-8);
maxo  = getf(opts, 'max_outer', 200);                % ALIGNED with both solvers' default
eso   = getf(opts, 'emt_static', struct());          % static-closure opts, threaded (P1-F)
eso.warn = false;   % node loop gates on so.converged; suppress the per-node console flood
integral_mode = string(getf(opts, 'hmf_integral_mode', 'full_profile'));
missing_area_approx = integral_mode == "missing_area_approx";
strict_profile = integral_mode == "full_profile";
approximation_only = missing_area_approx;
if ~(strict_profile || missing_area_approx)
    error('invz:hmfIntegralMode', ...
        ['hmf_integral_mode must be ''full_profile'' or ' ...
         '''missing_area_approx''.']);
end
missing_area_factor = getf(opts,'hmf_missing_area_factor',NaN);
approximation_branch = string(getf(opts,'hmf_approx_branch', ...
    'picard_attracting_contiguous_high_h_component'));
supported_approximation_branch = ...
    "picard_attracting_contiguous_high_h_component";
if missing_area_approx
    if ~(isnumeric(missing_area_factor) && isreal(missing_area_factor) && ...
            isscalar(missing_area_factor) && isfinite(missing_area_factor) && ...
            missing_area_factor >= 0.5)
        error('invz:hmfMissingAreaFactor', ...
            ['missing_area_approx requires hmf_missing_area_factor to be ' ...
             'a finite scalar >= 0.5 for a nonnegative linear completion.']);
    end
    if approximation_branch ~= supported_approximation_branch
        error('invz:hmfApproxBranch', ...
            ['The only implemented approximation branch is ' ...
             '''picard_attracting_contiguous_high_h_component''.']);
    end
else
    missing_area_factor = NaN;
    approximation_branch = "not_applicable";
end

% Fixed-field nodes do not re-apply the ordering update.
sibase = struct('hyp', hyp, 'Jxx0', Jxx0);
hmin_abs = getf(opts, 'hmin_abs', NaN);              % resolved after hmax below (P1-C)

prof = struct('hgrid', [], 'r', [], 'h0', [], 'm', [], 'Sigma0', [], 'K0', [], ...
              'D_uni', [], 'G0bare', [], 'Gstat', [], 'node_conv', [], 'F', [], ...
              'static_status', strings(1,0), 'predictor_static_status', "not_evaluated", ...
              'static_status_star', "not_evaluated", ...
              'refinement_failure_status', "not_evaluated", ...
              'refinement_primary_failure_status', "not_evaluated", ...
              'target_carrier_fallback_used', false, ...
              'refinement_fallback_used', false, ...
              'refinement_fallback_kind', "not_used", ...
              'slope0', NaN, 'Sigma0_pm0', NaN, 'K0_pm0', NaN, 'J0eff', J0eff, ...
              'n_extend', 0, 'hmin_initial', NaN, 'status', 'no_bare_order', ...
              'redensified', false, ...
              'predictor_converged', false, 'converged_node_count', 0, ...
              'm_star', NaN, 'D_uni_star', NaN, 'r_star', NaN, 'Gstat_star', NaN, ...
              'Sigma_star', [], 'K0_star', NaN, ...
              'integral_mode', char(integral_mode), ...
              'approximation_only',approximation_only, ...
              'production_approximation',missing_area_approx, ...
              'missing_area_factor',missing_area_factor, ...
              'missing_area',NaN,'missing_area_model',"not_applicable", ...
              'missing_area_integral',struct(), ...
              'approximation_branch',char(approximation_branch), ...
              'branch_selection_status',"not_thermodynamically_selected", ...
              'bare_mf_phase', 0, 'bare_mf_mass_pm', NaN, ...
              'bare_mf_method', "not_evaluated", ...
              'bare_mf_root_residual', NaN, ...
              'root_bracket_indices', [NaN NaN], ...
              'root_bracket_count',0, ...
              'root_bracket_bridged', false);
hmf_star = NaN;

% Bracket ceiling from the selected BARE MF state. The PM mass decides the
% phase; on its ordered side a bracketed fixed-hz solve avoids the critical
% slowing of direct ordered Picard iteration near the bare-MF boundary.
[sib, bare_mf_phase, bare_mf_mass, bare_mf] = invz_bare_mf_state( ...
    ion, T, Bx, J0eff, Jxx0, hyp);
prof.bare_mf_phase = bare_mf_phase;
prof.bare_mf_mass_pm = bare_mf_mass;
prof.bare_mf_method = bare_mf.method;
prof.bare_mf_root_residual = bare_mf.root_residual;
if bare_mf_phase ~= 1, return; end                    % bare state is paramagnetic
hmax = hfac * abs(sib.hz);
if isnan(hmin_abs), hmin_abs = 1e-10*hmax; end

% --- Matsubara grid, weights, beta: MIRROR invz_solve_point_ordered's setup block
% verbatim (wn, wts, beta, eopts from opts -- honors Ecut and EMT options, P1-6).
Ecut  = getf(opts, 'Ecut', 40);
eopts = getf(opts, 'emt', struct());
[wn, wts, beta] = invz_matsubara(T, Ecut);
[profile_state_seed, profile_state_seeded] = ...
    ordered_profile_state_seed(opts,numel(wn));
profile_sweep_direction = string(getf(opts, ...
    'hmf_profile_sweep_direction','ascending'));
if ~isscalar(profile_sweep_direction) || ...
        ~ismember(profile_sweep_direction, ...
        ["ascending" "descending" "descending_local"])
    error('invz:hmfProfileSweepDirection', ...
        ['hmf_profile_sweep_direction must be ''ascending'', ''descending'', ' ...
         'or ''descending_local''.']);
end
if profile_sweep_direction == "descending" && ~profile_state_seeded
    error('invz:hmfProfileStateSeed', ...
        'A descending profile sweep requires hmf_profile_state_seed.');
end
if profile_sweep_direction == "descending_local" && profile_state_seeded
    error('invz:hmfProfileStateSeed', ...
        ['A descending_local profile is target-local and does not accept ' ...
         'hmf_profile_state_seed; use descending for a cross-field retry.']);
end
profile_lambda_seed = zeros(3,1);
if profile_state_seeded
    profile_lambda_seed = profile_state_seed.lambda;
end
prof.profile_state_seeded = profile_state_seeded;
prof.profile_sweep_direction = char(profile_sweep_direction);

% Independent h = 0 PM predictor node (round-3 P0-3; doubles as Gate 6b's comparator):
% ONE node solve at hz_fixed = 0 gives THIS machinery's PM fixed point. Its mass
%   slope_pred = r(0) + J0eff*G0bare(0) = 1 + Sigma0(0) - J0eff*chi_path(0)   (= crit, SS5)
% predicts root existence INDEPENDENTLY of any sampled profile value.
Sigma = [];  K0s = 0;                                % last accepted node state
if profile_state_seeded
    Sigma = profile_state_seed.Sigma;
    K0s = profile_state_seed.K0;
end
[r0n, ~, S0pm, K0pm, ~, Gb0, ~, ok0, Sigma_candidate, K0_candidate, st0] = ...
    eval_node(0, Sigma, K0s);
if ok0
    Sigma = Sigma_candidate;
    K0s = K0_candidate;
end
prof.predictor_converged = ok0;
prof.predictor_static_status = st0;
predictor_usable = ok0;
% Under the default strict rule, a predictor-node convergence failure is NOT one of the three enumerated
% 'node_failed' triggers (round-5 P2: profile/bisection/final-evaluation), and it is
% NOT 'unresolved' either (that label presupposes a computed slope_pred). It is its
% own case: h0(h) = int_0^h r seeds on r(0), so without a converged predictor the
% cumulative integral (hence F, hence the root search) is categorically undefined --
% but the per-node grid quantities below remain well-defined direct
% diagonalizations. The overall verdict is forced to node_failed/NaN.
if predictor_usable
    slope_pred = r0n + J0eff*Gb0;
    prof.Sigma0_pm0 = S0pm;  prof.K0_pm0 = K0pm;  prof.slope0 = slope_pred;
else
    slope_pred = NaN;
end

ratio = hfrac^(1/(nH-1));
hgrid = hmax * ratio.^((nH-1):-1:0);                 % geometric, clustered at 0 (P1-4)
prof.hmin_initial = hgrid(1);

[rv, mv, S0v, K0v, Dv, Gbv, Gsv, cnv, Sigma, K0s, stv, ...
    target_high_Sigma, target_high_K0] = run_sweep(hgrid, Sigma, K0s);

if missing_area_approx
    integral_eligible = cnv & isfinite(rv) & isfinite(mv);
    [~, component_meta] = invz_missing_area_integral( ...
        hgrid,rv,integral_eligible,1);
    if component_meta.node_count >= 2 && ...
            component_meta.component_edge > 0 && component_meta.edge_r > 0
        missing_area = missing_area_factor* ...
            component_meta.component_edge*component_meta.edge_r;
        [h0, missing_meta] = invz_missing_area_integral( ...
            hgrid,rv,integral_eligible,missing_area);
        missing_meta.branch_selector = approximation_branch;
        prof.missing_area = missing_area;
        prof.missing_area_model = "linear_completion_shape_factor";
    else
        h0 = nan(size(hgrid));
        missing_meta = component_meta;
        missing_meta.missing_area = NaN;
        missing_meta.branch_selector = approximation_branch;
    end
    prof.missing_area_integral = missing_meta;
elseif predictor_usable
    h0 = cumtrapz([0 hgrid], [r0n rv]);  h0 = h0(2:end); % first panel seeded with r(0)
else
    h0 = nan(1, nH);                                     % undefined without a real r(0)
end
F  = h0 - J0eff*mv;

% ADAPTIVE lower extension (round-3 P0-3): predictor-driven, NOT self-referential.
% slope_pred < 0 predicts an ordered root; extend geometrically downward until a
% negative F sample appears or the absolute floor is reached.
n_extend = 0;
while strict_profile && predictor_usable && slope_pred < 0 && ...
        all(F >= 0) && hgrid(1) > hmin_abs
    n_extend = n_extend + 1;
    hext = hgrid(1) * ratio.^(3:-1:1);                % three more decades-fraction nodes
    [re, me, S0e, K0e, De, Gbe, Gse] = deal(nan(1, 3));  ce = false(1, 3);
    ste = strings(1,3);
    for k = 1:3
        [re(k), me(k), S0e(k), K0e(k), De(k), Gbe(k), Gse(k), ce(k), ...
            Sigma_candidate, K0_candidate, ste(k)] = eval_node(hext(k), Sigma, K0s);
        if ce(k)
            Sigma = Sigma_candidate;
            K0s = K0_candidate;
        end
    end
    hgrid = [hext hgrid];  rv = [re rv];  mv = [me mv];  cnv = [ce cnv]; %#ok<AGROW>
    stv = [ste stv]; %#ok<AGROW>
    S0v = [S0e S0v];  K0v = [K0e K0v];  Dv = [De Dv];  Gbv = [Gbe Gbv];  Gsv = [Gse Gsv]; %#ok<AGROW>
    h0 = cumtrapz([0 hgrid], [r0n rv]);  h0 = h0(2:end);
    F  = h0 - J0eff*mv;
end

if strict_profile && n_extend > 0 && any(F < 0)
    % RE-DENSIFY (execution amendment 3, 2026-07-22): the extension's sparse geometric
    % panels feed O(coarse-grid) quadrature error into h0 exactly where F is a small
    % difference of large terms (measured: 11% root error at Bc_1z - 0.01 on a
    % deliberately coarse grid vs the fine default). Rebuild the profile at FULL nH
    % resolution anchored to the discovered bracket scale, so adaptive-path roots match
    % default-path quality. Cost: one extra nH-sweep, only when extension fired.
    idx0 = find(F < 0, 1, 'first');
    hfrac_eff = max(hmin_abs/hmax, 0.25*hgrid(idx0)/hmax);
    ratio2 = hfrac_eff^(1/(nH-1));
    hgrid = hmax * ratio2.^((nH-1):-1:0);
    [rv, mv, S0v, K0v, Dv, Gbv, Gsv, cnv, Sigma, K0s, stv, ...
        target_high_Sigma, target_high_K0] = run_sweep(hgrid, Sigma, K0s);
    h0 = cumtrapz([0 hgrid], [r0n rv]);  h0 = h0(2:end);
    F  = h0 - J0eff*mv;
    prof.redensified = true;
end
prof.n_extend = n_extend;

prof.hgrid = hgrid;  prof.r = rv;  prof.h0 = h0;  prof.m = mv;
prof.Sigma0 = S0v;   prof.K0 = K0v;  prof.D_uni = Dv;  prof.node_conv = cnv;  prof.F = F;
prof.G0bare = Gbv;   prof.Gstat = Gsv;  prof.static_status = stv;
prof.converged_node_count = nnz(cnv);

if missing_area_approx && prof.missing_area_integral.node_count < 2
    prof.status = 'missing_area_no_certified_component';
    return;
elseif missing_area_approx && ~isfinite(prof.missing_area)
    prof.status = 'missing_area_invalid_edge';
    return;
elseif strict_profile && ~predictor_usable % predictor never produced a usable finite value
    prof.status = 'node_failed'; % above -- report the honest verdict now that the grid's
    return;                      % own (convergence-independent) diagnostics are exported;
end                              % NEVER fall through to the F-based search on NaN data
if strict_profile && slope_pred < 0 && all(F >= 0) % floor hit without a bracket:
    prof.status = 'unresolved';                       % NEVER silently PM (round-3 P0-3)
    warning('invz:hmfUnresolved', ...
        'ordering predicted (slope_pred = %.3g) but no negative F above hmin_abs = %.3g', ...
        slope_pred, hmin_abs);
    return;                                           % hmf_star stays NaN; the jensen solver
end                                                   % must return converged = false here
if strict_profile && (any(~cnv) || any(~isfinite([rv, mv, F])))
    prof.status = 'node_failed';                      % on node failure -- never 'ok'
    return;
end
if missing_area_approx
    valid = prof.missing_area_integral.selected_mask & ...
        isfinite(mv) & isfinite(F);
    pair = valid(1:end-1) & valid(2:end) & ...
        F(1:end-1) < 0 & F(2:end) >= 0;
    brackets = find(pair);
    prof.root_bracket_count = numel(brackets);
    prof.status = 'ok_missing_area_approx';
    if isempty(brackets)
        prof.status = 'missing_area_no_bracket';
        return;
    elseif numel(brackets) > 1
        prof.status = 'missing_area_multiple_brackets';
        return;
    end
    idx = brackets(1);
    ia = idx;
    ib = idx+1;
else
    prof.status = 'ok';
    s = sign(F);
    idx = find(s(1:end-1) < 0 & s(2:end) >= 0, 1, 'last');
    if isempty(idx), return; end                      % no nonzero root: PM side
    ia = idx;
    ib = idx+1;
end
prof.root_bracket_indices = [ia ib];
prof.root_bracket_bridged = ib > ia+1;

% --- Root refinement by DIRECT evaluation (P1-4): bisection on F between the
% bracketing nodes, fresh node solve per iterate, cumulative h0 via local trapezoid
% panel from the bracket's left node.
% Preserve the established source-seed refinement as the primary path. A
% descending retry also owns a certified target-field carrier at the low edge
% of its accepted component (Sigma/K0s returned by run_sweep). Use that carrier
% only if the legacy source seed loses the static root during direct refinement.
% A target-local descending sweep uses its component-edge carrier first. At the
% exact zero-field degeneracy that carrier can lose the static root even though
% the independently certified high-h endpoint remains valid; retain that
% endpoint as a same-field fallback. Prior successful paths remain untouched.
target_edge_Sigma = Sigma;
target_edge_K0 = K0s;
ncarrier = 1;
if profile_state_seeded && profile_sweep_direction == "descending"
    ncarrier = 2;
elseif profile_sweep_direction == "descending_local" && ...
        ~isempty(target_high_Sigma) && all(isfinite(target_high_Sigma)) && ...
        isfinite(target_high_K0)
    ncarrier = 2;
end
refined = false;
for carrier_attempt = 1:ncarrier
    seeded_descending = profile_state_seeded && ...
        profile_sweep_direction == "descending";
    local_descending = profile_sweep_direction == "descending_local";
    if carrier_attempt == 1 && seeded_descending
        Sigma = profile_state_seed.Sigma;
        K0s = profile_state_seed.K0;
    elseif carrier_attempt == 2 && local_descending
        Sigma = target_high_Sigma;
        K0s = target_high_K0;
    else
        Sigma = target_edge_Sigma;
        K0s = target_edge_K0;
    end
    a = hgrid(ia);  b = hgrid(ib);  Fa = F(ia);  h0a = h0(ia);  ra = rv(ia);
    failure_stage = "none";
    failure_status = "not_evaluated";
    for it = 1:12
        c = 0.5*(a + b);
        [rc, mc, ~, ~, ~, ~, ~, okc, Sigma_candidate, K0_candidate, stc] = ...
            eval_node(c, Sigma, K0s);
        if ~okc
            failure_stage = "refinement";
            failure_status = stc;
            break;
        end
        Sigma = Sigma_candidate;
        K0s = K0_candidate;
        h0c = h0a + 0.5*(ra + rc)*(c - a);
        Fc  = h0c - J0eff*mc;
        if sign(Fc) == sign(Fa), a = c; Fa = Fc; h0a = h0c; ra = rc; else, b = c; end
        if (b - a) < trt*b, break; end
    end
    if failure_stage == "none" && (b - a) >= trt*b
        prof.status = 'unresolved';  hmf_star = NaN;
        warning('invz:hmfUnresolved', ...
            'root bracket not refined to tol_root: (b-a)/b = %.3g', (b-a)/b);
        return;
    end
    if failure_stage == "none"
        hmf_candidate = 0.5*(a + b);
        [r_s, m_s, ~, ~, D_s, ~, Gs_s, ok_s, Sigma_star, K0_star, st_s] = ...
            eval_node(hmf_candidate, Sigma, K0s);
        if ~ok_s
            failure_stage = "final";
            failure_status = st_s;
        end
    end
    if failure_stage ~= "none"
        if carrier_attempt == 1 && ncarrier == 2
            prof.refinement_primary_failure_status = failure_status;
            continue;
        end
        prof.refinement_failure_status = failure_status;
        if missing_area_approx && failure_stage == "refinement"
            prof.status = 'missing_area_refinement_failed';
        elseif missing_area_approx
            prof.status = 'missing_area_final_failed';
        else
            prof.status = 'node_failed';
        end
        hmf_star = NaN;
        return;
    end
    hmf_star = hmf_candidate;
    prof.static_status_star = st_s;
    prof.target_carrier_fallback_used = ...
        seeded_descending && carrier_attempt == 2;
    prof.refinement_fallback_used = carrier_attempt == 2;
    if prof.refinement_fallback_used && local_descending
        prof.refinement_fallback_kind = "target_high_node";
    elseif prof.refinement_fallback_used
        prof.refinement_fallback_kind = "target_component_edge";
    end
    refined = true;
    break;
end
if ~refined
    error('invz:hmfRefinementInvariant', ...
        'H_MF refinement exited without a success or classified failure.');
end
prof.m_star = m_s;  prof.D_uni_star = D_s;  prof.r_star = r_s;  prof.Gstat_star = Gs_s;
prof.Sigma_star = Sigma_star;
prof.K0_star = K0_star;

    function [rv, mv, S0v, K0v, Dv, Gbv, Gsv, cnv, Sigma, K0s, stv, ...
            Sigma_high, K0_high] = run_sweep(hgrid, Sigma, K0s)
    % Evaluate in the declared diagnostic direction, warm-starting only from
    % the last accepted node. A failed candidate is retained in the diagnostic
    % arrays but is never committed as the next node's seed.
    n = numel(hgrid);
    [rv, mv, S0v, K0v, Dv, Gbv, Gsv] = deal(nan(1, n));  cnv = false(1, n);
    stv = strings(1,n);
    Sigma_high = [];
    K0_high = NaN;
    if profile_sweep_direction == "ascending"
        sweep_indices = 1:n;
    else
        sweep_indices = n:-1:1;
    end
    for is = sweep_indices
        [rv(is), mv(is), S0v(is), K0v(is), Dv(is), Gbv(is), Gsv(is), cnv(is), ...
            Sigma_candidate, K0_candidate, stv(is)] = eval_node(hgrid(is), Sigma, K0s);
        if cnv(is)
            Sigma = Sigma_candidate;
            K0s = K0_candidate;
            if is == n
                Sigma_high = Sigma_candidate;
                K0_high = K0_candidate;
            end
        end
    end
    end

    function [rk, mk, S0k, K0k, Dk, Gbk, Gsk, ok, Sigma, K0s, node_status] = eval_node(hp, Sigma, K0s)
    % One fixed-field node: si (hz_fixed, NO order), tl, c0/G0, then the ordered
    % Sigma<->EMT loop WITH the static-sector closure each pass (Interfaces bullet).
    % This function computes a candidate from value-input seeds. Its caller
    % commits the returned Sigma/K0s only if ok=true, so a failed node cannot
    % contaminate later profile nodes while its last iterate remains inspectable.
    sio = sibase;  sio.hz_fixed = hp;
    si = invz_single_ion(ion, T, Bvec, sio);
    if abs(si.hz - hp) > 1e-12
        error('invz:hzFixed', 'hz_fixed not held: si.hz = %.6g vs %.6g', si.hz, hp);
    end
    mk = si.Jexp(3);
    % A near-degenerate electronic doublet is outside the retained vertex
    % domain. Report that node as rejected without throwing or mutating the
    % accepted continuation carrier. This is load-bearing for the opt-in
    % missing-area mode, which may still possess a certified terminal high-h
    % component; strict full_profile remains fail-closed because its predictor
    % and complete-path gates still require every node.
    tl = invz_twolevel_ordered(ion, T, Bx, hp, ...
        struct('Jxx0', Jxx0, 'domain_policy', 'return'));
    if ~tl.valid
        [rk,S0k,K0k,Dk,Gbk,Gsk] = deal(NaN);
        ok = false;
        node_status = "twolevel_domain_invalid";
        return;
    end
    c0 = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
    G0 = -real(squeeze(c0(3,3,:)));
    c0i = invz_chi0z(si, T, 1i*wn(1), struct('elastic', false));   % static inelastic only
    G0inel0 = -real(c0i(3,3,1));                                   % fixed-Hamiltonian slot
    % Path-consistent bare static response -dm/dh along the single-axis
    % transverse mean-field path.
    X = real(c0(:, :, 1));                                         % static chi tensor (chi = -G)
    fb = X(3, 1) * (Jxx0 / (1 - Jxx0*X(1, 1))) * X(1, 3);
    G0bare0 = -(X(3, 3) + fb);
    G0el0   = G0bare0 - G0inel0;                                   % elastic + feedback (SS4a)
    g  = real(invz_g(tl, 1i*wn));
    if isempty(Sigma), Sigma = zeros(size(wn)); end
    K = zeros(size(wn));  lam = profile_lambda_seed;  ok = false;
    for outer = 1:maxo
        % (1) dynamic sector -- MIRROR invz_solve_point_ordered's emt call verbatim
        eopts.K0 = K;
        med = invz_emt_scalar(G0, Sigma, Jnu_flat, eopts);
        K   = med.K;
        % (2) static sector (P0-2/P0-A), threaded opts (P1-F):
        [K0s, ~, sout] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), Jnu_flat, K0s, ...
                                                 beta, J0eff, G0inel0, G0el0, eso);
        if ~sout.converged
            break; % never propagate a failed/pseudo-root static value into lambdas or Sigma
        end
        K(1) = K0s;
        % (3)-(5) lambdas, ordered Sigma, damped mix -- MIRROR the solver's statements
        lam = invz_lambdas(K, g, wts, beta, [1 2 3]);
        sg  = invz_sigma_ordered(tl, lam, K, g, beta);
        dS  = max(abs(sg.Sigma - Sigma));
        Sigma = Sigma + mixo*(sg.Sigma - Sigma);
        if dS < tolo && sout.converged, ok = true; break; end
    end
    [K0s, Gsk, so] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), Jnu_flat, K0s, ...
                                             beta, J0eff, G0inel0, G0el0, eso);
    ctol = getf(eso, 'resid_tol', 1e-10);              % documented closure tolerance (meV^-1),
    % matching invz_emt_static_ordered's own default so the outer gate never disagrees
    % with the inner closure's own converged flag.
    ok = ok && so.converged && isfinite(so.resid) && so.resid < ctol;
    if ~so.converged
        node_status = string(so.status);
    elseif ~ok
        node_status = "outer_iteration_failed";
    else
        node_status = "ok";
    end
    % round-5 P1-B: the final refresh must ITSELF converge and close -- an unconverged
    % refresh makes this node not-ok (callers then mark node_failed), never silent export.
    % round-4 P1-B: the final refresh runs on the newly mixed Sigma(1), so its closed K0
    % differs from the seed -- KEEP it, and report the SAME value the returned
    % Gstat/r/D_uni were computed with (exported below as K0k = K0s).
    rk = so.r;  S0k = Sigma(1);  K0k = K0s;  Dk = so.D_uni;  Gbk = G0bare0;
    end
end

function [seed, seeded] = ordered_profile_state_seed(opts,nwn)
raw = getf(opts,'hmf_profile_state_seed',[]);
seed = struct('Sigma',zeros(nwn,1),'lambda',zeros(3,1),'K0',0);
seeded = ~isempty(raw);
if ~seeded
    return;
end
if ~(isstruct(raw) && isscalar(raw) && ...
        all(isfield(raw,{'Sigma','lambda'})))
    error('invz:hmfProfileStateSeed', ...
        'hmf_profile_state_seed must contain Sigma and lambda.');
end
sigma = raw.Sigma(:);
lambda = raw.lambda(:);
if ~(isnumeric(sigma) && isreal(sigma) && numel(sigma) == nwn && ...
        all(isfinite(sigma)) && isnumeric(lambda) && isreal(lambda) && ...
        numel(lambda) >= 3 && all(isfinite(lambda)))
    error('invz:hmfProfileStateSeed', ...
        'Seed Sigma/lambda must be finite real vectors matching the solve.');
end
K0 = 0;
if isfield(raw,'K0')
    K0 = raw.K0;
    if ~(isnumeric(K0) && isreal(K0) && isscalar(K0) && isfinite(K0))
        error('invz:hmfProfileStateSeed', ...
            'Seed K0 must be a finite real scalar.');
    end
end
seed = struct('Sigma',sigma,'lambda',lambda(1:3),'K0',K0);
end
