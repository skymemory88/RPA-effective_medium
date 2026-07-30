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

hopts = opts;
hopts.J0eff = J0eff;
for f = {'ordered_mode', 'forced_moment', 'transverse_mf', 'bz_tol'}
    if isfield(hopts, f{1}), hopts = rmfield(hopts, f{1}); end
end
[hstar, hprof] = invz_hmf_ordered(ion, T, Bx, Jnu_flat, hopts);
if ~isfinite(hstar)
    si = invz_single_ion(ion, T, Bvec, struct('hyp', hyp, 'hz_fixed', 0, 'Jxx0', Jxx0));
    pt = early_return(si, hprof);
    return;
end

si = invz_single_ion(ion, T, Bvec, ...
    struct('hyp', hyp, 'hz_fixed', hstar, 'Jxx0', Jxx0));
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
pt.D_uni = 1 + (J0eff - K(1))*med.G(1);
pt.static_status = so.status;
if oddOn
    pt.odd = struct('d', d, 'Xp', Xp);
end
if abs(pt.si.hz - hstar) > 1e-12
    error('invz:hzFixed', 'Jensen final state did not hold hmf: %.6g vs %.6g.', pt.si.hz, hstar);
end
end

function pt = early_return(si, hprof)
pt = struct('m0', 0, 'is_ordered', false, 'converged', false, ...
    'Sigma0', NaN, 'crit', NaN, 'si', si, 'tl', [], ...
    'hmf', NaN, 'hmf_prof', hprof, 'hmf_status', hprof.status, ...
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
endpoint_visual = integral_mode == "endpoint_trapezoid_visual";
filtered_visual = integral_mode == "filtered_profile_visual";
strict_profile = integral_mode == "full_profile";
visual_only = endpoint_visual || filtered_visual;
if ~(strict_profile || endpoint_visual || filtered_visual)
    error('invz:hmfIntegralMode', ...
        ['hmf_integral_mode must be ''full_profile'', ' ...
         '''endpoint_trapezoid_visual'', or ''filtered_profile_visual''.']);
end
filtered_anchor_label = "not_applicable";
if filtered_visual, filtered_anchor_label = "pm_endpoint"; end
allow_unconverged_pm = getf(opts, 'hmf_endpoint_allow_unconverged_pm', false);
needs_pm_endpoint = endpoint_visual || filtered_visual;
filtered_include_unconverged = getf(opts, ...
    'hmf_filtered_include_unconverged',false);
if ~(isscalar(filtered_include_unconverged) && ...
        (islogical(filtered_include_unconverged) || ...
         (isnumeric(filtered_include_unconverged) && ...
          isfinite(filtered_include_unconverged) && ...
          any(filtered_include_unconverged == [0 1]))))
    error('invz:hmfFilteredOption', ...
        'hmf_filtered_include_unconverged must be a finite logical/scalar.');
end
filtered_include_unconverged = logical(filtered_include_unconverged);

% Fixed-field nodes do not re-apply the ordering update.
sibase = struct('hyp', hyp, 'Jxx0', Jxx0);
hmin_abs = getf(opts, 'hmin_abs', NaN);              % resolved after hmax below (P1-C)

prof = struct('hgrid', [], 'r', [], 'h0', [], 'm', [], 'Sigma0', [], 'K0', [], ...
              'D_uni', [], 'G0bare', [], 'Gstat', [], 'node_conv', [], 'F', [], ...
              'static_status', strings(1,0), 'predictor_static_status', "not_evaluated", ...
              'static_status_star', "not_evaluated", ...
              'refinement_failure_status', "not_evaluated", ...
              'slope0', NaN, 'Sigma0_pm0', NaN, 'K0_pm0', NaN, 'J0eff', J0eff, ...
              'n_extend', 0, 'hmin_initial', NaN, 'status', 'no_bare_order', ...
              'redensified', false, ...
              'predictor_converged', false, 'converged_node_count', 0, ...
              'm_star', NaN, 'D_uni_star', NaN, 'r_star', NaN, 'Gstat_star', NaN, ...
              'integral_mode', char(integral_mode), 'visual_only', visual_only, ...
              'endpoint_r0', NaN, 'endpoint_r0_source', "not_evaluated", ...
              'endpoint_pm_converged', false, 'endpoint_pm_iters', NaN, ...
              'endpoint_pm_crit', NaN, 'filtered_integral', struct(), ...
              'filtered_anchor', char(filtered_anchor_label), ...
              'filtered_include_unconverged', filtered_include_unconverged, ...
              'filtered_unconverged_used_count', 0, ...
              'root_bracket_indices', [NaN NaN], ...
              'root_bracket_bridged', false);
hmf_star = NaN;

% Bracket ceiling from the BARE ordered fixed point: SAME MF option base plus
% order=true and J0z (P1-F -- the bracket runs under the caller's MF knobs too).
sibo = sibase;  sibo.order = true;  sibo.J0z = J0eff;
sib = invz_single_ion(ion, T, Bvec, sibo);
if ~(sib.mf_converged && abs(sib.Jexp(3)) > 1e-6), return; end     % bare does not order
hmax = hfac * abs(sib.hz);
if isnan(hmin_abs), hmin_abs = 1e-10*hmax; end

% --- Matsubara grid, weights, beta: MIRROR invz_solve_point_ordered's setup block
% verbatim (wn, wts, beta, eopts from opts -- honors Ecut and EMT options, P1-6).
Ecut  = getf(opts, 'Ecut', 40);
eopts = getf(opts, 'emt', struct());
[wn, wts, beta] = invz_matsubara(T, Ecut);

% Independent h = 0 PM predictor node (round-3 P0-3; doubles as Gate 6b's comparator):
% ONE node solve at hz_fixed = 0 gives THIS machinery's PM fixed point. Its mass
%   slope_pred = r(0) + J0eff*G0bare(0) = 1 + Sigma0(0) - J0eff*chi_path(0)   (= crit, SS5)
% predicts root existence INDEPENDENTLY of any sampled profile value.
Sigma = [];  K0s = 0;                                % warm-start carriers across nodes
[r0n, ~, S0pm, K0pm, ~, Gb0, ~, ok0, Sigma, K0s, st0] = eval_node(0, Sigma, K0s);
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

% VISUAL-ONLY PM endpoint used by the two-endpoint approximation and,
% optionally, by the filtered positive-h profile:
%   h0(h) = h*[r_PM(0)+r(h)]/2              (two endpoints)
%   h0(h_k) = trapz([0,h_valid<=h_k],r)     (filtered profile)
% At m=0, r_PM(0)=1+Sigma0 exactly. The PM fixed point is evaluated
% independently so failure of the strict ordered h=0 static closure does not
% contaminate this endpoint. An unconverged but finite PM last iterate is
% permitted only behind the explicit visual-only option; it is exported as
% such and never changes the default full-profile contract.
endpoint_usable = false;
endpoint_r0 = NaN;
if needs_pm_endpoint
    pmopts = opts;
    remove = {'odd','odd_blocks','odd_retarded','odd_retarded_exact', ...
        'odd_rn_override','nH','hmax_fac','hmin_frac','hmin_abs','tol_root', ...
        'emt_static','hmf_integral_mode','hmf_endpoint_allow_unconverged_pm', ...
        'hmf_endpoint_pm_max_outer','hmf_filtered_include_unconverged'};
    for k = 1:numel(remove)
        if isfield(pmopts,remove{k}), pmopts = rmfield(pmopts,remove{k}); end
    end
    pmopts.J0eff = J0eff;
    pmopts.Jxx0 = Jxx0;
    pmopts.max_outer = getf(opts,'hmf_endpoint_pm_max_outer',maxo);
    pm = invz_solve_point(ion,T,Bx,Jnu_flat,pmopts);
    endpoint_r0 = 1+pm.Sigma0;
    endpoint_usable = isfinite(endpoint_r0) && ...
        (pm.converged || allow_unconverged_pm);
    prof.endpoint_r0 = endpoint_r0;
    prof.endpoint_pm_converged = pm.converged;
    prof.endpoint_pm_iters = pm.outer_iters;
    prof.endpoint_pm_crit = pm.crit;
    if pm.converged
        prof.endpoint_r0_source = "converged_pm_limit";
    elseif endpoint_usable
        prof.endpoint_r0_source = "unconverged_pm_last_iterate_visual";
    else
        prof.endpoint_r0_source = "pm_endpoint_unusable";
    end
end

ratio = hfrac^(1/(nH-1));
hgrid = hmax * ratio.^((nH-1):-1:0);                 % geometric, clustered at 0 (P1-4)
prof.hmin_initial = hgrid(1);

[rv, mv, S0v, K0v, Dv, Gbv, Gsv, cnv, Sigma, K0s, stv] = run_sweep(hgrid, Sigma, K0s);

if endpoint_visual && endpoint_usable
    h0 = 0.5*hgrid.*(endpoint_r0+rv);
elseif filtered_visual
    if filtered_include_unconverged
        integral_eligible = isfinite(mv);
    else
        integral_eligible = cnv & isfinite(mv);
    end
    if endpoint_usable
        [h0, filtered_meta] = invz_filtered_profile_integral( ...
            hgrid,rv,integral_eligible,endpoint_r0);
    else
        [h0, filtered_meta] = invz_filtered_profile_integral( ...
            hgrid,rv,integral_eligible,0);
        h0(:) = NaN;
        filtered_meta.anchor = "pm_endpoint_unusable";
    end
    prof.filtered_integral = filtered_meta;
    prof.filtered_unconverged_used_count = ...
        nnz(filtered_meta.used_mask & ~cnv);
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
        [re(k), me(k), S0e(k), K0e(k), De(k), Gbe(k), Gse(k), ce(k), Sigma, K0s, ste(k)] = ...
            eval_node(hext(k), Sigma, K0s);
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
    [rv, mv, S0v, K0v, Dv, Gbv, Gsv, cnv, Sigma, K0s, stv] = run_sweep(hgrid, Sigma, K0s);
    h0 = cumtrapz([0 hgrid], [r0n rv]);  h0 = h0(2:end);
    F  = h0 - J0eff*mv;
    prof.redensified = true;
end
prof.n_extend = n_extend;

prof.hgrid = hgrid;  prof.r = rv;  prof.h0 = h0;  prof.m = mv;
prof.Sigma0 = S0v;   prof.K0 = K0v;  prof.D_uni = Dv;  prof.node_conv = cnv;  prof.F = F;
prof.G0bare = Gbv;   prof.Gstat = Gsv;  prof.static_status = stv;
prof.converged_node_count = nnz(cnv);

if endpoint_visual && ~endpoint_usable
    prof.status = 'endpoint_lower_failed';
    return;
elseif filtered_visual && ~endpoint_usable
    prof.status = 'filtered_anchor_failed';
    return;
elseif filtered_visual && prof.filtered_integral.node_count < 2
    prof.status = 'filtered_too_few_nodes';
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
if endpoint_visual
    valid = cnv & isfinite(rv) & isfinite(mv) & isfinite(F);
    pair = valid(1:end-1) & valid(2:end) & ...
        F(1:end-1) < 0 & F(2:end) >= 0;
    idx = find(pair,1,'last');
    prof.status = 'ok_endpoint_visual';
    if isempty(idx)
        prof.status = 'endpoint_no_bracket';
        return;
    end
    ia = idx;
    ib = idx+1;
elseif filtered_visual
    valid = prof.filtered_integral.used_mask & isfinite(mv) & isfinite(F);
    retained = find(valid);
    pair = F(retained(1:end-1)) < 0 & F(retained(2:end)) >= 0;
    j = find(pair,1,'last');
    prof.status = 'ok_filtered_profile_visual';
    if isempty(j)
        prof.status = 'filtered_no_bracket';
        return;
    end
    ia = retained(j);
    ib = retained(j+1);
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
a = hgrid(ia);  b = hgrid(ib);  Fa = F(ia);  h0a = h0(ia);  ra = rv(ia);
for it = 1:12
    c = 0.5*(a + b);
    [rc, mc, ~, ~, ~, ~, ~, okc, Sigma, K0s, stc] = eval_node(c, Sigma, K0s);
    if ~okc
        prof.refinement_failure_status = stc;
        if filtered_visual
            prof.status = 'filtered_refinement_failed';
        else
            prof.status = 'node_failed';
        end
        hmf_star = NaN;                              % TERMINATES the solve -- never a root
        return;                                       % from a partial bracket
    end
    if endpoint_visual
        h0c = 0.5*c*(endpoint_r0+rc);
    else
        h0c = h0a + 0.5*(ra + rc)*(c - a);
    end
    Fc  = h0c - J0eff*mc;
    if sign(Fc) == sign(Fa), a = c; Fa = Fc; h0a = h0c; ra = rc; else, b = c; end
    if (b - a) < trt*b, break; end
end
if (b - a) >= trt*b                                   % round-5 P1-A: tol_root not reached --
    prof.status = 'unresolved';  hmf_star = NaN;      % a distinct refinement failure
    warning('invz:hmfUnresolved', 'root bracket not refined to tol_root: (b-a)/b = %.3g', (b-a)/b);
    return;
end
hmf_star = 0.5*(a + b);
[r_s, m_s, ~, ~, D_s, ~, Gs_s, ok_s, ~, ~, st_s] = eval_node(hmf_star, Sigma, K0s);
prof.static_status_star = st_s;
if ~ok_s
    if filtered_visual
        prof.status = 'filtered_final_failed';
    else
        prof.status = 'node_failed';
    end
    hmf_star = NaN;
    return;
end
prof.m_star = m_s;  prof.D_uni_star = D_s;  prof.r_star = r_s;  prof.Gstat_star = Gs_s;

    function [rv, mv, S0v, K0v, Dv, Gbv, Gsv, cnv, Sigma, K0s, stv] = run_sweep(hgrid, Sigma, K0s)
    % Behavior-neutral factoring of the per-node nH sweep (execution amendment 3):
    % evaluate eval_node at every point of hgrid in ascending order, warm-starting
    % Sigma/K0s across calls. Shared by the initial profile sweep and the re-densify
    % pass; the extension's 3-node prepends stay inline (unchanged).
    n = numel(hgrid);
    [rv, mv, S0v, K0v, Dv, Gbv, Gsv] = deal(nan(1, n));  cnv = false(1, n);
    stv = strings(1,n);
    for is = 1:n
        [rv(is), mv(is), S0v(is), K0v(is), Dv(is), Gbv(is), Gsv(is), cnv(is), Sigma, K0s, stv(is)] = ...
            eval_node(hgrid(is), Sigma, K0s);
    end
    end

    function [rk, mk, S0k, K0k, Dk, Gbk, Gsk, ok, Sigma, K0s, node_status] = eval_node(hp, Sigma, K0s)
    % One fixed-field node: si (hz_fixed, NO order), tl, c0/G0, then the ordered
    % Sigma<->EMT loop WITH the static-sector closure each pass (Interfaces bullet).
    sio = sibase;  sio.hz_fixed = hp;
    si = invz_single_ion(ion, T, Bvec, sio);
    if abs(si.hz - hp) > 1e-12
        error('invz:hzFixed', 'hz_fixed not held: si.hz = %.6g vs %.6g', si.hz, hp);
    end
    tl = invz_twolevel_ordered(ion, T, Bx, hp, struct('Jxx0', Jxx0));
    mk = si.Jexp(3);
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
    K = zeros(size(wn));  lam = [0; 0; 0];  ok = false;
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
