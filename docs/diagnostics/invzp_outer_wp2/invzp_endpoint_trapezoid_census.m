function result = invzp_endpoint_trapezoid_census()
%INVZP_ENDPOINT_TRAPEZOID_CENSUS Two-endpoint approximation to Jensen Eq. 45.
% For every candidate upper endpoint h on the production 33-node ordered
% profile, replace the full integral by
%   h0_two(h) = h*(r0 + r(h))/2.
% The strict ordered h=0 endpoint and the independently converged PM-limit
% endpoint r0=1+Sigma0 are reported separately. A finite but unconverged PM
% last iterate is retained only to show what the uncontrolled substitution
% would produce; it is never classified as certified.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

fixture_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat');
F = load(fixture_path);
ion = invz_ion();
T = F.provenance.T;
fields = (0.5:0.5:9).';
nfield = numel(fields);

static_opts = struct('Jsup',F.info.Jcc0,'warn',false, ...
    'scan_points',4097,'endpoint_margin',1e-10,'resid_tol',1e-10);
ordered_opts = struct('mix_outer',0.35,'max_outer',200, ...
    'tol_outer',1e-8,'Ecut',40,'J0eff',F.info.Jcc0, ...
    'Jxx0',F.info.Jaa0,'nH',33,'emt_static',static_opts);
pm_opts = ordered_opts;
pm_opts = rmfield(pm_opts,{'nH','emt_static'});
pm_opts.max_outer = 1000;

pm_converged = false(nfield,1);
pm_outer_iters = nan(nfield,1);
pm_sigma0 = nan(nfield,1);
pm_crit = nan(nfield,1);
r0_pm_limit = nan(nfield,1);
pm_Gtilde0 = nan(nfield,1);
pm_uniform_mass = nan(nfield,1);

% Follow the PM fixed point from the stable high-field side. A rejected trial
% never replaces the last certified continuation seed.
pm_seed = [];
for k = nfield:-1:1
    opts = pm_opts;
    if ~isempty(pm_seed), opts.Sigma_seed = pm_seed; end
    pm = invz_solve_point(ion,T,fields(k),F.J,opts);
    pm_converged(k) = pm.converged;
    pm_outer_iters(k) = pm.outer_iters;
    pm_sigma0(k) = pm.Sigma0;
    pm_crit(k) = pm.crit;
    if pm.converged
        pm_seed = pm.Sigma;
        % Exact m=0 identity of invz_gstat_ordered:
        % Gtilde0=G0/(1+Sigma0), hence r0=1+Sigma0.
        r0_pm_limit(k) = 1+pm.Sigma0;
        pm_Gtilde0(k) = -pm.chi0cc0/r0_pm_limit(k);
        pm_uniform_mass(k) = 1+F.info.Jcc0*pm_Gtilde0(k);
    end
end

bare_orders = false(nfield,1);
strict_r0_converged = false(nfield,1);
strict_r0_status = strings(nfield,1);
n_upper_valid = zeros(nfield,1);
upper_valid_high_block = false(nfield,1);
pm_limit_bracket = false(nfield,1);
pm_limit_h_linear = nan(nfield,1);
exploratory_bracket = false(nfield,1);
exploratory_h_linear = nan(nfield,1);
classification = strings(nfield,1);
details = cell(nfield,1);

for k = 1:nfield
    Bx = fields(k);
    pt = invz_solve_point_ordered(ion,T,Bx,F.J,ordered_opts);
    p = pt.hmf_prof;
    bare_orders(k) = ~isempty(p.hgrid);
    strict_r0_converged(k) = p.predictor_converged;
    strict_r0_status(k) = string(p.predictor_static_status);

    if ~bare_orders(k)
        if pm_converged(k) && pm_crit(k) > 0
            classification(k) = "accepted_pm_side_no_order_integral";
        else
            classification(k) = "no_order_profile_and_no_stable_pm_endpoint";
        end
        details{k} = struct('profile_status',p.status);
        fprintf('endpoint census B=%.1f T: PM region, no ordered integral\n',Bx);
        continue;
    end

    Gtilde = p.G0bare./p.r;
    valid = p.node_conv & isfinite(p.hgrid) & isfinite(p.r) & ...
        isfinite(p.m) & isfinite(Gtilde) & Gtilde < 0 & ...
        isfinite(p.D_uni) & p.D_uni > 0;
    n_upper_valid(k) = nnz(valid);
    upper_valid_high_block(k) = top_true_count(valid) == nnz(valid);

    Fpm = nan(size(p.hgrid));
    if pm_converged(k)
        Fpm(valid) = 0.5*p.hgrid(valid).* ...
            (r0_pm_limit(k)+p.r(valid))-F.info.Jcc0*p.m(valid);
        [pm_limit_bracket(k),pm_limit_h_linear(k),pm_idx] = ...
            last_linear_bracket(p.hgrid,Fpm,valid);
    else
        pm_idx = NaN;
    end

    % Deliberately uncontrolled comparator: use the finite PM last iterate
    % even when it did not satisfy the PM fixed-point tolerance.
    Flast = nan(size(p.hgrid));
    r0_last = 1+pm_sigma0(k);
    if isfinite(r0_last)
        Flast(valid) = 0.5*p.hgrid(valid).* ...
            (r0_last+p.r(valid))-F.info.Jcc0*p.m(valid);
        [exploratory_bracket(k),exploratory_h_linear(k),expl_idx] = ...
            last_linear_bracket(p.hgrid,Flast,valid);
    else
        expl_idx = NaN;
    end

    classification(k) = classify_result(pm_converged(k),valid,Fpm,Flast, ...
        pm_limit_bracket(k),exploratory_bracket(k), ...
        strict_r0_converged(k));
    details{k} = struct('hgrid',p.hgrid,'r_upper',p.r,'moment',p.m, ...
        'G0bare_upper',p.G0bare,'Gtilde_upper',Gtilde, ...
        'D_uni_upper',p.D_uni,'Sigma0_upper',p.Sigma0,'K0_upper',p.K0, ...
        'upper_valid',valid,'static_status',p.static_status, ...
        'F_pm_limit',Fpm,'F_last_iterate',Flast, ...
        'pm_limit_bracket_index',pm_idx, ...
        'exploratory_bracket_index',expl_idx, ...
        'profile_status',p.status);
    fprintf(['endpoint census B=%.1f T: strict0=%s PM0=%d upper=%d/%d ' ...
        'pm_limit_h=%g exploratory_h=%g class=%s\n'],Bx,strict_r0_status(k), ...
        pm_converged(k),n_upper_valid(k),numel(valid), ...
        pm_limit_h_linear(k),exploratory_h_linear(k),classification(k));
end

tab = table(fields,pm_converged,pm_outer_iters,pm_sigma0,pm_crit, ...
    r0_pm_limit,pm_Gtilde0,pm_uniform_mass,bare_orders, ...
    strict_r0_converged,strict_r0_status,n_upper_valid, ...
    upper_valid_high_block,pm_limit_bracket,pm_limit_h_linear, ...
    exploratory_bracket,exploratory_h_linear, ...
    classification,'VariableNames',{'Bx','pm_converged','pm_outer_iters', ...
    'pm_sigma0','pm_crit','r0_pm_limit','pm_Gtilde0', ...
    'pm_uniform_mass','bare_orders','strict_r0_converged', ...
    'strict_r0_status','n_upper_valid','upper_valid_high_block', ...
    'pm_limit_bracket','pm_limit_h_linear', ...
    'exploratory_bracket','exploratory_h_linear','classification'});
result = struct('table',tab,'details',{details}, ...
    'ordered_opts',ordered_opts,'pm_opts',pm_opts, ...
    'formula','h0_two(h)=0.5*h*(r0+r(h)); F_two=h0_two-J0eff*m(h)', ...
    'provenance',struct('fixture',fixture_path, ...
                        'coupling_opts',F.provenance.coupling_opts), ...
    'interpretation',['A pm_limit bracket requires a converged, finite PM ' ...
        'fixed-point lower limit and two physically gated ordered upper-endpoint ' ...
        'evaluations. The PM limit can be unstable (negative uniform mass); ' ...
        'the bracket does not certify the strict ordered lower endpoint, the ' ...
        'skipped interior, or thermodynamic selection. exploratory_h_linear ' ...
        'uses an unconverged PM last iterate when needed.']);
save(fullfile(here,'wp2_endpoint_trapezoid_census.mat'),'result','-v7');
disp(tab);
end

function n = top_true_count(v)
v = logical(v(:));
idx = find(~v,1,'last');
if isempty(idx), n = numel(v);
else,            n = numel(v)-idx;
end
end

function [found,hlin,idx] = last_linear_bracket(h,F,valid)
pair = valid(1:end-1) & valid(2:end) & ...
       F(1:end-1) <= 0 & F(2:end) >= 0;
idx = find(pair,1,'last');
found = ~isempty(idx);
hlin = NaN;
if ~found
    idx = NaN;
    return;
end
dF = F(idx+1)-F(idx);
if ~(isfinite(dF) && dF > 0)
    found = false;
    idx = NaN;
    return;
end
hlin = h(idx)-F(idx)*(h(idx+1)-h(idx))/dF;
end

function label = classify_result(pm_ok,valid,Fpm,Flast,pm_br,expl_br,strict_ok)
if ~strict_ok
    prefix = "strict_ordered_lower_failed:";
else
    prefix = "strict_ordered_lower_ok:";
end
if ~any(valid)
    label = prefix+"no_valid_upper_endpoint";
elseif pm_br
    label = prefix+"converged_pm_limit_coarse_bracket";
elseif ~pm_ok && expl_br
    label = prefix+"uncertified_pm_last_iterate_bracket";
elseif pm_ok && all(Fpm(valid) < 0)
    label = prefix+"coarse_root_above_profile";
elseif pm_ok && all(Fpm(valid) >= 0)
    label = prefix+"coarse_root_below_valid_block_or_none";
elseif ~pm_ok && all(Flast(valid) < 0)
    label = prefix+"uncertified_root_above_profile";
elseif ~pm_ok && all(Flast(valid) >= 0)
    label = prefix+"uncertified_root_below_valid_block_or_none";
elseif ~pm_ok
    label = prefix+"pm_lower_endpoint_unconverged";
else
    label = prefix+"no_contiguous_coarse_bracket";
end
end
