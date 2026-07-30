function result = invzp_saturation_tail_census()
%INVZP_SATURATION_TAIL_CENSUS Test prerequisites for an upper-anchored H0.
% A certified high-h outer-map component is continued downward from a nearly
% saturated state. The diagnostic measures r-1 and connected fluctuations;
% it neither assumes delta_h(infinity)=0 nor solves equation (45).
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

fixture_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat');
control_path = fullfile(repo,'docs','diagnostics','invzp_outer_wp2', ...
    'wp2_hmf_node_transaction_census.mat');
F = load(fixture_path);
C = load(control_path);
ion = invz_ion();
T = F.provenance.T;
opts = C.result.base_opts;
[wn,wts,beta] = invz_matsubara(T,opts.Ecut);
B = [1;4];
factor = logspace(0,log10(128),25).';
nB = numel(B);
details = cell(nB,1);
tail_exponent = nan(nB,1);
tail_amplitude = nan(nB,1);
tail_integral_beyond_max = nan(nB,1);
tail_fit_count = zeros(nB,1);
tail_sign_constant = false(nB,1);
variance_exponent = nan(nB,1);
closed_tail_exponent = nan(nB,1);
closed_variance_exponent = nan(nB,1);
closed_accepted_count = zeros(nB,1);
validity_probe_h = 50*ones(nB,1);
validity_probe_Jz = nan(nB,1);
validity_probe_variance = nan(nB,1);
validity_probe_energy_span = nan(nB,1);
hmax_profile = nan(nB,1);
accepted_count = zeros(nB,1);
lowest_accepted_h = nan(nB,1);

for ib = 1:nB
    sib = invz_single_ion(ion,T,[B(ib) 0 0], ...
        struct('hyp',true,'order',true,'J0z',opts.J0eff,'Jxx0',opts.Jxx0));
    hmax_profile(ib) = 1.25*abs(sib.hz);
    h = hmax_profile(ib)*factor;
    n = numel(h);
    accepted = false(n,1);
    probe_status = strings(n,1);
    map_status = strings(n,1);
    outer_iterations = nan(n,1);
    Sigma0 = nan(n,1);
    K0 = nan(n,1);
    r = nan(n,1);
    G0bare = nan(n,1);
    Gtil0 = nan(n,1);
    D_uni = nan(n,1);
    supremum_mass = nan(n,1);
    dynamic_min_abs = nan(n,1);
    Jz_variance = nan(n,1);
    M2 = nan(n,1);
    used_fresh_fallback = false(n,1);
    Sigma_seed = zeros(size(wn));

    for j = n:-1:1
        [o,Sigma_candidate] = hybrid_node(ion,T,B(ib),h(j),F.J,opts, ...
            wn,wts,beta,Sigma_seed);
        accepted(j) = o.converged;
        probe_status(j) = o.probe_status;
        map_status(j) = o.map_status;
        outer_iterations(j) = o.iterations;
        Sigma0(j) = o.Sigma0;
        K0(j) = o.K0;
        r(j) = o.r;
        G0bare(j) = o.G0bare;
        Gtil0(j) = o.Gtil0;
        D_uni(j) = o.D_uni;
        supremum_mass(j) = o.supremum_mass;
        dynamic_min_abs(j) = o.dynamic_min_abs;
        Jz_variance(j) = o.Jz_variance;
        M2(j) = o.M2;
        used_fresh_fallback(j) = o.used_fresh_fallback;
        if o.converged, Sigma_seed = Sigma_candidate; end
    end
    closed_accepted = false(n,1);
    closed_r_minus_1 = nan(n,1);
    closed_variance = nan(n,1);
    closed_seed = zeros(size(wn));
    for j = n:-1:1
        [oc,closed_candidate] = closed_node(ion,T,B(ib),h(j),F.J,opts, ...
            wn,wts,beta,closed_seed);
        closed_accepted(j) = oc.converged;
        closed_r_minus_1(j) = oc.r-1;
        closed_variance(j) = oc.variance;
        if oc.converged, closed_seed = closed_candidate; end
    end
    r_minus_1 = r-1;
    [tail_exponent(ib),tail_amplitude(ib),tail_integral_beyond_max(ib), ...
        tail_fit_count(ib),tail_sign_constant(ib)] = ...
        power_tail(h,r_minus_1,accepted,factor);
    [variance_exponent(ib),~,~,~,~] = ...
        power_tail(h,Jz_variance,accepted,factor);
    [closed_tail_exponent(ib),~,~,~,~] = ...
        power_tail(h,closed_r_minus_1,closed_accepted,factor);
    [closed_variance_exponent(ib),~,~,~,~] = ...
        power_tail(h,closed_variance,closed_accepted,factor);
    closed_accepted_count(ib) = nnz(closed_accepted);
    si_validity = invz_single_ion(ion,T,[B(ib) 0 0], ...
        struct('hyp',true,'Jxx0',opts.Jxx0,'hz_fixed',validity_probe_h(ib)));
    validity_probe_Jz(ib) = si_validity.Jexp(3);
    validity_probe_variance(ib) = si_validity.JzJz_fluct;
    validity_probe_energy_span(ib) = max(si_validity.E)-min(si_validity.E);
    accepted_count(ib) = nnz(accepted);
    highest_accepted_h = NaN;
    if any(accepted)
        lowest_accepted_h(ib) = min(h(accepted));
        highest_accepted_h = max(h(accepted));
    end
    details{ib} = table(factor,h,accepted,probe_status,map_status, ...
        outer_iterations,Sigma0,K0,r,r_minus_1,G0bare,Gtil0,D_uni, ...
        supremum_mass,dynamic_min_abs,Jz_variance,M2,used_fresh_fallback, ...
        closed_accepted,closed_r_minus_1,closed_variance);
    fprintf(['B=%.1f: %d/%d accepted, h %.5g..%.5g, p(r-1)=%.4g ' ...
        'from %d points, p(var)=%.4g; closed p=%.4g p(var)=%.4g\n'], ...
        B(ib),accepted_count(ib),n, ...
        lowest_accepted_h(ib),highest_accepted_h,tail_exponent(ib), ...
        tail_fit_count(ib),variance_exponent(ib),closed_tail_exponent(ib), ...
        closed_variance_exponent(ib));
end

summary = table(B,hmax_profile,accepted_count,lowest_accepted_h, ...
    tail_exponent,tail_amplitude,tail_integral_beyond_max,tail_fit_count, ...
    tail_sign_constant,variance_exponent,closed_accepted_count, ...
    closed_tail_exponent,closed_variance_exponent,validity_probe_h, ...
    validity_probe_Jz,validity_probe_variance,validity_probe_energy_span);
result = struct('summary',summary,'details',{details},'base_opts',opts, ...
    'provenance',struct('fixture',fixture_path,'control',control_path, ...
    'commit',git_commit(repo)), ...
    'note',['Tail fits use accepted points with factor>=8 and |quantity| ' ...
    'above a roundoff-scaled floor. A fitted exponent >1 supports ' ...
    'integrability over the sampled asymptotic window but does not by itself ' ...
    'prove delta_h(infinity)=0 or validate extrapolation through level changes.']);
save(fullfile(here,'wp5_saturation_tail_census.mat'),'result','-v7');
disp(summary);
end

function [out,Sigma_candidate] = hybrid_node( ...
        ion,T,B,h,J,opts,wn,wts,beta,Sigma_seed)
si = invz_single_ion(ion,T,[B 0 0], ...
    struct('hyp',true,'Jxx0',opts.Jxx0,'hz_fixed',h));
tl = invz_twolevel_ordered(ion,T,B,h,struct('Jxx0',opts.Jxx0));
c0 = invz_chi0z(si,T,1i*wn,struct('elastic',true));
G0 = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si,T,0,struct('elastic',false));
G0inel0 = -real(c0i(3,3,1));
X = real(c0(:,:,1));
feedback = X(3,1)*(opts.Jxx0/(1-opts.Jxx0*X(1,1)))*X(1,3);
G0bare0 = -(X(3,3)+feedback);
G0el0 = G0bare0-G0inel0;
g = real(invz_g(tl,1i*wn));
ctx = struct('G0',G0,'g',g,'tl',tl,'wts',wts,'beta',beta, ...
    'Jnu_flat',J,'J0eff',opts.J0eff,'G0inel0',G0inel0,'G0el0',G0el0);
mapopts = struct('emt_static',opts.emt_static,'dynamic_diagnostics',true);
mapfun = @(s) invz_ordered_outer_map(s,ctx,mapopts);
probe_opts = struct('mix',opts.mix_outer,'tol',opts.tol_outer, ...
    'maxit',opts.max_outer);
probe = invz_outer_picard_diagnostic(mapfun,Sigma_seed,probe_opts);
used_fresh = false;
if ~probe.converged && any(Sigma_seed ~= 0)
    fresh = invz_outer_picard_diagnostic(mapfun,zeros(size(wn)),probe_opts);
    used_fresh = true;
    if fresh.converged, probe = fresh; end
end
out = struct('probe_status',string(probe.status),'map_status',"", ...
    'converged',probe.converged,'iterations',probe.iterations, ...
    'Sigma0',NaN,'K0',NaN,'r',NaN,'G0bare',G0bare0,'Gtil0',NaN, ...
    'D_uni',NaN,'supremum_mass',NaN,'dynamic_min_abs',NaN, ...
    'Jz_variance',si.JzJz_fluct,'M2',tl.M2, ...
    'used_fresh_fallback',used_fresh);
Sigma_candidate = Sigma_seed;
if isstruct(probe.last_map) && isfield(probe.last_map,'status')
    out.map_status = string(probe.last_map.status);
end
if probe.converged
    final = probe.last_map;
    Sigma_candidate = probe.Sigma;
    out.Sigma0 = probe.Sigma(1);
    out.K0 = final.K(1);
    out.r = final.static.r;
    out.Gtil0 = final.static.Gtil0;
    out.D_uni = final.static.D_uni;
    out.supremum_mass = ...
        final.static.root_table.supremum_mass(final.static.selected_index);
    out.dynamic_min_abs = final.dynamic_min_abs;
end
end

function [out,Sigma_candidate] = closed_node( ...
        ion,T,B,h,J,opts,wn,wts,beta,Sigma_seed)
tl = invz_twolevel_ordered(ion,T,B,h,struct('Jxx0',opts.Jxx0));
g = real(invz_g(tl,1i*wn));
G0 = -tl.M2*g;
G0(1) = G0(1)-tl.m^2*tl.h0;
ctx = struct('G0',G0,'g',g,'tl',tl,'wts',wts,'beta',beta, ...
    'Jnu_flat',J,'J0eff',opts.J0eff, ...
    'G0inel0',-tl.M2*tl.g0,'G0el0',-tl.m^2*tl.h0);
mapopts = struct('emt_static',opts.emt_static);
mapfun = @(s) invz_ordered_outer_map(s,ctx,mapopts);
probe_opts = struct('mix',opts.mix_outer,'tol',opts.tol_outer, ...
    'maxit',opts.max_outer);
probe = invz_outer_picard_diagnostic(mapfun,Sigma_seed,probe_opts);
if ~probe.converged && any(Sigma_seed ~= 0)
    fresh = invz_outer_picard_diagnostic(mapfun,zeros(size(wn)),probe_opts);
    if fresh.converged, probe = fresh; end
end
out = struct('converged',probe.converged,'r',NaN, ...
    'variance',tl.M2+tl.m^2*(1-tl.n01^2));
Sigma_candidate = Sigma_seed;
if probe.converged
    Sigma_candidate = probe.Sigma;
    out.r = probe.last_map.static.r;
end
end

function [p,A,tail,nfit,sign_constant] = power_tail(h,y,accepted,factor)
p = NaN; A = NaN; tail = NaN; nfit = 0; sign_constant = false;
values = abs(y(accepted & isfinite(y)));
if isempty(values), return; end
scale = max(values);
floor_y = max(1e-14,256*eps(max(1,scale)));
use = accepted & factor >= 8 & isfinite(y) & abs(y) > floor_y;
idx = find(use);
nfit = numel(idx);
if nfit < 4, return; end
sy = sign(y(idx));
sign_constant = all(sy == sy(1));
if ~sign_constant, return; end
c = polyfit(log(h(idx)),log(abs(y(idx))),1);
p = -c(1);
A = sy(1)*exp(c(2));
if p > 1
    tail = A/(p-1)*max(h(idx))^(1-p);
end
end

function commit = git_commit(repo)
[status,text] = system(sprintf('git -C "%s" rev-parse HEAD',repo));
if status == 0, commit = strtrim(text); else, commit = 'unknown'; end
end
