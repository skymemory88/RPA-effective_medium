function result = invzp_closed_twolevel_landmarks()
%INVZP_CLOSED_TWOLEVEL_LANDMARKS Compare hybrid nodes with a closed model.
% The bare response, ordered static weights, and Jensen vertex all come from
% the same electronic two-level system. This tests representation consistency
% at four landmarks; it does not rebuild equation (45) or prove root absence
% when the bounded Picard probe fails.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

source_path = fullfile(repo,'docs','diagnostics','invzp_outer_wp2', ...
    'wp2_low_field_m2_census.mat');
fixture_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat');
S = load(source_path);
F = load(fixture_path);
ion = invz_ion();
T = 0.1;
opts = S.result.base_opts;
[wn,wts,beta] = invz_matsubara(T,opts.Ecut);
Bkeep = [0.5 1.0 1.8 2.2];
labels = ["high_endpoint";"first_accepted";"last_failed";"predictor"];
nrow = numel(Bkeep)*numel(labels);

B = nan(nrow,1);
sample = strings(nrow,1);
h = nan(nrow,1);
hybrid_accepted = false(nrow,1);
hybrid_status = strings(nrow,1);
zero_map_status = strings(nrow,1);
probe_status = strings(nrow,1);
probe_last_map_status = strings(nrow,1);
closed_converged = false(nrow,1);
outer_iterations = nan(nrow,1);
Sigma0 = nan(nrow,1);
static_D_uni = nan(nrow,1);
static_supremum_mass = nan(nrow,1);
dynamic_min_abs = nan(nrow,1);
used_fresh_fallback = false(nrow,1);
details = cell(numel(Bkeep),1);

row = 0;
for ib = 1:numel(Bkeep)
    is = find(abs(S.result.summary.B-Bkeep(ib)) < 1e-12,1);
    t = S.result.details{is};
    ilastfail = find(~t.accepted,1,'last');
    ifirstok = find(t.accepted,1);
    h_all = [0;t.h];
    hybrid_accepted_all = [false;t.accepted];
    hybrid_status_all = [S.result.summary.predictor_status(is);t.status];
    n = numel(h_all);
    zero_status_all = strings(n,1);
    probe_status_all = strings(n,1);
    last_status_all = strings(n,1);
    converged_all = false(n,1);
    iterations_all = nan(n,1);
    Sigma0_all = nan(n,1);
    Duni_all = nan(n,1);
    Dsup_all = nan(n,1);
    Ddyn_all = nan(n,1);
    fallback_all = false(n,1);
    Sigma_seed = zeros(size(wn));

    % Use every nested profile node to avoid mistaking a coarse continuation
    % jump for a component endpoint. Rejected states never replace the carrier.
    for j = n:-1:1
        [o,Sigma_candidate] = closed_node(ion,T,Bkeep(ib),h_all(j),F.J, ...
            opts,wn,wts,beta,Sigma_seed);
        zero_status_all(j) = o.zero_map_status;
        probe_status_all(j) = o.probe_status;
        last_status_all(j) = o.last_map_status;
        converged_all(j) = o.converged;
        iterations_all(j) = o.iterations;
        Sigma0_all(j) = o.Sigma0;
        Duni_all(j) = o.D_uni;
        Dsup_all(j) = o.supremum_mass;
        Ddyn_all(j) = o.dynamic_min_abs;
        fallback_all(j) = o.used_fresh_fallback;
        if o.converged, Sigma_seed = Sigma_candidate; end
    end
    d = table(h_all,hybrid_accepted_all,hybrid_status_all,zero_status_all, ...
        probe_status_all,last_status_all,converged_all,iterations_all, ...
        Sigma0_all,Duni_all,Dsup_all,Ddyn_all,fallback_all, ...
        'VariableNames',{'h','hybrid_accepted','hybrid_status', ...
        'zero_map_status','probe_status','probe_last_map_status', ...
        'closed_converged','outer_iterations','Sigma0','static_D_uni', ...
        'static_supremum_mass','dynamic_min_abs','used_fresh_fallback'});
    details{ib} = d;
    indices = [n ifirstok+1 ilastfail+1 1];

    for ic = 1:numel(labels)
        row = row+1;
        B(row) = Bkeep(ib);
        sample(row) = labels(ic);
        j = indices(ic);
        h(row) = d.h(j);
        hybrid_accepted(row) = d.hybrid_accepted(j);
        hybrid_status(row) = d.hybrid_status(j);
        zero_map_status(row) = d.zero_map_status(j);
        probe_status(row) = d.probe_status(j);
        probe_last_map_status(row) = d.probe_last_map_status(j);
        closed_converged(row) = d.closed_converged(j);
        outer_iterations(row) = d.outer_iterations(j);
        Sigma0(row) = d.Sigma0(j);
        static_D_uni(row) = d.static_D_uni(j);
        static_supremum_mass(row) = d.static_supremum_mass(j);
        dynamic_min_abs(row) = d.dynamic_min_abs(j);
        used_fresh_fallback(row) = d.used_fresh_fallback(j);
        fprintf(['B=%.1f %-14s h=%.7g hybrid=%d/%s closed=%d/%s ' ...
            'map0=%s\n'],B(row),sample(row),h(row),hybrid_accepted(row), ...
            hybrid_status(row),closed_converged(row),probe_status(row), ...
            zero_map_status(row));
    end
end

tab = table(B,sample,h,hybrid_accepted,hybrid_status,zero_map_status, ...
    probe_status,probe_last_map_status,closed_converged,outer_iterations, ...
    Sigma0,static_D_uni,static_supremum_mass,dynamic_min_abs, ...
    used_fresh_fallback);
result = struct('table',tab,'details',{details},'base_opts',opts,'source',source_path, ...
    'provenance',struct('fixture',fixture_path,'commit',git_commit(repo)), ...
    'note',['Closed two-level means G0=-M2*g-m2*h, with the elastic term ' ...
    'only at zero Matsubara frequency; the same tl supplies Sigma and the ' ...
    'ordered static closure. Probe failure is algorithmic/domain evidence, ' ...
    'not a completeness proof. Every 33-node source grid plus h=0 is ' ...
    'traversed high-to-low per field with accepted-state rollback and a ' ...
    'fresh fallback; the displayed table retains four landmarks.']);
save(fullfile(here,'wp3_closed_twolevel_landmarks.mat'),'result','-v7');
disp(tab);
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
mapopts = struct('emt_static',opts.emt_static,'dynamic_diagnostics',true);
mapfun = @(s) invz_ordered_outer_map(s,ctx,mapopts);
zmap = mapfun(zeros(size(wn)));
probe_opts = struct('mix',opts.mix_outer,'tol',opts.tol_outer, ...
    'maxit',opts.max_outer);
probe = invz_outer_picard_diagnostic(mapfun,Sigma_seed,probe_opts);
used_fresh = false;
if ~probe.converged && any(Sigma_seed ~= 0)
    fresh = invz_outer_picard_diagnostic(mapfun,zeros(size(wn)),probe_opts);
    used_fresh = true;
    if fresh.converged, probe = fresh; end
end
out = struct('zero_map_status',string(zmap.status), ...
    'probe_status',string(probe.status),'last_map_status',"", ...
    'converged',probe.converged,'iterations',probe.iterations, ...
    'Sigma0',NaN,'D_uni',NaN,'supremum_mass',NaN, ...
    'dynamic_min_abs',NaN,'used_fresh_fallback',used_fresh);
Sigma_candidate = Sigma_seed;
if isstruct(probe.last_map) && isfield(probe.last_map,'status')
    out.last_map_status = string(probe.last_map.status);
end
if probe.converged
    final = probe.last_map;
    Sigma_candidate = probe.Sigma;
    out.Sigma0 = probe.Sigma(1);
    out.D_uni = final.static.D_uni;
    out.supremum_mass = ...
        final.static.root_table.supremum_mass(final.static.selected_index);
    out.dynamic_min_abs = final.dynamic_min_abs;
end
end

function commit = git_commit(repo)
[status,text] = system(sprintf('git -C "%s" rev-parse HEAD',repo));
if status == 0, commit = strtrim(text); else, commit = 'unknown'; end
end
