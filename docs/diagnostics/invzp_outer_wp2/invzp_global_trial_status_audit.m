function result = invzp_global_trial_status_audit()
%INVZP_GLOBAL_TRIAL_STATUS_AUDIT Classify undefined global-search trials.
% Adversarial check for roots hidden by multiple/unresolved dynamic scalar
% equations rather than ordinary domain failure.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
A = load(fullfile(here,'wp2_reduced_global_audit.mat'));
C = load(fullfile(repo,'docs','diagnostics','invzp_representation_wp3', ...
    'wp3_closed_global_audit.mat'));
T = A.result.T;
Bx = A.result.Bx;

labels = ["endpoint_h0";"node22";"closed_h0"];
contexts = {make_hybrid_context(Bx,0,T,F), ...
    make_hybrid_context(Bx,A.result.cases(2).h,T,F), ...
    make_closed_context(Bx,0,T,F)};
search_records = {A.result.records{1},A.result.records{2},C.result};
details = cell(numel(labels),1);
summary = table('Size',[numel(labels),7], ...
    'VariableTypes',{'string','double','double','double','double','double','double'}, ...
    'VariableNames',{'case_name','trial_count','undefined_count', ...
    'multiple_root_count','unresolved_count','no_domain_count','no_root_count'});

for ic = 1:numel(labels)
    R = search_records{ic};
    trials = R.global_trials;
    idx = find(trials.Fval > 999);
    status = strings(numel(idx),1);
    max_root_count = nan(numel(idx),1);
    unresolved_frequencies = nan(numel(idx),1);
    for j = 1:numel(idx)
        z = trials.X(idx(j),:).';
        y = decode(z,R,F.info.Jcc0);
        q = invz_ordered_reduced_residual(y,contexts{ic}, ...
            struct('Jsup',F.info.Jcc0, ...
            'dynamic',struct('resid_tol',1e-12)));
        status(j) = string(q.status);
        if ~isempty(q.dynamic)
            max_root_count(j) = max(q.dynamic.root_count(2:end));
            unresolved_frequencies(j) = ...
                sum(q.dynamic.fallback_unresolved(2:end));
        end
    end
    details{ic} = table(idx,status,max_root_count,unresolved_frequencies);
    summary.case_name(ic) = labels(ic);
    summary.trial_count(ic) = numel(trials.Fval);
    summary.undefined_count(ic) = numel(idx);
    summary.multiple_root_count(ic) = sum(status == "dynamic_multiple_roots");
    summary.unresolved_count(ic) = sum(status == "dynamic_search_unresolved");
    summary.no_domain_count(ic) = sum(status == "dynamic_no_physical_domain");
    summary.no_root_count(ic) = sum(status == "dynamic_no_root");
end

result = struct('summary',summary,'details',{details}, ...
    'sources',struct('hybrid_global', ...
    fullfile(here,'wp2_reduced_global_audit.mat'),'closed_global', ...
    fullfile(repo,'docs','diagnostics','invzp_representation_wp3', ...
    'wp3_closed_global_audit.mat')));
save(fullfile(here,'wp2_global_trial_status_audit.mat'),'result','-v7');
disp(summary);
end

function y = decode(z,R,Jsup)
mid = (R.lambda_lo+R.lambda_hi)/2;
half = (R.lambda_hi-R.lambda_lo)/2;
y = [mid+half.*z(1:3);-z(4)/Jsup];
end

function ctx = make_hybrid_context(Bx,h,T,F)
ion = invz_ion();
[wn,wts,beta] = invz_matsubara(T,40);
si = invz_single_ion(ion,T,[Bx 0 0], ...
    struct('hyp',true,'Jxx0',F.info.Jaa0,'hz_fixed',h));
tl = invz_twolevel_ordered(ion,T,Bx,h,struct('Jxx0',F.info.Jaa0));
c0 = invz_chi0z(si,T,1i*wn,struct('elastic',true));
G0 = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si,T,1i*wn(1),struct('elastic',false));
G0i = -real(c0i(3,3,1));
X = real(c0(:,:,1));
feedback = X(3,1)*(F.info.Jaa0/(1-F.info.Jaa0*X(1,1)))*X(1,3);
G0e = -(X(3,3)+feedback)-G0i;
g = real(invz_g(tl,1i*wn));
ctx = struct('G0',G0,'g',g,'tl',tl,'wts',wts,'beta',beta, ...
    'Jnu_flat',F.J,'J0eff',F.info.Jcc0, ...
    'G0inel0',G0i,'G0el0',G0e);
end

function ctx = make_closed_context(Bx,h,T,F)
ion = invz_ion();
[wn,wts,beta] = invz_matsubara(T,40);
tl = invz_twolevel_ordered(ion,T,Bx,h,struct('Jxx0',F.info.Jaa0));
g = real(invz_g(tl,1i*wn));
G0 = -tl.M2*g;
G0(1) = G0(1)-tl.m^2*tl.h0;
ctx = struct('G0',G0,'g',g,'tl',tl,'wts',wts,'beta',beta, ...
    'Jnu_flat',F.J,'J0eff',F.info.Jcc0, ...
    'G0inel0',-tl.M2*tl.g0,'G0el0',-tl.m^2*tl.h0);
end
