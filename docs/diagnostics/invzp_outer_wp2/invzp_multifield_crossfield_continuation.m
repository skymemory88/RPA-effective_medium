function result = invzp_multifield_crossfield_continuation()
%INVZP_MULTIFIELD_CROSSFIELD_CONTINUATION Continue the high-s root in Bx.
% Diagnostic only.  Start from the independently certified high-s root at
% Bx=1.2 T and solve the exact four reduced equations while changing Bx.
% The diagnostic h(Bx) is piecewise-linear through the section-audit probes.
% Adaptive step reduction distinguishes a corrector failure from a robust
% termination of this graph over Bx; it does not cross a fold in Bx.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
S = load(fullfile(here,'wp2_multifield_sheet_section_audit.mat'));
sections = S.result.sections;
Bx_anchor = [sections.Bx];
h_anchor = [sections.h];
T = S.result.T;
Jsup = F.info.Jcc0;
start_index = find(abs(Bx_anchor-1.2) <= 1e-12,1);
assert(~isempty(start_index) && sections(start_index).root_count == 2, ...
    'The continuation requires the certified two-root 1.2 T section.');
[~,high_index] = max([sections(start_index).roots.s]);
y_start = sections(start_index).roots(high_index).y(:);

ls4 = optimoptions('lsqnonlin','Display','off', ...
    'FunctionTolerance',1e-20,'StepTolerance',2e-13, ...
    'OptimalityTolerance',2e-11,'MaxIterations',35, ...
    'MaxFunctionEvaluations',260,'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',3e-6);
config = struct('initial_step',0.05,'maximum_step',0.05, ...
    'minimum_step',0.00078125,'residual_tolerance',1e-8, ...
    'boundary_s',0.9999,'boundary_uniform_mass',1e-4);

start = certify_seed(1.2,h_anchor(start_index),y_start,T,F, ...
    Jsup,ls4,config);
assert(start.accepted,'The source high-s root failed recertification.');
[down,down_status] = advance_path(start,min(Bx_anchor), ...
    Bx_anchor,h_anchor,T,F,Jsup,ls4,config);
[up,up_status] = advance_path(start,max(Bx_anchor), ...
    Bx_anchor,h_anchor,T,F,Jsup,ls4,config);

anchor_summary = repmat(empty_anchor(),numel(Bx_anchor),1);
all_records = [flipud(down(2:end));start;up(2:end)];
for k = 1:numel(Bx_anchor)
    [distance,ix] = min(abs([all_records.Bx]-Bx_anchor(k)));
    if distance <= 5e-10
        high = all_records(ix);
        reached = true;
    else
        high = empty_record();
        reached = false;
    end
    known = sections(k).roots;
    [~,low_index] = min([known.s]);
    low = known(low_index);
    if reached
        r_difference = high.r-low.r;
        distinct = abs(high.s-low.s) > 1e-6;
        signature_pass = distinct && r_difference < 0;
    else
        r_difference = NaN;
        distinct = false;
        signature_pass = false;
    end
    anchor_summary(k) = struct('Bx',Bx_anchor(k),'h',h_anchor(k), ...
        'reached',reached,'high',high,'known_low_s',low.s, ...
        'known_low_r',low.r,'distinct_from_low',distinct, ...
        'r_high_minus_low',r_difference, ...
        'local_signature_pass',signature_pass);
end

result = struct('T',T,'Jsup',Jsup,'Bx_anchor',Bx_anchor, ...
    'h_anchor',h_anchor,'start',start,'down',down,'up',up, ...
    'down_status',down_status,'up_status',up_status, ...
    'anchor_summary',anchor_summary,'config',config, ...
    'source',fullfile(here,'wp2_multifield_sheet_section_audit.mat'), ...
    'interpretation',['Exact cross-field correction avoids the inner ' ...
    'Picard map and fixed-s gaps. Reaching an anchor certifies a distinct ' ...
    'high-s root there. Failure after adaptive step reduction localizes a ' ...
    'termination as a graph over Bx but does not by itself exclude a fold ' ...
    'that requires pseudo-arclength continuation in Bx.']);
save(fullfile(here,'wp2_multifield_crossfield_continuation.mat'), ...
    'result','-v7');

fprintf('cross-field down: %s at Bx %.9g, s %.9g\n', ...
    down_status.reason,down(end).Bx,down(end).s);
fprintf('cross-field up: %s at Bx %.9g, s %.9g\n', ...
    up_status.reason,up(end).Bx,up(end).s);
for k = 1:numel(anchor_summary)
    a = anchor_summary(k);
    fprintf(['cross-field anchor %.3g: reached %d, distinct %d, ' ...
        'r_high-r_low %.9g, signature %d\n'],a.Bx,a.reached, ...
        a.distinct_from_low,a.r_high_minus_low,a.local_signature_pass);
end
end

function [records,status] = advance_path(start,target,Bx_anchor,h_anchor, ...
    T,F,Jsup,ls4,config)
records = start;
direction = sign(target-start.Bx);
step = direction*config.initial_step;
failures = repmat(empty_failure(),0,1);
while direction*(target-records(end).Bx) > 5e-13
    remaining = target-records(end).Bx;
    trial_step = direction*min(abs(step),abs(remaining));
    Bx_trial = records(end).Bx+trial_step;
    h_trial = interp1(Bx_anchor,h_anchor,Bx_trial,'linear');
    trial = correct_seed(Bx_trial,h_trial,records(end).y, ...
        records(end).s,T,F,Jsup,ls4,config);
    if trial.accepted
        records(end+1,1) = trial; %#ok<AGROW>
        step = direction*min(config.maximum_step,1.35*abs(trial_step));
        fprintf(['cross-field Bx %.6g h %.8g s %.9g r %.9g ' ...
            'R %.3g m_uni %.3g\n'],trial.Bx,trial.h,trial.s, ...
            trial.r,trial.residual_norm,trial.uniform_mass);
        if trial.s >= config.boundary_s || ...
                trial.uniform_mass <= config.boundary_uniform_mass
            status = struct('reached_target',false, ...
                'reason',"uniform_boundary_reached", ...
                'failures',failures,'last_Bx',trial.Bx, ...
                'target_Bx',target);
            return;
        end
    else
        failures(end+1,1) = struct('from_Bx',records(end).Bx, ...
            'trial_Bx',Bx_trial,'step',trial_step, ...
            'best_residual_norm',trial.residual_norm, ...
            'best_s',trial.s,'best_status',trial.status); %#ok<AGROW>
        step = trial_step/2;
        if abs(step) < config.minimum_step
            status = struct('reached_target',false, ...
                'reason',"minimum_step_failure",'failures',failures, ...
                'last_Bx',records(end).Bx,'target_Bx',target);
            return;
        end
    end
end
status = struct('reached_target',true,'reason',"target_reached", ...
    'failures',failures,'last_Bx',records(end).Bx,'target_Bx',target);
end

function record = certify_seed(Bx,h,yseed,T,F,Jsup,ls4,config)
ctx = make_context(Bx,h,T,F);
[lo,hi] = lambda_bounds(ctx);
mid = (lo+hi)/2;
half = (hi-lo)/2;
zseed = [(yseed(1:3)-mid)./half;-Jsup*yseed(4)];
record = solve_scaled(Bx,h,zseed,[1e-8 1-1e-8],ctx,mid,half, ...
    Jsup,ls4,config);
end

function record = correct_seed(Bx,h,yseed,s_previous,T,F,Jsup,ls4,config)
ctx = make_context(Bx,h,T,F);
[lo,hi] = lambda_bounds(ctx);
mid = (lo+hi)/2;
half = (hi-lo)/2;
zseed = [(yseed(1:3)-mid)./half;-Jsup*yseed(4)];
windows = [0.04 0.04;0.10 0.10;0.20 0.20];
record = empty_record();
for k = 1:size(windows,1)
    sbounds = [max(1e-8,s_previous-windows(k,1)), ...
        min(1-1e-8,s_previous+windows(k,2))];
    trial = solve_scaled(Bx,h,zseed,sbounds,ctx,mid,half, ...
        Jsup,ls4,config);
    if trial.accepted
        record = trial;
        return;
    end
    if ~isfinite(record.residual_norm) || ...
            trial.residual_norm < record.residual_norm
        record = trial;
    end
end
end

function record = solve_scaled(Bx,h,zseed,sbounds,ctx,mid,half, ...
    Jsup,ls4,config)
ropts = struct('Jsup',Jsup,'dynamic',struct('resid_tol',1e-12));
lb = [-ones(3,1)+1e-9;sbounds(1)];
ub = [ ones(3,1)-1e-9;sbounds(2)];
zseed = min(max(zseed,lb),ub);
[z,~,rscaled,exitflag,output] = lsqnonlin( ...
    @(zz) full_vector(zz,ctx,mid,half,Jsup,ropts), ...
    zseed,lb,ub,ls4);
y = [mid+half.*z(1:3);-z(4)/Jsup];
q = invz_ordered_reduced_residual(y,ctx,ropts);
if q.defined
    outer = invz_ordered_outer_map(q.Sigma,ctx, ...
        struct('emt_static',struct('Jsup',Jsup,'warn',false)));
    outer_norm = outer.residual_norm;
else
    outer_norm = Inf;
end
accepted = exitflag > 0 && q.defined && q.trial_admissible && ...
    q.residual_norm <= config.residual_tolerance && ...
    max(abs(rscaled)) <= config.residual_tolerance && ...
    outer_norm <= config.residual_tolerance;
record = struct('Bx',Bx,'h',h,'y',y,'s',z(4),'r',q.r, ...
    'residual_norm',q.residual_norm, ...
    'scaled_residual_norm',max(abs(rscaled)), ...
    'outer_residual_norm',outer_norm,'accepted',accepted, ...
    'exitflag',exitflag,'status',string(q.status),'output',output, ...
    'supremum_mass',q.supremum_mass,'uniform_mass',q.D_uni, ...
    'static_mesh_mass',q.static_mesh_min_mass, ...
    'dynamic_lattice_mass',q.dynamic.dynamic_min_lattice_mass, ...
    'dynamic_medium_mass',q.dynamic.dynamic_min_medium_mass);
end

function rr = full_vector(z,ctx,mid,half,Jsup,ropts)
q = invz_ordered_reduced_residual( ...
    [mid+half.*z(1:3);-z(4)/Jsup],ctx,ropts);
if q.defined && all(isfinite(q.residual))
    rr = [q.lambda_residual./half;Jsup*q.static_residual];
else
    rr = 1e3*(1+abs(z));
end
end

function [lo,hi] = lambda_bounds(ctx)
nw = numel(ctx.g);
if isvector(ctx.Jnu_flat)
    jlo = repmat(min(ctx.Jnu_flat),nw,1);
    jhi = repmat(max(ctx.Jnu_flat),nw,1);
else
    jlo = min(ctx.Jnu_flat,[],1).';
    jhi = max(ctx.Jnu_flat,[],1).';
end
lo = nan(3,1);
hi = nan(3,1);
for p = 1:3
    a = ctx.wts(:).*ctx.g(:).^p/ctx.beta;
    lo(p) = sum(min(a.*jlo,a.*jhi));
    hi(p) = sum(max(a.*jlo,a.*jhi));
end
end

function ctx = make_context(Bx,h,T,F)
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

function r = empty_record()
r = struct('Bx',NaN,'h',NaN,'y',nan(4,1),'s',NaN,'r',NaN, ...
    'residual_norm',Inf,'scaled_residual_norm',Inf, ...
    'outer_residual_norm',Inf,'accepted',false,'exitflag',NaN, ...
    'status',"",'output',struct(),'supremum_mass',NaN, ...
    'uniform_mass',NaN,'static_mesh_mass',NaN, ...
    'dynamic_lattice_mass',NaN,'dynamic_medium_mass',NaN);
end

function f = empty_failure()
f = struct('from_Bx',NaN,'trial_Bx',NaN,'step',NaN, ...
    'best_residual_norm',Inf,'best_s',NaN,'best_status',"");
end

function a = empty_anchor()
a = struct('Bx',NaN,'h',NaN,'reached',false, ...
    'high',empty_record(),'known_low_s',NaN,'known_low_r',NaN, ...
    'distinct_from_low',false,'r_high_minus_low',NaN, ...
    'local_signature_pass',false);
end
