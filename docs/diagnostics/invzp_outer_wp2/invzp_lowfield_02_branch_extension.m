function result = invzp_lowfield_02_branch_extension()
%INVZP_LOWFIELD_02_BRANCH_EXTENSION Extend both roots from 0.5 to 0.2 T.
% Diagnostic only.  Below the measured 0.5 T probe, extrapolate the local
% diagnostic h(Bx) linearly using the 0.5/0.8 T section choices, then
% independently correct both exact reduced roots.  The extrapolation chooses
% test sections; it is not a claim about the production acceptance boundary.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
S = load(fullfile(here,'wp2_multifield_sheet_section_audit.mat'));
C = load(fullfile(here,'wp2_multifield_crossfield_continuation.mat'));
T = S.result.T;
Jsup = F.info.Jcc0;
sections = S.result.sections;
i05 = find(abs([sections.Bx]-0.5) <= 1e-12,1);
i08 = find(abs([sections.Bx]-0.8) <= 1e-12,1);
assert(~isempty(i05) && ~isempty(i08));
h_slope = (sections(i08).h-sections(i05).h)/0.3;
h_of_B = @(B) sections(i05).h+h_slope*(B-0.5);

low_start = sections(i05).roots(1).y(:);
assert(abs(C.result.down(end).Bx-0.5) <= 1e-12);
high_start = C.result.down(end).y(:);
Bx_grid = [(0.5:-0.025:0.375),(0.37:-0.005:0.2)].';
ls4 = optimoptions('lsqnonlin','Display','off', ...
    'FunctionTolerance',1e-20,'StepTolerance',2e-13, ...
    'OptimalityTolerance',2e-11,'MaxIterations',35, ...
    'MaxFunctionEvaluations',260,'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',3e-6);
config = struct('residual_tolerance',1e-8,'s_window',0.15);

low = continue_branch(Bx_grid,h_of_B,low_start,T,F,Jsup,ls4,config);
high = continue_branch(Bx_grid,h_of_B,high_start,T,F,Jsup,ls4,config);
anchors = [0.5,0.4,0.3,0.2];
summary = repmat(empty_summary(),numel(anchors),1);
for k = 1:numel(anchors)
    il = find(abs([low.Bx]-anchors(k)) <= 5e-12,1);
    ih = find(abs([high.Bx]-anchors(k)) <= 5e-12,1);
    reached = ~isempty(il) && ~isempty(ih);
    if reached
        r_difference = high(ih).r-low(il).r;
        distinct = high(ih).s-low(il).s > 1e-6;
        signature = distinct && r_difference < 0;
        low_record = low(il);
        high_record = high(ih);
    else
        r_difference = NaN;
        distinct = false;
        signature = false;
        low_record = empty_record();
        high_record = empty_record();
    end
    summary(k) = struct('Bx',anchors(k),'h',h_of_B(anchors(k)), ...
        'reached',reached,'low',low_record,'high',high_record, ...
        'distinct',distinct,'r_high_minus_low',r_difference, ...
        'local_signature_pass',signature);
end

result = struct('T',T,'Jsup',Jsup,'Bx_grid',Bx_grid, ...
    'h_slope',h_slope,'h_at_05',sections(i05).h, ...
    'low',low,'high',high,'summary',summary,'config',config, ...
    'source_crossfield', ...
    fullfile(here,'wp2_multifield_crossfield_continuation.mat'), ...
    'interpretation',['The two exact roots are tested below 0.5 T on ' ...
    'explicit linearly extrapolated diagnostic h(Bx) sections. This ' ...
    'validates branch persistence and the local r ordering but does not ' ...
    'validate the extrapolated h choice as a production threshold.']);
save(fullfile(here,'wp2_lowfield_02_branch_extension.mat'), ...
    'result','-v7');
for k = 1:numel(summary)
    a = summary(k);
    fprintf(['lowfield anchor %.3g h %.9g reached %d: ' ...
        's low/high %.9g %.9g, r_high-r_low %.9g, signature %d\n'], ...
        a.Bx,a.h,a.reached,a.low.s,a.high.s, ...
        a.r_high_minus_low,a.local_signature_pass);
end
end

function records = continue_branch(Bx_grid,h_of_B,y_start,T,F,Jsup, ...
    ls4,config)
records = repmat(empty_record(),numel(Bx_grid),1);
yseed = y_start;
s_previous = -Jsup*y_start(4);
for k = 1:numel(Bx_grid)
    Bx = Bx_grid(k);
    h = h_of_B(Bx);
    ctx = make_context(Bx,h,T,F);
    [lo,hi] = lambda_bounds(ctx);
    mid = (lo+hi)/2;
    half = (hi-lo)/2;
    zseed = [(yseed(1:3)-mid)./half;s_previous];
    sbounds = [max(1e-8,s_previous-config.s_window), ...
        min(1-1e-8,s_previous+config.s_window)];
    rec = solve_scaled(Bx,h,zseed,sbounds,ctx,mid,half, ...
        Jsup,ls4,config);
    if ~rec.accepted
        records = records(1:k-1);
        return;
    end
    records(k) = rec;
    yseed = rec.y;
    s_previous = rec.s;
    fprintf('lowfield branch Bx %.4g h %.8g s %.9g r %.9g R %.3g\n', ...
        rec.Bx,rec.h,rec.s,rec.r,rec.residual_norm);
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

function a = empty_summary()
a = struct('Bx',NaN,'h',NaN,'reached',false, ...
    'low',empty_record(),'high',empty_record(),'distinct',false, ...
    'r_high_minus_low',NaN,'local_signature_pass',false);
end
