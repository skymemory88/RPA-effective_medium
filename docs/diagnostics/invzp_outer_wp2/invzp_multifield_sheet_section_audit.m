function result = invzp_multifield_sheet_section_audit()
%INVZP_MULTIFIELD_SHEET_SECTION_AUDIT Independent multi-Bx two-sheet check.
% Diagnostic only.  At five transverse fields requested for validation,
% solve the original full-Sigma equations at prescribed static coordinate s,
% bracket every sampled static crossing, and compare r on the low/high roots.
% This tests the necessary local-potential signature r_high<r_low without
% importing the 1 T reduced continuation or its fixed-h root corrector.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
L = load(fullfile(here,'wp2_low_field_m2_census.mat'));
TX = load(fullfile(here,'wp2_hmf_node_transaction_census.mat'));
T = TX.result.T;
Bx_values = [0.5,0.8,0.9,1.1,1.2,2.5,2.7,2.9].';
% Interpolate only the choice of fixed-h probe from measured production
% first-accepted nodes, including the independent 3 T census value. The
% original-equation root search below does not use production states.
B_anchor = [L.result.summary.B;3.0];
h_anchor = [L.result.summary.first_accepted_h;0.0115294];
min_h = interp1(B_anchor,h_anchor,Bx_values,'linear');
h_values = 1.10*min_h;
sgrid = [0.05:0.05:0.95,0.975,0.99,0.995,0.999].';
opts = struct('mix',0.5,'tol',1e-11,'maxit',180);
sections = repmat(empty_section(),numel(Bx_values),1);

for ib = 1:numel(Bx_values)
    Bx = Bx_values(ib);
    h = h_values(ib);
    ctx = make_context(Bx,h,T,F);
    Sigma = zeros(size(ctx.G0));
    profile = repmat(empty_direct(numel(ctx.G0)),numel(sgrid),1);
    for k = 1:numel(sgrid)
        rec = direct_solve(ctx,sgrid(k),F.info.Jcc0,Sigma,opts);
        if ~rec.converged
            rec = direct_solve(ctx,sgrid(k),F.info.Jcc0, ...
                zeros(size(ctx.G0)),opts);
        end
        profile(k) = rec;
        if rec.converged, Sigma = rec.Sigma; end
        fprintf('multifield Bx %.3g h %.9g s %.4g R %.6g conv %d\n', ...
            Bx,h,sgrid(k),rec.static_residual,rec.converged);
    end
    brackets = zeros(0,2);
    for k = 1:numel(sgrid)-1
        if profile(k).converged && profile(k+1).converged && ...
                isfinite(profile(k).static_residual) && ...
                isfinite(profile(k+1).static_residual) && ...
                sign(profile(k).static_residual) ~= ...
                sign(profile(k+1).static_residual)
            brackets(end+1,:) = sgrid([k,k+1]); %#ok<AGROW>
        end
    end
    roots = repmat(empty_root(),size(brackets,1),1);
    for k = 1:size(brackets,1)
        rec = bisect_static(ctx,F.info.Jcc0,brackets(k,:), ...
            profile(find(sgrid == brackets(k,1),1)).Sigma,opts);
        q = invz_ordered_reduced_residual([rec.lambda;rec.x],ctx, ...
            struct('Jsup',F.info.Jcc0, ...
            'dynamic',struct('resid_tol',1e-12)));
        outer = invz_ordered_outer_map(rec.Sigma,ctx, ...
            struct('emt_static',struct('Jsup',F.info.Jcc0,'warn',false)));
        roots(k) = struct('s',rec.s,'y',[rec.lambda;rec.x], ...
            'r',q.r,'original_static_residual',rec.static_residual, ...
            'reduced_residual_norm',q.residual_norm, ...
            'outer_residual_norm',outer.residual_norm, ...
            'admissible',q.trial_admissible, ...
            'supremum_mass',q.supremum_mass,'uniform_mass',q.D_uni, ...
            'static_mesh_mass',q.static_mesh_min_mass, ...
            'dynamic_lattice_mass',q.dynamic.dynamic_min_lattice_mass, ...
            'dynamic_medium_mass',q.dynamic.dynamic_min_medium_mass);
    end
    if numel(roots) >= 2
        [~,ix] = sort([roots.s]);
        roots = roots(ix);
        r_difference = roots(end).r-roots(1).r;
        signature_pass = numel(roots) == 2 && r_difference < 0 && ...
            all([roots.admissible]) && ...
            max([roots.reduced_residual_norm]) <= 1e-8 && ...
            max([roots.outer_residual_norm]) <= 1e-8;
    else
        r_difference = NaN;
        signature_pass = false;
    end
    sections(ib) = struct('Bx',Bx,'h',h,'min_accepted_h_estimate', ...
        min_h(ib),'sgrid',sgrid,'profile',profile,'brackets',brackets, ...
        'roots',roots,'root_count',numel(roots), ...
        'r_high_minus_low',r_difference, ...
        'local_signature_pass',signature_pass);
    fprintf(['multifield section Bx %.3g: roots %d, ' ...
        'r_high-r_low %.9g, signature %d\n'], ...
        Bx,numel(roots),r_difference,signature_pass);
end

result = struct('T',T,'Bx_values',Bx_values,'h_values',h_values, ...
    'sgrid',sgrid,'sections',sections, ...
    'all_sections_two_roots',all([sections.root_count] == 2), ...
    'all_local_signatures_pass',all([sections.local_signature_pass]), ...
    'source_boundary',fullfile(here,'wp2_low_field_m2_census.mat'), ...
    'interpretation',['Each passing section verifies only the necessary ' ...
        'local sign r_high<r_low for the fold-anchored selector. A full ' ...
        'potential comparison at that Bx still requires locating and ' ...
        'integrating from its own fold.']);
save(fullfile(here,'wp2_multifield_sheet_section_audit.mat'), ...
    'result','-v7');
fprintf('multifield summary: two roots all %d, signatures all %d\n', ...
    result.all_sections_two_roots,result.all_local_signatures_pass);
end

function rec = bisect_static(ctx,Jsup,bracket,Sigma,opts)
a = bracket(1);
b = bracket(2);
ra = direct_solve(ctx,a,Jsup,Sigma,opts);
rb = direct_solve(ctx,b,Jsup,ra.Sigma,opts);
assert(ra.converged && rb.converged && ...
    sign(ra.static_residual) ~= sign(rb.static_residual));
rec = ra;
for it = 1:100
    m = (a+b)/2;
    rm = direct_solve(ctx,m,Jsup,rec.Sigma,opts);
    assert(rm.converged);
    if abs(rm.static_residual) < abs(rec.static_residual), rec = rm; end
    if sign(ra.static_residual) ~= sign(rm.static_residual)
        b = m;
        rb = rm; %#ok<NASGU>
    else
        a = m;
        ra = rm;
    end
    if abs(rec.static_residual) <= 2e-10 || b-a <= 2e-12, break; end
end
end

function rec = direct_solve(ctx,s,Jsup,Sigma,opts)
J = ctx.Jnu_flat(:);
x = -s/Jsup;
dx = 1+J*x;
M = mean(1./dx);
Phi = x*M;
K0 = mean(J./dx)/M;
map_residual = Inf;
for it = 1:opts.maxit
    med = invz_emt_scalar(ctx.G0,Sigma,ctx.Jnu_flat,struct());
    if ~med.converged, break; end
    K = med.K(:);
    K(1) = K0;
    lambda = invz_lambdas(K,ctx.g,ctx.wts,ctx.beta,[1 2 3]);
    sig = invz_sigma_ordered(ctx.tl,lambda,K,ctx.g,ctx.beta);
    map_residual = max(abs(sig.Sigma-Sigma));
    if map_residual <= opts.tol, break; end
    Sigma = (1-opts.mix)*Sigma+opts.mix*sig.Sigma;
end
if ~(isfinite(map_residual) && map_residual <= opts.tol)
    rec = empty_direct(numel(ctx.G0));
    rec.s = s;
    rec.x = x;
    rec.iterations = it;
    rec.map_residual = map_residual;
    return;
end
med = invz_emt_scalar(ctx.G0,Sigma,ctx.Jnu_flat,struct());
K = med.K(:);
K(1) = K0;
lambda = invz_lambdas(K,ctx.g,ctx.wts,ctx.beta,[1 2 3]);
sig = invz_sigma_ordered(ctx.tl,lambda,K,ctx.g,ctx.beta);
Sigma = sig.Sigma(:);
[Gstat,go] = invz_gstat_ordered(ctx.tl,lambda(1:2),K0,Sigma(1), ...
    ctx.beta,ctx.G0inel0,ctx.G0el0,struct('stable_form',true));
rec = struct('s',s,'x',x,'converged',true,'iterations',it, ...
    'map_residual',map_residual,'Sigma',Sigma,'K',K,'lambda',lambda, ...
    'K0',K0,'Phi',Phi,'Gstat',Gstat,'Gtil0',go.Gtil0, ...
    'static_residual',x-go.Gtil0);
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

function s = empty_section()
s = struct('Bx',NaN,'h',NaN,'min_accepted_h_estimate',NaN, ...
    'sgrid',nan(0,1),'profile',repmat(empty_direct(0),0,1), ...
    'brackets',zeros(0,2),'roots',repmat(empty_root(),0,1), ...
    'root_count',0,'r_high_minus_low',NaN, ...
    'local_signature_pass',false);
end

function r = empty_direct(nw)
r = struct('s',NaN,'x',NaN,'converged',false,'iterations',NaN, ...
    'map_residual',NaN,'Sigma',nan(nw,1),'K',nan(nw,1), ...
    'lambda',nan(3,1),'K0',NaN,'Phi',NaN,'Gstat',NaN, ...
    'Gtil0',NaN,'static_residual',NaN);
end

function r = empty_root()
r = struct('s',NaN,'y',nan(4,1),'r',NaN, ...
    'original_static_residual',NaN,'reduced_residual_norm',NaN, ...
    'outer_residual_norm',NaN,'admissible',false, ...
    'supremum_mass',NaN,'uniform_mass',NaN,'static_mesh_mass',NaN, ...
    'dynamic_lattice_mass',NaN,'dynamic_medium_mass',NaN);
end
