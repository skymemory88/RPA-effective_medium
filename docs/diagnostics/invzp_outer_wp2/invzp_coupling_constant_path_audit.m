function result = invzp_coupling_constant_path_audit()
%INVZP_COUPLING_CONSTANT_PATH_AUDIT Fixed-state coupling-path sheet audit.
% Diagnostic only.  At Bx=1 T, h=0.006 meV, and the production temperature,
% continue both independently verified ordered roots under
%
%   J_cc(q) -> rho J_cc(q),  J_cc(0) -> rho J_cc(0),  0 < rho <= 1,
%
% while holding the single-ion state (Bx,h,T,J_aa(0)) fixed.  This is the
% fixed-moment coupling path: it scales Jensen's fluctuation Hamiltonian H1,
% not the transverse mean-field interaction that defines the reference ion.
%
% For an exact fixed-moment constrained potential,
%
%   d Phi/d rho = <H1(rho)>/rho,
%   <H1(rho)>/N = (1/(2 beta)) sum_n w_n K_n(rho) G_n(rho)       (J 2.21),
%
% so a branch that reaches the common rho=0 anchor supplies a candidate
% fluctuation-potential difference.  In this boundary-preserving hybrid
% implementation that identity is an audit hypothesis, not proof that the
% closure is Phi-derivable; agreement with an independent field/S3 route is
% still required before using the result as a production sheet selector.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
A = load(fullfile(here,'wp2_pseudoarclength_original_equations_audit.mat'));
T = A.result.T;
Bx = A.result.Bx;
h = A.result.hsection;
Jbase = F.J(:);
J0base = F.info.Jcc0;
Jsupbase = A.result.Jsup;
base = make_context(Bx,h,T,F);
[lambda_lo_base,lambda_hi_base] = lambda_bounds(base);
lambda_mid_base = (lambda_lo_base+lambda_hi_base)/2;
lambda_half_base = (lambda_hi_base-lambda_lo_base)/2;
assert(all(lambda_half_base > 0));

% A moderate uniform grid resolves the coupling integral; logarithmic tail
% nodes test the noninteracting limit.  Failure is refined separately and is
% not silently skipped.
rho_grid = [1:-0.05:0.1,0.075,0.05,0.035,0.025,0.015,0.01, ...
    0.007,0.005,0.003,0.002,0.001].';
rho_grid = unique(rho_grid,'stable');
lsopts = optimoptions('lsqnonlin','Display','off', ...
    'FunctionTolerance',1e-22,'StepTolerance',2e-13, ...
    'OptimalityTolerance',2e-12,'MaxIterations',12, ...
    'MaxFunctionEvaluations',60,'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',3e-6);
solve_opts = struct('Jsup',NaN,'dynamic',struct('resid_tol',1e-12));
branches = repmat(empty_branch(),numel(A.result.roots),1);
checkpoint_file = fullfile(here, ...
    'wp2_coupling_constant_path_checkpoint.mat');
checkpoint_schema = 2;
first_branch = 1;
if exist(checkpoint_file,'file')
    C = load(checkpoint_file,'checkpoint');
    cp = C.checkpoint;
    compatible = isfield(cp,'schema') && cp.schema == checkpoint_schema && ...
        isequal(cp.rho_grid,rho_grid) && cp.T == T && ...
        cp.Bx == Bx && cp.h == h && cp.J0base == J0base;
    if compatible && cp.completed_branch >= 1 && ...
            cp.completed_branch < numel(branches)
        branches = cp.branches;
        first_branch = cp.completed_branch+1;
        fprintf('resuming after completed coupling branch %d\n', ...
            cp.completed_branch);
    end
end

for ib = first_branch:numel(A.result.roots)
    source = A.result.roots(ib);
    z0 = encode_state([source.lambda(:);source.x],1);
    records = repmat(empty_record(),numel(rho_grid),1);
    records(1) = evaluate_record(1,z0,z0,0,struct(),false);
    records(1) = add_jacobian(records(1));
    assert(records(1).accepted, ...
        'Source root %d does not satisfy the coupling-path acceptance gate.',ib);
    naccepted = 1;
    failure = empty_failure();
    fprintf(['coupling branch %d: rho %.6g, ||R||inf %.3g, ' ...
        'I %.9g, masses [%.4g %.4g %.4g %.4g %.4g]\n'], ...
        ib,records(1).rho,records(1).physical_residual_norm, ...
        records(1).coupling_integrand,records(1).supremum_mass, ...
        records(1).uniform_mass,records(1).static_mesh_mass, ...
        records(1).dynamic_lattice_mass,records(1).dynamic_medium_mass);

    for ir = 2:numel(rho_grid)
        rho = rho_grid(ir);
        if naccepted >= 2
            r1 = records(naccepted);
            r0 = records(naccepted-1);
            zpred = r1.z+(rho-r1.rho)*(r1.z-r0.z)/(r1.rho-r0.rho);
        else
            zpred = records(naccepted).z;
        end
        zprev = records(naccepted).z;
        starts = [zpred,zprev];
        rec = empty_record();
        found = false;
        for is = 1:size(starts,2)
            candidate = solve_at_rho(rho,starts(:,is),zpred);
            if candidate.accepted
                rec = candidate;
                found = true;
                break;
            end
        end
        if ~found
            failure = refine_failure(records(naccepted),rho);
            break;
        end

        % A large discontinuous jump is evidence of branch loss, not a
        % license to relabel another root as the continued sheet.
        if rec.predictor_distance > 0.35
            failure = refine_failure(records(naccepted),rho);
            failure.status = "candidate_jump";
            failure.candidate = rec;
            break;
        end
        if mod(ir-1,4) == 0 || rho <= 0.01
            rec = add_jacobian(rec);
        end
        naccepted = naccepted+1;
        records(naccepted) = rec;
        fprintf(['coupling branch %d: rho %.6g, ||R||inf %.3g, ' ...
            'I %.9g, svmin %.3g, pred %.3g\n'],ib,rho, ...
            rec.physical_residual_norm,rec.coupling_integrand, ...
            rec.jacobian_svmin,rec.predictor_distance);
    end
    records = records(1:naccepted);
    summary = summarize_branch(records,failure);
    branches(ib) = struct('source_root',source,'records',records, ...
        'failure',failure,'summary',summary);
    checkpoint = struct('schema',checkpoint_schema, ...
        'T',T,'Bx',Bx,'h',h,'J0base',J0base, ...
        'rho_grid',rho_grid,'completed_branch',ib,'branches',branches);
    save(checkpoint_file,'checkpoint','-v7');
end

common_anchor = branches(1).summary.reached_tail && ...
    branches(2).summary.reached_tail && ...
    norm(branches(1).records(end).y-branches(2).records(end).y,Inf) ...
        <= 1e-5;
summaries = [branches.summary];
candidate_phi = [summaries.candidate_delta_phi];
potential_comparison_defined = common_anchor && all(isfinite(candidate_phi));
if potential_comparison_defined
    delta_phi = candidate_phi;
    [~,preferred_branch] = min(delta_phi);
else
    delta_phi = [NaN NaN];
    preferred_branch = NaN;
end

result = struct('T',T,'Bx',Bx,'h',h,'J0base',J0base, ...
    'Jsupbase',Jsupbase,'rho_grid',rho_grid,'branches',branches, ...
    'common_anchor',common_anchor, ...
    'potential_comparison_defined',potential_comparison_defined, ...
    'candidate_delta_phi',delta_phi,'candidate_preferred_branch', ...
    preferred_branch,'lambda_lo_base',lambda_lo_base, ...
    'lambda_hi_base',lambda_hi_base, ...
    'source',fullfile(here, ...
        'wp2_pseudoarclength_original_equations_audit.mat'), ...
    'convention',['Only longitudinal fluctuation couplings Jcc and Jcc(0) ' ...
        'are scaled; Bx, h, T, and transverse Jaa(0) are fixed.'], ...
    'normalization',['Per-site interaction energy is (2 beta)^-1 times ' ...
        'the doubled nonnegative-Matsubara sum of K*G; the coupling ' ...
        'integrand divides that energy by rho.'], ...
    'interpretation',['A finite integral is only a candidate constrained ' ...
        'fluctuation potential for this hybrid closure. It becomes a sheet ' ...
        'selector only after a common rho=0 anchor and an independent ' ...
        'thermodynamic consistency check are established.']);
save(fullfile(here,'wp2_coupling_constant_path_audit.mat'),'result','-v7');

fprintf(['coupling-path summary: reached tail %d/%d, common anchor %d, ' ...
    'comparison defined %d, candidate dPhi [%.12g %.12g]\n'], ...
    branches(1).summary.reached_tail,branches(2).summary.reached_tail, ...
    common_anchor,potential_comparison_defined,delta_phi);

    function rec = solve_at_rho(rho,zstart,zpred)
        [lb,ub] = z_bounds(rho);
        zstart = min(max(zstart,lb+1e-10),ub-1e-10);
        tic;
        [z,resnorm,rscaled,exitflag,output] = lsqnonlin( ...
            @(zz) scaled_vector(zz,rho),zstart,lb,ub,lsopts);
        elapsed = toc;
        rec = evaluate_record(rho,z,zpred,exitflag,output,false);
        rec.start = zstart;
        rec.resnorm = resnorm;
        rec.scaled_residual = rscaled(:);
        rec.scaled_residual_norm = max(abs(rscaled));
        rec.elapsed_s = elapsed;
        rec.accepted = rec.accepted && exitflag > 0 && ...
            rec.scaled_residual_norm <= 2e-8;
    end

    function rec = evaluate_record(rho,z,zpred,exitflag,output,need_jacobian)
        ctxeval = scaled_context(rho);
        solve_opts.Jsup = rho*Jsupbase;
        y = decode_state(z,rho);
        qeval = invz_ordered_reduced_residual(y,ctxeval,solve_opts);
        if qeval.defined && all(isfinite(qeval.residual))
            rscaled = scale_residual(qeval,rho);
            if need_jacobian
                Jscaled = finite_jacobian(@(zz) scaled_vector(zz,rho),z);
                sv = svd(Jscaled);
            else
                sv = nan(4,1);
            end
            eint = real(sum(base.wts(:).*qeval.K(:).*qeval.G(:))) ...
                /(2*base.beta);
            cint = eint/rho;
            dynL = qeval.dynamic.dynamic_min_lattice_mass;
            dynM = qeval.dynamic.dynamic_min_medium_mass;
            unresolved = sum(qeval.dynamic.fallback_unresolved(2:end));
            xi_denom = qeval.static.xi_denom;
            gstat_local_denom = qeval.static.gstat_local_denom;
        else
            rscaled = nan(4,1);
            sv = nan(4,1);
            eint = NaN;
            cint = NaN;
            dynL = NaN;
            dynM = NaN;
            unresolved = NaN;
            xi_denom = NaN;
            gstat_local_denom = NaN;
        end
        accepted = qeval.defined && qeval.trial_admissible && ...
            qeval.residual_norm <= 2e-8 && max(abs(rscaled)) <= 2e-8 && ...
            all(isfinite([eint,cint])) && unresolved == 0;
        rec = struct('rho',rho,'start',nan(4,1),'z',z(:),'y',y, ...
            'predictor',zpred(:),'predictor_distance',norm(z-zpred,Inf), ...
            'exitflag',exitflag,'output',output,'resnorm',NaN, ...
            'scaled_residual',rscaled,'scaled_residual_norm', ...
            max(abs(rscaled)),'physical_residual',qeval.residual, ...
            'physical_residual_norm',qeval.residual_norm,'status', ...
            string(qeval.status),'defined',qeval.defined,'admissible', ...
            qeval.trial_admissible,'accepted',accepted,'K',qeval.K, ...
            'G',qeval.G,'Sigma',qeval.Sigma,'interaction_energy',eint, ...
            'coupling_integrand',cint,'jacobian_singular_values',sv, ...
            'jacobian_svmin',min(sv),'supremum_mass',qeval.supremum_mass, ...
            'uniform_mass',qeval.D_uni,'static_mesh_mass', ...
            qeval.static_mesh_min_mass,'dynamic_lattice_mass',dynL, ...
            'dynamic_medium_mass',dynM,'xi_denom',xi_denom, ...
            'gstat_local_denom',gstat_local_denom, ...
            'dynamic_unresolved_frequencies',unresolved,'elapsed_s',0);
    end

    function rec = add_jacobian(rec)
        Jscaled = finite_jacobian( ...
            @(zz) scaled_vector(zz,rec.rho),rec.z);
        sv = svd(Jscaled);
        rec.jacobian_singular_values = sv;
        rec.jacobian_svmin = min(sv);
    end

    function failure = refine_failure(last,rho_failed)
        lo = rho_failed;
        hi = last.rho;
        best = last;
        attempts = repmat(struct('rho',NaN,'record',empty_record()),0,1);
        for it = 1:6
            mid = (lo+hi)/2;
            rec = solve_at_rho(mid,best.z,best.z);
            attempts(end+1) = struct('rho',mid,'record',rec); %#ok<AGROW>
            if rec.accepted && rec.predictor_distance <= 0.35
                hi = mid;
                best = rec;
            else
                lo = mid;
            end
            if hi-lo <= 5e-4, break; end
        end
        best = add_jacobian(best);
        failure = struct('status',"corrector_failed", ...
            'last_accepted_rho',best.rho,'first_failed_rho',lo, ...
            'last_accepted',best,'candidate',empty_record(), ...
            'attempts',attempts);
    end

    function summary = summarize_branch(records,failure)
        rho = [records.rho].';
        cint = [records.coupling_integrand].';
        reached_tail = rho(end) <= min(rho_grid)+10*eps;
        if reached_tail
            [rho_asc,ix] = sort(rho,'ascend');
            cint_asc = cint(ix);
            p = tail_power(rho_asc,cint_asc);
            if isfinite(p) && p > -0.9
                tail = cint_asc(1)*rho_asc(1)/(p+1);
            else
                tail = 0.5*cint_asc(1)*rho_asc(1);
            end
            candidate = trapz(rho_asc,cint_asc)+tail;
            coarse = coarse_integral(rho_asc,cint_asc)+tail;
            quadrature_change = abs(candidate-coarse);
        else
            p = NaN;
            tail = NaN;
            candidate = NaN;
            quadrature_change = NaN;
        end
        last = records(end);
        % The ordered hybrid's zero-coupling static reference includes the
        % transverse-feedback term used to construct G0el0; it is therefore
        % G0inel0+G0el0, not the ordinary-Dyson ctx.G0(1).
        yanchor = [zeros(3,1);base.G0inel0+base.G0el0];
        summary = struct('reached_tail',reached_tail, ...
            'last_rho',last.rho,'last_y',last.y, ...
            'noninteracting_expected_y',yanchor, ...
            'anchor_y_error',norm(last.y-yanchor,Inf), ...
            'last_scaled_moments',last.y(1:3)/last.rho, ...
            'last_sigma_norm',norm(last.Sigma,Inf), ...
            'tail_power',p,'tail_integral',tail, ...
            'candidate_delta_phi',candidate, ...
            'quadrature_coarsening_change',quadrature_change, ...
            'minimum_masses',[min([records.supremum_mass]), ...
                min([records.uniform_mass]), ...
                min([records.static_mesh_mass]), ...
                min([records.dynamic_lattice_mass]), ...
                min([records.dynamic_medium_mass])], ...
            'minimum_jacobian_sv', ...
                min([records.jacobian_svmin],[],'omitnan'), ...
            'failure_status',failure.status);
    end

    function ctxout = scaled_context(rho)
        ctxout = base;
        ctxout.Jnu_flat = rho*Jbase;
        ctxout.J0eff = rho*J0base;
    end

    function z = encode_state(y,rho)
        z = [(y(1:3)/rho-lambda_mid_base)./lambda_half_base; ...
            -J0base*y(4)];
    end

    function y = decode_state(z,rho)
        y = [rho*(lambda_mid_base+lambda_half_base.*z(1:3)); ...
            -z(4)/J0base];
    end

    function [lb,ub] = z_bounds(rho)
        lb = [-ones(3,1)+1e-9;1e-9];
        ub = [ones(3,1)-1e-9;min(3,(1-1e-8)/rho)];
    end

    function r = scaled_vector(z,rho)
        ctxscale = scaled_context(rho);
        solve_opts.Jsup = rho*Jsupbase;
        qscale = invz_ordered_reduced_residual( ...
            decode_state(z,rho),ctxscale,solve_opts);
        if qscale.defined && all(isfinite(qscale.residual))
            r = scale_residual(qscale,rho);
        else
            r = 1e3*(1+abs(z));
        end
    end

    function r = scale_residual(qin,rho)
        r = [qin.lambda_residual./(rho*lambda_half_base); ...
            J0base*qin.static_residual];
    end
end

function J = finite_jacobian(fun,z)
f0 = fun(z);
J = nan(numel(f0),numel(z));
for k = 1:numel(z)
    dz = 2e-5*max(1,abs(z(k)));
    zp = z;
    zp(k) = zp(k)+dz;
    J(:,k) = (fun(zp)-f0)/dz;
end
end

function p = tail_power(rho,cint)
n = min(4,numel(rho));
rr = rho(1:n);
ii = cint(1:n);
if n >= 2 && all(rr > 0) && all(isfinite(ii)) && ...
        all(abs(ii) > realmin) && all(sign(ii) == sign(ii(1)))
    c = polyfit(log(rr),log(abs(ii)),1);
    p = c(1);
else
    p = NaN;
end
end

function v = coarse_integral(rho,cint)
if numel(rho) <= 3
    v = trapz(rho,cint);
    return;
end
ix = unique([1:2:numel(rho),numel(rho)]);
v = trapz(rho(ix),cint(ix));
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

function b = empty_branch()
b = struct('source_root',struct(),'records',repmat(empty_record(),0,1), ...
    'failure',empty_failure(),'summary',struct());
end

function f = empty_failure()
f = struct('status',"not_applicable",'last_accepted_rho',NaN, ...
    'first_failed_rho',NaN,'last_accepted',empty_record(), ...
    'candidate',empty_record(), ...
    'attempts',repmat(struct('rho',NaN,'record',empty_record()),0,1));
end

function r = empty_record()
r = struct('rho',NaN,'start',nan(4,1),'z',nan(4,1),'y',nan(4,1), ...
    'predictor',nan(4,1),'predictor_distance',NaN,'exitflag',NaN, ...
    'output',struct(),'resnorm',NaN,'scaled_residual',nan(4,1), ...
    'scaled_residual_norm',NaN,'physical_residual',nan(4,1), ...
    'physical_residual_norm',NaN,'status',"",'defined',false, ...
    'admissible',false,'accepted',false,'K',nan(0,1),'G',nan(0,1), ...
    'Sigma',nan(0,1),'interaction_energy',NaN, ...
    'coupling_integrand',NaN,'jacobian_singular_values',nan(4,1), ...
    'jacobian_svmin',NaN,'supremum_mass',NaN,'uniform_mass',NaN, ...
    'static_mesh_mass',NaN,'dynamic_lattice_mass',NaN, ...
    'dynamic_medium_mass',NaN,'xi_denom',NaN, ...
    'gstat_local_denom',NaN,'dynamic_unresolved_frequencies',NaN, ...
    'elapsed_s',NaN);
end
