function tests = test_invzt_convergence_diagnostics
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
repo = fullfile(here, '..', '..');
addpath(repo, fullfile(repo, 'invz_common'), fullfile(repo, 'invz_projected'), ...
    fullfile(repo, 'invz_tensor'), ...
    fullfile(repo, 'docs', 'diagnostics', 'invzt_convergence'));
tc.TestData.cfg = struct( ...
    'fields_T', [3 6], ...
    'grid_specs', struct('N', 4, 'conv', 'halfopen'), ...
    'lattice_specs', {{struct('name', 'test_brute10', ...
        'dipole', 'bruteforce', 'dpRng', 10)}}, ...
    'repeat_fields_T', 6, ...
    'output_mat', '');
end

function test_wp0_cold_capture_and_repeat(tc)
r = invzt_wp0_baseline_measure(tc.TestData.cfg);
verifyTrue(tc, r.acceptance.passed);
verifyEqual(tc, numel(r.repeat_checks), 1);
verifyTrue(tc, r.repeat_checks.passed);

pm3 = find([r.records.B_T] == 3 & strcmp({r.records.leg}, 'pm'), 1);
ord3 = find([r.records.B_T] == 3 & strcmp({r.records.leg}, 'ordered'), 1);
pm6 = find([r.records.B_T] == 6 & strcmp({r.records.leg}, 'pm'), 1);
verifyTrue(tc, r.records(pm3).converged);
verifyEqual(tc, r.records(pm3).reason, 'unstable');
verifyEqual(tc, r.records(ord3).reason, 'accepted_ordered');
verifyEqual(tc, r.records(pm6).reason, 'accepted_pm');
verifyGreaterThan(tc, r.records(pm3).raw_index, 0);
verifyGreaterThan(tc, r.records(ord3).raw_index, 0);
verifyGreaterThan(tc, r.records(pm6).raw_index, 0);
end

function test_wp05_frozen_candidates_share_inputs(tc)
r = invzt_wp0_baseline_measure(tc.TestData.cfg);
q = invzt_wp05_frozen_compare(r, struct('grid_N', 4, 'dpRng', 10, ...
    'output_mat', ''));
verifyEqual(tc, numel(q.anchors), 3);
verifyTrue(tc, all([q.anchors.inputs_identical]));
verifyFalse(tc, q.interpretation.formal_equivalence_proved);
verifyEqual(tc, q.interpretation.preferred_representation, 'unselected');
for k = 1:numel(q.anchors)
    verifyLessThanOrEqual(tc, q.anchors(k).decomposition_error, 1e-12);
    verifyEqual(tc, q.anchors(k).split.min_one_plus_sigma, ...
        q.anchors(k).whole.min_one_plus_sigma, 'AbsTol', 0);
    verifyEqual(tc, q.anchors(k).split.crit, ...
        q.anchors(k).split.gamma_mass_eigenvalues(1), 'AbsTol', 1e-12);
end
end

function test_auto_attempt_capture_is_opt_in(tc)
ion = invz_ion();
lat = invzt_jq_tensor(ion, invzt_qgrid(4, 'halfopen'), ...
    struct('dpRng', 10, 'cache', false));
[~, ~, d0] = invzt_solve_auto(ion, 0.1, [3 0 0], lat, struct());
verifyEmpty(tc, d0.para.point);
verifyEmpty(tc, d0.ordered.point);
[~, phase, d1] = invzt_solve_auto(ion, 0.1, [3 0 0], lat, ...
    struct('capture_attempts', true));
verifyEqual(tc, phase, 1);
verifyNotEmpty(tc, d1.para.point);
verifyNotEmpty(tc, d1.ordered.point);
verifyEqual(tc, d1.ordered.hmf_J0z, d1.ordered.point.J0z_mf, ...
    'RelTol', 1e-14);
end

function test_ordered_map_reproduces_current_whole_costate(tc)
ion = invz_ion();
lat = invzt_jq_tensor(ion, invzt_qgrid(4, 'halfopen'), ...
    struct('dpRng', 10, 'cache', false));
[~, phase, di] = invzt_solve_auto(ion, 0.1, [3 0 0], lat, ...
    struct('capture_attempts', true));
verifyEqual(tc, phase, 1);
p = di.ordered.point;
whole_ctx = invzt_a1_ordered_context(p.si, p.tl, 0.1, lat, ...
    struct('dominant_count', numel(p.si.E), ...
    'selector_source', 'full_rank_current_whole_reference'));
q = invzt_a1_ordered_map(whole_ctx, p.Sigma);
verifyEqual(tc, q.status, 'evaluated');
verifyEqual(tc, q.representation, 'full_rank_reference');
verifyEqual(tc, q.K, p.K, 'AbsTol', 1e-12);
verifyEqual(tc, q.G, p.G, 'AbsTol', 1e-12);
verifyEqual(tc, q.lambda, p.lambda, 'AbsTol', 1e-12);
verifyEqual(tc, q.residual_inf, p.outer_residual, 'AbsTol', 1e-12);
verifyEqual(tc, q.crit, p.crit, 'AbsTol', 1e-12);

split_ctx = invzt_a1_ordered_context(p.si, p.tl, 0.1, lat, ...
    struct('dominant_count', di.para.point.mspec.ndom, ...
    'selector_source', 'same_point_pm_fixed_rank'));
qs = invzt_a1_ordered_map(split_ctx, p.Sigma);
verifyEqual(tc, qs.status, 'evaluated');
verifyTrue(tc, qs.valid);
verifyEqual(tc, qs.representation, 'fixed_rank_split_candidate');
verifyEqual(tc, split_ctx.mspec.ndom, di.para.point.mspec.ndom);
bad = invzt_a1_ordered_map(split_ctx, -2*ones(size(p.Sigma)));
verifyEqual(tc, bad.status, 'sigma_domain');
verifyEqual(tc, bad.failure, 'sigma_domain');
verifyEmpty(tc, bad.residual);
end

function test_ordered_split_map_reduces_to_pm_boundary_map(tc)
ion = invz_ion();
lat = invzt_jq_tensor(ion, invzt_qgrid(4, 'halfopen'), ...
    struct('dpRng', 10, 'cache', false));
[~, phase, di] = invzt_solve_auto(ion, 0.1, [3 0 0], lat, ...
    struct('capture_attempts', true));
verifyEqual(tc, phase, 1);
p = di.para.point;
ctx = invzt_a1_ordered_context(p.si, p.tl, 0.1, lat, ...
    struct('dominant_count', p.mspec.ndom, ...
    'selector_source', 'same_point_pm_fixed_rank'));
q = invzt_a1_ordered_map(ctx, p.Sigma);
verifyEqual(tc, q.status, 'evaluated');
verifyEqual(tc, q.G, p.G, 'AbsTol', 1e-12);
verifyEqual(tc, q.K, p.K, 'AbsTol', 1e-12);
verifyEqual(tc, q.lambda(1:2), p.lambda, 'AbsTol', 1e-12);
verifyEqual(tc, q.residual_inf, p.outer_residual, 'AbsTol', 1e-12);

ratio = real(q.chi_tilde(3,3,1)) / real(ctx.cfull(3,3,1));
verifyEqual(tc, ratio, di.ordered.handoff_ratio, 'AbsTol', 1e-14);
verifyEqual(tc, lat.info.Jcc0*ratio, di.ordered.hmf_J0z, 'AbsTol', 1e-14);
end

function test_diagnostic_iterator_reproduces_captured_whole_root(tc)
ion = invz_ion();
lat = invzt_jq_tensor(ion, invzt_qgrid(4, 'halfopen'), ...
    struct('dpRng', 10, 'cache', false));
[~, phase, di] = invzt_solve_auto(ion, 0.1, [3 0 0], lat, ...
    struct('capture_attempts', true));
verifyEqual(tc, phase, 1);
p = di.ordered.point;
ctx = invzt_a1_ordered_context(p.si, p.tl, 0.1, lat, ...
    struct('dominant_count', numel(p.si.E), ...
    'selector_source', 'full_rank_current_map_reference'));
r = invzt_a1_map_iterate(ctx, p.Sigma, struct());
verifyTrue(tc, r.converged);
verifyEqual(tc, r.status, 'converged');
verifyEqual(tc, r.iterations, 1);
verifyEqual(tc, r.terminal.residual_inf, p.outer_residual, 'AbsTol', 1e-12);
verifyEqual(tc, r.terminal.Sigma, p.Sigma, 'AbsTol', 0);
end

function test_cartesian_c_only_map_reduces_to_scalar_medium(tc)
T = 1.0; C = invz_const(); beta = 1/(C.kB*T);
Delta = 0.5; m = 0.4; M = 0.9;
E = [0; Delta; 100]; P = exp(-beta*E); P = P/sum(P);
Mz = [m M 0; M -m 0; 0 0 0];
si = struct('E',E,'P',P,'Mx',zeros(3),'My',zeros(3),'Mz',Mz, ...
    'Jexp',[0;0;real(diag(Mz)).'*P], ...
    'mf_converged',true,'hz',0.01);
n01 = tanh(Delta/(2*C.kB*T));
tl = struct('Delta',Delta,'M2',M^2,'m',m,'n01',n01,'g0',2*n01/Delta);

j = 1e-3; u = ones(4,1)/2; Jpage = zeros(12); cc = 3:3:12;
Jpage(cc,cc) = j*(u*u.');
lat = struct('Jt',reshape(Jpage,12,12,1),'JtGamma',Jpage, ...
    'w',1,'conv','synthetic_c_only', ...
    'info',struct('Jcc0',j,'Jaa0',0));
ctx = invzt_a1_ordered_context(si,tl,T,lat, ...
    struct('Ecut',1,'dominant_count',2, ...
    'selector_source','synthetic_c_only_fixed_rank'));
Sigma = 0.2*ones(numel(ctx.wn),1);
q = invzt_a1_ordered_map(ctx,Sigma);
verifyLessThanOrEqual(tc,max(abs(ctx.crest(:))),1e-14);
verifyTrue(tc,q.valid);

G0 = -ctx.cdom_cc;
med = invz_emt_scalar(G0,Sigma,[j;0;0;0],struct('debug',true));
lam = invz_lambdas(med.K,ctx.g,ctx.wts,ctx.beta,[1 2 3]);
sg = invz_sigma_ordered(ctx.tl,lam,med.K,ctx.g,ctx.beta);
verifyTrue(tc,med.converged);
verifyLessThanOrEqual(tc,med.closure,1e-12);
verifyEqual(tc,q.G,med.G,'AbsTol',5e-12);
verifyEqual(tc,q.K,med.K,'AbsTol',5e-12);
verifyEqual(tc,q.lambda,lam,'AbsTol',1e-11);
verifyEqual(tc,q.Sigma_next,sg.Sigma,'AbsTol',1e-11);
verifyEqual(tc,q.residual,sg.Sigma-Sigma,'AbsTol',1e-11);

c0 = real(ctx.cdom(3,3,1));
mass_projected = 1+Sigma(1)-j*c0;
verifyEqual(tc,q.crit,mass_projected/(1+Sigma(1)),'AbsTol',1e-12);
verifyEqual(tc,q.crit_active_rank,4);
verifyEqual(tc,q.crit_clipped_mass,0,'AbsTol',1e-14);
end

function test_wp1_selfconsistent_runner_is_diagnostic_only(tc)
b = invzt_wp0_baseline_measure(tc.TestData.cfg);
r = invzt_wp1_selfconsistent_compare(b,struct('grid_N',4,'dpRng',10, ...
    'iterator_options',struct('max_iter',2,'tol',1e-8, ...
    'mix',0.7,'anderson_depth',2),'output_mat',''));
verifyEqual(tc,numel(r.runs),12);
verifyTrue(tc,r.acceptance.all_declared_runs_executed);
verifyTrue(tc,r.acceptance.no_production_integration);
verifyFalse(tc,r.acceptance.physical_representation_selected);
verifyEqual(tc,sort(unique(string({r.runs.seed_label}))), ...
    sort(["captured_state" "cold_zero"]));
verifyEqual(tc,sort(unique(string({r.runs.representation}))), ...
    sort(["full_rank_reference" "split_fixed_rank"]));
end

function test_ordered_realaxis_full_rank_reproduces_current_path(tc)
ion=invz_ion(); T=0.1;
lat=invzt_jq_tensor(ion,invzt_qgrid(4,'halfopen'), ...
    struct('dpRng',10,'cache',false));
[~,phase,di]=invzt_solve_auto(ion,T,[3 0 0],lat, ...
    struct('capture_attempts',true));
verifyEqual(tc,phase,1);
p=di.ordered.point;
ctx=invzt_a1_ordered_context(p.si,p.tl,T,lat, ...
    struct('dominant_count',numel(p.si.E), ...
    'selector_source','full_rank_current_map_reference'));
q=invzt_a1_ordered_map(ctx,p.Sigma);
w=[0.005;0.02;0.05];
a=invzt_a1_ordered_realaxis_diagnostic(ion,ctx,q,w, ...
    struct('qsel','gamma_uniform'));
b=invzt_chi_realaxis(ion,T,[3 0 0],p,w, ...
    struct('qsel','gamma_uniform'));
verifyEqual(tc,a.Sigma_w,b.Sigma_w,'AbsTol',1e-12);
verifyEqual(tc,a.chi_uniform,b.chi_uniform,'AbsTol',1e-11);
verifyLessThanOrEqual(tc,max(abs(a.crest(:))),1e-12);
verifyLessThanOrEqual(tc,a.decomposition_error,1e-12);

qvec=[0.1 0 0;0.03 -0.02 0.01];
ae=invzt_a1_ordered_realaxis_diagnostic(ion,ctx,q,w, ...
    struct('qsel',qvec,'cache',false));
be=invzt_chi_realaxis(ion,T,[3 0 0],p,w, ...
    struct('qsel',qvec,'cache',false));
verifyEqual(tc,ae.Sigma_w,be.Sigma_w,'AbsTol',1e-12);
verifyEqual(tc,ae.chi_cc_q,be.chi_cc_q,'AbsTol',1e-11);
verifyEqual(tc,ae.explicit_q_dipole,be.explicit_q_dipole);
end

function test_split_realaxis_reduces_to_pm_path_at_m_zero(tc)
ion=invz_ion(); T=0.1;
lat=invzt_jq_tensor(ion,invzt_qgrid(4,'halfopen'), ...
    struct('dpRng',10,'cache',false));
[~,phase,di]=invzt_solve_auto(ion,T,[3 0 0],lat, ...
    struct('capture_attempts',true));
verifyEqual(tc,phase,1);
p=di.para.point;
ctx=invzt_a1_ordered_context(p.si,p.tl,T,lat, ...
    struct('dominant_count',p.mspec.ndom, ...
    'selector_source','same_point_pm_fixed_rank'));
q=invzt_a1_ordered_map(ctx,p.Sigma);
w=[0.005;0.02;0.05];
a=invzt_a1_ordered_realaxis_diagnostic(ion,ctx,q,w, ...
    struct('qsel','gamma_uniform'));
b=invzt_chi_realaxis(ion,T,[3 0 0],p,w, ...
    struct('qsel','gamma_uniform'));
verifyEqual(tc,p.tl.m,0,'AbsTol',1e-12);
verifyEqual(tc,a.Sigma_w,b.Sigma_w,'AbsTol',1e-12);
verifyEqual(tc,a.chi_uniform,b.chi_uniform,'AbsTol',1e-11);

ao=invzt_a1_ordered_realaxis_diagnostic(ion,ctx,q,w, ...
    struct('qsel','gamma_uniform','force_sigma0',true));
fullctx=invzt_a1_ordered_context(p.si,p.tl,T,lat, ...
    struct('dominant_count',numel(p.si.E), ...
    'selector_source','full_rank_bare_identity'));
fq=invzt_a1_ordered_map(fullctx,p.Sigma);
fr=invzt_a1_map_iterate(fullctx,p.Sigma,struct());
verifyTrue(tc,fr.converged);
af=invzt_a1_ordered_realaxis_diagnostic(ion,fullctx,fr.terminal,w, ...
    struct('qsel','gamma_uniform','force_sigma0',true));
verifyEqual(tc,ao.chi_tilde,af.chi_tilde,'AbsTol',1e-12);
verifyEqual(tc,ao.chi_uniform,af.chi_uniform,'AbsTol',1e-11);
verifyGreaterThan(tc,fq.residual_inf,0);
end

function test_relaxed_longitudinal_response_matches_single_ion_difference(tc)
ion=invz_ion(); T=0.1; B=[3 0 0];
lat=invzt_jq_tensor(ion,invzt_qgrid(4,'halfopen'), ...
    struct('dpRng',10,'cache',false));
[~,phase,di]=invzt_solve_auto(ion,T,B,lat,struct('capture_attempts',true));
verifyEqual(tc,phase,1);
si=di.ordered.point.si;
chi=real(invz_chi0z(si,T,0,struct('elastic',true)));
a=invzt_relaxed_longitudinal_response(chi,lat.info.Jaa0,1);
dh=1e-6;
base=struct('hyp',true,'Jxx0',lat.info.Jaa0, ...
    'transverse_mf','legacy_x');
op=base; op.hz_fixed=si.hz+dh;
om=base; om.hz_fixed=si.hz-dh;
sp=invz_single_ion(ion,T,B,op); sm=invz_single_ion(ion,T,B,om);
dm=(sp.Jexp-sm.Jexp)/(2*dh);
dhx=(sp.hx-sm.hx)/(2*dh);
verifyTrue(tc,sp.mf_converged&&sm.mf_converged);
verifyEqual(tc,a.moment_slope,dm,'RelTol',2e-6,'AbsTol',2e-8);
verifyEqual(tc,a.field_slope,dhx,'RelTol',2e-6,'AbsTol',2e-8);
verifyLessThanOrEqual(tc,a.closure_residual_inf,1e-14);
end

function test_field_differential_reduces_to_decoupled_scalar_ratio(tc)
Cbare=diag([0.3 0.2 2.0]); Ctilde=diag([0.25 0.15 1.25]);
J=diag([0.01 0.02 0.003]); q=[0;0;1];
r=invzt_field_differential_path(Cbare,Ctilde,J,q,struct());
verifyEqual(tc,r.longitudinal_integrand,2/1.25,'AbsTol',1e-14);
verifyEqual(tc,r.moment_slope,[0;0;2],'AbsTol',0);
verifyEqual(tc,r.transverse_compatibility_inf,0,'AbsTol',0);
verifyEqual(tc,r.equation_residual_inf,0,'AbsTol',1e-14);
verifyEqual(tc,r.active_rank,3);
verifyTrue(tc,r.fixed_applied_transverse_path.defined);
verifyEqual(tc,r.fixed_applied_transverse_path.qpath,q,'AbsTol',0);
verifyEqual(tc,r.fixed_applied_transverse_path.longitudinal_integrand, ...
    2/1.25,'AbsTol',1e-14);
end

function test_sigma0_reduction_reproduces_known_fixed_point(tc)
ion=invz_ion();lat=invzt_jq_tensor(ion,invzt_qgrid(4,'halfopen'), ...
    struct('dpRng',10,'cache',false));
[~,phase,di]=invzt_solve_auto(ion,0.1,[3 0 0],lat, ...
    struct('capture_attempts',true));
verifyEqual(tc,phase,1);p=di.ordered.point;
ctx=invzt_a1_ordered_context(p.si,p.tl,0.1,lat, ...
    struct('dominant_count',numel(p.si.E), ...
    'selector_source','full_rank_current_map_reference'));
r=invzt_a1_sigma0_reduced(ctx,p.Sigma(1),p.Sigma,struct());
verifyTrue(tc,r.tail_converged);
verifyEqual(tc,r.iterations,1);
verifyLessThanOrEqual(tc,abs(r.R0),p.outer_residual+1e-12);
verifyLessThanOrEqual(tc,r.tail_residual_inf,p.outer_residual+1e-12);

s=invzt_a1_sigma0_scan(ctx,p.Sigma,struct('output_mat',''));
verifyEqual(tc,s.certified_root_points,p.Sigma(1),'AbsTol',0);
verifyEmpty(tc,s.sign_change_intervals);
verifyEqual(tc,s.records.iterations,1);

c=invzt_a1_sigma0_tail_census(ctx,p.Sigma(1),[p.Sigma p.Sigma], ...
    ["known_root_1" "known_root_2"],struct('output_mat',''));
verifyEqual(tc,c.point_clusters.root_count,1);
verifyTrue(tc,all([c.runs.tail_converged]));
verifyEqual(tc,[c.runs.cluster],[1 1]);
end
