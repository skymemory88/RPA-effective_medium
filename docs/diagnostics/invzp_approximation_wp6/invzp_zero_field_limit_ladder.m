function result = invzp_zero_field_limit_ladder()
%INVZP_ZERO_FIELD_LIMIT_LADDER Positive-field and basis-invariant zero audit.
% Diagnostic only. The scalar factor-one missing-area solve uses the current
% production interaction, node count, damping, iteration budget, and real-axis
% grid. The raw single-ion table is evaluated independently of the two-level
% domain floor to distinguish a subspace limit from solver availability.

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));
F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
ion = invz_ion();
C = invz_const();
T = 0.1;
J = F.J;
J0eff = F.info.Jcc0;
Jxx0 = F.info.Jaa0;
fields = unique([0 logspace(-6,-2,5) 0.02 0.03 0.04 0.042 0.043 ...
    0.0435 0.044 0.045 0.05 0.06 0.075 0.09 0.12]);
w_GHz = (0:0.01:6).';
w = w_GHz*C.Gh2mV;
eta = 5e-5;

electronic = electronic_subspace_table(ion,T,fields,Jxx0);
full_ion = full_ion_table(ion,T,fields,Jxx0);

opts = struct('J0eff',J0eff,'Jxx0',Jxx0,'mix_outer',0.35, ...
    'max_outer',2000,'nH',129, ...
    'hmf_integral_mode','missing_area_approx', ...
    'hmf_missing_area_factor',1, ...
    'hmf_approx_branch', ...
    'picard_attracting_contiguous_high_h_component');
points = cell(size(fields));
errors = strings(size(fields));
parfor k = 1:numel(fields)
    try
        points{k} = invz_solve_point_ordered(ion,T,fields(k),J,opts);
    catch err
        if ~strncmp(err.identifier,'invz:',5), rethrow(err); end
        errors(k) = string(err.identifier);
        points{k} = [];
    end
end

spectra = nan(numel(w),numel(fields));
records = repmat(empty_record(),numel(fields),1);
for k = 1:numel(fields)
    records(k) = summarize(fields(k),points{k},errors(k),J,J0eff);
    if records(k).converged
        copts = struct('Jsel',J0eff,'eta',eta,'Jxx0',Jxx0, ...
            'Jshape',0,'hyp',true,'si',points{k}.si);
        try
            q = invz_chi_realaxis(ion,T,fields(k),points{k},w,copts);
            spectra(:,k) = imag(q.chi_cc_q(1,:)).';
            records(k).spectrum_finite = all(isfinite(spectra(:,k)));
            records(k).peak_meV = invz_peak_energy(spectra(:,k),w,0);
        catch err
            if ~strncmp(err.identifier,'invz:',5), rethrow(err); end
            records(k).spectrum_error = string(err.identifier);
        end
    end
end
summary = struct2table(records);

domain_floor = 1e-4;
above = find(electronic.Delta_meV >= domain_floor,1);
first_domain_field = NaN;
if ~isempty(above), first_domain_field = fields(above); end
valid_solve = summary.converged;
first_converged_field = NaN;
if any(valid_solve), first_converged_field = fields(find(valid_solve,1)); end

assert(all(isfinite(electronic.projector_overlap_B0)) && ...
    all(electronic.projector_min_singular_B0 > 0), ...
    'The electronic lowest-doublet projector comparison is undefined.');
assert(all(electronic.gap23_meV > 0.01), ...
    'The electronic lowest doublet is not spectrally isolated.');
assert(all(isfinite(full_ion.projector_overlap_B0)) && ...
    all(full_ion.projector_min_singular_B0 > 0), ...
    'The full-ion lowest-doublet projector comparison is undefined.');
assert(all(summary.final_resid(valid_solve) < 1e-8) && ...
    all(summary.D_uni(valid_solve) > 0) && ...
    all(summary.supremum_mass(valid_solve) > 0) && ...
    all(summary.min_mesh_x_mass(valid_solve) > 0) && ...
    all(summary.min_mesh_medium_mass(valid_solve) > 0), ...
    'Every accepted positive-field point must pass the retained gates.');

result = struct('T',T,'fields',fields,'nH',opts.nH, ...
    'area_factor',1,'mix_outer',opts.mix_outer, ...
    'max_outer',opts.max_outer,'twolevel_domain_floor_meV',domain_floor, ...
    'first_sample_above_domain_floor_T',first_domain_field, ...
    'first_converged_sample_T',first_converged_field, ...
    'electronic_subspace',electronic,'full_ion',full_ion, ...
    'summary',summary,'w_meV',w,'w_GHz',w_GHz, ...
    'eta_meV',eta,'spectra',spectra, ...
    'production_changed',false);
save(fullfile(here,'wp6_zero_field_limit_ladder.mat'),'result','-v7');
disp(electronic(:,{'B_T','Delta_meV','gap23_meV','m_basis', ...
    'M2_offdiag','half_trace_PJzPJz','static_weight','projector_overlap_B0'}));
disp(full_ion(:,{'B_T','Delta_meV','gap23_meV','Jz_variance', ...
    'chi_path','projector_overlap_B0'}));
disp(summary);
end

function tab = electronic_subspace_table(ion,T,fields,Jxx0)
C = invz_const();
oJ = stevens_ops(ion.J);
n = numel(fields);
Delta = nan(n,1);
gap23 = nan(n,1);
m = nan(n,1);
M2 = nan(n,1);
halftrace = nan(n,1);
static_weight = nan(n,1);
overlap = nan(n,1);
minsv = nan(n,1);
projectors = cell(n,1);
bases = cell(n,1);
for k = 1:n
    si = invz_single_ion(ion,T,[fields(k) 0 0], ...
        struct('hyp',false,'hz_fixed',0,'Jxx0',Jxx0));
    V2 = si.V(:,1:2);
    M = V2'*oJ.Jz*V2;
    Delta(k) = si.E(2)-si.E(1);
    gap23(k) = si.E(3)-si.E(2);
    m(k) = real(M(1,1));
    M2(k) = abs(M(1,2))^2;
    halftrace(k) = real(trace(M*M))/2;
    n01 = tanh(Delta(k)/(2*C.kB*T));
    g0 = 2*n01/Delta(k);
    h0 = (1-n01^2)/(C.kB*T);
    static_weight(k) = M2(k)*g0+m(k)^2*h0;
    projectors{k} = V2*V2';
    bases{k} = V2;
end
for k = 1:n
    overlap(k) = real(trace(projectors{1}*projectors{k}))/2;
    minsv(k) = min(svd(bases{1}'*bases{k}));
end
tab = table(fields(:),Delta,gap23,m,M2,halftrace,static_weight, ...
    overlap,minsv,'VariableNames',{'B_T','Delta_meV','gap23_meV', ...
    'm_basis','M2_offdiag','half_trace_PJzPJz','static_weight', ...
    'projector_overlap_B0','projector_min_singular_B0'});
end

function tab = full_ion_table(ion,T,fields,Jxx0)
n = numel(fields);
Delta = nan(n,1);
gap23 = nan(n,1);
variance = nan(n,1);
chi_path = nan(n,1);
overlap = nan(n,1);
minsv = nan(n,1);
bases = cell(n,1);
for k = 1:n
    si = invz_single_ion(ion,T,[fields(k) 0 0], ...
        struct('hyp',true,'hz_fixed',0,'Jxx0',Jxx0));
    Delta(k) = si.E(2)-si.E(1);
    gap23(k) = si.E(3)-si.E(2);
    variance(k) = si.JzJz_fluct;
    c0 = invz_chi0z(si,T,0,struct('elastic',true));
    X = real(c0(:,:,1));
    feedback = X(3,1)*(Jxx0/(1-Jxx0*X(1,1)))*X(1,3);
    chi_path(k) = X(3,3)+feedback;
    bases{k} = si.V(:,1:2);
end
for k = 1:n
    overlap(k) = norm(bases{1}'*bases{k},'fro')^2/2;
    minsv(k) = min(svd(bases{1}'*bases{k}));
end
tab = table(fields(:),Delta,gap23,variance,chi_path,overlap,minsv, ...
    'VariableNames',{'B_T','Delta_meV','gap23_meV','Jz_variance', ...
    'chi_path','projector_overlap_B0','projector_min_singular_B0'});
end

function out = summarize(B,point,error_id,J,J0eff)
out = empty_record();
out.B = B;
out.exception_id = error_id;
if isempty(point)
    out.status = "exception";
    return;
end
out.status = string(point.hmf_status);
out.converged = point.converged;
out.hmf = point.hmf;
out.Sigma0 = point.Sigma0;
out.D_uni = point.D_uni;
out.final_resid = point.final_resid;
p = point.hmf_prof;
out.predictor_converged = p.predictor_converged;
if isfield(p,'missing_area_integral') && ...
        isfield(p.missing_area_integral,'node_count')
    out.component_nodes = p.missing_area_integral.node_count;
    out.component_edge = p.missing_area_integral.component_edge;
end
if ~point.converged
    return;
end
out.final_Delta = point.tl.Delta;
out.final_m = point.tl.m;
out.final_M2 = point.tl.M2;
out.final_n01 = point.tl.n01;
Gstat = p.Gstat_star;
K0 = point.K(1);
Gtilde = Gstat/(1-K0*Gstat);
out.Fprime = p.r_star*(1+J0eff*Gtilde);
out.supremum_mass = 1+J0eff*Gtilde;
out.min_mesh_x_mass = min(1+J(:)*Gtilde);
out.min_mesh_medium_mass = min(1+(J(:)-K0)*Gstat);
end

function out = empty_record()
out = struct('B',NaN,'status',"",'exception_id',"", ...
    'converged',false,'predictor_converged',false,'hmf',NaN, ...
    'Sigma0',NaN,'D_uni',NaN,'Fprime',NaN,'final_resid',NaN, ...
    'component_nodes',0,'component_edge',NaN,'final_Delta',NaN, ...
    'final_m',NaN,'final_M2',NaN,'final_n01',NaN, ...
    'supremum_mass',NaN,'min_mesh_x_mass',NaN, ...
    'min_mesh_medium_mass',NaN,'spectrum_finite',false, ...
    'spectrum_error',"",'peak_meV',NaN);
end
