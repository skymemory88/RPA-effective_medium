function result = invzp_low_field_m2_census()
%INVZP_LOW_FIELD_M2_CENSUS Relate low-field node status to projected weight.
% This is an observational diagnostic: it calls the strict production point
% solver and reconstructs scalar kinematic/static diagnostics from its
% exported profile. It does not alter the ordered self-energy or acceptance.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

baseline_path = fullfile(here,'wp2_hmf_node_transaction_census.mat');
fixture_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat');
B = [0.5 0.8 1.0 1.2 1.5 1.8 2.0 2.2].';
N = load(baseline_path);
F = load(fixture_path);
opts = N.result.base_opts;
opts.nH = 33;
ion = invz_ion();
T = F.provenance.T;

nB = numel(B);
elapsed_seconds = nan(nB,1);
profile_status = strings(nB,1);
predictor_status = strings(nB,1);
n_nodes = nan(nB,1);
n_accepted = nan(nB,1);
first_accepted_h = nan(nB,1);
M2_first_accepted = nan(nB,1);
min_M2 = nan(nB,1);
accepted_min_M2 = nan(nB,1);
failed_max_M2 = nan(nB,1);
max_m2_over_M2 = nan(nB,1);
dominant_share_at_first = nan(nB,1);
projected_inelastic_ratio_at_first = nan(nB,1);
predictor_M2 = nan(nB,1);
predictor_projected_inelastic_ratio = nan(nB,1);
max_prefactor_product_relerr = nan(nB,1);
details = cell(nB,1);

for ib = 1:nB
    tic;
    pt = invz_solve_point_ordered(ion,T,B(ib),F.J,opts);
    elapsed_seconds(ib) = toc;
    p = pt.hmf_prof;
    h = p.hgrid(:);
    n = numel(h);
    M2 = nan(n,1);
    m_electronic = nan(n,1);
    n01 = nan(n,1);
    Delta = nan(n,1);
    g0 = nan(n,1);
    G0inel = nan(n,1);

    for ih = 1:n
        tl = invz_twolevel_ordered(ion,T,B(ib),h(ih), ...
            struct('Jxx0',opts.Jxx0));
        si = invz_single_ion(ion,T,[B(ib) 0 0], ...
            struct('hyp',true,'Jxx0',opts.Jxx0,'hz_fixed',h(ih)));
        c0i = invz_chi0z(si,T,0,struct('elastic',false));
        M2(ih) = tl.M2;
        m_electronic(ih) = tl.m;
        n01(ih) = tl.n01;
        Delta(ih) = tl.Delta;
        g0(ih) = tl.g0;
        G0inel(ih) = -real(c0i(3,3,1));
    end

    Sigma0 = p.Sigma0(:);
    K0 = p.K0(:);
    Gstat = p.Gstat(:);
    G0bare = p.G0bare(:);
    G0el = G0bare-G0inel;
    U = G0inel./(1+Sigma0+K0.*G0inel);
    V = Gstat-U;
    M2_over_n01sq = M2./n01.^2;
    m2_over_M2 = m_electronic.^2./M2;
    projected_static_weight = M2.*g0;
    full_path_static_weight = -G0bare;
    dominant_share = projected_static_weight./full_path_static_weight;
    projected_inelastic_ratio = projected_static_weight./(-G0inel);

    tl0 = invz_twolevel_ordered(ion,T,B(ib),0, ...
        struct('Jxx0',opts.Jxx0));
    si0 = invz_single_ion(ion,T,[B(ib) 0 0], ...
        struct('hyp',true,'Jxx0',opts.Jxx0,'hz_fixed',0));
    c0i0 = invz_chi0z(si0,T,0,struct('elastic',false));
    predictor_M2(ib) = tl0.M2;
    predictor_projected_inelastic_ratio(ib) = ...
        (tl0.M2*tl0.g0)/real(c0i0(3,3,1));

    % The M2-sensitive part of (2m^2/M2)*gamma0 is the product below.
    % Its exact reassociation is 2m^2/n01^2, before multiplying by the
    % common bracket [lambda1-(1-n01^2)K0].
    naive_prefactor_product = 2*m2_over_M2.*M2_over_n01sq;
    stable_prefactor_product = 2*m_electronic.^2./n01.^2;
    prefactor_product_relerr = abs(naive_prefactor_product- ...
        stable_prefactor_product)./max(abs(stable_prefactor_product),realmin);

    accepted = p.node_conv(:);
    first = find(accepted,1);
    profile_status(ib) = string(p.status);
    predictor_status(ib) = string(p.predictor_static_status);
    n_nodes(ib) = n;
    n_accepted(ib) = nnz(accepted);
    min_M2(ib) = min(M2);
    max_m2_over_M2(ib) = max(m2_over_M2);
    max_prefactor_product_relerr(ib) = max(prefactor_product_relerr);
    if ~isempty(first)
        first_accepted_h(ib) = h(first);
        M2_first_accepted(ib) = M2(first);
        dominant_share_at_first(ib) = dominant_share(first);
        projected_inelastic_ratio_at_first(ib) = ...
            projected_inelastic_ratio(first);
        accepted_min_M2(ib) = min(M2(accepted));
    end
    if any(~accepted)
        failed_max_M2(ib) = max(M2(~accepted));
    end

    details{ib} = table(h,accepted,p.static_status(:),Delta,M2, ...
        m_electronic,n01,g0,m2_over_M2,M2_over_n01sq, ...
        projected_static_weight,full_path_static_weight,dominant_share, ...
        projected_inelastic_ratio, ...
        G0inel,G0el,U,V,Sigma0,K0,Gstat,p.D_uni(:), ...
        prefactor_product_relerr, ...
        'VariableNames',{'h','accepted','status','Delta','M2', ...
        'm_electronic','n01','g0','m2_over_M2','M2_over_n01sq', ...
        'projected_static_weight','full_path_static_weight','dominant_share', ...
        'projected_inelastic_ratio', ...
        'G0inel','G0el','U','V','Sigma0','K0','Gstat','D_uni', ...
        'prefactor_product_relerr'});
    fprintf('B=%.1f T: %d/%d accepted, first h=%.7g, %.2f s\n', ...
        B(ib),n_accepted(ib),n,first_accepted_h(ib),elapsed_seconds(ib));
end

summary = table(B,elapsed_seconds,profile_status,predictor_status,n_nodes, ...
    n_accepted,first_accepted_h,M2_first_accepted,min_M2, ...
    accepted_min_M2,failed_max_M2,max_m2_over_M2, ...
    dominant_share_at_first,projected_inelastic_ratio_at_first, ...
    predictor_M2,predictor_projected_inelastic_ratio, ...
    max_prefactor_product_relerr);
result = struct('summary',summary,'details',{details},'base_opts',opts, ...
    'provenance',struct('transaction_census',baseline_path, ...
    'fixture',fixture_path,'commit',git_commit(repo)), ...
    'note',['Strict-profile observational census. projected_static_weight=' ...
    'M2*g0 is a bare two-level static weight, not an exported real-axis ' ...
    'susceptibility. dominant_share uses the full path response, whereas ' ...
    'projected_inelastic_ratio compares with the full electronuclear ' ...
    'inelastic spectral response. U and V are reconstructed from the exact ' ...
    'production ' ...
    'split and exported final candidate values; interpret them only where ' ...
    'accepted=true.']);
save(fullfile(here,'wp2_low_field_m2_census.mat'),'result','-v7');
disp(summary);
end

function commit = git_commit(repo)
[status,text] = system(sprintf('git -C "%s" rev-parse HEAD',repo));
if status == 0, commit = strtrim(text); else, commit = 'unknown'; end
end
