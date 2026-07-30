function result = invzp_dominant_manifold_census()
%INVZP_DOMINANT_MANIFOLD_CENSUS Test the fixed-rank electronuclear partition.
% This census does not dress either sector. It asks whether the framework's
% lowest-16 field-adapted manifold is isolated, smooth across the observed
% node-status boundary, and converged enough to justify a later bridge test.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'), ...
    fullfile(repo,'invz_tensor'));

source_path = fullfile(repo,'docs','diagnostics','invzp_outer_wp2', ...
    'wp2_low_field_m2_census.mat');
S = load(source_path);
ion = invz_ion();
T = 0.1;
opts = S.result.base_opts;
Bkeep = [0.5 1.0 1.8 2.2];
labels = ["predictor";"last_failed";"first_accepted";"high_endpoint"];
nrank = [16 24 32 48 64 96 136];
nrow = numel(Bkeep)*numel(labels);

B = nan(nrow,1);
sample = strings(nrow,1);
h = nan(nrow,1);
profile_accepted = false(nrow,1);
profile_status = strings(nrow,1);
M2 = nan(nrow,1);
electronic_inelastic_weight = nan(nrow,1);
full_static = nan(nrow,1);
full_inelastic = nan(nrow,1);
dom16_static = nan(nrow,1);
dom16_inelastic = nan(nrow,1);
dom16_static_share = nan(nrow,1);
dom16_inelastic_share = nan(nrow,1);
electronic_to_dom16_inelastic = nan(nrow,1);
elastic_convention_share = nan(nrow,1);
p_mass = nan(nrow,numel(nrank));
chi_share = nan(nrow,numel(nrank));
var_share = nan(nrow,numel(nrank));
edge_gap = nan(nrow,numel(nrank));
gap_kBT = nan(nrow,numel(nrank));
previous_overlap16 = nan(nrow,1);
previous_projector_distance16 = nan(nrow,1);

row = 0;
for ib = 1:numel(Bkeep)
    is = find(abs(S.result.summary.B-Bkeep(ib)) < 1e-12,1);
    if isempty(is), error('invzp:representationField','Missing B=%.3g T.',Bkeep(ib)); end
    t = S.result.details{is};
    ilastfail = find(~t.accepted,1,'last');
    ifirstok = find(t.accepted,1);
    indices = [0 ilastfail ifirstok height(t)];
    vb_previous = [];

    for ic = 1:numel(labels)
        row = row+1;
        B(row) = Bkeep(ib);
        sample(row) = labels(ic);
        if indices(ic) == 0
            h(row) = 0;
            profile_accepted(row) = false;
            profile_status(row) = S.result.summary.predictor_status(is);
        else
            j = indices(ic);
            h(row) = t.h(j);
            profile_accepted(row) = t.accepted(j);
            profile_status(row) = t.status(j);
        end

        si = invz_single_ion(ion,T,[B(row) 0 0], ...
            struct('hyp',true,'Jxx0',opts.Jxx0,'hz_fixed',h(row)));
        tl = invz_twolevel_ordered(ion,T,B(row),h(row), ...
            struct('Jxx0',opts.Jxx0));
        M2(row) = tl.M2;
        electronic_inelastic_weight(row) = tl.M2*tl.g0;

        [cdom,crest,mspec] = invzt_chi0_split(si,T,0, ...
            struct('dominant_count',16,'elastic',true));
        [cdom_i,crest_i] = invzt_chi0_split(si,T,0, ...
            struct('dominant_count',16,'elastic',false));
        full_static(row) = real(cdom(3,3,1)+crest(3,3,1));
        full_inelastic(row) = real(cdom_i(3,3,1)+crest_i(3,3,1));
        dom16_static(row) = real(cdom(3,3,1));
        dom16_inelastic(row) = real(cdom_i(3,3,1));
        dom16_static_share(row) = dom16_static(row)/full_static(row);
        dom16_inelastic_share(row) = dom16_inelastic(row)/full_inelastic(row);
        electronic_to_dom16_inelastic(row) = ...
            electronic_inelastic_weight(row)/dom16_inelastic(row);
        elastic_convention_share(row) = mspec.elastic_conv_share;

        for ir = 1:numel(nrank)
            vo = struct('n_vertex',nrank(ir));
            if ir == 1 && ~isempty(vb_previous), vo.vb_prev = vb_previous; end
            vb = invzt_ordered_vertex_basis(ion,T,si,vo);
            p_mass(row,ir) = vb.p_mass;
            chi_share(row,ir) = vb.chi_share;
            var_share(row,ir) = vb.var_share;
            edge_gap(row,ir) = vb.gap_16_17;
            gap_kBT(row,ir) = vb.gap_kBT;
            if ir == 1
                previous_overlap16(row) = vb.min_subspace_overlap;
                previous_projector_distance16(row) = vb.projector_distance;
                vb_previous = vb;
            end
        end
        fprintf('B=%.1f %-14s h=%.7g: chi16 %.4f var16 %.4f\n', ...
            B(row),sample(row),h(row),chi_share(row,1),var_share(row,1));
    end
end

min_rank_chi90 = nan(nrow,1);
min_rank_var90 = nan(nrow,1);
for k = 1:nrow
    jchi = find(chi_share(k,:) >= 0.9,1);
    jvar = find(var_share(k,:) >= 0.9,1);
    if ~isempty(jchi), min_rank_chi90(k) = nrank(jchi); end
    if ~isempty(jvar), min_rank_var90(k) = nrank(jvar); end
end
tab = table(B,sample,h,profile_accepted,profile_status,M2, ...
    electronic_inelastic_weight,full_static,full_inelastic,dom16_static, ...
    dom16_inelastic,dom16_static_share,dom16_inelastic_share, ...
    electronic_to_dom16_inelastic,elastic_convention_share, ...
    p_mass(:,1),chi_share(:,1),var_share(:,1),edge_gap(:,1),gap_kBT(:,1), ...
    previous_overlap16,previous_projector_distance16, ...
    chi_share(:,2),var_share(:,2),chi_share(:,3),var_share(:,3), ...
    min_rank_chi90,min_rank_var90, ...
    'VariableNames',{'B','sample','h','profile_accepted','profile_status','M2', ...
    'electronic_inelastic_weight','full_static','full_inelastic', ...
    'dom16_static','dom16_inelastic','dom16_static_share', ...
    'dom16_inelastic_share','electronic_to_dom16_inelastic', ...
    'elastic_convention_share','p_mass16','chi_share16','var_share16', ...
    'edge_gap16','gap_kBT16','previous_overlap16', ...
    'previous_projector_distance16','chi_share24','var_share24', ...
    'chi_share32','var_share32','min_rank_chi90','min_rank_var90'});
rank_diagnostics = struct('p_mass',p_mass,'chi_share',chi_share, ...
    'var_share',var_share,'edge_gap',edge_gap,'gap_kBT',gap_kBT);
result = struct('table',tab,'ranks',nrank, ...
    'rank_diagnostics',rank_diagnostics,'source',source_path, ...
    'provenance',struct('commit',git_commit(repo)), ...
    'note',['No dressing is performed. invzt_chi0_split gives an exact ' ...
    'full=dominant+rest response partition with the locked elastic convention; ' ...
    'invzt_ordered_vertex_basis normalizes the selected subspace for vertex ' ...
    'coverage diagnostics. Their static shares can differ slightly if p_mass<1.']);
save(fullfile(here,'wp3_dominant_manifold_census.mat'),'result','-v7');
disp(tab(:,{'B','sample','h','profile_accepted','M2','dom16_inelastic_share', ...
    'electronic_to_dom16_inelastic','chi_share16','var_share16', ...
    'chi_share24','var_share24','chi_share32','var_share32', ...
    'min_rank_chi90','min_rank_var90'}));
end

function commit = git_commit(repo)
[status,text] = system(sprintf('git -C "%s" rev-parse HEAD',repo));
if status == 0, commit = strtrim(text); else, commit = 'unknown'; end
end
