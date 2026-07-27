% Task 17 follow-up: G11 trace-node detail table only (fixes an fprintf type-mismatch bug in
% the first characterisation pass -- %.6g was fed a num2str'd STRING for nd.h instead of the
% numeric value, which desynced the whole row via fprintf's char-array flattening). Histograms
% from the first pass already stand; this just adds per-node h-locations.
REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(fullfile(REPO, 'invz_projected'));
addpath(REPO);
addpath(fullfile(REPO, 'invz_common'));

ion = invz_ion();
prov = struct('grid', [16 16 16], 'dpRng', 30, 'dipole', 'bruteforce', 'cache', false);
[Jnu11, info11] = invz_bz_couplings(ion, prov);
mom11 = invz_coupling_moments(Jnu11(:));
o11 = struct('J0eff', info11.Jcc0, 'Jxx0', ion.Jxx0, 'hyp', true, ...
             'static_medium', 'strict_1z_dyson_ref', 'Jmom', mom11, 'trace', true);
[hstar11, p11, trc11] = invz_hmf_ordered(ion, 0.1, [1 0 0], Jnu11(:), o11);
fprintf('p.status=%s hstar=%.6g crit_star=%.6g  numel(trc.nodes)=%d\n', ...
    p11.status, hstar11, p11.crit_star, numel(trc11.nodes));
fprintf('fieldnames(trc.nodes): %s\n', strjoin(fieldnames(trc11.nodes), ', '));
has_m = isfield(trc11.nodes, 'm');
fprintf('\n%4s  %4s  %14s  %8s  %10s  %-24s  %s\n', 'idx', 'id', 'h', 'accepted', 'm', 'medium_status', 'term_reason');
for i = 1:numel(trc11.nodes)
    nd = trc11.nodes(i);
    if has_m, mval = nd.m; else, mval = NaN; end
    line = sprintf('%4d  %4d  %14.8g  %8d  %10.6g  %-24s  %s', ...
        i, nd.id, nd.h, double(nd.accepted), mval, nd.medium_status, nd.term_reason);
    fprintf('%s\n', line);
end
fprintf('\n-- non-ok nodes only --\n');
for i = 1:numel(trc11.nodes)
    nd = trc11.nodes(i);
    if has_m, mval = nd.m; else, mval = NaN; end
    if ~strcmp(nd.medium_status, 'ok')
        fprintf('  idx=%3d id=%4d h=%.8g m=%.6g medium_status=%s term_reason=%s\n', ...
            i, nd.id, nd.h, mval, nd.medium_status, nd.term_reason);
    end
end
fprintf('\n### done ###\n');
