ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT); addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff',6.42e-3,'Jxx0',ion.Jxx0,'hyp',true,'static_medium','strict_1z_dyson_ref');
G = cell(1,3); nHs = [33 65 129];
for k = 1:3
    ok_ = o; ok_.nH = nHs(k);
    [~, p] = invz_hmf_ordered(ion, 0.31, [2.85 0 0], Jnu, ok_);
    G{k} = p.hgrid;
    fprintf('nH=%3d  status=%-22s n_extend=%d  redensified=%d  numel(hgrid)=%d\n', ...
        nHs(k), p.status, p.n_extend, p.redensified, numel(p.hgrid));
end
fprintf('\n-- nesting of the ACTUAL adapted grids --\n');
for a = 1:2
    for b = a+1:3
        coarse = G{a}; fine = G{b};
        hit = 0;
        for i = 1:numel(coarse)
            if any(abs(fine - coarse(i)) <= 1e-14*max(1,abs(coarse(i)))), hit = hit + 1; end
        end
        fprintf('  nH=%3d nodes found in nH=%3d grid: %d / %d  -> %s\n', ...
            nHs(a), nHs(b), hit, numel(coarse), ...
            string(hit == numel(coarse)));
    end
end
fprintf('\n### done ###\n');
