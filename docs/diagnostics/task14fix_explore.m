ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT, 'invz_common'));
addpath(genpath(fullfile(ROOT, 'invz_projected')));
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o   = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true);

cfg = {};
cfg{end+1} = {'dyson', struct('static_medium','strict_1z_dyson_ref')};
cfg{end+1} = {'bare_ref', struct('static_medium','strict_1z_bare_ref')};
cfg{end+1} = {'dyson_m2', struct('static_medium','strict_1z_dyson_ref','ref_margin',2)};
cfg{end+1} = {'bare_ref_m2', struct('static_medium','strict_1z_bare_ref','ref_margin',2)};
cfg{end+1} = {'dyson_nohyp', struct('static_medium','strict_1z_dyson_ref','hyp',false)};
cfg{end+1} = {'bare_ref_nohyp', struct('static_medium','strict_1z_bare_ref','hyp',false)};

for k = 1:numel(cfg)
    nm = cfg{k}{1};  ex = cfg{k}{2};
    oo = o;  f = fieldnames(ex);
    for i = 1:numel(f), oo.(f{i}) = ex.(f{i}); end
    for Bv = [0.5 2.0 2.85]
        try
            p = invz_solve_point_ordered(ion, 0.31, [Bv 0 0], Jnu, oo);
            fprintf('%-16s B=%.2f ord=%d conv=%d iters=%3d status=%-16s denom=%.6g\n', ...
                nm, Bv, p.is_ordered, p.converged, p.outer_iters, p.medium_status, p.medium_denom);
        catch ME
            fprintf('%-16s B=%.2f ERROR %s\n', nm, Bv, ME.identifier);
        end
    end
end
fprintf('EXPLORE DONE\n');
