REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(genpath(fullfile(REPO, 'invz_projected')));
addpath(fullfile(REPO, 'invz_common'));
ion = invz_ion();
prov = struct('grid', [16 16 16], 'dpRng', 30, 'dipole', 'bruteforce', 'cache', false);
[Jnu, info] = invz_bz_couplings(ion, prov);
fprintf('TEST B (my genpath-style addpath, fresh, no part-1): OK, numel(Jnu)=%d backend=%s\n', numel(Jnu), info.dipole.backend);
