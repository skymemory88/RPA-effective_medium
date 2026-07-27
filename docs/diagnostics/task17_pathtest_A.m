REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
here = fullfile(REPO, 'invz_projected', 'tests');
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
ion = invz_ion();
prov = struct('grid', [16 16 16], 'dpRng', 30, 'dipole', 'bruteforce', 'cache', false);
[Jnu, info] = invz_bz_couplings(ion, prov);
fprintf('TEST A (setupOnce-style flat addpath, fresh, no part-1): OK, numel(Jnu)=%d backend=%s\n', numel(Jnu), info.dipole.backend);
