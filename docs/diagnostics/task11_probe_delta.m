ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(genpath(fullfile(ROOT, 'invz_projected')));
addpath(fullfile(ROOT, 'invz_common'));

ion = invz_ion();
T = 0.31;
Bx = [0 0 0];
hz = 0;
Jxx0 = ion.Jxx0;
tmf = 'legacy_x';
si = invz_single_ion(ion, T, invz_field_vec(Bx), struct('hyp', false, 'hz_fixed', hz, 'Jxx0', Jxx0, 'transverse_mf', tmf));
Delta = si.E(2) - si.E(1);
fprintf('PROBE Delta (Bx=[0 0 0], hz=0) = %.15e meV\n', Delta);
fprintf('PROBE Delta < 1e-4 ? %d\n', Delta < 1e-4);
