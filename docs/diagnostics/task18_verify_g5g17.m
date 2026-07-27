ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
SDD  = fullfile(ROOT, '.superpowers', 'sdd');
addpath(fullfile(ROOT, 'invz_projected'));
addpath(fullfile(ROOT, 'invz_common'));
addpath(ROOT);

ion = invz_ion();
L = load(fullfile(SDD, 'task18_couplings.mat'));
rep = invz_gate0_report(ion, 0.10, 2.0, [], L.Jnu_flat, struct('J0eff', L.info.Jcc0));

fprintf('g5.pass_Sigma0=%d pass_r_minus_1=%d\n', rep.g5.pass_Sigma0, rep.g5.pass_r_minus_1);
fprintf('g5.d_65_129_Sigma0=%.6g tol=%.6g\n', rep.g5.d_65_129_Sigma0, rep.g5.tol_65_129_Sigma0);
fprintf('g5.d_65_129_r_minus_1=%.6g tol=%.6g\n', rep.g5.d_65_129_r_minus_1, rep.g5.tol_65_129_r_minus_1);
fprintf('g17.pass=%d\n', rep.g17.pass);
fprintf('g17.g17_r_65.status=%s max_jump=%.6g\n', rep.g17.g17_r_65.status, rep.g17.g17_r_65.max_jump);
fprintf('g17.g17_r_129.status=%s max_jump=%.6g\n', rep.g17.g17_r_129.status, rep.g17.g17_r_129.max_jump);
fprintf('g17.g17_crit_65.status=%s max_jump=%.6g\n', rep.g17.g17_crit_65.status, rep.g17.g17_crit_65.max_jump);
fprintf('g17.g17_crit_129.status=%s max_jump=%.6g\n', rep.g17.g17_crit_129.status, rep.g17.g17_crit_129.max_jump);

save(fullfile(SDD, 'task18_chunk_B2p0.mat'), 'rep');
fprintf('RE-SAVED task18_chunk_B2p0.mat with g5/g17 fields\n');
