% Task 18 Step 4: G7 scheme-jump measurement (brief verbatim; prereg SS6, non-gating). Put in a
% script because the command is too long for one -batch argument (controller notes trap #9).
% Absolute addpath below: `run` executes with THIS script's own directory as its context, not
% the repo root (this is what bit Task 16's implementer).
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
addpath(genpath(fullfile(ROOT, 'invz_projected')));
addpath(fullfile(ROOT, 'invz_common'));
addpath(ROOT);

ion = invz_ion();
[Jnu, ~] = invz_bz_couplings(ion, struct('grid', [16 16 16], 'dpRng', 30, ...
    'dipole', 'bruteforce', 'cache', false));
Jf = Jnu(:);  mom = invz_coupling_moments(Jf);
for T = [0.05 0.1 0.31 1.0]
    [wn, ~, ~] = invz_matsubara(T, 40);
    si = invz_single_ion(ion, T, [6 0 0], struct('hyp', true, 'Jxx0', ion.Jxx0, ...
        'transverse_mf', 'legacy_x'));
    c0 = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
    G0 = -real(squeeze(c0(3,3,:)));  S = zeros(size(wn));
    a = invz_emt_scalar(G0, S, Jf, struct());
    b = invz_emt_scalar(G0, S, Jf, struct('static_medium', 'strict_1z_dyson_ref', 'Jmom', mom));
    fprintf('T=%-5.2f jump|K1s-K1r|=%.4g  dispersion|K2-K1|=%.4g  ratio=%.4g\n', T, ...
        abs(b.K(1)-a.K(1)), abs(a.K(2)-a.K(1)), abs(b.K(1)-a.K(1))/abs(a.K(2)-a.K(1)));
end
