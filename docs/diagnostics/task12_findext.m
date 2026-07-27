root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root); addpath(fullfile(root,'invz_projected')); addpath(fullfile(root,'invz_common'));
ion = invz_ion(); T = 0.31; Jnu = linspace(-2e-3,6.0e-3,24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
for B = [3.00 3.05 3.10 3.15]
    [h_ref, p_ref] = invz_hmf_ordered(ion, T, [B 0 0], Jnu, o);
    hmax_ref = max(p_ref.hgrid);
    frac = min(0.5, 4*h_ref/hmax_ref);
    od = o; od.hmin_frac = frac; od.nH = 17;
    [h, p] = invz_hmf_ordered(ion, T, [B 0 0], Jnu, od);
    fprintf('B=%.2f h_ref=%.6g hmax=%.6g frac=%.6g -> h=%.6g status=%s n_extend=%d redens=%d hmin_init=%.4g\n', ...
        B, h_ref, hmax_ref, frac, h, p.status, p.n_extend, p.redensified, p.hmin_initial);
end
