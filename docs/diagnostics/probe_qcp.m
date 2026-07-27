root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(fullfile(root, 'invz_projected'));  addpath(fullfile(root, 'invz_common'));  addpath(root);
ion = invz_ion(); T = 0.31;
info = struct('Jcc0', 6.4e-3); Jnu = linspace(-2e-3, 6.0e-3, 24).';
sopts = struct('hyp', true, 'J0eff', info.Jcc0, 'Jxx0', ion.Jxx0, 'bz_tol', 1e-9);
for B = 2.5:0.05:5.5
    [pt, ph] = invz_solve_auto(ion, T, [B 0 0], Jnu, sopts);
    c = NaN; cv = false;
    try
        ptp = invz_solve_point(ion, T, [B 0 0], Jnu, sopts); c = ptp.crit; cv = ptp.converged;
    catch err
        if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    end
    fprintf('B=%.2f  auto_phase=%d  pm_conv=%d  pm_crit=%.4g\n', B, ph, cv, c);
end
