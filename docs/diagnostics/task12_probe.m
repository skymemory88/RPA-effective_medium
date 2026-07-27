% Task 12 fixture-physics probe (controller, pre-dispatch).
% The diagnostics test makes several UNVERIFIED numeric/physics claims about the fixture:
%   (a) prof.status == 'ok' and hstar finite            (test_path_integrals..., test_crit_star...)
%   (b) crit_star = r_star + J0eff*G0bare_star > 0      (test_crit_star_at_the_root)
%   (c) |r-1 - Sigma0| > 1e-9 at a finite-m node        (test_r_minus_1_is_not_sigma0...)
% Strict mode is not wired into invz_hmf_ordered until Task 12 itself, so measure the RESUMMED
% path at the same (ion,T,Bx,Jnu) to establish whether the profile converges and the physics
% claims are even reachable.
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root);
addpath(fullfile(root, 'invz_projected'));
addpath(fullfile(root, 'invz_common'));

ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
J0eff = 6.42e-3;
o = struct('J0eff', J0eff, 'Jxx0', ion.Jxx0, 'hyp', true);

[hstar, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, o);

fprintf('PROBE status        = %s\n', prof.status);
fprintf('PROBE hstar         = %.17g   finite=%d\n', hstar, isfinite(hstar));
fprintf('PROBE nNodes        = %d   n_extend=%d   redensified=%d\n', ...
        numel(prof.hgrid), prof.n_extend, prof.redensified);
fprintf('PROBE hgrid(1)      = %.6g    hgrid(end) = %.6g\n', prof.hgrid(1), prof.hgrid(end));
fprintf('PROBE slope0        = %.17g\n', prof.slope0);
fprintf('PROBE r_star        = %.17g\n', prof.r_star);
fprintf('PROBE m_star        = %.17g\n', prof.m_star);
fprintf('PROBE D_uni_star    = %.17g\n', prof.D_uni_star);

% (b) crit over the grid, and at the node nearest the root
crit = prof.r + J0eff*prof.G0bare;
fprintf('PROBE crit: min=%.6g  max=%.6g  nNeg=%d  nNaN=%d\n', ...
        min(crit), max(crit), nnz(crit < 0), nnz(isnan(crit)));
if isfinite(hstar)
    [~, inear] = min(abs(prof.hgrid - hstar));
    fprintf('PROBE crit nearest root (h=%.6g): %.17g   SIGN=%+d\n', ...
            prof.hgrid(inear), crit(inear), sign(crit(inear)));
    fprintf('PROBE G0bare nearest root = %.17g\n', prof.G0bare(inear));
    fprintf('PROBE crit_star ESTIMATE (r_star + J0eff*G0bare_near) = %.17g\n', ...
            prof.r_star + J0eff*prof.G0bare(inear));
end

% (c) r-1 vs Sigma0 at finite moment
k = find(abs(prof.m) > 1e-3 & isfinite(prof.Sigma0) & isfinite(prof.r), 1, 'last');
if isempty(k)
    fprintf('PROBE *** NO finite-moment converged node -- claim (c) UNREACHABLE ***\n');
else
    d = abs((prof.r(k) - 1) - prof.Sigma0(k));
    fprintf('PROBE finite-m node k=%d  h=%.6g  m=%.6g\n', k, prof.hgrid(k), prof.m(k));
    fprintf('PROBE   r-1    = %.17g\n', prof.r(k) - 1);
    fprintf('PROBE   Sigma0 = %.17g\n', prof.Sigma0(k));
    fprintf('PROBE   |diff| = %.6g   clears 1e-9 = %d\n', d, d > 1e-9);
end

% node convergence census
fprintf('PROBE node_conv: nTrue=%d / %d\n', nnz(prof.node_conv), numel(prof.node_conv));
fprintf('PROBE Sigma0_pm0 = %.17g   K0_pm0 = %.17g\n', prof.Sigma0_pm0, prof.K0_pm0);
fprintf('PROBE existing int fields present: int_Sigma0=%d int_r_minus_1=%d\n', ...
        isfield(prof, 'int_Sigma0'), isfield(prof, 'int_r_minus_1'));
fprintf('\nPROBE DONE\n');
