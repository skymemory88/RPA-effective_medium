% Task 14 probe P6: is local_sigma_loop's strict domain halt REACHABLE?
% denom (dyson ref) = 1 + Sigma0, and Tier 2 re-enters the inner loop WARM from the converged
% Tier-1 Sigma. P5 measured denom_tier1 = 1.2106255568070807 > denom_tier2 = 1.2081157623207905,
% so a ref_margin strictly between them trips Tier 2 and not Tier 1 -- PROVIDED Tier 1 is
% seeded above the floor too (its cold Sigma = 0 start would otherwise trip on iteration 1).
% opts.Sigma_seed is the existing public hook for exactly that.
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT, 'invz_common'));
addpath(genpath(fullfile(ROOT, 'invz_projected')));

ion = invz_ion();
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 10, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
base = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S, ...
              'static_medium', 'strict_1z_dyson_ref');

t1 = invz_solve_point(ion, 1.80, 0.05, [], base);          % converged Tier-1 reference
fprintf('P6 tier1 (cold, default margin): converged=%d denom=%.17g outer=%d\n', ...
    t1.converged, t1.medium_denom, t1.outer_iters);

o = base;  o.odd_tier2 = true;  o.Sigma_seed = t1.Sigma;   % warm PM seed keeps Tier 1 above
o.ref_margin = 1.2095;                                     % strictly between the two denoms
p = invz_solve_point(ion, 1.80, 0.05, [], o);
fprintf('P6 tier2 with ref_margin=%.6g: medium_status=%s converged=%d outer_iters=%d ...\n', ...
    o.ref_margin, p.medium_status, p.converged, p.outer_iters);
fprintf('P6   tier2_iters=%d denom=%.17g margin=%.17g K1=%g G1=%g alpha=%g crit=%g\n', ...
    p.tier2_iters, p.medium_denom, p.medium_margin, p.K(1), p.G(1), p.alpha, p.crit);
fprintf('P6   has C=%d tier2_iters field=%d tla=%d tier2_resid=%d\n', ...
    isfield(p,'C'), isfield(p,'tier2_iters'), isfield(p,'tla'), isfield(p,'tier2_resid'));
ref = invz_solve_point(ion, 1.80, 0.05, [], setfield(base, 'Sigma_seed', t1.Sigma));  %#ok<SFLD>
fa = fieldnames(ref);  fb = fieldnames(p);
fprintf('P6 FIELDSET (vs tier1-only, which has no Tier-2 members): only_in_tier1={%s} only_in_tier2={%s}\n', ...
    strjoin(setdiff(fa,fb)', ','), strjoin(setdiff(fb,fa)', ','));
