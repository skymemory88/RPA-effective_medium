ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));

nw = 8;
G0    = -(0.5 ./ (1:nw)').^0 * 0.5;
G0(1) = -0.9;
Sigma = 0.05*ones(nw, 1);
Jnu   = linspace(-2e-3, 6e-3, 24).';

mom = invz_coupling_moments(Jnu);
fprintf('Jbar=%.10g mu2=%.10g mu3=%.10g mu4=%.10g\n', mom.Jbar, mom.mu2, mom.mu3, mom.mu4);

[Gref_d, refd] = invz_static_medium_reference(G0(1), Sigma(1), 'strict_1z_dyson_ref');
[Gref_b, refb] = invz_static_medium_reference(G0(1), Sigma(1), 'strict_1z_bare_ref');
fprintf('Gref_dyson=%.10g  Gref_bare=%.10g\n', Gref_d, Gref_b);

[K0_d, cld] = invz_medium_moment_closure(Gref_d, mom, 'strict_1z_dyson_ref');
[K0_b, clb] = invz_medium_moment_closure(Gref_b, mom, 'strict_1z_bare_ref');
fprintf('K0_dyson=%.10g  K0_bare=%.10g\n', K0_d, K0_b);
ratio = abs(K0_d - K0_b)/abs(K0_b);
fprintf('ratio (direct primitives) = %.10g\n', ratio);

d = invz_emt_scalar(G0, Sigma, Jnu, struct('static_medium', 'strict_1z_dyson_ref'));
b = invz_emt_scalar(G0, Sigma, Jnu, struct('static_medium', 'strict_1z_bare_ref'));
fprintf('d.K(1)=%.10g  b.K(1)=%.10g\n', d.K(1), b.K(1));
ratio2 = abs(d.K(1)-b.K(1))/abs(b.K(1));
fprintf('ratio (via invz_emt_scalar) = %.10g\n', ratio2);

fprintf('\n--- scanning Sigma(1) to find where ratio crosses 1e-3 ---\n');
for s0 = [0.05 0.1 0.2 0.3 0.383 0.4 0.5 0.8 1.0 1.5 2.0]
    [Gd,~] = invz_static_medium_reference(G0(1), s0, 'strict_1z_dyson_ref');
    [Gb,~] = invz_static_medium_reference(G0(1), s0, 'strict_1z_bare_ref');
    [K0d,~] = invz_medium_moment_closure(Gd, mom, 'strict_1z_dyson_ref');
    [K0b,~] = invz_medium_moment_closure(Gb, mom, 'strict_1z_bare_ref');
    fprintf('Sigma0=%.4g  ratio=%.6g\n', s0, abs(K0d-K0b)/abs(K0b));
end
