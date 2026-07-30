function results = test_invzp_static_domain_resolution()
%TEST_INVZP_STATIC_DOMAIN_RESOLUTION Adversarial scan/margin sensitivity gate.
here = fileparts(mfilename('fullpath'));
projected = fileparts(here);
repo = fileparts(projected);
addpath(projected, fullfile(repo,'invz_common'), repo);

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
D = load(fullfile(repo,'diag_rev3_check.mat'));
[~,~,beta] = invz_matsubara(F.provenance.T,F.provenance.solve_opts.Ecut);
[G0i4,G0e4] = static_split(F.pt.si,F.provenance.T,F.info.Jaa0);

scan_points = [2049;4097;8193];
endpoint_margin = [1e-8;1e-10;1e-12];
healthy_roots = nan(3,1);
node9_roots = nan(3,1);
threeT_roots = nan(3,7);
healthy_x = nan(3,1);
for ic = 1:3
    so = struct('Jsup',F.info.Jcc0,'warn',false, ...
        'scan_points',scan_points(ic),'endpoint_margin',endpoint_margin(ic), ...
        'resid_tol',1e-10);
    [~,~,o4] = invz_emt_static_ordered(F.pt.tl,F.pt.lambda(1:2), ...
        F.pt.Sigma0,F.J,0,beta,F.info.Jcc0,G0i4,G0e4,so);
    healthy_roots(ic) = o4.n_admissible_roots;
    if o4.n_admissible_roots == 1, healthy_x(ic) = o4.root_table.x(1); end

    n9 = D.out(1).nodes{10};
    [~,~,o9] = invz_emt_static_ordered(n9.tl,n9.lam(1:2),n9.Sigma0, ...
        F.J,0,beta,F.info.Jcc0,n9.G0inel0,n9.G0el0,so);
    node9_roots(ic) = o9.n_admissible_roots;

    for inode = 21:27
        nd = D.out(2).nodes{inode+1};
        [~,~,oo] = invz_emt_static_ordered(nd.tl,nd.lam(1:2),nd.Sigma0, ...
            F.J,0,beta,F.info.Jcc0,nd.G0inel0,nd.G0el0,so);
        threeT_roots(ic,inode-20) = oo.n_admissible_roots;
    end
end

classification_ok = all(healthy_roots == 1) && all(node9_roots == 0) && ...
                    all(threeT_roots(:) == 0);
x_spread = max(healthy_x)-min(healthy_x);
if ~classification_ok || x_spread > 1e-7
    error('invz:testStaticResolution', ...
        'static classification changed under scan/margin refinement (x spread %.3g).',x_spread);
end

results = struct();
results.config = table(scan_points,endpoint_margin);
results.healthy_roots = healthy_roots;
results.healthy_x = healthy_x;
results.healthy_x_spread = x_spread;
results.node9_roots = node9_roots;
results.threeT_roots = threeT_roots;
results.classification_ok = classification_ok;
diag_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'wp1_static_resolution_gate.mat');
save(diag_path,'results','-v7');
fprintf(['test_invzp_static_domain_resolution: stable across 2049/4097/8193 ' ...
    'points and margins 1e-8/1e-10/1e-12; x spread %.3g\n'],x_spread);
end

function [G0i,G0e] = static_split(si,T,Jxx0)
c0i = invz_chi0z(si,T,0,struct('elastic',false));
G0i = -real(c0i(3,3,1));
c0 = invz_chi0z(si,T,0,struct('elastic',true));
X = real(c0(:,:,1));
feedback = X(3,1)*(Jxx0/(1-Jxx0*X(1,1)))*X(1,3);
G0e = -(X(3,3)+feedback)-G0i;
end
