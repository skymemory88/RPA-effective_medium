% Task 17 characterisation script (controller-notes.md §7, and the G1/G11 "measure the nH
% scaling" / "measure the medium_status histogram + trc.nodes term reasons" requirements).
% READ-ONLY measurement: reproduces test bodies to print intermediate numbers the failure
% diagnostics did not show. Does not modify any committed file.
REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
% Mirror the committed tests' own setupOnce EXACTLY (flat, not genpath, with the repo root
% explicitly added) -- confirmed by controlled A/B test that a genpath(invz_projected)-only
% recipe run via run() (which changes cwd to this script's own folder) fails to resolve
% qVec_generator, because repo root then has no explicit path entry. Not a defect in the
% committed tests, whose setupOnce always adds the repo root explicitly.
addpath(fullfile(REPO, 'invz_projected'));
addpath(REPO);
addpath(fullfile(REPO, 'invz_common'));

fprintf('##################################################################\n');
fprintf('### (1) G1d nH-scaling: e(33), e(65), e(129) and ratios\n');
fprintf('##################################################################\n');
ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, 'static_medium', 'strict_1z_dyson_ref');
e = nan(1, 3);  nHs = [33 65 129];
statuses1d = cell(1,3);
for k = 1:3
    ok_ = o;  ok_.nH = nHs(k);
    [~, p] = invz_hmf_ordered(ion, T, Bx, Jnu, ok_);
    statuses1d{k} = p.status;
    dF = diff(p.F)./diff(p.hgrid);
    cm = 0.5*(p.crit(1:end-1) + p.crit(2:end));
    m_ = isfinite(dF) & isfinite(cm);
    e(k) = max(abs(dF(m_) - cm(m_)));
end
for k = 1:3
    fprintf('nH=%4d  status=%-6s  e=%.15e\n', nHs(k), statuses1d{k}, e(k));
end
fprintf('ratio e(1)/e(2) = %.6f  (2nd order expects ~4)\n', e(1)/e(2));
fprintf('ratio e(2)/e(3) = %.6f  (2nd order expects ~4)\n', e(2)/e(3));
fprintf('log2 ratios (empirical order): p12=%.4f  p23=%.4f\n', log2(e(1)/e(2)), log2(e(2)/e(3)));

fprintf('\n##################################################################\n');
fprintf('### (2) G11: full status/medium_status/term_reason/trace characterisation\n');
fprintf('##################################################################\n');
prov = struct('grid', [16 16 16], 'dpRng', 30, 'dipole', 'bruteforce', 'cache', false);
[Jnu11, info11] = invz_bz_couplings(ion, prov);
mom11 = invz_coupling_moments(Jnu11(:));
fprintf('digest = %s\n', invz_exact_numeric_digest(Jnu11(:)));
fprintf('digest matches frozen  = %d\n', strcmp(invz_exact_numeric_digest(Jnu11(:)), ...
    'ddb9532d11326458554b93b2ce09c80a3299cca9baa69202dc311f138b4fae17'));
fprintf('info.dipole.backend = %s ; isfield(info,''grid'') = %d\n', info11.dipole.backend, isfield(info11,'grid'));
fprintf('mom.n = %d ; mom.Jbar = %.10e ; mom.mu2 = %.10e ; info.Jcc0 = %.10e\n', mom11.n, mom11.Jbar, mom11.mu2, info11.Jcc0);
o11 = struct('J0eff', info11.Jcc0, 'Jxx0', ion.Jxx0, 'hyp', true, ...
             'static_medium', 'strict_1z_dyson_ref', 'Jmom', mom11, 'trace', true);
[hstar11, p11, trc11] = invz_hmf_ordered(ion, 0.1, [1 0 0], Jnu11(:), o11);
fprintf('p.status = %s\n', p11.status);
fprintf('hstar = %.15g  isfinite(hstar)=%d  hstar>0=%d\n', hstar11, isfinite(hstar11), hstar11 > 0);
fprintf('p.crit_star = %.15g  isfinite=%d\n', p11.crit_star, isfinite(p11.crit_star));
fprintf('p.slope0 = %.15g\n', p11.slope0);
fprintf('numel(p.hgrid) = %d ; numel(p.medium_status) = %d\n', numel(p11.hgrid), numel(p11.medium_status));
print_hist('p11.medium_status', p11.medium_status);
print_hist('p11.node_term_reason', p11.node_term_reason);
fprintf('trc.enabled = %d ; numel(trc.nodes) = %d\n', trc11.enabled, numel(trc11.nodes));
print_hist('trc11.nodes.medium_status', {trc11.nodes.medium_status});
print_hist('trc11.nodes.term_reason', {trc11.nodes.term_reason});
fprintf('\n-- trc.nodes detail (id,h,accepted,medium_status,term_reason) for every non-ok node --\n');
nbad = 0;
for i = 1:numel(trc11.nodes)
    nd = trc11.nodes(i);
    if ~strcmp(nd.medium_status, 'ok') || ~strcmp(nd.term_reason, 'accepted')
        nbad = nbad + 1;
        fprintf('  node %3d: id=%s h=%.6g accepted=%d medium_status=%-20s term_reason=%s\n', ...
            i, safe_field(nd,'id'), safe_field(nd,'h'), nd.accepted, nd.medium_status, nd.term_reason);
    end
end
fprintf('total non-(ok & accepted) trace nodes: %d / %d\n', nbad, numel(trc11.nodes));

fprintf('\n##################################################################\n');
fprintf('### (3) G13: state.K(1) vs med.K(1), and the leak-detectability margin\n');
fprintf('##################################################################\n');
T13 = 0.31;  Bx13 = [2.85 0 0];  hz13 = 0.02;
Jnu13 = linspace(-2e-3, 6.0e-3, 24).';  J0eff13 = 6.42e-3;
[wn13, wts13, beta13] = invz_matsubara(T13, 40);
si13 = invz_single_ion(ion, T13, Bx13, struct('hyp', true, 'hz_fixed', hz13, 'Jxx0', ion.Jxx0, ...
                                        'transverse_mf', 'legacy_x'));
tl13 = invz_twolevel_ordered(ion, T13, Bx13, hz13, struct('Jxx0', ion.Jxx0, ...
                                                  'transverse_mf', 'legacy_x'));
c013  = invz_chi0z(si13, T13, 1i*wn13, struct('elastic', true));
G013  = -real(squeeze(c013(3,3,:)));
c0i13 = invz_chi0z(si13, T13, 1i*wn13(1), struct('elastic', false));
G0inel013 = -real(c0i13(3,3,1));
X13  = real(c013(:, :, 1));
fb13 = X13(3,1) * (ion.Jxx0/(1 - ion.Jxx0*X13(1,1))) * X13(1,3);
G0bare013 = -(X13(3,3) + fb13);
node13 = struct('tl', tl13, 'G0', G013, 'g', real(invz_g(tl13, 1i*wn13)), 'wts', wts13, 'wn', wn13, ...
    'beta', beta13, 'J0eff', J0eff13, 'G0inel0', G0inel013, 'G0el0', G0bare013 - G0inel013, ...
    'G0bare0', G0bare013, 'eso', struct('warn', false, 'static_medium', 'strict_1z_dyson_ref'), ...
    'eopts', struct('static_medium', 'strict_1z_dyson_ref'), 'Jnu_flat', Jnu13, ...
    'Jmom', invz_coupling_moments(Jnu13));
[state13, ~] = invz_ordered_node_solve(node13, [], struct('trace', false));
med13 = invz_emt_scalar(node13.G0, state13.Sigma, node13.Jnu_flat, node13.eopts);
fprintf('state.K(1) = %.15g\n', state13.K(1));
fprintf('med.K(1)   = %.15g\n', med13.K(1));
fprintf('state.K(1) - med.K(1) = %.6e  (relative: %.6e)\n', state13.K(1)-med13.K(1), ...
    (state13.K(1)-med13.K(1))/med13.K(1));
lam_from_exported13 = invz_lambdas(state13.K, node13.g, node13.wts, node13.beta, [1 2 3]);
diff_lam = state13.lam - lam_from_exported13;
relerr_lam = diff_lam ./ lam_from_exported13;
fprintf('state.lam vs lam_from_exported: abs diff = [%s]\n', num2str(diff_lam(:).', '%.6e  '));
fprintf('                                 rel diff = [%s]\n', num2str(relerr_lam(:).', '%.6e  '));
fprintf('max abs relative error = %.6e  (RelTol asserted = 1e-10)\n', max(abs(relerr_lam)));
Kleak13 = state13.K;  Kleak13(1) = med13.K(1);
lam_leaked13 = invz_lambdas(Kleak13, node13.g, node13.wts, node13.beta, [1 2 3]);
margin13 = max(abs(lam_leaked13 - state13.lam));
fprintf('leak-detectability margin max(abs(lam_leaked - state.lam)) = %.6e  (must exceed 1e-12)\n', margin13);

fprintf('\n##################################################################\n');
fprintf('### (4) G2: all four comparisons, values and headroom\n');
fprintf('##################################################################\n');
o2 = o;  % same fx() fixture as G1/G2 (J0eff=6.42e-3, hyp=true, strict scheme, no nH override)
Jscale2 = max(abs(Jnu));  K_atol2 = 1e-14;  K_rtol2 = 1e-12;
ptp2 = invz_solve_point(ion, T, Bx, Jnu, o2);
[~, p2] = invz_hmf_ordered(ion, T, Bx, Jnu, o2);
fprintf('p.slope0 = %.15g ; ptp.crit = %.15g\n', p2.slope0, ptp2.crit);
fprintf('|p.slope0 - ptp.crit| = %.6e  (AbsTol asserted = 1e-6)\n', abs(p2.slope0 - ptp2.crit));
G0pm2 = -ptp2.chi0cc0;
GrefP2 = invz_static_medium_reference(G0pm2, ptp2.Sigma0, 'strict_1z_dyson_ref');
GrefO2 = invz_static_medium_reference(p2.G0bare_pm0, p2.Sigma0_pm0, 'strict_1z_dyson_ref');
fprintf('GrefO = %.15g ; GrefP = %.15g\n', GrefO2, GrefP2);
fprintf('|GrefO-GrefP|/|GrefP| = %.6e  (RelTol asserted = 1e-6)\n', abs(GrefO2-GrefP2)/abs(GrefP2));
mom2 = invz_coupling_moments(Jnu);
KP2 = invz_medium_moment_closure(GrefP2, mom2, 'strict_1z_dyson_ref');
KO2 = invz_medium_moment_closure(GrefO2, mom2, 'strict_1z_dyson_ref');
Kgate2 = K_atol2 + K_rtol2*max([abs(KO2), abs(KP2), Jscale2]);
fprintf('KO = %.15g ; KP = %.15g\n', KO2, KP2);
fprintf('|KO-KP| = %.6e\n', abs(KO2-KP2));
fprintf('p.K0_pm0 = %.15g ; ptp.K(1) = %.15g\n', p2.K0_pm0, ptp2.K(1));
fprintf('|p.K0_pm0-ptp.K(1)| = %.6e\n', abs(p2.K0_pm0-ptp2.K(1)));
fprintf('Kgate = %.6e  (= %.3e + %.3e*max([%.6e,%.6e,%.6e]))\n', Kgate2, K_atol2, K_rtol2, abs(KO2), abs(KP2), Jscale2);
fprintf('headroom: |KO-KP|/Kgate = %.6f  (<=1 passes)\n', abs(KO2-KP2)/Kgate2);
fprintf('headroom: |p.K0_pm0-ptp.K(1)|/Kgate = %.6f  (<=1 passes)\n', abs(p2.K0_pm0-ptp2.K(1))/Kgate2);

fprintf('\n##################################################################\n');
fprintf('### (5) G15b: status + which term_reason/medium_status path fired\n');
fprintf('##################################################################\n');
o15b = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
           'static_medium', 'strict_1z_dyson_ref', 'ref_margin', 1e9, 'trace', true);
[hstar15b, p15b, trc15b] = invz_hmf_ordered(ion, 0.31, [2.85 0 0], Jnu, o15b);
fprintf('p.status = %s\n', p15b.status);
print_hist('p15b.medium_status', p15b.medium_status);
print_hist('p15b.node_term_reason', p15b.node_term_reason);
print_hist('trc15b.nodes.medium_status', {trc15b.nodes.medium_status});

fprintf('\n### done ###\n');

function print_hist(label, c)
u = unique(c);
fprintf('%s histogram (n=%d):\n', label, numel(c));
for i = 1:numel(u)
    fprintf('    %-24s : %d\n', u{i}, sum(strcmp(c, u{i})));
end
end

function s = safe_field(nd, name)
if isfield(nd, name)
    v = nd.(name);
    if ischar(v), s = v; else, s = num2str(v); end
else
    s = '<absent>';
end
end
