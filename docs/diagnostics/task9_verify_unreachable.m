% Task 9 self-verification (controller-requested): confirm directly, on the committed
% invz_ordered_residual.m, whether Step 4's strict `else rB = NaN; passB = false;
% convB = false;` branch is reached on the accepted synthetic-control fixture, i.e. whether
% the independent preflight's 'ok' verdict and local_blockB's own internal medium_status ever
% disagree on this fixture.
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root);
addpath(fullfile(root, 'invz_projected'));
addpath(fullfile(root, 'invz_common'));

ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];  hz = 0.02;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
J0eff = 6.42e-3;  Jxx0 = ion.Jxx0;  tmf = 'legacy_x';  Ecut = 40;
[wn, wts, beta] = invz_matsubara(T, Ecut);
si = invz_single_ion(ion, T, Bx, struct('hyp', true, 'hz_fixed', hz, 'Jxx0', Jxx0, ...
                                        'transverse_mf', tmf));
tl = invz_twolevel_ordered(ion, T, Bx, hz, struct('Jxx0', Jxx0, 'transverse_mf', tmf));
c0  = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
G0  = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si, T, 1i*wn(1), struct('elastic', false));
G0inel0 = -real(c0i(3,3,1));
X  = real(c0(:, :, 1));
fb = X(3,1) * (Jxx0 / (1 - Jxx0*X(1,1))) * X(1,3);
G0bare0 = -(X(3,3) + fb);
G0el0   = G0bare0 - G0inel0;
g = real(invz_g(tl, 1i*wn));
mom = invz_coupling_moments(Jnu);
eso = struct('warn', false, 'static_medium', 'strict_1z_dyson_ref', 'Jmom', mom);
eopts = struct('static_medium', 'strict_1z_dyson_ref', 'Jmom', mom);
node = struct('tl', tl, 'G0', G0, 'g', g, 'wts', wts, 'wn', wn, 'beta', beta, ...
              'J0eff', J0eff, 'G0inel0', G0inel0, 'G0el0', G0el0, 'G0bare0', G0bare0, ...
              'eso', eso, 'eopts', eopts, 'Jnu_flat', Jnu, 'Jmom', mom);
[state, info] = invz_ordered_node_solve(node, [], struct('trace', false));
fprintf('info.accepted = %d\n', info.accepted);

% Independent preflight recomputation, standalone (mirrors invz_ordered_residual's own).
[Gref_pf, ref_pf] = invz_static_medium_reference(node.G0bare0, state.Sigma(1), 'strict_1z_dyson_ref', eso);
[~, clo_pf] = invz_medium_moment_closure(Gref_pf, node.Jmom, 'strict_1z_dyson_ref');
fprintf('PREFLIGHT ref.status=%s  clo.status=%s\n', ref_pf.status, clo_pf.status);

res = invz_ordered_residual(node, state);
fprintf('res.blockB.status    = %s\n', res.blockB.status);
fprintf('res.blockB.resid     = %.17g   (exactly 0 => the else rB=NaN branch was NOT taken)\n', res.blockB.resid);
fprintf('res.blockB.pass      = %d\n', res.blockB.pass);
fprintf('res.blockB.converged = %d\n', res.blockB.converged);
fprintf('res.accepted         = %d\n', res.accepted);

% Direct proof of bitwise equality between the two G0bare sources this task's preflight and
% local_blockB's internal invz_emt_static_ordered call each independently consume.
fprintf('\nnode.G0bare0            = %.17g\n', node.G0bare0);
fprintf('node.G0inel0+node.G0el0 = %.17g\n', node.G0inel0 + node.G0el0);
fprintf('bitwise equal?          = %d\n', node.G0bare0 == (node.G0inel0 + node.G0el0));
fprintf('isequaln(node.Jmom, eso.Jmom) = %d\n', isequaln(node.Jmom, eso.Jmom));

if strcmp(res.blockB.status, 'ok') && res.blockB.resid == 0 && res.blockB.converged
    fprintf('\nCONCLUSION: else rB=NaN;passB=false;convB=false branch UNREACHABLE on this fixture.\n');
else
    fprintf('\nCONCLUSION: branch WAS reached or resid nonzero -- investigate.\n');
end
