% Task-9 review probe: three things source alone cannot settle.
%  (1) fieldnames parity of res and res.blockB between the preflight path and the normal path,
%      with and without opts.debug_resummed.
%  (2) whether medD.converged == medD.dynamic_converged at the strict fixture (i.e. whether
%      the Step-5 Block-D switch is exercised by any committed test).
%  (3) res.stall on both paths for the same opts.dS.
ROOT = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(ROOT, 'invz_projected'))); addpath(fullfile(ROOT, 'invz_common'));
addpath(ROOT);

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

% ---- (1) fieldnames parity -----------------------------------------------------------------
resN = invz_ordered_residual(node, state, struct('dS', 1e-10));
bad = state;  bad.Sigma(1) = -1;                        % forces the preflight short circuit
resP = invz_ordered_residual(node, bad, struct('dS', 1e-10));
fprintf('res fieldnames identical (normal vs preflight): %d\n', ...
        isequal(fieldnames(resN), fieldnames(resP)));
fprintf('blockB fieldnames identical: %d\n', ...
        isequal(fieldnames(resN.blockB), fieldnames(resP.blockB)));
fprintf('blockA fieldnames identical: %d   blockC: %d   blockD: %d\n', ...
        isequal(fieldnames(resN.blockA), fieldnames(resP.blockA)), ...
        isequal(fieldnames(resN.blockC), fieldnames(resP.blockC)), ...
        isequal(fieldnames(resN.blockD), fieldnames(resP.blockD)));
fprintf('stability fieldnames identical: %d\n', ...
        isequal(fieldnames(resN.stability), fieldnames(resP.stability)));
fprintf('blockA scaleN=[%g %g] scaleP=[%g %g]\n', resN.blockA.scale_abs, ...
        resN.blockA.scale_rel, resP.blockA.scale_abs, resP.blockA.scale_rel);
fprintf('blockD scaleN=[%g %g] scaleP=[%g %g]\n', resN.blockD.scale_abs, ...
        resN.blockD.scale_rel, resP.blockD.scale_abs, resP.blockD.scale_rel);
fprintf('blockB scaleN=[%g %g] scaleP=[%g %g]\n', resN.blockB.scale_abs, ...
        resN.blockB.scale_rel, resP.blockB.scale_abs, resP.blockB.scale_rel);

dN = invz_ordered_residual(node, state, struct('dS', 1e-10, 'debug_resummed', true));
dP = invz_ordered_residual(node, bad,   struct('dS', 1e-10, 'debug_resummed', true));
fprintf('DEBUG_RESUMMED: blockB fieldnames identical: %d  (normal has resid_resummed=%d, preflight=%d)\n', ...
        isequal(fieldnames(dN.blockB), fieldnames(dP.blockB)), ...
        isfield(dN.blockB, 'resid_resummed'), isfield(dP.blockB, 'resid_resummed'));

% ---- (2) is the Block-D switch exercised anywhere? ------------------------------------------
medS = invz_emt_scalar(node.G0, state.Sigma, node.Jnu_flat, node.eopts);
fprintf('STRICT fixture: medD.converged=%d  medD.dynamic_converged=%d  (differ=%d)\n', ...
        medS.converged, medS.dynamic_converged, medS.converged ~= medS.dynamic_converged);
legopts = rmfield(node.eopts, 'static_medium');
medL = invz_emt_scalar(node.G0, state.Sigma, node.Jnu_flat, legopts);
fprintf('RESUMMED eopts: medD.converged=%d  medD.dynamic_converged=%d  (differ=%d)\n', ...
        medL.converged, medL.dynamic_converged, medL.converged ~= medL.dynamic_converged);

% ---- (3) stall + status on both paths -------------------------------------------------------
fprintf('normal:    accepted=%d finite=%d stall=%d  blockB.status=%s conv=%d refden=%.17g omit3=%g omit4=%g\n', ...
        resN.accepted, resN.finite, resN.stall, resN.blockB.status, resN.blockB.converged, ...
        resN.blockB.ref_denom, resN.blockB.omit_mu3, resN.blockB.omit_cubic);
fprintf('preflight: accepted=%d finite=%d stall=%d  blockB.status=%s conv=%d refden=%.17g omit3=%g omit4=%g\n', ...
        resP.accepted, resP.finite, resP.stall, resP.blockB.status, resP.blockB.converged, ...
        resP.blockB.ref_denom, resP.blockB.omit_mu3, resP.blockB.omit_cubic);
fprintf('normal stability: class=%s crit=%.17g pass=%d\n', resN.stability.class, ...
        resN.stability.crit, resN.stability.pass);
disp('PROBE DONE');
