% Task-9 re-review probe (fix round 1). Targeted measurements only, no suite runs.
%  (A) F5: fieldnames(res.blockB) parity under opts.debug_resummed, both paths, incl. order.
%  (B) F3: does node.G0(1)=NaN split medD.converged from medD.dynamic_converged, leave
%      slots 2:end finite, and leave all four block residuals BIT-identical?
%  (C) F3 discrimination: recompute passD both ways from the same medD, off-implementation.
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
fprintf('info.accepted = %d   nw = %d\n', info.accepted, numel(wn));

% ---- (A) F5 fieldnames parity under debug_resummed ------------------------------------------
bad = state;  bad.Sigma(1) = -1;
dN = invz_ordered_residual(node, state, struct('debug_resummed', true));
dP = invz_ordered_residual(node, bad,   struct('debug_resummed', true));
fnN = fieldnames(dN.blockB);  fnP = fieldnames(dP.blockB);
fprintf('F5 blockB fieldnames identical incl ORDER: %d   (last field: %s / %s)\n', ...
        isequal(fnN, fnP), fnN{end}, fnP{end});
fprintf('F5 res top-level fieldnames identical: %d\n', isequal(fieldnames(dN), fieldnames(dP)));
fprintf('F5 resid_resummed: normal=%.6g  preflight=%.6g (isnan=%d)\n', ...
        dN.blockB.resid_resummed, dP.blockB.resid_resummed, isnan(dP.blockB.resid_resummed));
% and WITHOUT the flag, the field must be absent on BOTH paths
nN = invz_ordered_residual(node, state);  nP = invz_ordered_residual(node, bad);
fprintf('F5 no-flag: resid_resummed absent on normal=%d preflight=%d, fieldnames equal=%d\n', ...
        ~isfield(nN.blockB,'resid_resummed'), ~isfield(nP.blockB,'resid_resummed'), ...
        isequal(fieldnames(nN.blockB), fieldnames(nP.blockB)));

% ---- (B)+(C) F3 discrimination --------------------------------------------------------------
res0 = invz_ordered_residual(node, state);
node2 = node;  node2.G0(1) = NaN;
medD = invz_emt_scalar(node2.G0, state.Sigma, node2.Jnu_flat, node2.eopts);
fprintf('F3 medD.converged=%d  medD.dynamic_converged=%d  K(2:end) finite=%d  G(2:end) finite=%d\n', ...
        medD.converged, medD.dynamic_converged, all(isfinite(medD.K(2:end))), ...
        all(isfinite(medD.G(2:end))));
fprintf('F3 medD.K(1) isnan=%d  medD.G(1) isnan=%d  medD.medium_status=%s\n', ...
        isnan(medD.K(1)), isnan(medD.G(1)), medD.medium_status);
med0 = invz_emt_scalar(node.G0, state.Sigma, node.Jnu_flat, node.eopts);
fprintf('F3 medD.K(2:end) BITWISE equal to unperturbed: %d\n', ...
        isequaln(medD.K(2:end), med0.K(2:end)));
res = invz_ordered_residual(node2, state);
fprintf('F3 bit-identical resid A=%d B=%d C=%d D=%d\n', ...
        isequaln(res.blockA.resid,res0.blockA.resid), isequaln(res.blockB.resid,res0.blockB.resid), ...
        isequaln(res.blockC.resid,res0.blockC.resid), isequaln(res.blockD.resid,res0.blockD.resid));
fprintf('F3 blockD.pass=%d accepted=%d  finite=%d  (res0: pass=%d accepted=%d)\n', ...
        res.blockD.pass, res.accepted, res.finite, res0.blockD.pass, res0.accepted);
% off-implementation recomputation of passD under BOTH gate variants
rD = max(abs(state.K(2:end) - medD.K(2:end)));
refD = max(abs(medD.K(2:end)));
sD_abs = 1e-8*max(abs(node.Jnu_flat));  sD_rel = 1e-8;
finite_D = all(isfinite(state.K)) && all(isfinite(state.lam)) && all(isfinite(state.Sigma));
base = isfinite(rD) && (rD < sD_abs + sD_rel*refD) && finite_D;
fprintf('F3 OFF-IMPL passD: with dynamic_converged=%d  with whole-PM converged=%d\n', ...
        base && medD.dynamic_converged, base && medD.converged);
fprintf('F3 other blocks pass under corruption: A=%d B=%d C=%d, blockB.status=%s\n', ...
        res.blockA.pass, res.blockB.pass, res.blockC.pass, res.blockB.status);
disp('PROBE DONE');
