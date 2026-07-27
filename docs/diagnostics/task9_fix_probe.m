% Probe for task-9 review fix F3: verify that corrupting node.G0(1) alone (leaving eso/eopts,
% state untouched) discriminates invz_ordered_residual.m:273/276's medD.dynamic_converged
% switch from the reverted medD.converged. Foreground, one-shot, not committed.
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(genpath(fullfile(root, 'invz_projected')));
addpath(fullfile(root, 'invz_common'));

% --- build_strict_fixture, copied verbatim from test_invz_ordered_residual_strict.m ---------
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
fprintf('fixture accepted = %d\n', info.accepted);

% --- baseline residual (uncorrupted) ---------------------------------------------------------
res0 = invz_ordered_residual(node, state);
fprintf('baseline: accepted=%d  blockD.pass=%d  blockD.resid=%.6e\n', ...
    res0.accepted, res0.blockD.pass, res0.blockD.resid);

% --- corrupt ONLY node.G0(1); eso/eopts/state untouched --------------------------------------
node2 = node;
node2.G0(1) = NaN;

medD_direct = invz_emt_scalar(node2.G0, state.Sigma, node2.Jnu_flat, node2.eopts);
fprintf('\ndirect invz_emt_scalar on corrupted G0: converged=%d dynamic_converged=%d\n', ...
    medD_direct.converged, medD_direct.dynamic_converged);
fprintf('K(1)=%.6g (expect NaN)  K(2)=%.6g  K(end)=%.6g  (both expect finite, unaffected)\n', ...
    medD_direct.K(1), medD_direct.K(2), medD_direct.K(end));
fprintf('G(1)=%.6g (expect NaN)  G(2)=%.6g  G(end)=%.6g\n', ...
    medD_direct.G(1), medD_direct.G(2), medD_direct.G(end));

res2 = invz_ordered_residual(node2, state);
fprintf('\nresidual w/ corrupted node.G0(1), CURRENT code (dynamic_converged):\n');
fprintf('  blockA.pass=%d resid=%.6e (expect == baseline %.6e)\n', res2.blockA.pass, res2.blockA.resid, res0.blockA.resid);
fprintf('  blockB.pass=%d resid=%.6e (expect == baseline %.6e)\n', res2.blockB.pass, res2.blockB.resid, res0.blockB.resid);
fprintf('  blockC.pass=%d resid=%.6e (expect == baseline %.6e)\n', res2.blockC.pass, res2.blockC.resid, res0.blockC.resid);
fprintf('  blockD.pass=%d resid=%.6e (expect == baseline %.6e)\n', res2.blockD.pass, res2.blockD.resid, res0.blockD.resid);
fprintf('  finite=%d  accepted=%d  (expect finite=1, accepted=1: dynamic_converged ignores slot 1)\n', ...
    res2.finite, res2.accepted);

% Sanity: also confirm the SAME discrimination holds if node2.eso/eopts are switched to
% 'resummed' instead (i.e. this is not an artifact of the strict scheme) -- informational only,
% not required for the fix, just to corroborate the reviewer's own measurement that both
% schemes show medD.converged==medD.dynamic_converged==1 UNPERTURBED, and to see whether the
% resummed leg even reaches Block D the same way (it does -- no strict-only preflight).
node3 = node2;
node3.eso.static_medium = 'resummed';
node3.eopts.static_medium = 'resummed';
node3 = rmfield(node3, 'Jmom');   % harmless-absent under resummed (test_jmom_required_only_under_strict)
try
    res3 = invz_ordered_residual(node3, state);
    fprintf('\n[info-only, resummed both legs] blockD.pass=%d accepted=%d\n', res3.blockD.pass, res3.accepted);
catch ME
    fprintf('\n[info-only, resummed both legs] threw: %s (%s)\n', ME.identifier, ME.message);
end

fprintf('\nDONE\n');
