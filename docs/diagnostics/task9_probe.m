% Task 9 fixture-viability probe (controller, pre-dispatch).
% The brief's build_strict_fixture() asserts info.accepted from invz_ordered_node_solve, which
% itself gates on invz_ordered_residual -- the function Task 9 rewrites. Measure whether the
% fixture reaches an accepted strict node at CURRENT HEAD before any code is written.
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

fprintf('PROBE G0bare0        = %.17g\n', G0bare0);
fprintf('PROBE G0inel0+G0el0  = %.17g\n', G0inel0 + G0el0);
fprintf('PROBE bitwise equal  = %d\n', G0bare0 == (G0inel0 + G0el0));

for which = {'strict_1z_dyson_ref', 'resummed'}
    sm = which{1};
    eso   = struct('warn', false, 'static_medium', sm, 'Jmom', mom);
    eopts = struct('static_medium', sm, 'Jmom', mom);
    node = struct('tl', tl, 'G0', G0, 'g', g, 'wts', wts, 'wn', wn, 'beta', beta, ...
                  'J0eff', J0eff, 'G0inel0', G0inel0, 'G0el0', G0el0, 'G0bare0', G0bare0, ...
                  'eso', eso, 'eopts', eopts, 'Jnu_flat', Jnu, 'Jmom', mom);
    fprintf('\n===== %s =====\n', sm);
    try
        [state, info] = invz_ordered_node_solve(node, [], struct('trace', false));
        r = info.res;
        fprintf('%s accepted     = %d\n', sm, info.accepted);
        fprintf('%s term_reason  = %s\n', sm, info.term_reason);
        fprintf('%s A resid=%.6g scale=%.3g pass=%d err=%s\n', sm, r.blockA.resid, ...
                r.blockA.scale_abs, r.blockA.pass, r.blockA.err);
        fprintf('%s B resid=%.6g scale=%.3g pass=%d conv=%d err=%s\n', sm, r.blockB.resid, ...
                r.blockB.scale_abs, r.blockB.pass, r.blockB.converged, r.blockB.err);
        fprintf('%s C resid=%.6g scale=%.3g pass=%d err=%s\n', sm, r.blockC.resid, ...
                r.blockC.scale_abs, r.blockC.pass, r.blockC.err);
        fprintf('%s D resid=%.6g scale=%.3g pass=%d err=%s\n', sm, r.blockD.resid, ...
                r.blockD.scale_abs, r.blockD.pass, r.blockD.err);
        fprintf('%s finite=%d D_uni=%.6g Dq_min=%.6g Dq_max=%.6g\n', sm, r.finite, ...
                r.D_uni, r.Dq_min, r.Dq_max);
        fprintf('%s K0s=%.17g Sigma(1)=%.17g\n', sm, state.K0s, state.Sigma(1));
        % the strict algebraic check Task 9's new block B will use
        if ~strcmp(sm, 'resummed')
            [Gref, refi] = invz_static_medium_reference(G0bare0, state.Sigma(1), sm, ...
                struct('ref_margin', 1e-6));
            [Kstrict, clo] = invz_medium_moment_closure(Gref, mom, sm);
            fprintf('%s PREFLIGHT ref.status=%s clo.status=%s ref.denom=%.17g\n', sm, ...
                    refi.status, clo.status, refi.denom);
            fprintf('%s |K0s-Kstrict| = %.6g   gate(1e-14+1e-12*max) = %.6g\n', sm, ...
                    abs(state.K0s - Kstrict), ...
                    1e-14 + 1e-12*max([abs(state.K0s), abs(Kstrict), max(abs(Jnu))]));
            [~, ~, so] = invz_emt_static_ordered(tl, state.lam(1:2), state.Sigma(1), Jnu, ...
                state.K0s, beta, J0eff, G0inel0, G0el0, eso);
            fprintf('%s crit = so.r + J0eff*G0bare0 = %.17g\n', sm, so.r + J0eff*G0bare0);
            fprintf('%s so.medium_status=%s so.resid=%.6g omit_max=%.6g\n', sm, ...
                    so.medium_status, so.resid, so.omit_max);
            % what class the Task 9 classifier would assign
            Dq_tol = 1e-6*max(1, abs(so.G0bare)*max(abs(Jnu)));
            fprintf('%s D_tol = %.6g\n', sm, Dq_tol);
        end
    catch err
        fprintf('%s THREW %s : %s\n', sm, err.identifier, err.message);
    end
end
fprintf('\nPROBE DONE\n');
