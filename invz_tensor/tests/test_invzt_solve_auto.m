function tests = test_invzt_solve_auto
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
ion = invz_ion();
g   = invzt_qgrid(16, 'halfopen');
tc.TestData.ion = ion;
tc.TestData.lat = invzt_jq_tensor(ion, g, struct('dpRng', 30, 'cache', true));
end

function test_below_Bc_is_phase1(tc)
[pt, phase, di] = invzt_solve_auto(tc.TestData.ion, 0.1, [3.0 0 0], tc.TestData.lat, struct());
verifyEqual(tc, phase, 1);
verifyTrue(tc, pt.is_ordered && pt.converged);
verifyTrue(tc, di.para.attempted && ~di.para.accepted);   % PM leg was consulted and rejected
end

function test_at_4p8_is_phase2_HARD(tc)
% THE P0-1 gate: 4.8 T is on the tensor-PM side (crit = +0.023 measured) even
% though the bare-MF ordered leg still has m0 = 1.17 there. Stability-based
% selection MUST return PM. An ordered-first implementation fails here.
[pt, phase, di] = invzt_solve_auto(tc.TestData.ion, 0.1, [4.8 0 0], tc.TestData.lat, struct());
verifyEqual(tc, phase, 2);
verifyTrue(tc, pt.converged && pt.crit > 0);
verifyFalse(tc, di.ordered.attempted);    % PM valid -> ordered leg never consulted
end

function test_above_Bc_is_phase2(tc)
[pt, phase] = invzt_solve_auto(tc.TestData.ion, 0.1, [5.5 0 0], tc.TestData.lat, struct());
verifyEqual(tc, phase, 2);
verifyTrue(tc, pt.converged && pt.crit > 0);
end

function test_spurious_pm_rejected_deterministic(tc)
% Review P1-2: force the ordered leg to early-return DETERMINISTICALLY
% (m_tol = Inf -> |m0| <= m_tol always -> is_ordered = false), leaving the PM
% leg untouched. At 4.0 T the PM leg CONVERGES with crit = -0.186 (measured):
% an implementation missing the crit > 0 rule returns phase 2 here and fails.
[~, phase, di] = invzt_solve_auto(tc.TestData.ion, 0.1, [4.0 0 0], tc.TestData.lat, ...
    struct('m_tol', Inf));
verifyEqual(tc, phase, 0);
verifyTrue(tc, di.para.attempted && ~di.para.accepted);
verifyTrue(tc, di.para.converged);                   % the spurious point DID converge...
verifyLessThan(tc, di.para.crit, 0);                 % ...and was rejected on stability
verifyTrue(tc, di.ordered.attempted && ~di.ordered.accepted);
end

function test_boundary_bracket(tc)
% Phase boundary consistency (P0-1): scanning across the measured crit bracket,
% no PM label at or below 4.65 T, PM at and above 4.75 T, and no re-entrant
% ordered label above the first PM field.
ion = tc.TestData.ion;  lat = tc.TestData.lat;
Bs = 4.55:0.05:4.85;  ph = zeros(size(Bs));
for k = 1:numel(Bs)
    [~, ph(k)] = invzt_solve_auto(ion, 0.1, [Bs(k) 0 0], lat, struct());
end
fprintf('boundary scan: B = [%s], phase = [%s]\n', num2str(Bs, '%.2f '), num2str(ph, '%d '));
verifyTrue(tc, ~any(ph(Bs <= 4.651) == 2));          % no PM at/below the lower bracket edge
verifyTrue(tc, all(ph(Bs >= 4.749) == 2));           % PM at/above the upper bracket edge
i2 = find(ph == 2, 1);
verifyTrue(tc, isempty(i2) || ~any(ph(i2:end) == 1));  % monotone: no FM re-entry
end

function test_longitudinal_rejected_at_entry(tc)
% Round-2 P1-5: the Bz guard is at DISPATCHER ENTRY -- it must fire even where
% the PM leg would be valid (5.5 T), not only on the ordered side. The round-1
% draft guarded only inside the ordered leg, so a PM-valid point silently
% accepted a longitudinal field.
verifyError(tc, @() invzt_solve_auto(tc.TestData.ion, 0.1, [5.5 0 1e-6], ...
    tc.TestData.lat, struct()), 'invzt:orderedLongitudinal');
verifyError(tc, @() invzt_solve_auto(tc.TestData.ion, 0.1, [3.0 0 0.1], ...
    tc.TestData.lat, struct()), 'invzt:orderedLongitudinal');
end

function test_invalid_mode_rejected(tc)
% a3d/a3 never classify the phase; a non-a1 mode is a caller error, not phase 0.
verifyError(tc, @() invzt_solve_auto(tc.TestData.ion, 0.1, [3.0 0 0], ...
    tc.TestData.lat, struct('mode', 'a3d')), 'invzt:autoMode');
end

function test_config_errors_rethrown_not_masked(tc)
% Round-2 P1-5: the catch allowlist must NOT convert configuration errors into
% phase 0. A malformed lattice struct raises a non-allowlisted identifier,
% which must propagate out of the dispatcher.
badlat = struct('nothing', 1);           % no Jt/JtGamma/info fields
try
    invzt_solve_auto(tc.TestData.ion, 0.1, [3.0 0 0], badlat, struct());
    verifyTrue(tc, false, 'malformed lattice must throw, not return phase 0');
catch err
    verifyFalse(tc, ismember(err.identifier, ...
        {'invz:degenerateDoublet', 'invzt:a1ZeroField', 'invz:orderedPhase'}));
end
end
