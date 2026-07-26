function tests = test_invz_gstat_removable_pole
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% A two-level struct whose hybrid denominator 1 + Sigma0 + K0*G0inel0 can be driven to zero,
% so Gstat sweeps through its own pole. m ~= 0 so the elastic xi term is live (the m = 0 case
% cannot pole this way).
function [tl, beta] = fixture_tl()
beta = 1/(0.0862*0.31);                     % kB*T at T = 0.31 K, meV
tl = struct('Delta', 0.02, 'M2', 0.8, 'm', 0.6, 'n01', tanh(0.02*beta/2), 'g0', 0.5);
tl.g0 = 2*tl.n01/tl.Delta;
end

% G17: r, Gtil0 and crit are finite and continuous through the Gstat crossing, and match the
% analytic limits. Sigma0 is chosen so the hybrid denominator crosses zero.
function test_r_and_crit_finite_through_the_crossing(testCase)
[tl, beta] = fixture_tl();
K0 = 5e-3;  G0inel0 = -300;  G0el0 = -20;  J0eff = 6.42444e-3;
lam = [0.01; 0.02; 0.005];
% denominator 1 + Sigma0 + K0*G0inel0 = 0 at Sigma0 = -1 - K0*G0inel0
S0 = -1 - K0*G0inel0;
eps_list = [-1e-6 -1e-9 0 1e-9 1e-6];
rv = nan(size(eps_list));  cv = nan(size(eps_list));
for k = 1:numel(eps_list)
[Gs, out] = invz_gstat_ordered(tl, lam, K0, S0 + eps_list(k), beta, G0inel0, G0el0, ...
                               struct('stable_form', true));
    verifyTrue(testCase, isfinite(out.r), ...
        sprintf('r must be finite at eps=%g (Gstat=%g)', eps_list(k), Gs));
    verifyTrue(testCase, isfinite(out.Gtil0), sprintf('Gtil0 finite at eps=%g', eps_list(k)));
    rv(k) = out.r;
    cv(k) = out.r + J0eff*out.G0bare;                                  % crit
end
% continuity: the spread across the crossing is small relative to |r| itself
verifyLessThan(testCase, (max(rv) - min(rv))/max(abs(rv)), 1e-3);
verifyLessThan(testCase, (max(cv) - min(cv))/max(abs(cv)), 1e-3);
end

% The analytic limits, hit exactly by driving Gstat to +-Inf. Under the OLD arrangement these
% are NaN, so this test is the direct guard on the reassociation.
function test_analytic_limits_at_infinite_gstat(testCase)
[tl, beta] = fixture_tl();
K0 = 5e-3;  G0inel0 = -300;  G0el0 = 0;   % G0el0 = 0 => Gstat is exactly the hybrid term
J0eff = 6.42444e-3;                       % same constant as the sibling crossing test; the
                                           % crit assertion below is J0eff-invariant (cancels
                                           % against the r = -G0bare*K0 identity above it)
lam = [0.01; 0.02; 0.005];
S0inf = -1 - K0*G0inel0;                  % denominator exactly 0 => Gstat = +-Inf
[Gs, out] = invz_gstat_ordered(tl, lam, K0, S0inf, beta, G0inel0, G0el0, ...
                               struct('stable_form', true));
verifyTrue(testCase, isinf(Gs), 'fixture must actually drive Gstat to Inf');
verifyEqual(testCase, out.Gtil0, -1/K0, 'RelTol', 1e-12);
verifyEqual(testCase, out.r, -out.G0bare*K0, 'RelTol', 1e-12);
verifyEqual(testCase, out.r + J0eff*out.G0bare, ...
            out.G0bare*(J0eff-K0), 'RelTol', 1e-12);
verifyTrue(testCase, isfinite(out.r));
end

% Away from any pole the new arrangement must agree with the old expression to float noise.
function test_agrees_with_old_expression_away_from_pole(testCase)
[tl, beta] = fixture_tl();
K0 = 5e-3;  G0inel0 = -300;  G0el0 = -20;  lam = [0.01; 0.02; 0.005];
for S0 = [0, 0.1, -0.3]
    [Gs, out] = invz_gstat_ordered(tl, lam, K0, S0, beta, G0inel0, G0el0);
    old_Gtil0 = Gs/(1 - K0*Gs);
    verifyEqual(testCase, out.Gtil0, old_Gtil0, 'RelTol', 1e-12);
    verifyEqual(testCase, out.r, out.G0bare/old_Gtil0, 'RelTol', 1e-12);
end
end

% The m = 0 pinned identity is untouched by the arrangement (it holds for ANY K0).
function test_m_zero_identity_still_r_equals_one_plus_sigma(testCase)
[tl, beta] = fixture_tl();
tl.m = 0;
lam = [0.01; 0.02; 0.005];
for K0 = [0, 1e-3, 5e-3]
    [~, out] = invz_gstat_ordered(tl, lam, K0, 0.25, beta, -300, 0);
    verifyEqual(testCase, out.r, 1 + 0.25, 'RelTol', 1e-12);
end
end

% G9 (bitwise-identity requirement): a seven-argument call (opts omitted entirely, not just
% stable_form=false) must reproduce the PRE-TASK-6 arithmetic bit-for-bit. Unlike
% test_agrees_with_old_expression_away_from_pole above (RelTol 1e-12, and its "old" value is
% recomputed FROM the function's own returned Gstat), this reconstructs the historical
% four-line computation independently from the SAME raw caller inputs -- mirroring
% invz_gstat_ordered.m exactly as it stood before this task, same operation order throughout
% -- so it is a genuine before/after check, not self-consistency, and asserted with AbsTol 0
% (exact equality), not a tolerance.
function test_g9_seven_arg_call_bit_identical_to_historical_arithmetic(testCase)
[tl, beta] = fixture_tl();
K0 = 5e-3;  Sigma0 = 0.2;  G0inel0 = -300;  G0el0 = -20;  lam = [0.01; 0.02; 0.005];
[Gstat, out] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0);   % 7-arg, no opts

% Historical formula (invz_gstat_ordered.m pre-task-6), transcribed independently here rather
% than reused from the current source, same operation order as the original inline form:
xi_ref = (1 + tanh(tl.m^2*tl.n01^2*beta*K0 - tl.M2*beta*lam(1))) / ...
         (1 + (4*tl.n01^2*K0*tl.g0 + 2*lam(2) + tl.g0*lam(1))*tl.M2/tl.n01^2);
Gstat_ref  = G0inel0/(1 + Sigma0 + K0*G0inel0) + xi_ref*G0el0;
G0bare_ref = G0inel0 + G0el0;
Gtil0_ref  = Gstat_ref/(1 - K0*Gstat_ref);
r_ref      = G0bare_ref/Gtil0_ref;

verifyEqual(testCase, Gstat, Gstat_ref, 'AbsTol', 0);
verifyEqual(testCase, out.Gtil0, Gtil0_ref, 'AbsTol', 0);
verifyEqual(testCase, out.r, r_ref, 'AbsTol', 0);
end
