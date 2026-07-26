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
% NOTE: this is an algebraic restatement of the r-limit assertion immediately above, not an
% independent third analytic-limit check. crit is DEFINED as r + J0eff*G0bare, so substituting
% the just-asserted out.r = -out.G0bare*K0 gives out.G0bare*(J0eff-K0) for ANY J0eff -- this
% cannot fail unless the out.r assertion above it already failed. Kept because it documents the
% derived `crit` limit used elsewhere in the solver; genuine continuity of `crit` THROUGH the
% crossing (as opposed to only at this Inf limit) is covered separately by
% test_r_and_crit_finite_through_the_crossing.
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

% Review fix (task-6 round 2): the two out.Gtil0/out.r ARRANGEMENTS -- legacy seven-argument
% (Gtil0 = Gstat/(1-K0*Gstat); r = G0bare/Gtil0) and opts.stable_form=true's exact reassociation
% (invGtil0 = 1/Gstat - K0; Gtil0 = 1/invGtil0; r = G0bare*invGtil0) -- must agree to float noise
% on IDENTICAL caller inputs, away from the Gstat pole. None of the tests above call both
% branches on the same inputs and compare them directly: the crossing/infinite-limit tests above
% use stable_form=true only, and test_agrees_with_old_expression_away_from_pole /
% test_g9_seven_arg_call_bit_identical_to_historical_arithmetic use the legacy 7-arg call only,
% each checked against a hand-recomputed reference rather than against the OTHER branch's actual
% output. This test is the permanent regression guard for the moment Tasks 8/10/12 begin wiring
% stable_form=true into production: it proves the two branches are still the same expression,
% not merely that each is separately finite/continuous.
%
% Tolerance: an explicit ULP bound, abs(a-b) <= 4*eps(abs(a)), not a bare RelTol. 4 ulp is the
% implementation plan's own stop-and-report threshold -- "If any value shows more than ~4 ulp,
% stop and report -- that would mean the two forms are not the same expression"
% (docs/superpowers/plans/2026-07-25-invzp-strict-static-medium.md, Task 6 Step 1). It is NOT
% stated in invz_gstat_ordered.m's docstring, which records only the measured value.
% The bound is not an arbitrary loosening: a separate out-of-band 16-point sweep spanning
% near-pole and sign-crossing Gstat values found a MAX of 1.0 ulp (6 of those 16 differing from
% bitwise-exact at all), and the 6-point subset committed below reproduces that same 1.0 ulp
% maximum -- both comfortably inside 4 ulp. A control run with K0 perturbed by 1e-7 relative
% produced ~2.3e6 ulp, proving the comparison discriminates real divergence rather than being
% vacuously loose. A looser tolerance (a RelTol of 1e-9 or more) would stop catching that.
function test_stable_form_matches_legacy_branch_on_identical_inputs(testCase)
[tl, beta] = fixture_tl();
K0 = 5e-3;  G0inel0 = -300;  G0el0 = -20;  lam = [0.01; 0.02; 0.005];
% Sigma0 values chosen (via a probe against this exact fixture) to span SMALL- and LARGE-
% magnitude Gstat of BOTH signs, all strictly away from the pole at
% Sigma0 = -1 - K0*G0inel0 = 0.5 (see test_agrees_with_old_expression_away_from_pole above) --
% never AT the pole itself, where the two branches are EXPECTED to differ (that is the entire
% point of the task, and is covered separately by the crossing/infinite-limit tests above):
%   Sigma0=-1e4  -> Gstat ~ -3.5   (small magnitude, negative)
%   Sigma0=-50   -> Gstat ~ +2.4   (small magnitude, positive)
%   Sigma0=50    -> Gstat ~ -9.6   (small magnitude, negative)
%   Sigma0=0.2   -> Gstat ~ +996   (medium/large magnitude, positive)
%   Sigma0=0.499 -> Gstat ~ +3e5   (large magnitude, positive; near-pole from below)
%   Sigma0=0.501 -> Gstat ~ -3e5   (large magnitude, negative; near-pole from above)
Sigma0_list = [-1e4, -50, 50, 0.2, 0.499, 0.501];
for k = 1:numel(Sigma0_list)
    S0 = Sigma0_list(k);
    [Gs_legacy, out_legacy] = invz_gstat_ordered(tl, lam, K0, S0, beta, G0inel0, G0el0);
    [Gs_stable, out_stable] = invz_gstat_ordered(tl, lam, K0, S0, beta, G0inel0, G0el0, ...
                                                  struct('stable_form', true));
    verifyTrue(testCase, isfinite(Gs_legacy) && Gs_legacy ~= 0, ...
        sprintf('fixture must be finite and non-pole at Sigma0=%g (Gstat=%g)', S0, Gs_legacy));

    d_Gstat = abs(Gs_stable - Gs_legacy);
    d_Gtil0 = abs(out_stable.Gtil0 - out_legacy.Gtil0);
    d_r     = abs(out_stable.r - out_legacy.r);
    u_Gstat = d_Gstat / eps(abs(Gs_legacy));
    u_Gtil0 = d_Gtil0 / eps(abs(out_legacy.Gtil0));
    u_r     = d_r / eps(abs(out_legacy.r));
    % No unconditional printing: the suite's output must stay pristine. Every one of these
    % figures is reproduced in the verifyTrue failure messages below if a bound is ever
    % breached, which is the only time they are worth reading.

    verifyTrue(testCase, d_Gstat <= 4*eps(abs(Gs_legacy)), ...
        sprintf(['Gstat itself differs by %.3f ulp at Sigma0=%g (both branches compute Gstat ' ...
                 'via the identical unbranched code, so this should be exactly 0)'], u_Gstat, S0));
    verifyTrue(testCase, d_Gtil0 <= 4*eps(abs(out_legacy.Gtil0)), ...
        sprintf('Gtil0 branches disagree by %.3f ulp (> 4) at Sigma0=%g', u_Gtil0, S0));
    verifyTrue(testCase, d_r <= 4*eps(abs(out_legacy.r)), ...
        sprintf('r branches disagree by %.3f ulp (> 4) at Sigma0=%g', u_r, S0));
end
end
