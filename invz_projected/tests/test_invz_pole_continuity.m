function tests = test_invz_pole_continuity
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% ---- 1. smooth removable crossing: y perfectly linear through a sign-change bracket. ----
function test_smooth_removable_crossing_is_ok(testCase)
h = 1:10;
d = h - 5.5;                     % sign change between h=5 (d=-0.5) and h=6 (d=+0.5)
y = h;                            % perfectly linear -> both extrapolants exact
res = invz_pole_continuity(h, d, y, 1e-3);
verifyEqual(testCase, res.status, 'ok');
verifyEqual(testCase, res.max_jump, 0, 'AbsTol', 1e-12);
verifyEqual(testCase, numel(res.crossings), 1);
verifyEqual(testCase, res.crossings(1).kind, 'sign_change');
verifyTrue(testCase, res.crossings(1).resolved);
verifyEqual(testCase, res.crossings(1).h_cross, 5.5, 'AbsTol', 1e-12);
end

% ---- 2. injected jump: y is genuinely discontinuous right at the crossing -> must NOT be 'ok'. ----
function test_injected_jump_is_detected(testCase)
h = 1:10;
d = h - 5.5;
y = h;  y(6:end) = y(6:end) + 100;      % a real discontinuity straddling the crossing
res = invz_pole_continuity(h, d, y, 1e-3);
verifyEqual(testCase, res.status, 'jump_exceeded');
verifyGreaterThan(testCase, res.max_jump, 1e-3);
end

% ---- 3. no crossing: d never changes sign and is never exactly zero. ----
function test_no_crossing(testCase)
h = 1:10;
d = h;                      % strictly positive throughout
y = sin(h);                 % irrelevant -- never evaluated
res = invz_pole_continuity(h, d, y, 1e-3);
verifyEqual(testCase, res.status, 'no_crossing');
verifyEmpty(testCase, res.crossings);
verifyTrue(testCase, isnan(res.max_jump));
end

% ---- 4. exact-zero node, value agrees with both extrapolants -> 'ok'. ----
function test_exact_zero_node_agrees(testCase)
h = 1:9;
d = h - 5;                  % d(5) == 0 exactly
y = h;                       % linear -> the node value agrees exactly with both extrapolants
res = invz_pole_continuity(h, d, y, 1e-3);
verifyEqual(testCase, res.status, 'ok');
verifyEqual(testCase, numel(res.crossings), 1);
c = res.crossings(1);
verifyEqual(testCase, c.kind, 'exact_zero');
verifyEqual(testCase, c.h_cross, 5);
verifyTrue(testCase, c.resolved);
verifyEqual(testCase, c.jump_lr, 0, 'AbsTol', 1e-12);
verifyEqual(testCase, c.jump_node_L, 0, 'AbsTol', 1e-12);
verifyEqual(testCase, c.jump_node_R, 0, 'AbsTol', 1e-12);
end

% ---- 5. exact-zero node, FINITE but DISPLACED value -- the case a finiteness-only rule would
% have passed (prereg SS5). yL/yR stay mutually consistent (jump_lr ~ 0); only the node-vs-
% extrapolant checks catch the displacement. THIS TEST MUST NOT REPORT 'ok'. ----
function test_exact_zero_node_displaced_must_fail(testCase)
h = 1:9;
d = h - 5;                  % d(5) == 0 exactly
y = h;
y(5) = y(5) + 100;          % finite, but arbitrarily wrong -- NOT what a finiteness-only
                             % check would catch (isfinite(105) is true)
res = invz_pole_continuity(h, d, y, 1e-3);
verifyNotEqual(testCase, res.status, 'ok');
verifyEqual(testCase, res.status, 'jump_exceeded');
c = res.crossings(1);
verifyTrue(testCase, c.resolved);
verifyEqual(testCase, c.jump_lr, 0, 'AbsTol', 1e-12);          % extrapolants alone agree...
verifyGreaterThan(testCase, c.jump_node_L, 1e-3);               % ...but the node disagrees with both
verifyGreaterThan(testCase, c.jump_node_R, 1e-3);
verifyGreaterThan(testCase, res.max_jump, 1e-3);
end

% ---- 6. insufficient side coverage -> 'unresolved', for BOTH a sign-change crossing and an
% exact-zero crossing missing their left pair. Never silently omitted (prereg SS5). ----
function test_insufficient_coverage_is_unresolved(testCase)
% sign-change crossing at the very first bracket: only 1 (not 2) nodes to its left.
h1 = 1:5;  d1 = [-1 0.4 1 2 3];  y1 = h1;
res1 = invz_pole_continuity(h1, d1, y1, 1e-3);
verifyEqual(testCase, res1.status, 'unresolved');
verifyEqual(testCase, numel(res1.crossings), 1);
verifyFalse(testCase, res1.crossings(1).resolved);

% exact-zero crossing at the second node: only 1 (not 2) nodes to its left.
h2 = 1:6;  d2 = [-2 0 1 2 3 4];  y2 = h2;
res2 = invz_pole_continuity(h2, d2, y2, 1e-3);
verifyEqual(testCase, res2.status, 'unresolved');
verifyEqual(testCase, numel(res2.crossings), 1, ...
    'the zero neighbour of an exact-zero node must not ALSO be counted as a sign-change crossing');
verifyFalse(testCase, res2.crossings(1).resolved);

% insufficient coverage on the RIGHT instead of the left (symmetric edge case).
h3 = 1:3;  d3 = [-2 -1 0.5];  y3 = h3;
res3 = invz_pole_continuity(h3, d3, y3, 1e-3);
verifyEqual(testCase, res3.status, 'unresolved');
verifyFalse(testCase, res3.crossings(1).resolved);
end

% ---- Multiple independent crossings: one 'ok', one exceeding -> overall must NOT be 'ok'.
% d has EXACTLY two sign changes (at i=4 and i=11, both resolved: verified by hand, every
% other adjacent pair is same-sign). y is linear (hence crossing-A-exact) except for a +50
% shift from index 12 onward, which lands only inside crossing B's right-side extrapolation. ----
function test_multiple_crossings_worst_case_wins(testCase)
h = 1:16;
d = [-4 -3 -2 -1 1 2 3 4 3 2 1 -1 -2 -3 -4 -5];
y = h;  y(12:end) = y(12:end) + 50;
res = invz_pole_continuity(h, d, y, 1e-3);
verifyEqual(testCase, numel(res.crossings), 2);
verifyTrue(testCase, all([res.crossings.resolved]));
verifyEqual(testCase, res.crossings(1).jump_lr, 0, 'AbsTol', 1e-12, 'crossing A is untouched by the shift');
verifyGreaterThan(testCase, res.crossings(2).jump_lr, 1e-3, 'crossing B straddles the shift');
verifyEqual(testCase, res.status, 'jump_exceeded');
end

% ---- Input-contract guards (finite/aligned precondition; no interp1-across-duplicates path). ----
function test_rejects_nonfinite_input(testCase)
h = 1:5;  d = [-1 -0.5 0.5 1 NaN];  y = h;
verifyError(testCase, @() invz_pole_continuity(h, d, y, 1e-3), 'invz:poleContinuityInput');
end

function test_rejects_non_increasing_h(testCase)
h = [1 2 2 3 4];  d = [-1 -0.5 0.5 1 2];  y = h;
verifyError(testCase, @() invz_pole_continuity(h, d, y, 1e-3), 'invz:poleContinuityInput');
end

function test_rejects_mismatched_lengths(testCase)
h = 1:5;  d = [-1 -0.5 0.5 1];  y = 1:5;
verifyError(testCase, @() invz_pole_continuity(h, d, y, 1e-3), 'invz:poleContinuityInput');
end

function test_rejects_nonpositive_tol(testCase)
h = 1:5;  d = [-1 -0.5 0.5 1 2];  y = h;
verifyError(testCase, @() invz_pole_continuity(h, d, y, 0), 'invz:poleContinuityInput');
verifyError(testCase, @() invz_pole_continuity(h, d, y, -1e-3), 'invz:poleContinuityInput');
end
