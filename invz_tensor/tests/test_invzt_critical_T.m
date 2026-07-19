function tests = test_invzt_critical_T
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion engine
end

% ---------- fast: pure crossing/slide policy (synthetic votes, no solves) -----

function test_pick_crossing_policy(testCase)
% Simple bracket
[act, ka, kb, nc] = invzt_tc_pick([-1 1]);
verifyEqual(testCase, act, 'bracket');
verifyEqual(testCase, [ka kb], [1 2]);
verifyEqual(testCase, nc, 1);
% Multiple crossings: HIGHEST ordered->para pair wins, re-entrance counted
[act, ka, kb, nc] = invzt_tc_pick([1 -1 1]);
verifyEqual(testCase, act, 'bracket');
verifyEqual(testCase, [ka kb], [2 3]);
verifyEqual(testCase, nc, 2);
% All-PM window: boundary below
act = invzt_tc_pick([1 2 3]);
verifyEqual(testCase, act, 'down');
% All-ordered window: boundary above
act = invzt_tc_pick([-1 -2]);
verifyEqual(testCase, act, 'up');
% Re-entrant LOWER leg (para below, ordered above): the physical high-T side
% is paramagnetic, so the highest ordered->para crossing is ABOVE -- must
% move up, not give up (the projected classifier's inherited gap).
act = invzt_tc_pick([1 1 -1 -1]);
verifyEqual(testCase, act, 'up');
% Singletons
verifyEqual(testCase, invzt_tc_pick(1),  'down');
verifyEqual(testCase, invzt_tc_pick(-1), 'up');
end

function test_pick_exact_zero(testCase)
% A sampled crit == 0 IS the boundary (the projected classifier mis-reads
% [-1, 0, 1] as two sign changes with no returnable crossing).
[act, ka, ~, nc] = invzt_tc_pick([-1 0 1]);
verifyEqual(testCase, act, 'zero');
verifyEqual(testCase, ka, 2);
verifyEqual(testCase, nc, 1);                    % one boundary, not two
% Zero RUN counted once
[act, ka, ~, nc] = invzt_tc_pick([-1 0 0 1]);
verifyEqual(testCase, act, 'zero');
verifyEqual(testCase, ka, 3);
verifyEqual(testCase, nc, 1);
% Zero with ordered ABOVE it: the true highest boundary is above -> up
act = invzt_tc_pick([1 0 -1]);
verifyEqual(testCase, act, 'up');
% Zero ABOVE a strict crossing: the zero is the higher root
[act, ka, ~, nc] = invzt_tc_pick([-1 1 -1 0 1]);
verifyEqual(testCase, act, 'zero');
verifyEqual(testCase, ka, 4);
verifyEqual(testCase, nc, 3);
% Lone exactly-critical voter
[act, ka] = invzt_tc_pick(0);
verifyEqual(testCase, act, 'zero');
verifyEqual(testCase, ka, 1);
end

% ---------- fast: input-contract validation (guards fire pre-compute) --------

function test_anchor_and_window_validation(testCase)
% lat is a dummy struct throughout -- every guard fires before any compute.
% Malformed-opts fixtures are built by FIELD ASSIGNMENT: the struct('f',[a b])
% constructor would create a struct ARRAY, not a field holding a vector.
% First, the hardened safe formatter itself: raw mat2str throws on structs
% (masking the intended error id); invzt_str must return a placeholder.
verifyEqual(testCase, invzt_str(struct()), '<1x1 struct>');
ion = invz_ion();
% invzt:tcAnchor -- adaptive mode needs a USABLE Tc0
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), struct()), 'invzt:tcAnchor');
o = struct(); o.Tc0 = [1 2];      % vector: must NOT surface MATLAB:nonLogicalConditional
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcAnchor');
o = struct(); o.Tc0 = 1 + 1i;     % complex: isfinite(1+1i) is true -- needs isreal
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcAnchor');
o = struct(); o.Tc0 = -1;
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcAnchor');
o = struct(); o.Tc0 = 0.01;       % at/below the 0.02 K solve floor
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcAnchor');
% invzt:tcWindow -- explicit window malformed / nonnumeric / below-floor
o = struct(); o.window = [1.8 1.0];
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcWindow');
o = struct(); o.window = struct();  % nonnumeric: message building must not mat2str-throw
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcWindow');
o = struct(); o.window = [0.005 0.015];
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcWindow');
% invzt:tcOpts -- numerical controls must be finite positive real scalars
o = struct(); o.Tc0 = 1.5; o.gridstep = 0;
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcOpts');
o = struct(); o.Tc0 = 1.5; o.tol = -1;
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcOpts');
o = struct(); o.Tc0 = 1.5; o.width = Inf;
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcOpts');
% CHEAP bracket diagnostics (third review T1/T3): with B = 0 every sample is
% invalidated by the PRE-LATTICE invzt:a1ZeroField guard (absorbed by
% invzt_crit_at), so whole windows run in milliseconds on a dummy lat.
% (a) Hard window, one pass: the invzt:bracket MESSAGE must report EXACTLY
% the user's window -- never a slide-mutated pair. NB explicit try/catch:
% verifyError does NOT return the caught MException (its outputs are the
% evaluated function's outputs -- <missing> on a throw).
o = struct(); o.window = [1.1 1.4];
didThrow = false;
try
    invzt_critical_T(ion, 0, struct(), o);
catch ME
    didThrow = true;
    verifyEqual(testCase, ME.identifier, 'invzt:bracket');
    verifySubstring(testCase, ME.message, '[1.100, 1.400]');
end
verifyTrue(testCase, didThrow, 'expected invzt:bracket was not thrown');
% (b) Adaptive no-valid grow path terminates at the 0.02 K floor (R3) --
% Tc0 just above the floor makes the first window floor-clamped already, so
% the scan must stop instead of re-sampling the identical grid.
o = struct(); o.Tc0 = 0.05;
didThrow = false;
try
    invzt_critical_T(ion, 0, struct(), o);
catch ME
    didThrow = true;
    verifyEqual(testCase, ME.identifier, 'invzt:bracket');
end
verifyTrue(testCase, didThrow, 'expected invzt:bracket was not thrown');
end

% ---------- fast: sampler contract (pre-lattice guard paths, dummy lat) ------

function test_crit_at_contract(testCase)
% Physics signals absorbed as invalid samples; config errors rethrow. Both
% paths fire inside invzt_solve_point BEFORE any lattice access, so a dummy
% lat keeps these at millisecond cost.
ion = invz_ion();
[c, ok, pt0] = invzt_crit_at(ion, 1.6, [0 0 0], struct(), struct());  % zero transverse field
verifyFalse(testCase, ok);
verifyTrue(testCase, isnan(c));                                   % invzt:a1ZeroField absorbed
verifyEqual(testCase, pt0, struct());                             % absorbed sample: empty pt (E1)
o = struct(); o.mode = 'bogus';
verifyError(testCase, @() invzt_crit_at(ion, 1.6, [2 0 0], struct(), o), ...
    'invzt:mode');                                                % misconfiguration rethrows
end

% ---------- fast: the moved invz_common helper, directly ----------------------

function test_refine_crossing_helper(testCase)
% Regula-falsi between converged bracket ends; non-converged interior ->
% midpoint retry; total interior failure -> linear-interpolation fallback.
froot = 1.3;
f1 = @(x) deal(2*(x - froot), true);
bx = invz_refine_crossing(f1, 1.0, 2*(1.0 - froot), 1.8, 2*(1.8 - froot), 1e-4);
verifyEqual(testCase, bx, froot, 'AbsTol', 1e-3);
% Interior dead zone straddling the root: midpoint retries + the final
% linear interpolation must still land the root (f is linear, so the
% fallback is exact up to the surviving bracket).
f2 = @(x) deal(2*(x - froot), ~(x > 1.25 && x < 1.35));
bx = invz_refine_crossing(f2, 1.0, 2*(1.0 - froot), 1.8, 2*(1.8 - froot), 1e-4);
verifyEqual(testCase, bx, froot, 'AbsTol', 0.06);
end

% ---------- slow: integration gates (8^3, warm caches) ------------------------

function test_tcut_matches_field_cut_slow(testCase)
% Crossing consistency (the projected branch's own validation pattern): at
% T* = 1.4 K find B* with the validated field-cut finder, then the T-cut at
% B* must return T*. Also proves the Sigma_seed strip (a length-MATCHED,
% poisonous-VALUED seed is passed -- unstripped, the matching endpoint
% consumes the NaNs and turns invalid) and the out-struct contract.
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
g   = invzt_qgrid(8, 'halfopen');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
o   = struct('odd', false);                             % lighter; round-trip is self-consistent
Bstar = invzt_critical(ion, 1.4, lat, [0.05 6], o);
oT = o;
oT.window = [1.0 1.8];                                  % HARD window containing T* = 1.4
% Deliberately POISONOUS seed, length-MATCHED to the top endpoint's Matsubara
% count (Ecut default 40): a wrong-length seed would be silently ignored by
% invzt_solve_point's numel(Sigma_seed) == nwn guard and prove nothing about
% the strip (second review R1). Unstripped, the 1.8 K endpoint consumes the
% NaN seed and turns invalid; stripped, it starts from the zero seed and
% stays a converged PM voter.
wnHi = invz_matsubara(oT.window(2), 40);
oT.Sigma_seed = nan(numel(wnHi), 1);
[tc, out] = invzt_critical_T(ion, Bstar, lat, oT);
verifyEqual(testCase, tc, 1.4, 'AbsTol', 0.05);
verifyTrue(testCase, all(isfield(out, {'Tg', 'c', 'ok', 'window', 'ncross', 'B'})));
verifyEqual(testCase, numel(out.c), numel(out.Tg));
verifyEqual(testCase, numel(out.ok), numel(out.Tg));
verifyGreaterThanOrEqual(testCase, out.ncross, 1);
verifyTrue(testCase, out.ok(end), ...
    'top-endpoint sample invalid: the poisonous Sigma_seed was not stripped');
% sigma_floor threading + rejection (second review R6): an impossible floor
% must invalidate a converged point whose crit is still finite (one solve).
oF = o;  oF.sigma_floor = Inf;
[cF, okF] = invzt_crit_at(ion, 1.8, [Bstar 0 0], lat, oF);
verifyTrue(testCase, isfinite(cF));
verifyFalse(testCase, okF);
fprintf('T-cut round-trip: Bc(1.4 K) = %.4f T -> Tc = %.4f K (ncross = %d)\n', ...
        Bstar, tc, out.ncross);
end

function test_tcut_odd_on_slow(testCase)
% The capability the projected T-cut lacks: bracketing with ODD on (the
% tensor A1 solver converges metastable PM points inside the ordered phase,
% so valid votes exist on both sides of the crossing -- ODD-LOG T2 records
% that the projected finder cannot bracket there at all).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
g   = invzt_qgrid(8, 'halfopen');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
oT = struct('odd', true);
oT.window = [1.2 1.7];
tc = invzt_critical_T(ion, 1.5, lat, oT);
verifyGreaterThan(testCase, tc, 1.35);   % loose physical band at 8^3: the odd-on
verifyLessThan(testCase, tc, 1.60);      % boundary has Bc(1.4 K) = 1.916 T > 1.5 T
fprintf('odd-on T-cut: Tc(1.5 T) = %.4f K\n', tc);
% Hard window entirely inside the ordered phase: ONE pass, no sliding ->
% invzt:bracket with the widen-the-window message (F4's contract end-to-end).
oT2 = struct('odd', true);
oT2.window = [0.25 0.45];
% Physical deep-ordered hard window: one pass, no sliding, and the error
% reports the window the user asked for (R2). Explicit try/catch: verifyError
% does NOT return the caught MException (third review T1).
didThrow = false;
try
    invzt_critical_T(ion, 1.5, lat, oT2);
catch ME
    didThrow = true;
    verifyEqual(testCase, ME.identifier, 'invzt:bracket');
    verifySubstring(testCase, ME.message, '[0.250, 0.450]');
end
verifyTrue(testCase, didThrow, 'expected invzt:bracket was not thrown');
end
