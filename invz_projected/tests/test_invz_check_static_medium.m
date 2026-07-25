function tests = test_invz_check_static_medium
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_default_is_resummed_and_not_strict(testCase)
sm = invz_check_static_medium(struct());
verifyEqual(testCase, sm.scheme, 'resummed');
verifyFalse(testCase, sm.is_strict);
end

function test_accepts_the_three_canonical_ids(testCase)
for id = {'resummed', 'strict_1z_dyson_ref', 'strict_1z_bare_ref'}
    sm = invz_check_static_medium(struct('static_medium', id{1}));
    verifyEqual(testCase, sm.scheme, id{1});
end
verifyTrue(testCase, invz_check_static_medium(struct('static_medium', ...
    'strict_1z_dyson_ref')).is_strict);
verifyFalse(testCase, invz_check_static_medium(struct('static_medium', 'resummed')).is_strict);
end

function test_unknown_id_throws(testCase)
verifyError(testCase, @() invz_check_static_medium(struct('static_medium', 'strict')), ...
            'invz:staticMedium');
verifyError(testCase, @() invz_check_static_medium(struct('static_medium', 42)), ...
            'invz:staticMedium');
% the rejected self-consistent comparator is NOT a selectable production scheme (spec §B)
verifyError(testCase, @() invz_check_static_medium(struct('static_medium', ...
            'strict_1z_selfconsistent')), 'invz:staticMedium');
end

% The stamp is what makes the two legs unable to disagree (spec §4.2).
function test_stamps_both_leg_structs_identically(testCase)
[sm, eopts, eso] = invz_check_static_medium( ...
    struct('static_medium', 'strict_1z_dyson_ref'), struct('K0', 7), struct('warn', false));
verifyEqual(testCase, eopts.static_medium, sm.scheme);
verifyEqual(testCase, eso.static_medium,   sm.scheme);
verifyEqual(testCase, eopts.ref_margin,    sm.ref_margin);
verifyEqual(testCase, eso.ref_margin,      sm.ref_margin);
verifyEqual(testCase, eopts.K0, 7);              % pre-existing fields preserved
verifyFalse(testCase, eso.warn);
end

% A per-leg value that DISAGREES with the resolved scheme is a conflict: it would split the two
% sectors across different truncation orders.
function test_disagreeing_per_leg_values_are_conflicts(testCase)
verifyError(testCase, @() invz_check_static_medium( ...
    struct('emt', struct('static_medium', 'strict_1z_dyson_ref'))), 'invz:staticMedium');
verifyError(testCase, @() invz_check_static_medium(struct( ...
    'static_medium', 'strict_1z_dyson_ref', ...
    'emt_static', struct('static_medium', 'resummed'))), 'invz:staticMedium');
end

% IDEMPOTENCY: re-validating an ALREADY-STAMPED opts struct must not be a conflict.
% invz_solve_point_ordered forwards its full numerical context (including opts.emt /
% opts.emt_static) into invz_hmf_ordered, which validates again -- so a stamp that agrees with
% the resolved scheme has to be accepted, or the second pass would throw on its own output.
function test_restamping_an_already_stamped_opts_is_idempotent(testCase)
[sm1, eopts, eso] = invz_check_static_medium(struct('static_medium', 'strict_1z_dyson_ref'));
forwarded = struct('static_medium', sm1.scheme, 'emt', eopts, 'emt_static', eso);
[sm2, eopts2, eso2] = invz_check_static_medium(forwarded, eopts, eso);
verifyEqual(testCase, sm2.scheme, sm1.scheme);
verifyEqual(testCase, eopts2.static_medium, sm1.scheme);
verifyEqual(testCase, eso2.static_medium, sm1.scheme);
verifyEqual(testCase, sm2.ref_margin, sm1.ref_margin, 'AbsTol', 0);
end

function test_ref_margin_default_and_override(testCase)
verifyEqual(testCase, invz_check_static_medium(struct()).ref_margin, 1e-6, 'AbsTol', 0);
sm = invz_check_static_medium(struct('ref_margin', 1e-8));
verifyEqual(testCase, sm.ref_margin, 1e-8, 'AbsTol', 0);
end

% §7.5: strict mode + an [nJ,nw] ordered-retarded coupling matrix is rejected, not silently
% flattened (invz_emt_static_ordered.m:43 would average all frequency columns together).
function test_validator_does_not_guess_coupling_shape(testCase)
% Retarded compatibility is decided only after Jnu_flat has actually been resolved. The PM
% leg supports [nJ,nw]; the ordered solver rejects that matrix at its public entry.
sm = invz_check_static_medium(struct('static_medium', 'strict_1z_dyson_ref'));
verifyTrue(testCase, sm.is_strict);

% This validator never receives Jnu_flat/coupling data at all -- it structurally cannot guess
% coupling shape, because coupling shape never reaches it. Confirm: a strict scheme resolves
% with ZERO coupling data supplied (the call below completes without error using only
% opts.static_medium); the returned sm carries EXACTLY scheme/is_strict/ref_margin and no
% coupling-related field; the stamped eopts/eso likewise gain only static_medium/ref_margin and
% no coupling-related field.
% NOTE: the actual [nJ,nw] rejection under strict mode is enforced LATER, at the ordered
% solver's public entry (invz_solve_point_ordered.m / invz_emt_static_ordered.m:43), not here --
% a future reader should not expect this test to cover that.
[sm2, eopts2, eso2] = invz_check_static_medium(struct('static_medium', 'strict_1z_dyson_ref'));
verifyEqual(testCase, sort(fieldnames(sm2)), sort({'scheme'; 'is_strict'; 'ref_margin'}));
verifyEqual(testCase, sort(fieldnames(eopts2)), sort({'static_medium'; 'ref_margin'}));
verifyEqual(testCase, sort(fieldnames(eso2)), sort({'static_medium'; 'ref_margin'}));
end

% Fix 1: a non-struct eopts/eso must escape as this validator's OWN invz:staticMedium
% identifier, not an un-namespaced MATLAB error from the dot-indexed stamping at the bottom of
% the function -- a validator whose whole purpose is policing malformed option structs should
% not itself raise an un-namespaced error for exactly that class of mistake. Omitted/empty
% eopts/eso (exercised by every other test in this file, which pass only opts) must keep
% defaulting to struct() exactly as before; only a SUPPLIED non-struct is an error.
function test_nonstruct_eopts_or_eso_raises_static_medium_error(testCase)
verifyError(testCase, @() invz_check_static_medium(struct(), 5, 7), 'invz:staticMedium');
verifyError(testCase, @() invz_check_static_medium(struct(), 5), 'invz:staticMedium');
verifyError(testCase, @() invz_check_static_medium(struct(), struct(), 7), 'invz:staticMedium');
end

% Fix 3: conflict detection must also run on the EXPLICIT eopts/eso positional arguments, not
% just on nested opts.emt/opts.emt_static (test_disagreeing_per_leg_values_are_conflicts
% above) -- without this, a disagreeing eopts/eso would simply be clobbered by the final stamp
% and the conflict hidden.
function test_disagreeing_explicit_eopts_or_eso_are_conflicts(testCase)
% (1) eopts (2nd positional arg) disagrees with a resolved NON-default scheme.
verifyError(testCase, @() invz_check_static_medium( ...
    struct('static_medium', 'strict_1z_dyson_ref'), ...
    struct('static_medium', 'resummed')), 'invz:staticMedium');
% (2) eso (3rd positional arg) disagrees with a resolved NON-default scheme.
verifyError(testCase, @() invz_check_static_medium( ...
    struct('static_medium', 'strict_1z_dyson_ref'), struct(), ...
    struct('static_medium', 'strict_1z_bare_ref')), 'invz:staticMedium');
% (3) disagreement present when the scheme resolves to the 'resummed' DEFAULT (opts carries no
% static_medium field at all).
verifyError(testCase, @() invz_check_static_medium( ...
    struct(), struct('static_medium', 'strict_1z_dyson_ref')), 'invz:staticMedium');

% Idempotent counterpart: an eopts/eso that AGREES with the resolved scheme must be ACCEPTED,
% not thrown -- this is what lets invz_solve_point_ordered forward its context into
% invz_hmf_ordered for a second validation pass without the validator throwing on its own
% output.
[sm, eopts, eso] = invz_check_static_medium(struct('static_medium', 'strict_1z_dyson_ref'), ...
    struct('static_medium', 'strict_1z_dyson_ref'), struct('static_medium', 'strict_1z_dyson_ref'));
verifyEqual(testCase, eopts.static_medium, sm.scheme);
verifyEqual(testCase, eso.static_medium, sm.scheme);
end
