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
end
