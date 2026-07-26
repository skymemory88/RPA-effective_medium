function tests = test_invz_spectra_map_phase_reasons
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function o = base_opts()
o = struct('grid', [4 4 4], 'dpRng', 12, 'cache', false, 'verbose', false);
end
function w = wgrid()
w = (0.02:0.04:0.42).';
end

% Provenance is scalar and mandatory.
function test_scheme_provenance_present(testCase)
if isempty(getenv('INVZ_SLOW')), assumeFail(testCase, 'set INVZ_SLOW=1'); end
ion = invz_ion();
S = invz_spectra_map(ion, 0.31, 8, wgrid(), base_opts());
verifyEqual(testCase, S.static_medium, 'resummed');
end

% Every masked column carries a reason. This is the honesty gate: no phase_1z = 0 without one.
function test_no_masked_column_without_a_reason(testCase)
if isempty(getenv('INVZ_SLOW')), assumeFail(testCase, 'set INVZ_SLOW=1'); end
ion = invz_ion();
S = invz_spectra_map(ion, 0.31, [0 3 8], wgrid(), base_opts());
verifyEqual(testCase, numel(S.phase_1z_reason), numel(S.phase_1z));
for k = 1:numel(S.phase_1z)
    verifyFalse(testCase, isempty(S.phase_1z_reason{k}), sprintf('column %d has no reason', k));
    if S.phase_1z(k) == 0
        verifyTrue(testCase, any(strcmp(S.phase_1z_reason{k}, ...
            {'unstable_endpoint', 'medium_out_of_domain', 'degenerate_doublet', ...
             'solver_failed', 'pm_probe_unknown', 'boundary_indeterminate', ...
             'not_attempted_longitudinal', 'bare_not_ordered'})), S.phase_1z_reason{k});
    end
end
end

% B = 0 is labelled, never a silent proxy and never "solver_failed".
function test_zero_field_column_is_labelled_degenerate(testCase)
if isempty(getenv('INVZ_SLOW')), assumeFail(testCase, 'set INVZ_SLOW=1'); end
ion = invz_ion();
o = base_opts();  o.ordered_1z = 'jensen';
S = invz_spectra_map(ion, 0.31, [0 8], wgrid(), o);
verifyTrue(testCase, any(strcmp(S.phase_1z_reason{1}, ...
    {'degenerate_doublet', 'bare_not_ordered', 'pm'})));
verifyNotEqual(testCase, S.phase_1z_reason{1}, 'solver_failed');
end

% The three-way rule is GATED: under 'resummed' the historical dispatch is byte-preserved.
function test_dispatcher_is_gated_on_strict(testCase)
if isempty(getenv('INVZ_SLOW')), assumeFail(testCase, 'set INVZ_SLOW=1'); end
ion = invz_ion();
a = invz_spectra_map(ion, 0.31, [0 3 8], wgrid(), base_opts());
b = invz_spectra_map(ion, 0.31, [0 3 8], wgrid(), ...
    setfield(base_opts(), 'static_medium', 'resummed'));  %#ok<SFLD>
verifyEqual(testCase, b.phase_1z, a.phase_1z);
verifyEqual(testCase, b.Sigma0, a.Sigma0, 'AbsTol', 0);
end

% An unknown PM probe must never vote 'ordered' under strict.
function test_unknown_pm_probe_cannot_label_ordered(testCase)
S = struct('crit_pm', NaN, 'pm_probe_status', {{'nonconverged'}});
verdict = invz_pm_verdict(NaN, false, 1e-6);            % helper under test
verifyEqual(testCase, verdict, 'unknown');
verifyEqual(testCase, invz_pm_verdict(1e-3, true, 1e-6), 'pm');
verifyEqual(testCase, invz_pm_verdict(-1e-3, true, 1e-6), 'ordered_eligible');
verifyEqual(testCase, invz_pm_verdict(1e-9, true, 1e-6), 'unknown');   % inside the band
end

% ---------------------------------------------------------------------------------------------
% invz_boundary_interval: the strict-side Bc_1z reduction, unit-tested WITHOUT running a
% spectrum (Step 6). The historical scalar invz_boundary_field reduction is a SEPARATE helper
% and stays the 'resummed' source of S.Bc_1z -- these cases pin the new helper only.

function test_boundary_interval_valid(testCase)
f = [1 2 3 4];
[iv, st, Bc] = invz_boundary_interval(f, logical([1 1 0 0]), logical([0 0 1 1]), false(1, 4));
verifyEqual(testCase, st, 'valid');
verifyEqual(testCase, iv, [2 3]);
verifyEqual(testCase, Bc, 2.5);
end

function test_boundary_interval_widened_by_intervening_unknown(testCase)
% An unknown column BETWEEN the anchors genuinely widens the bracket.
f = [1 2 3 4];
[iv, st, Bc] = invz_boundary_interval(f, logical([1 0 0 0]), logical([0 0 0 1]), logical([0 1 1 0]));
verifyEqual(testCase, st, 'widened');
verifyEqual(testCase, iv, [1 4]);
verifyEqual(testCase, Bc, 2.5);
end

function test_boundary_interval_unknown_outside_bracket_does_not_widen(testCase)
% The defect the earlier min/max reduction had: an unrelated unknown column on the far side
% of the sweep must NOT widen the boundary.
f = [1 2 3 4];
[iv, st, Bc] = invz_boundary_interval(f, logical([0 1 0 0]), logical([0 0 1 0]), logical([1 0 0 1]));
verifyEqual(testCase, st, 'valid');
verifyEqual(testCase, iv, [2 3]);
verifyEqual(testCase, Bc, 2.5);
end

function test_boundary_interval_unbracketed(testCase)
% Labels present but the flip is not bracketed (first PM anchor precedes the last ordered one).
f = [1 2 3 4];
[iv, st, Bc] = invz_boundary_interval(f, logical([0 0 1 1]), logical([1 1 0 0]), false(1, 4));
verifyEqual(testCase, st, 'unbracketed');
verifyTrue(testCase, all(isnan(iv)));
verifyTrue(testCase, isnan(Bc));
% and with no PM anchor at all
[~, st2, Bc2] = invz_boundary_interval(f, logical([1 1 0 0]), false(1, 4), false(1, 4));
verifyEqual(testCase, st2, 'unbracketed');
verifyTrue(testCase, isnan(Bc2));
end

function test_boundary_interval_all_unknown(testCase)
% Nothing is labelled: no anchors, so no boundary -- never a fabricated midpoint.
f = [1 2 3 4];
[iv, st, Bc] = invz_boundary_interval(f, false(1, 4), false(1, 4), true(1, 4));
verifyEqual(testCase, st, 'unbracketed');
verifyTrue(testCase, all(isnan(iv)));
verifyTrue(testCase, isnan(Bc));
end

function test_boundary_interval_inconsistent_labels_are_invalid(testCase)
% A column carrying two mutually exclusive labels is a CALLER predicate defect, not a
% physics outcome: report it rather than silently anchoring on one of them.
f = [1 2 3 4];
[iv, st, Bc] = invz_boundary_interval(f, logical([1 1 0 0]), logical([0 1 1 0]), false(1, 4));
verifyEqual(testCase, st, 'invalid');
verifyTrue(testCase, all(isnan(iv)));
verifyTrue(testCase, isnan(Bc));
end

function test_boundary_interval_requires_increasing_fields(testCase)
% Hard input contract: "last ordered" / "first PM" must not depend on presentation order.
verifyError(testCase, @() invz_boundary_interval([1 3 2], logical([1 0 0]), logical([0 0 1]), ...
    false(1, 3)), 'invz:boundaryInterval');
verifyError(testCase, @() invz_boundary_interval([1 2 NaN], logical([1 0 0]), logical([0 0 1]), ...
    false(1, 3)), 'invz:boundaryInterval');
verifyError(testCase, @() invz_boundary_interval([1 2 3], logical([1 0]), logical([0 0 1]), ...
    false(1, 3)), 'invz:boundaryInterval');
end
