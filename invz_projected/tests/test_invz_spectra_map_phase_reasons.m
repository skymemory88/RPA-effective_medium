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
function o = syn_opts()
% Cheap synthetic couplings: no 16^3 lattice sum, so these cases need no INVZ_SLOW gate.
o = struct('Jnu', linspace(-2e-3, 6.0e-3, 24).', 'info', struct('Jcc0', 6.4e-3), ...
           'verbose', false);
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
             'not_attempted_longitudinal', 'bare_not_ordered', 'response_failed'})), ...
             S.phase_1z_reason{k});
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

% Review F2, STRICT counterpart of the honesty gate above (the resummed sibling accepts only
% {degenerate_doublet, bare_not_ordered, pm} and would FAIL under a strict scheme). At B = 0 the
% strict PM probe dies with invz:degenerateDoublet, so the PM-probe-dominant phase reason is
% 'pm_probe_unknown' BY DESIGN -- a failed PM probe must never vote 'ordered'. The KNOWABLE cause
% must therefore be recorded on the diagnostic field, and the sweep summary must attribute it:
% before this fix the summary said "0 degenerate-doublet" for a column that is one, which made
% the strict path LESS informative than resummed at exactly the knowable column.
%
% UNGATED (review M5). This is the ONLY regression test for F2, so it must run in the DEFAULT
% batch. The synthetic couplings of syn_opts() reproduce all six assertions in seconds at
% fields = [0 5.50] -- measured: reason 'pm_probe_unknown | pm', ordered_diag_reason
% 'degenerate_doublet | not_attempted', pm_probe_error_id 'invz:degenerateDoublet', and the same
% invz:spectraMapMasked text -- so the 4^3 lattice sum the INVZ_SLOW gate existed for is not
% needed. Reverting the F2 counter attribution (the diag_cause override) makes the two counter
% assertions fail, so the gate still discriminates.
function test_zero_field_column_is_labelled_degenerate_strict(testCase)
ion = invz_ion();
o = syn_opts();  o.ordered_1z = 'jensen';  o.static_medium = 'strict_1z_dyson_ref';
lastwarn('', 'invz:none');
S = invz_spectra_map(ion, 0.31, [0 5.50], wgrid(), o);
verifyEqual(testCase, S.phase_1z(1), 0);
verifyEqual(testCase, S.ordered_diag_reason{1}, 'degenerate_doublet');
verifyNotEqual(testCase, S.ordered_diag_reason{1}, 'solver_failed');
verifyEqual(testCase, S.pm_probe_error_id{1}, 'invz:degenerateDoublet');
% the counters must be factually true -- Stage-4 G16 reads them
[wmsg, wid] = lastwarn();
verifyEqual(testCase, wid, 'invz:spectraMapMasked');
verifyTrue(testCase, contains(wmsg, '1 degenerate-doublet'), wmsg);
verifyTrue(testCase, contains(wmsg, '0 unknown-PM-probe'), wmsg);
end

% Review F1 regression: an EMPTY field sweep must RETURN, not throw. all([]) is true and 0 < 2 is
% true, so an incomplete pre-screen let fields = [] reach invz_boundary_interval, whose (correct)
% nonempty precondition then raised invz:boundaryInterval -- on the DEFAULT 'resummed' path,
% which must stay behaviour-identical. The helper's contract is NOT loosened; the screen is.
function test_empty_field_sweep_returns(testCase)
ion = invz_ion();
w = wgrid();
S = invz_spectra_map(ion, 0.31, [], w, syn_opts());        % default scheme = 'resummed'
verifySize(testCase, S.chiz,   [numel(w) 0]);
verifySize(testCase, S.chirpa, [numel(w) 0]);
verifyEmpty(testCase, S.phase_1z);
verifyEmpty(testCase, S.phase_1z_reason);
verifyTrue(testCase, isnan(S.Bc_1z));
verifyEqual(testCase, S.Bc_1z_status, 'invalid');   % reduction never ran: no sweep to reduce
% and under a strict scheme, where the interval reduction actually feeds Bc_1z
o = syn_opts();  o.static_medium = 'strict_1z_dyson_ref';
Ss = invz_spectra_map(ion, 0.31, [], w, o);
verifySize(testCase, Ss.chiz, [numel(w) 0]);
verifyTrue(testCase, isnan(Ss.Bc_1z));
end

% Review M1, COMPLETING F1. The driver pre-screen must be a SUPERSET of
% invz_boundary_interval's precondition -- the call is pure and lands in new fields, so throwing
% is the only way it can change default-path behaviour. It screened ~isempty / isfinite /
% strictly-increasing but NOT isnumeric / isreal, so a char `fields` (which sweeps perfectly
% well, char*double being double) reached the helper and raised invz:boundaryInterval on the
% DEFAULT 'resummed' path, where b9082fd returned. Same defect class as F1.
function test_nonnumeric_field_sweep_returns(testCase)
ion = invz_ion();
w = wgrid();
S = invz_spectra_map(ion, 0.31, char([3 5]), w, syn_opts());    % default scheme = 'resummed'
verifySize(testCase, S.chiz, [numel(w) 2]);
verifyEqual(testCase, S.Bc_1z_status, 'invalid');    % screened out: the reduction never ran
verifyTrue(testCase, all(isnan(S.Bc_1z_interval)));
% and under strict, where that reduction is what feeds Bc_1z
o = syn_opts();  o.static_medium = 'strict_1z_dyson_ref';
Ss = invz_spectra_map(ion, 0.31, char([3 5]), w, o);
verifyEqual(testCase, Ss.Bc_1z_status, 'invalid');
verifyTrue(testCase, isnan(Ss.Bc_1z));
end

% Review I2. n_dom/n_deg/n_unk between them attribute only 4 of the 9 masked reason strings, so
% a column masked for one of the other five contributed ZERO to the reported total. With a
% FINITE Bc_1z the boundary_lost trigger is disarmed too (it requires ~isfinite(Bc_1z)), so such
% a sweep emitted NO WARNING AT ALL. MEASURED at 73dc55c on exactly this fixture: Bc_1z = 1.75
% ('valid'), phase_1z_reason 'solver_failed | ordered | pm | pm', warning identifier = none.
% The B = 0 column of these weakened synthetic couplings is masked 'solver_failed' because BOTH
% auto legs fail there -- too weak to order spontaneously, and a degenerate doublet in the PM
% leg -- while 1.5/2.0 T still bracket a boundary. The headline count must be
% nnz(phase_1z == 0), not the attributed subtotal.
function test_masked_column_with_uncounted_cause_is_not_silent(testCase)
ion = invz_ion();
w = wgrid();
o = syn_opts();
o.Jnu = 0.40 * o.Jnu;  o.info.Jcc0 = 0.40 * o.info.Jcc0;   % too weak to order at B = 0
o.static_medium = 'strict_1z_dyson_ref';
F = [0 1.5 2.0 5.50];
lastwarn('', 'invz:none');
S = invz_spectra_map(ion, 0.31, F, w, o);
n_msk = nnz(S.phase_1z == 0);
verifyGreaterThan(testCase, n_msk, 0);
verifyTrue(testCase, isfinite(S.Bc_1z));      % boundary_lost cannot be what fires the warning
uncounted = {'unstable_endpoint', 'solver_failed', 'not_attempted_longitudinal', ...
             'bare_not_ordered', 'response_failed'};
verifyTrue(testCase, any(ismember(S.phase_1z_reason(S.phase_1z == 0), uncounted)), ...
           strjoin(S.phase_1z_reason, ' | '));
[wmsg, wid] = lastwarn();
verifyEqual(testCase, wid, 'invz:spectraMapMasked');       % never silent
verifyTrue(testCase, contains(wmsg, ...
    sprintf('left %d of %d columns unlabelled', n_msk, numel(F))), wmsg);
verifyTrue(testCase, contains(wmsg, 'other'), wmsg);       % the remainder closes the sum
end

% Review F3 REACHABILITY PIN. No whitelisted recoverable identifier is reachable through
% invz_chi_realaxis today (task-15 report SS6), so every response-call boundary must absorb
% nothing and no column may be masked 'response_failed'. The strict branch now records the
% identifier AND masks the column if one ever fires; this pins the premise, so a change that
% makes it reachable fails here instead of silently producing an all-NaN spectrum on a column
% still labelled 'pm'/'ordered' -- which would anchor Bc_1z with no recorded cause.
function test_response_boundaries_absorb_nothing(testCase)
ion = invz_ion();
w = wgrid();
schemes = {'resummed', 'strict_1z_dyson_ref'};
for k = 1:numel(schemes)
    o = syn_opts();  o.static_medium = schemes{k};
    S = invz_spectra_map(ion, 0.31, [2.85 5.50], w, o);
    verifyNumElements(testCase, S.response_error_id, 2, schemes{k});
    verifyNumElements(testCase, S.ordered_diag_reason, 2, schemes{k});
    verifyTrue(testCase, all(cellfun(@isempty, S.response_error_id)), schemes{k});
    verifyFalse(testCase, any(strcmp(S.phase_1z_reason, 'response_failed')), schemes{k});
    % consequence: every LABELLED column really carries a spectrum, never an all-NaN one
    lab = S.phase_1z > 0;
    verifyTrue(testCase, all(any(isfinite(S.chiz(:, lab)), 1)), schemes{k});
end
end

% Review F7: under strict + ordered_1z = 'bare' every ordered column reports
% 'bare_escape_hatch', which by design cannot anchor a 1/z boundary, so `ord` is empty and Bc_1z
% comes back NaN -- while 'resummed' returns a FINITE one from the same run. With no masked
% columns the count-based warning trigger never fired, so that asymmetry was silent.
function test_lost_boundary_under_strict_bare_is_not_silent(testCase)
ion = invz_ion();
w = wgrid();
o = syn_opts();  o.ordered_1z = 'bare';  o.static_medium = 'strict_1z_dyson_ref';
S = verifyWarning(testCase, ...
    @() invz_spectra_map(ion, 0.31, [2.85 3.30 5.50], w, o), 'invz:spectraMapMasked');
verifyTrue(testCase, isnan(S.Bc_1z));
verifyTrue(testCase, any(strcmp(S.phase_1z_reason, 'bare_escape_hatch')));
verifyTrue(testCase, all(strcmp(S.static_medium_used(S.phase_1z == 1), 'n/a_bare_escape')));
% the same sweep under 'resummed' DOES reduce a finite boundary: that is the asymmetry the
% warning has to surface rather than leave to a reader comparing two runs by eye.
r = syn_opts();  r.ordered_1z = 'bare';
Sr = invz_spectra_map(ion, 0.31, [2.85 3.30 5.50], w, r);
verifyTrue(testCase, isfinite(Sr.Bc_1z));
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
