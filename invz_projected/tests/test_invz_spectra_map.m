function tests = test_invz_spectra_map
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_map_shape_and_phase_codes(testCase)
% Structural contract of the (omega, B) map: correct shape, valid per-field phase codes,
% a high-field column that is paramagnetic (phase 2) with a finite dissipative spectrum.
% A cheap synthetic coupling is injected so the solve is fast (no 16^3 lattice sum).
ion = invz_ion();
T = 0.31;
fields = [2.0 5.5];                       % 2 T: ordered candidate;  5.5 T: paramagnet
w = (0.02:0.02:0.6).';
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';     % stand-in branch eigenvalues (meV), all < Jcc0
S = invz_spectra_map(ion, T, fields, w, ...
                     struct('Jnu', Jnu, 'info', info, 'verbose', false));

verifySize(testCase, S.chiz,   [numel(w) numel(fields)]);
verifySize(testCase, S.chirpa, [numel(w) numel(fields)]);
verifyTrue(testCase, all(ismember(S.phase, [0 1 2])));      % 0 masked, 1 FM, 2 PM

% Bx = 5.5 T: paramagnet -> phase 2, finite dissipative RPA overlay
verifyEqual(testCase, S.phase(2), 2);
finite2 = isfinite(S.chirpa(:, 2));
verifyTrue(testCase, any(finite2));
verifyGreaterThanOrEqual(testCase, min(S.chirpa(finite2, 2)), -1e-9);
% and a genuine resonance peak in the 1/z spectrum somewhere on the grid
verifyTrue(testCase, any(isfinite(S.chiz(:))));
verifyGreaterThan(testCase, max(S.chiz(isfinite(S.chiz))), 0);

% per-field peak energy (invz_peak_energy, shared with invz_spectra_qpath): right shape,
% and any censored-survivor is a genuine in-window frequency, not a fabricated edge value.
verifySize(testCase, S.Epeak,     [1 numel(fields)]);
verifySize(testCase, S.Epeak_rpa, [1 numel(fields)]);
finiteEpeak = S.Epeak(isfinite(S.Epeak));
if ~isempty(finiteEpeak)
    verifyGreaterThanOrEqual(testCase, min(finiteEpeak), min(w));
    verifyLessThanOrEqual(testCase, max(finiteEpeak), max(w));
end
end

function test_split_phase_diagnostics(testCase)
% Stage-1 QCP split (docs/superpowers/plans/2026-07-21-invzp-qcp-stage1-split-overlays.md):
% the auto overlay keeps the auto-solve (bare-MF) selection, the 1/z curve carries its own
% stability-based phase label, and the per-field PM 1/z mass plus midpoint boundary
% estimates are returned. Same cheap synthetic coupling as the structural test above.
ion = invz_ion();
T = 0.31;
fields = [2.0 5.5];                       % 2 T: ordered candidate;  5.5 T: paramagnet
w = (0.02:0.02:0.6).';
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
S = invz_spectra_map(ion, T, fields, w, ...
                     struct('Jnu', Jnu, 'info', info, 'verbose', false));

verifySize(testCase, S.phase_1z, [1 numel(fields)]);
verifySize(testCase, S.crit_pm,  [1 numel(fields)]);
verifyTrue(testCase, all(ismember(S.phase_1z, [0 1 2])));
% the phase_1z = 2 label is GATED on a stable PM mass -- unconditional by construction
lab2 = S.phase_1z == 2;
verifyTrue(testCase, all(isfinite(S.crit_pm(lab2))));
if any(lab2), verifyGreaterThan(testCase, S.crit_pm(lab2), 0); end
% a PM auto phase with stable mass carries the PM 1/z label; spurious PM (crit <= 0)
% must land in S.suspect with phase_1z = 0, never in phase_1z = 2
verifyTrue(testCase, isfield(S, 'suspect') && islogical(S.suspect));
verifyEqual(testCase, S.suspect, (S.phase == 2) & isfinite(S.crit_pm) & (S.crit_pm <= 0));
verifyTrue(testCase, all(S.phase_1z(S.phase == 2 & ~S.suspect & isfinite(S.crit_pm)) == 2));
% validity masks reach the returned spectra (re-review findings 2, R3-1): no valid 1/z
% label -> chiz column all-NaN; NO accepted auto state (phase 0) OR suspect -> chirpa
% column all-NaN (the legacy phase-0 Sigma-zero fallback is no longer surfaced)
verifyTrue(testCase, all(all(isnan(S.chiz(:,  S.phase_1z == 0)))));
verifyTrue(testCase, all(all(isnan(S.chirpa(:, S.phase == 0 | S.suspect)))));
% boundary estimates are always present as fields (NaN allowed: two fields cannot bracket
% a boundary that is not crossed between them)
verifyTrue(testCase, isfield(S, 'Bc_auto') && isfield(S, 'Bc_1z'));
% Bc_auto predicate (re-review R3-2): a suspect PM point between the last ordered and the
% first VALID PM point is excluded and WIDENS the bracket -- direct unit check on the
% public helper with crafted labels (index 3 plays the suspect point)
fieldsU = [1 2 3 4];
verifyEqual(testCase, invz_boundary_field(fieldsU, logical([1 1 0 0]), logical([0 0 0 1])), 3.0);
verifyEqual(testCase, invz_boundary_field(fieldsU, logical([1 1 0 0]), logical([0 0 1 1])), 2.5);
verifyTrue(testCase, isnan(invz_boundary_field(fieldsU, logical([1 1 0 0]), false(1, 4))));
% the 1/z map itself must still be a well-formed spectrum
verifyTrue(testCase, any(isfinite(S.chiz(:))));
end
