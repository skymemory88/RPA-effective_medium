function tests = test_invz_spectra_map
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
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
end
