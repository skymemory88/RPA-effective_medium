function tests = test_invz_plot_spectra_qpath
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_line_xdata_matches_qaxis_not_color_scratch(testCase)
% Regression: a percentile/color-limit block reused the variable name `x` as scratch
% space for sorted chi values, clobbering the q-axis coordinate also named `x` before
% the peak-overlay plot() call consumed it. Pins that the plotted line's XData is
% exactly S.x (length nq), not the clobbered scratch array.
ion = invz_ion();
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
w = (0.02:0.02:0.6).';
qpath = [1 0 0; 1.5 0 0; 2 0 0];
S = invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, struct('Jnu', Jnu, 'info', info, 'dpRng', 10));

fig = figure('Visible', 'off');
cleanupObj = onCleanup(@() close(fig));
ax = axes(fig);
invz_plot_spectra_qpath(ax, S, S.chiz, S.Epeak, 'test', 'meV');

lines = findobj(ax, 'Type', 'line');
verifyNumElements(testCase, lines, 1);
verifyEqual(testCase, lines.XData, S.x);
verifyEqual(testCase, lines.YData, S.Epeak);
end
