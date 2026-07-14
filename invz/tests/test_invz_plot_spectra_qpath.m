function tests = test_invz_plot_spectra_qpath
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_line_xdata_matches_qaxis_not_color_scratch(testCase)
% Regression: the percentile/color-limit block used to reuse the variable name `x` as
% scratch space for sorting finite positive chi values (an array whose length is the
% number of finite chi'' samples, NOT the number of q-points). Once the q-axis coordinate
% was also bound to a persistent `x`, that scratch reuse silently clobbered it before the
% peak-overlay plot() call consumed it, so plot() saw mismatched vector lengths and
% errored (or, for accidentally-matching lengths, would have plotted the wrong x-data
% silently). Pin that the plotted line's XData is exactly S.x, not some other length.
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
