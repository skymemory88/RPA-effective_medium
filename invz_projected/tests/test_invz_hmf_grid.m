function tests = test_invz_hmf_grid
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% Byte-identical to the expression it replaces (invz_hmf_ordered.m:144-145).
function test_matches_the_inline_expression(testCase)
hmax = 0.37;  nH = 33;  hfrac = 1e-3;
ratio_ref = hfrac^(1/(nH-1));
hgrid_ref = hmax * ratio_ref.^((nH-1):-1:0);
[hgrid, ratio] = invz_hmf_grid(hmax, nH, hfrac);
verifyEqual(testCase, hgrid, hgrid_ref, 'AbsTol', 0);
verifyEqual(testCase, ratio, ratio_ref, 'AbsTol', 0);
end

function test_shape_and_endpoints(testCase)
[hgrid, ~] = invz_hmf_grid(0.5, 17, 1e-4);
verifyEqual(testCase, numel(hgrid), 17);
verifyEqual(testCase, hgrid(end), 0.5, 'RelTol', 1e-15);
verifyEqual(testCase, hgrid(1), 0.5*1e-4, 'RelTol', 1e-12);
verifyTrue(testCase, all(diff(hgrid) > 0));            % ascending, clustered at 0
end

function test_rejects_bad_arguments(testCase)
verifyError(testCase, @() invz_hmf_grid(0, 33, 1e-3), 'invz:hmfGrid');
verifyError(testCase, @() invz_hmf_grid(0.5, 1, 1e-3), 'invz:hmfGrid');
verifyError(testCase, @() invz_hmf_grid(0.5, 33, 0), 'invz:hmfGrid');
verifyError(testCase, @() invz_hmf_grid(0.5, 33, 1.5), 'invz:hmfGrid');
end
