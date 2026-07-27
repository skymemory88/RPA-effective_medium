function tests = test_invz_static_domain_scan
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function [ion, Jnu, o] = fx()
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'static_medium', 'strict_1z_dyson_ref');
end

function test_explicit_hgrid_is_honoured(testCase)
[ion, Jnu, o] = fx();
o.hgrid = [1e-4 1e-3 1e-2];
scan = invz_static_domain_scan(ion, 0.31, [2.85 0 0], Jnu, o);
verifyEqual(testCase, scan.hgrid, o.hgrid, 'AbsTol', 0);
verifyEqual(testCase, scan.grid_source, 'explicit');
verifyEqual(testCase, scan.n_nodes, 3);
verifyEqual(testCase, scan.n_required, 4);     % explicit profile nodes plus h=0 predictor
end

function test_default_grid_comes_from_the_shared_helper(testCase)
[ion, Jnu, o] = fx();
scan = invz_static_domain_scan(ion, 0.31, [2.85 0 0], Jnu, o);
verifyEqual(testCase, scan.grid_source, 'invz_hmf_grid');
verifyEqual(testCase, numel(scan.hgrid), 33);
verifyTrue(testCase, all(diff(scan.hgrid) > 0));
end

% Every array is aligned and every counter accounted for -- no silent caps (spec §5.5).
function test_arrays_aligned_and_counters_complete(testCase)
[ion, Jnu, o] = fx();
scan = invz_static_domain_scan(ion, 0.31, [2.85 0 0], Jnu, o);
n = scan.n_nodes;
for f = {'Delta', 'valid', 'G0bare', 'Gref', 'ref_denom', 'ref_margin', ...
         'omit_mu3', 'omit_cubic', 'omit_max'}
    verifyEqual(testCase, numel(scan.(f{1})), n, f{1});
end
verifyEqual(testCase, numel(scan.ref_status), n);
verifyEqual(testCase, scan.n_valid + scan.n_degenerate + ...
            scan.n_out_of_domain + scan.n_skipped, scan.n_required, ...
            'every predictor/profile node must be accounted for exactly once');
end

% Small transverse field: the lowest nodes are where Delta dips below the floor.
function test_degenerate_nodes_are_flagged_at_small_field(testCase)
[ion, Jnu, o] = fx();
scan = invz_static_domain_scan(ion, 0.31, [0 0 0], Jnu, o);
verifyGreaterThan(testCase, scan.n_degenerate, 0);
verifyEqual(testCase, scan.predictor.status, 'degenerate_doublet');
end

function test_resummed_scheme_is_rejected(testCase)
[ion, Jnu, o] = fx();
o.static_medium = 'resummed';
verifyError(testCase, @() invz_static_domain_scan(ion, 0.31, [2.85 0 0], Jnu, o), ...
            'invz:staticMedium');
end
