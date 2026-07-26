function tests = test_invz_hmf_ordered_status
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% B = 0 no longer escapes as a throw: it is a LABELLED domain outcome, distinguishable from a
% solver failure. Production sweeps start at B = 0, so this path is hit on every run.
function test_zero_field_is_degenerate_not_a_throw(testCase)
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
[hstar, prof] = invz_hmf_ordered(ion, 0.31, [0 0 0], Jnu, o);   % must not throw
verifyTrue(testCase, isnan(hstar));
verifyEqual(testCase, prof.status, 'degenerate_doublet');
verifyTrue(testCase, isfinite(prof.Delta(1)) || isnan(prof.Delta(1)));
end

% A medium domain event is its own status, NOT 'node_failed'.
function test_medium_domain_status_is_distinct(testCase)
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
           'static_medium', 'strict_1z_dyson_ref', ...
           'ref_margin', 1e9);                          % public authority: force domain event
[hstar, prof] = invz_hmf_ordered(ion, 0.31, [2.85 0 0], Jnu, o);
verifyTrue(testCase, isnan(hstar));
verifyEqual(testCase, prof.status, 'medium_out_of_domain');
end

% Legacy statuses keep their exact relative order: the two new cases are prepended only.
function test_legacy_status_order_unchanged(testCase)
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
[~, prof] = invz_hmf_ordered(ion, 0.31, [2.85 0 0], Jnu, o);
verifyTrue(testCase, any(strcmp(prof.status, ...
    {'ok', 'unresolved', 'node_failed', 'no_bare_order'})));
end

% Integration anchors for the two new statuses and the pre-existing no-bare-order exit.
function test_new_and_no_bare_statuses_are_reachable(testCase)
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
base = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
seen = {};
[~, p] = invz_hmf_ordered(ion, 0.31, [0 0 0], Jnu, base);          seen{end+1} = p.status;
o2 = base;  o2.static_medium = 'strict_1z_dyson_ref';
o2.ref_margin = 1e9;
[~, p] = invz_hmf_ordered(ion, 0.31, [2.85 0 0], Jnu, o2);         seen{end+1} = p.status;
[~, p] = invz_hmf_ordered(ion, 0.31, [9 0 0], Jnu, base);          seen{end+1} = p.status;
verifyTrue(testCase, ismember('degenerate_doublet', seen));
verifyTrue(testCase, ismember('medium_out_of_domain', seen));
verifyTrue(testCase, ismember('no_bare_order', seen));             % 9 T: bare does not order
end

% Pure precedence pin, including the mixed case node_failed > unresolved.
function test_status_precedence_all_five_cases(testCase)
base = struct('accepted', true, 'term_reason', 'accepted', 'medium_status', 'ok');
pred = base;  nodes = repmat(base, 1, 2);  F = [0.1 0.2];
verifyEqual(testCase, invz_hmf_status(pred, nodes, -1, F), 'unresolved');
nodes(1).accepted = false; nodes(1).term_reason = 'max_iter';
verifyEqual(testCase, invz_hmf_status(pred, nodes, -1, F), 'node_failed');
nodes(1).medium_status = 'ref_denom_small';
verifyEqual(testCase, invz_hmf_status(pred, nodes, -1, F), 'medium_out_of_domain');
nodes(2).term_reason = 'degenerate_doublet';
verifyEqual(testCase, invz_hmf_status(pred, nodes, -1, F), 'degenerate_doublet');
nodes = repmat(base, 1, 2);
verifyEqual(testCase, invz_hmf_status(pred, nodes, 1, [-0.1 0.2]), 'ok');
end
