% F2 evidence: does test_invz_ordered_trace.m's ismember assertion (line ~72-73) FAIL against
% the OLD 4-element allowed set under a genuine strict domain-event run, and PASS against the
% NEW 5-element set? Read-only measurement; no committed file touched by this script.
ROOT = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(ROOT); addpath(fullfile(ROOT,'invz_projected')); addpath(fullfile(ROOT,'invz_common'));

ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
solve_opts = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
                     'static_medium', 'strict_1z_dyson_ref', 'ref_margin', 1e9);

% Route through invz_ordered_trace itself (the exact wrapper test_invz_ordered_trace.m calls),
% synthetic (no qc/Jnu_unflat) since only trc.nodes.term_reason is under test here.
trc = invz_ordered_trace(ion, T, Bx, Jnu, struct('solve', solve_opts));

OLD_SET = {'converged', 'max_iter', 'refresh_failed', 'bare_shortcut'};
NEW_SET = {'converged', 'max_iter', 'refresh_failed', 'bare_shortcut', 'medium_out_of_domain'};

reasons = {trc.nodes.term_reason};
fprintf('numel(trc.nodes) = %d\n', numel(trc.nodes));
fprintf('term_reason set observed = {%s}\n', strjoin(unique(reasons, 'stable'), ', '));

old_ok = all(ismember(reasons, OLD_SET));
new_ok = all(ismember(reasons, NEW_SET));
fprintf('all(ismember(reasons, OLD_SET)) = %d   <-- would verifyTrue PASS with the OLD set?\n', old_ok);
fprintf('all(ismember(reasons, NEW_SET)) = %d   <-- would verifyTrue PASS with the NEW set?\n', new_ok);
n_not_in_old = nnz(~ismember(reasons, OLD_SET));
fprintf('count of nodes whose term_reason is NOT in the OLD set = %d / %d\n', n_not_in_old, numel(reasons));
if ~old_ok && new_ok
    fprintf('\nCONCLUSION: the OLD assertion (ismember against OLD_SET) WOULD HAVE FAILED here;\n');
    fprintf('the NEW assertion (ismember against NEW_SET) PASSES. This is a genuine latent\n');
    fprintf('failure under a strict domain-event run, not a hypothetical.\n');
else
    fprintf('\nCONCLUSION: unexpected -- old_ok=%d new_ok=%d (investigate further before reporting)\n', old_ok, new_ok);
end
