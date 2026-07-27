% Post-fix verification: F9 bitwise identity, F11 K0 provenance, F12 schema_version, on the
% CURRENT (edited) codebase. Read-only measurement.
ROOT = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(ROOT); addpath(fullfile(ROOT,'invz_projected')); addpath(fullfile(ROOT,'invz_common'));

ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
           'static_medium', 'strict_1z_dyson_ref');

fprintf('=== F9: bitwise identity of prof.slope0 (now rec0.crit) vs r_pm0+J0eff*G0bare_pm0 ===\n');
[~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
lhs = prof.slope0;  rhs = prof.r_pm0 + o.J0eff*prof.G0bare_pm0;
fprintf('slope0 = %.17g\nr_pm0+J0eff*G0bare_pm0 = %.17g\nbitwise isequal = %d\n', lhs, rhs, isequal(lhs, rhs));

% Also on the resummed default (a second fixture, matching the reviewer's own two-fixture check)
orr = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
[~, profr] = invz_hmf_ordered(ion, T, Bx, Jnu, orr);
lhs2 = profr.slope0;  rhs2 = profr.r_pm0 + orr.J0eff*profr.G0bare_pm0;
fprintf('[resummed] slope0 = %.17g ; bitwise isequal = %d\n', lhs2, isequal(lhs2, rhs2));

fprintf('\n=== F11: K0 provenance on a TRACED force_bare run (post-fix) ===\n');
ob = o;  ob.force_bare = true;  ob.trace = true;
[~, pb, tb] = invz_hmf_ordered(ion, T, Bx, Jnu, ob);
fprintf('trc.nodes numel=%d ; K0 set = %s ; term_reason set = {%s}\n', numel(tb.nodes), ...
    mat2str(unique([tb.nodes.K0])), strjoin(unique({tb.nodes.term_reason},'stable'), ','));
fprintf('prof.K0 set on this path = %s ; prof.Sigma0 set = %s ; prof.r set = %s\n', ...
    mat2str(unique(pb.K0)), mat2str(unique(pb.Sigma0)), mat2str(unique(pb.r)));

fprintf('\n=== F12: trc.schema_version (post-fix) ===\n');
fprintf('trc.schema_version = %d\n', tb.schema_version);
[~, ~, tplain] = invz_hmf_ordered(ion, T, Bx, Jnu, struct('J0eff',6.42e-3,'Jxx0',ion.Jxx0,'hyp',true,'trace',true));
fprintf('trc.schema_version (untraced-default fixture, traced) = %d ; enabled=%d\n', ...
    tplain.schema_version, tplain.enabled);

fprintf('\nPOSTFIX VERIFY DONE\n');
