% Task 12 bit-identity capture. Run BEFORE and AFTER the invz_hmf_ordered refactor with
% TASK12_TAG set to 'pre' / 'post'; the compare script then diffs the two .mat files.
% Cases:
%   A  resummed default profile (the diagnostics fixture's (ion,T,Bx,Jnu), resummed scheme)
%   B  force_bare shortcut profile (exercises eval_node's fbare exit)
%   C  traced run of case A  -> numel(trc.nodes) invariant (report item 2)
%   D  coarse hmin_frac so the ADAPTIVE EXTENSION + REDENSIFY sites actually run.
%      PROBE ONLY -- not a committed test, and the brief's fixture is untouched.
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root);
addpath(fullfile(root, 'invz_projected'));
addpath(fullfile(root, 'invz_common'));

tag = getenv('TASK12_TAG');
if isempty(tag), tag = 'pre'; end

ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
J0eff = 6.42e-3;
o = struct('J0eff', J0eff, 'Jxx0', ion.Jxx0, 'hyp', true);

C = struct();

% --- A: resummed default -------------------------------------------------------------------
[hA, pA] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
C.A = struct('h', hA, 'prof', pA);
fprintf('CAP A status=%-12s h=%.17g nNodes=%d n_extend=%d redens=%d\n', ...
        pA.status, hA, numel(pA.hgrid), pA.n_extend, pA.redensified);

% --- B: force_bare -------------------------------------------------------------------------
ob = o;  ob.force_bare = true;
[hB, pB] = invz_hmf_ordered(ion, T, Bx, Jnu, ob);
C.B = struct('h', hB, 'prof', pB);
fprintf('CAP B status=%-12s h=%.17g nNodes=%d\n', pB.status, hB, numel(pB.hgrid));

% --- C: traced run of A -> trc.nodes count -------------------------------------------------
ot = o;  ot.trace = true;
[hC, pC, tC] = invz_hmf_ordered(ion, T, Bx, Jnu, ot);
C.C = struct('h', hC, 'prof', pC, 'n_nodes', numel(tC.nodes), 'n_iters', numel(tC.iters), ...
             'node_h', [tC.nodes.h], 'node_phase', {{tC.nodes.phase}}, ...
             'node_term', {{tC.nodes.term_reason}}, 'node_ok', [tC.nodes.ok_final], ...
             'node_id', [tC.nodes.id], 'node_K0', [tC.nodes.K0], ...
             'node_Duni', [tC.nodes.D_uni], 'node_resid', [tC.nodes.resid_static], ...
             'node_seedkind', {{tC.nodes.seed_kind}}, 'node_seedfrom', [tC.nodes.seed_from], ...
             'node_outer', [tC.nodes.outer_iters], 'node_hitmax', [tC.nodes.outer_hit_max], ...
             'node_dSbreak', [tC.nodes.dS_break]);
fprintf('CAP C numel(trc.nodes)=%d  numel(trc.iters)=%d  phases={%s}\n', ...
        numel(tC.nodes), numel(tC.iters), strjoin(unique({tC.nodes.phase}, 'stable'), ','));

% --- D: coarse hmin_frac -> extension + redensification ------------------------------------
% Same deterministic construction test_invz_hmf_ordered/test_near_boundary_root_not_missed
% uses (fine-grid reference root first, then a coarse lower limit provably above it), at the
% fixed field B = 3.00 T where it actually fires (measured: n_extend = 7, redensified = 1).
Bd = [3.00 0 0];
[hDref, pDref] = invz_hmf_ordered(ion, T, Bd, Jnu, o);
fracD = min(0.5, 4*hDref/max(pDref.hgrid));
od = o;  od.hmin_frac = fracD;  od.nH = 17;  od.trace = true;
[hD, pD, tD] = invz_hmf_ordered(ion, T, Bd, Jnu, od);
C.Dref = struct('h', hDref, 'frac', fracD);
C.D = struct('h', hD, 'prof', pD, 'n_nodes', numel(tD.nodes), 'n_iters', numel(tD.iters), ...
             'node_h', [tD.nodes.h], 'node_phase', {{tD.nodes.phase}}, ...
             'node_term', {{tD.nodes.term_reason}}, 'node_ok', [tD.nodes.ok_final], ...
             'node_id', [tD.nodes.id], 'node_K0', [tD.nodes.K0], ...
             'node_Duni', [tD.nodes.D_uni], 'node_resid', [tD.nodes.resid_static], ...
             'node_seedkind', {{tD.nodes.seed_kind}}, 'node_seedfrom', [tD.nodes.seed_from], ...
             'node_outer', [tD.nodes.outer_iters], 'node_hitmax', [tD.nodes.outer_hit_max], ...
             'node_dSbreak', [tD.nodes.dS_break]);
fprintf('CAP D status=%-12s h=%.17g nNodes=%d n_extend=%d redens=%d trcNodes=%d phases={%s}\n', ...
        pD.status, hD, numel(pD.hgrid), pD.n_extend, pD.redensified, numel(tD.nodes), ...
        strjoin(unique({tD.nodes.phase}, 'stable'), ','));

save(fullfile(root, '.superpowers', 'sdd', ['task12_cap_' tag '.mat']), 'C');
fprintf('CAP DONE tag=%s\n', tag);
