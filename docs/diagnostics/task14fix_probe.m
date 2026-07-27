function task14fix_probe(PHASE)
% Task-14 fix probes (review finding 1): the ordered BARE loop's strict-scheme domain halt.
% Run identically before and after the fix; PHASE is 'pre' or 'post'.
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT, 'invz_common'));
addpath(genpath(fullfile(ROOT, 'invz_projected')));

if nargin < 1, PHASE = 'pre'; end

ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o   = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
ob  = o;  ob.static_medium = 'strict_1z_dyson_ref';     % strict + DEFAULT ordered_mode 'bare'
o2  = ob; o2.ref_margin = 2;                            % same knob task-14 probe P2 used

sweep(PHASE, 'ref_margin=2 B=2.85', ion, Jnu, o2, 0.31, 2.85, 4);
sweep(PHASE, 'default margin B=2.0', ion, Jnu, ob, 0.31, 2.0, 40);

% ---- A: strict + bare, ref_margin above the denominator (fires on iteration 1) -------------
A = invz_solve_point_ordered(ion, 0.31, [2.85 0 0], Jnu, o2);
fprintf('A %s bare strict ref_margin=2 B=2.85: is_ordered=%d conv=%d iters=%d status=%s\n', ...
    PHASE, A.is_ordered, A.converged, A.outer_iters, A.medium_status);
fprintf('A %s   denom=%.17g margin=%.17g alpha=%.6g alpha_m=%.6g K1=%.6g G1=%.6g\n', ...
    PHASE, A.medium_denom, A.medium_margin, A.alpha, A.alpha_m, A.K(1), A.G(1));
fprintf('A %s   lambda=[%s] Sigma0=%.6g crit=%.6g sumrule_rel=%.6g\n', ...
    PHASE, num2str(A.lambda(:).'), A.Sigma0, A.crit, A.sumrule_rel);

% ---- A2: strict + bare at B = 2.0 T, the NATURAL domain event on this fixture --------------
A2 = invz_solve_point_ordered(ion, 0.31, [2.0 0 0], Jnu, ob);
fprintf('A2 %s bare strict default-margin B=2.0: conv=%d iters=%d status=%s\n', ...
    PHASE, A2.converged, A2.outer_iters, A2.medium_status);

% ---- B: a strict + bare point whose medium stays OK, for the FIELD-SET comparison ----------
B = invz_solve_point_ordered(ion, 0.31, [2.85 0 0], Jnu, ob);
fprintf('B %s bare strict B=2.85: is_ordered=%d conv=%d iters=%d status=%s denom=%.6g\n', ...
    PHASE, B.is_ordered, B.converged, B.outer_iters, B.medium_status, B.medium_denom);
fa = sort(fieldnames(A));  fb = sort(fieldnames(B));
onlyA = setdiff(fa, fb);  onlyB = setdiff(fb, fa);
fprintf('AB %s FIELDSET identical=%d  only_in_halt={%s} only_in_ok={%s}\n', ...
    PHASE, isequal(fa, fb), strjoin(onlyA.', ','), strjoin(onlyB.', ','));

% ---- C/D: the DEFAULT resummed bare path -- must be untouched, bitwise --------------------
C = invz_solve_point_ordered(ion, 0.31, [2.0 0 0], Jnu, o);
fprintf('C %s bare resummed B=2.0: conv=%d iters=%d status=%s crit=%.17g Sigma0=%.17g\n', ...
    PHASE, C.converged, C.outer_iters, C.medium_status, C.crit, C.Sigma0);
D = invz_solve_point_ordered(ion, 0.31, [2.85 0 0], Jnu, o);
fprintf('D %s bare resummed B=2.85: conv=%d iters=%d status=%s crit=%.17g Sigma0=%.17g\n', ...
    PHASE, D.converged, D.outer_iters, D.medium_status, D.crit, D.Sigma0);

% ---- E: inertness under 'resummed', MEASURED not asserted ---------------------------------
% Every medium a resummed bare solve produces must report 'not_applicable', one of the two
% strings the guard whitelists -- so the guard can never fire on the default path.
stats = {};
for Bv = [0.1 0.5 1.0 1.5 2.0 2.5 2.85 3.0]
    p = invz_solve_point_ordered(ion, 0.31, [Bv 0 0], Jnu, o);
    stats{end+1} = p.medium_status;  %#ok<AGROW>
end
fprintf('E %s resummed bare status set over B sweep = {%s}\n', PHASE, strjoin(unique(stats), ','));

save(fullfile(ROOT, '.superpowers', 'sdd', ['task14fix_' PHASE '.mat']), 'A', 'A2', 'B', 'C', 'D');
fprintf('%s PROBE DONE\n', upper(PHASE));
end

function sweep(PHASE, tag, ion, Jnu, oo, T, Bv, kmax)
% TRUE per-iteration status sequence, measured non-invasively: max_outer = k runs exactly the
% first k iterations of the SAME trajectory, and pt.medium_status is the last invz_emt_scalar
% verdict, so sweeping k reads the true sequence straight off the loop.
seq = cell(1, kmax);
for k = 1:kmax
    ok = oo;  ok.max_outer = k;
    p = invz_solve_point_ordered(ion, T, [Bv 0 0], Jnu, ok);
    seq{k} = p.medium_status;
end
fb = find(~ismember(seq, {'ok', 'not_applicable'}), 1);
if isempty(fb)
    fprintf('A0 %s [%s] no non-ok status in k=1..%d\n', PHASE, tag, kmax);
else
    fprintf('A0 %s [%s] FIRST non-ok iteration = %d, TRUE cause = %s (seq: %s)\n', ...
        PHASE, tag, fb, seq{fb}, strjoin(seq(1:min(fb+1, kmax)), ' '));
end
end
