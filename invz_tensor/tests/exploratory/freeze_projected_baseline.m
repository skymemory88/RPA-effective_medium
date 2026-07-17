% freeze_projected_baseline.m  (A0 preflight, Task 2 of
% docs/superpowers/plans/2026-07-17-invz-tensor-full.md)
%
% DEV-TIME generator script (projected ON path). NEVER run by any test suite:
% it writes invz_tensor/tests/fixtures/projected_baseline.json, a FROZEN
% snapshot of projected reference numbers, so the invz_tensor CORE suite (which
% runs with invz_projected absent from the path) can compare against them
% without ever calling into invz_projected at test time.
%
% CLEAN-STATE SEMANTICS (v3, LOCKED -- plan Task 2 note; do not "simplify"):
%   This generator and its JSON output are themselves uncommitted while this
%   script produces the JSON, so a raw `git status --porcelain` dirty flag is
%   permanently true and meaningless (the JSON's own presence on disk, or this
%   very script before its Step-3a commit, always dirties a raw check). Instead:
%     (i)  this file is committed FIRST, ALONE, in a clean intermediate commit
%          (Step 3a) -- so the git hash recorded below corresponds to a tree
%          whose physics sources AND this generator are already committed.
%     (ii) this script records a FILTERED dirty flag = `git status --porcelain`
%          restricted to PHYSICS-SOURCE paths only: invz_projected/*.m
%          (top-level only, not tests/ or cache/), invz_common/*.m (top-level
%          only), and exactly {MF_dipole.m, exchange.m, qVec_generator.m} at
%          repo root -- EXCLUDING invz_tensor/ entirely (this script's own
%          fixtures/exploratory/cache paths) and every non-.m path.
%   A nonempty filtered status means the physics inputs changed under the
%   recorded git hash, so the baseline would not be reproducible: this script
%   STOPS (errors) and does NOT write the JSON. The JSON stores BOTH the raw
%   hash and the filtered_clean boolean so a later reader can detect physics
%   drift while tolerating the always-dirty scaffolding.
%
% GRID / CALL CONVENTIONS
%   "legacy_inclusive 8^3 / dpRng-15" grid: qVec_generator(ion.a, 'mode',
%   'grid', 'grid', [8 8 8], 'range', [-0.5 0.5]) + Gamma-exclusion filter --
%   the SAME generator + range + filter as the published Sigma_c benchmark
%   (test_invz_sigma_crit.m) and invz_odd_zero_field, just at a small n=8 /
%   dpRng=15 pair chosen for SPEED (Global Constraints: fast tests n <= 8,
%   dpRng <= 20). This is a regression-comparison anchor, not a physics-
%   converged claim -- do not expect Tc/Sigma_c numbers here to match the
%   production 12^3/24^3 @ dpRng=30 Richardson benchmarks.
%
% Run:  matlab -batch "run('.../invz_tensor/tests/exploratory/freeze_projected_baseline.m')"

format long g

here = fileparts(mfilename('fullpath'));                        % invz_tensor/tests/exploratory
repoRoot = fullfile(here, '..', '..', '..');
addpath(repoRoot);                                               % MF_dipole, exchange, qVec_generator
addpath(fullfile(repoRoot, 'invz_projected'));                   % invz_solve_point, invz_critical, invz_odd_zero_field, invz_jq_modes (dev-time only)
addpath(fullfile(repoRoot, 'invz_common'));                      % invz_ion, invz_const

fprintf('=== A0 freeze_projected_baseline  (%s) ===\n', char(datetime('now')));

% --- clean-state gate: FILTERED dirty flag over physics-source paths only --
[filteredClean, filteredLines] = filtered_physics_dirty();
[gstat, ghashRaw] = system('git rev-parse HEAD');
if gstat ~= 0
    error('invzt:gitHash', 'git rev-parse HEAD failed (exit %d): %s', gstat, ghashRaw);
end
gitHash = strtrim(ghashRaw);
fprintf('git hash = %s\n', gitHash);
fprintf('filtered_clean = %d\n', filteredClean);
if ~filteredClean
    fprintf('Filtered dirty lines (physics-source paths only):\n');
    for k = 1:numel(filteredLines), fprintf('  %s\n', filteredLines{k}); end
    error('invzt:physicsDirty', ['Physics-source paths (invz_projected/*.m, invz_common/*.m, ' ...
        'repo-root MF_dipole.m/exchange.m/qVec_generator.m) are uncommitted under the recorded ' ...
        'git hash %s: the frozen baseline would not be reproducible. STOP -- see filtered dirty ' ...
        'lines above. Commit or stash the physics inputs, then re-run.'], gitHash);
end

ion = invz_ion();
gridN  = 8;
dpRng  = 15;
gridRange = [-0.5 0.5];

% --- legacy_inclusive 8^3 grid, Gamma-excluded -------------------------------
[qvec, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', [gridN gridN gridN], ...
    'range', gridRange, 'verbose', false);
qvec = qvec(any(abs(qvec) > 1e-12, 2), :);
nq = size(qvec, 1);
fprintf('grid: legacy_inclusive, n=%d, range=[%.1f %.1f], Gamma-excluded -> nq=%d\n', ...
    gridN, gridRange(1), gridRange(2), nq);

% --- Jnu_flat + Gamma info (plain, pre-ODD path), dpRng 15 ------------------
[Jnu, info] = invz_jq_modes(ion, qvec, struct('dpRng', dpRng, 'cache', true));
Jnu_flat = Jnu(:);
fprintf('\n-- Jnu (dpRng %d, n=%d) --\n', dpRng, gridN);
fprintf('nq=%d nbranch=%d numel(Jnu_flat)=%d\n', size(Jnu,1), size(Jnu,2), numel(Jnu_flat));
fprintf('ANCHOR Jnu_flat min/max/mean = %.17g / %.17g / %.17g\n', ...
    min(Jnu_flat), max(Jnu_flat), mean(Jnu_flat));
fprintf('ANCHOR info.Jcc0 = %.17g   info.Jaa0 = %.17g\n', info.Jcc0, info.Jaa0);
fprintf('ANCHOR info.Jcc0_dipole = %.17g   info.Jaa0_dipole = %.17g\n', info.Jcc0_dipole, info.Jaa0_dipole);

% --- Sigma0/crit of invz_solve_point at (1.6 K, 0.5 T); Jxx0 = info.Jaa0 ----
Tpt = 1.6; Bpt = 0.5;
pt = invz_solve_point(ion, Tpt, Bpt, Jnu_flat, ...
    struct('hyp', true, 'J0eff', info.Jcc0, 'Jxx0', info.Jaa0));
fprintf('\n-- invz_solve_point(T=%.2gK, Bx=%.2gT), J0eff=info.Jcc0, Jxx0=info.Jaa0 --\n', Tpt, Bpt);
fprintf('ANCHOR pt.Sigma0 = %.17g\n', pt.Sigma0);
fprintf('ANCHOR pt.crit   = %.17g\n', pt.crit);
fprintf('pt.converged = %d   pt.outer_iters = %d   pt.chi0cc0 = %.17g   pt.sumrule_rel = %.6g\n', ...
    pt.converged, pt.outer_iters, pt.chi0cc0, pt.sumrule_rel);

% --- Bc(1.2 K) from invz_critical, same grid/info ---------------------------
Tbc = 1.2;
Bc = invz_critical(ion, Tbc, Jnu_flat, ...
    struct('hyp', true, 'J0eff', info.Jcc0, 'Jxx0', info.Jaa0, 'window', [0.05 7.0]));
fprintf('\n-- invz_critical(T=%.2gK) --\n', Tbc);
fprintf('ANCHOR Bc(%.2gK) = %.17g\n', Tbc, Bc);

% --- invz_odd_zero_field 'full' + 'off' on the SAME 8^3/dpRng-15 grid -------
zopts = struct('grids', gridN, 'dpRng', dpRng, 'cache', true);
[TcFull, outFull] = invz_odd_zero_field(ion, setfield(zopts, 'mode', 'full')); %#ok<SFLD>
[TcOff,  outOff ] = invz_odd_zero_field(ion, setfield(zopts, 'mode', 'off'));  %#ok<SFLD>
fprintf('\n-- invz_odd_zero_field, mode=full, n=%d, dpRng=%d --\n', gridN, dpRng);
fprintf('ANCHOR full.Tc_rich  = %.17g\n', outFull.Tc_rich);
fprintf('ANCHOR full.Sc_rich  = %.17g\n', outFull.Sc_rich);
fprintf('ANCHOR full.d_at_Tc  = %.17g\n', outFull.d_at_Tc);
fprintf('-- invz_odd_zero_field, mode=off,  n=%d, dpRng=%d --\n', gridN, dpRng);
fprintf('ANCHOR off.Tc_rich   = %.17g\n', outOff.Tc_rich);
fprintf('ANCHOR off.Sc_rich   = %.17g\n', outOff.Sc_rich);
fprintf('ANCHOR off.d_at_Tc   = %.17g\n', outOff.d_at_Tc);

% --- assemble + write the frozen JSON ---------------------------------------
S = struct();
S.schema = 'invzt_projected_baseline_v1';
S.generated_by = 'invz_tensor/tests/exploratory/freeze_projected_baseline.m';
S.date = char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss'));
S.git_hash = gitHash;
S.git_hash_note = ['HEAD at generation time; corresponds to a tree where physics sources ' ...
    'AND this generator are committed (Step 3a happens before this run, per the v3 clean-state note).'];
S.filtered_clean = filteredClean;
S.filtered_dirty_lines = {filteredLines{:}}; %#ok<CCAT1>
S.filtered_scope_note = ['git status --porcelain restricted to invz_projected/*.m (top-level), ' ...
    'invz_common/*.m (top-level), and repo-root MF_dipole.m/exchange.m/qVec_generator.m only; ' ...
    'excludes invz_tensor/ (this fixture''s own scaffolding) and every non-.m path.'];

S.grid = struct('conv', 'legacy_inclusive', 'n', gridN, 'dpRng', dpRng, ...
    'range', gridRange, 'gamma_excluded', true, 'nq', nq);

S.jq = struct('Jnu_flat', Jnu_flat(:).', 'nq', size(Jnu,1), 'nbranch', size(Jnu,2), ...
    'Jcc0', info.Jcc0, 'Jaa0', info.Jaa0, ...
    'Jcc0_dipole', info.Jcc0_dipole, 'Jaa0_dipole', info.Jaa0_dipole);

S.solve_point_1p6K_0p5T = struct('T', Tpt, 'Bx', Bpt, ...
    'J0eff', info.Jcc0, 'Jxx0', info.Jaa0, 'Jxx0_note', 'Jxx0 = info.Jaa0 (this grid/dpRng)', ...
    'Sigma0', pt.Sigma0, 'crit', pt.crit, 'converged', pt.converged, ...
    'outer_iters', pt.outer_iters, 'chi0cc0', pt.chi0cc0, 'sumrule_rel', pt.sumrule_rel);

S.critical_field_1p2K = struct('T', Tbc, 'Bc', Bc, 'window', [0.05 7.0]);

S.odd_zero_field = struct( ...
    'full', struct('Tc_rich', outFull.Tc_rich, 'Sc_rich', outFull.Sc_rich, 'd_at_Tc', outFull.d_at_Tc), ...
    'off',  struct('Tc_rich', outOff.Tc_rich,  'Sc_rich', outOff.Sc_rich,  'd_at_Tc', outOff.d_at_Tc) );

jsonText = jsonencode(S, 'PrettyPrint', true);
outDir = fullfile(here, '..', 'fixtures');
if ~exist(outDir, 'dir'), mkdir(outDir); end
outFile = fullfile(outDir, 'projected_baseline.json');
fid = fopen(outFile, 'w');
if fid < 0, error('invzt:jsonWrite', 'Could not open %s for writing.', outFile); end
fwrite(fid, jsonText, 'char');
fclose(fid);
fprintf('\nWrote %s (%d bytes)\n', outFile, numel(jsonText));
fprintf('=== freeze_projected_baseline done ===\n');

% ============================================================================
function [isClean, lines] = filtered_physics_dirty()
%FILTERED_PHYSICS_DIRTY Git status restricted to physics-source paths only.
% Physics-source = invz_projected/*.m (top-level only, NOT tests/ or cache/),
% invz_common/*.m (top-level only), and exactly {MF_dipole.m, exchange.m,
% qVec_generator.m} at repo root. Excludes invz_tensor/ (fixtures, exploratory,
% cache -- this generator's own scaffolding) entirely and every non-.m path,
% per the v3 clean-state semantics documented in this file's header.
[status, raw] = system('git status --porcelain');
if status ~= 0
    error('invzt:gitStatus', 'git status --porcelain failed (exit %d): %s', status, raw);
end
rawLines = strsplit(raw, sprintf('\n'));
lines = {};
for k = 1:numel(rawLines)
    ln = rawLines{k};
    if isempty(strtrim(ln)) || numel(ln) < 4, continue; end
    p = strtrim(ln(4:end));                          % porcelain v1: "XY path" or "XY old -> new"
    arrowIdx = strfind(p, ' -> ');
    if ~isempty(arrowIdx), p = p(arrowIdx(1)+4:end); end   % renames: keep the NEW path
    isPhysics = ~isempty(regexp(p, '^invz_projected/[^/]+\.m$', 'once')) || ...
                ~isempty(regexp(p, '^invz_common/[^/]+\.m$', 'once')) || ...
                any(strcmp(p, {'MF_dipole.m', 'exchange.m', 'qVec_generator.m'}));
    if isPhysics, lines{end+1} = ln; end %#ok<AGROW>
end
isClean = isempty(lines);
end
