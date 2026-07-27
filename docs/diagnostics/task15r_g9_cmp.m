% G9 comparison: worktree vs a54de5a (and vs b9082fd for the empty-sweep regression).
SP = ['/private/tmp/claude-503/-Users-yikaiyang-Library-CloudStorage-OneDrive-Nexus365-' ...
      'Programming-scripts-Matlab-Simulation-invZ-expansion/' ...
      'e72c63b7-dc11-4663-9d4b-867feace2d26/scratchpad'];
A = load(fullfile(SP, 'g9_worktree.mat'));  A = A.out;
B = load(fullfile(SP, 'g9_a54de5a.mat'));   B = B.out;
C = load(fullfile(SP, 'g9_b9082fd.mat'));   C = C.out;

% every numeric/label field a54de5a already had -- the new ones cannot be compared, they are new
flds = {'chiz','chirpa','Sigma0','phase','phase_1z','crit_pm','m_1z','D_ord','suspect', ...
        'Bc_auto','Bc_1z','Epeak','Epeak_rpa','fields','phase_1z_reason','stability_1z', ...
        'pm_probe_status','pm_probe_error_id','static_medium_used','static_medium', ...
        'Bc_1z_interval','Bc_1z_status'};
names = fieldnames(A);
verdict = true;
for k = 1:numel(names)
    n = names{k};
    a = A.(n);  b = B.(n);
    if a.threw ~= b.threw
        if a.threw, fprintf('%-14s DIFF worktree THREW %s / a54de5a ok\n', n, a.id);
        else,       fprintf('%-14s DIFF worktree ok / a54de5a THREW %s\n', n, b.id);  end
        % the ONLY intended difference: the F1 regression fix
        if strcmp(n, 'deg_empty') && ~a.threw && b.threw && ...
                strcmp(b.id, 'invz:boundaryInterval')
            fprintf('%-14s ^^^ INTENDED (F1 fix). b9082fd threw = %d -> restored\n', ...
                    n, C.(n).threw);
        else
            verdict = false;
        end
        continue;
    end
    if a.threw
        fprintf('%-14s both threw %s / %s\n', n, a.id, b.id);
        if ~strcmp(a.id, b.id), verdict = false; end
        continue;
    end
    bad = {};
    for j = 1:numel(flds)
        f = flds{j};
        if ~isfield(a.S, f) || ~isfield(b.S, f)
            if isfield(a.S, f) ~= isfield(b.S, f), bad{end+1} = [f '(presence)']; end %#ok<AGROW>
            continue;
        end
        if ~isequaln(a.S.(f), b.S.(f)), bad{end+1} = f; end %#ok<AGROW>
    end
    if isempty(bad)
        fprintf('%-14s BITWISE IDENTICAL (%d fields via isequaln)\n', n, numel(flds));
    else
        fprintf('%-14s DIFFERS IN: %s\n', n, strjoin(bad, ', '));
        verdict = false;
    end
end
fprintf('G9VERDICT %d\n', verdict);

% and the new schema fields are present with the right per-column shape everywhere
for k = 1:numel(names)
    n = names{k};
    if A.(n).threw, continue; end
    S = A.(n).S;
    ok = isfield(S,'ordered_diag_reason') && isfield(S,'response_error_id') && ...
         numel(S.ordered_diag_reason) == numel(S.phase_1z) && ...
         numel(S.response_error_id)   == numel(S.phase_1z);
    fprintf('NEWFIELDS %-14s ok=%d  diag = %s\n', n, ok, ...
            strjoin(S.ordered_diag_reason, ' | '));
end
fprintf('G9CMPDONE\n');
