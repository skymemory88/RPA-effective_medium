function S = invzp_exec_s2_column_effect(save_path)
%INVZP_EXEC_S2_COLUMN_EFFECT Column-level effect of the defactored closure coordinate.
%
% Execution packet S2 of docs/execution/invzp_plan_execution_diary.md.
% Runs the 3.825 T node census under matched (mix_outer, max_outer) pairs with the
% closure coordinate factored (production default) and defactored (opt-in), so the
% mask count attributable to the pole-crossing arithmetic can be read off directly.
% The first row also serves as the BEHAVIOUR-NEUTRALITY regression: with the default
% coordinate the failure set must be exactly {23}, as recorded before the change.
if nargin < 1, save_path = ''; end
cfgs = {struct('mix_outer',0.40,'max_outer',1000), ...
        struct('mix_outer',0.70,'max_outer',1000), ...
        struct('mix_outer',0.70,'max_outer',200), ...
        struct('mix_outer',0.30,'max_outer',1000)};
coords = {'factored','defactored'};
rows = {};
for a = 1:numel(cfgs)
    for b = 1:numel(coords)
        c = cfgs{a};
        c.emt_static = struct('closure_coordinate', coords{b});
        txt = evalc('Rk = invzp_exec_s1_failure_census(3.825, '''', c);'); %#ok<NASGU>
        fid = find(~Rk.tab.accepted).';
        rows{end+1} = struct('mix_outer',c.mix_outer,'max_outer',c.max_outer, ...
            'coord',coords{b},'status',Rk.meta.hmf_status,'n_failed',Rk.meta.n_failed, ...
            'n_nodes',Rk.meta.n_nodes,'hstar',Rk.meta.hstar,'failed_ids',fid); %#ok<AGROW>
        fprintf('mix=%.2f max=%5d coord=%-11s -> %-12s failed=%2d/%d hstar=%.8g ids=[%s]\n', ...
            c.mix_outer, c.max_outer, coords{b}, Rk.meta.hmf_status, ...
            Rk.meta.n_failed, Rk.meta.n_nodes, Rk.meta.hstar, num2str(fid));
    end
end
S = [rows{:}];
if ~isempty(save_path), save(save_path,'S'); end
end
