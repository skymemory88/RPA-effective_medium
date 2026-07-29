function S = invzp_exec_s1_config_sweep(Bx, save_path)
%INVZP_EXEC_S1_CONFIG_SWEEP Iteration-configuration dependence of the node census.
%
% Execution packet S1 of docs/execution/invzp_plan_execution_diary.md.
% MEASUREMENT ONLY. Re-runs INVZP_EXEC_S1_FAILURE_CENSUS at one field over a set
% of damping / iteration-budget / acceleration configurations, so the reported
% failure count can be separated into (i) an iteration-dynamics component that
% moves with mix_outer/max_outer and (ii) a residual hard core that does not.
% No solver behaviour is changed and nothing here is an acceptance gate.
if nargin < 1 || isempty(Bx), Bx = 3.825; end
if nargin < 2, save_path = ''; end

cfg = {struct('mix_outer',0.20,'max_outer',1000), ...
       struct('mix_outer',0.30,'max_outer',1000), ...
       struct('mix_outer',0.40,'max_outer',1000), ...
       struct('mix_outer',0.50,'max_outer',1000), ...
       struct('mix_outer',0.70,'max_outer',200), ...
       struct('mix_outer',0.70,'max_outer',1000), ...
       struct('mix_outer',0.40,'max_outer',5000), ...
       struct('mix_outer',0.40,'max_outer',1000,'cold_acceleration','signed_aitken1')};

S = struct('Bx', Bx, 'cfg', {cfg}, 'rows', []);
rows = cell(numel(cfg),1);
for k = 1:numel(cfg)
    c = cfg{k};
    txt = evalc('Rk = invzp_exec_s1_failure_census(Bx, '''', c);'); %#ok<NASGU>
    fid = find(~Rk.tab.accepted).';
    accel = 'none';
    if isfield(c,'cold_acceleration'), accel = c.cold_acceleration; end
    rows{k} = struct('mix_outer', c.mix_outer, 'max_outer', c.max_outer, ...
        'accel', accel, 'status', Rk.meta.hmf_status, ...
        'n_failed', Rk.meta.n_failed, 'n_nodes', Rk.meta.n_nodes, ...
        'failed_ids', fid, 'hstar', Rk.meta.hstar, 'wall_s', Rk.meta.wall_s);
    fprintf('CFG %d: mix=%.2f max=%5d accel=%-14s -> status=%-14s failed=%2d/%d hstar=%.6g ids=[%s]\n', ...
        k, c.mix_outer, c.max_outer, accel, Rk.meta.hmf_status, ...
        Rk.meta.n_failed, Rk.meta.n_nodes, Rk.meta.hstar, num2str(fid));
end
S.rows = [rows{:}];

allfail = {};
for k = 1:numel(S.rows), allfail{end+1} = S.rows(k).failed_ids; end %#ok<AGROW>
core = allfail{1};
for k = 2:numel(allfail), core = intersect(core, allfail{k}); end
S.hard_core_ids = core;
fprintf('\nHARD CORE (failed under EVERY configuration above): [%s]\n', num2str(core));

if ~isempty(save_path), save(save_path, 'S'); fprintf('saved: %s\n', save_path); end
end
