function E = invzp_exec_e3_root_sensitivity(save_path)
%INVZP_EXEC_E3_ROOT_SENSITIVITY Is the published 3.825 T root numerically determined?
%
% Execution packet E3 of docs/execution/invzp_plan_execution_diary.md.
% MEASUREMENT ONLY.
%
% Two independent perturbation axes are applied to a configuration that DOES close
% the 3.825 T column (mix_outer = 0.30), and the published root h* is compared:
%
%   axis A -- ARITHMETIC. closure_coordinate factored vs defactored. These are exact
%     algebraic reassociations of the same closure; the gate
%     invzp_exec_s2_defactor_gate.m measures their disagreement at <= 8e-13 relative
%     on a real node fixture. Any h* difference far above that is amplification, not
%     a different equation.
%
%   axis B -- ITERATION DAMPING. mix_outer over a window around 0.30. Damping cannot
%     move the fixed-point set; if h* moves with it, the solve is selecting among
%     several fixed points rather than converging to one.
%
% For axis A the full per-node tables are retained so the divergence can be located:
% a SMOOTH growth of the node-by-node difference indicates accumulated quadrature in
% h0 = int r dh', while a JUMP at one node indicates the two runs took different
% basins there.
if nargin < 1, save_path = ''; end

Bx = 3.825;
base = struct('mix_outer', 0.30, 'max_outer', 1000);

% --- axis A: arithmetic coordinate, tables retained ----------------------------------------
cfgF = base;  cfgF.emt_static = struct('closure_coordinate', 'factored');
cfgD = base;  cfgD.emt_static = struct('closure_coordinate', 'defactored');
txt = evalc('RF = invzp_exec_s1_failure_census(Bx, '''', cfgF);'); %#ok<NASGU>
txt = evalc('RD = invzp_exec_s1_failure_census(Bx, '''', cfgD);'); %#ok<NASGU>

hF = RF.meta.hstar;  hD = RD.meta.hstar;
relA = abs(hD - hF)/abs(hF);

nmin = min(height(RF.tab), height(RD.tab));
cmp = table((1:nmin).', RF.tab.h(1:nmin), RD.tab.h(1:nmin), ...
    abs(RD.tab.h(1:nmin) - RF.tab.h(1:nmin))./max(abs(RF.tab.h(1:nmin)), realmin), ...
    abs(RD.tab.rr(1:nmin) - RF.tab.rr(1:nmin))./max(abs(RF.tab.rr(1:nmin)), realmin), ...
    abs(RD.tab.mm(1:nmin) - RF.tab.mm(1:nmin))./max(abs(RF.tab.mm(1:nmin)), realmin), ...
    abs(RD.tab.K0(1:nmin) - RF.tab.K0(1:nmin))./max(abs(RF.tab.K0(1:nmin)), realmin), ...
    RF.tab.iters(1:nmin), RD.tab.iters(1:nmin), ...
    'VariableNames', {'node','h_fac','h_def','rel_dh','rel_dr','rel_dm','rel_dK0', ...
                      'iters_fac','iters_def'});

% --- axis B: damping window ----------------------------------------------------------------
mixes = [0.26 0.28 0.30 0.32 0.34 0.36];
hB = nan(size(mixes));  stB = cell(size(mixes));  nfB = nan(size(mixes));
for k = 1:numel(mixes)
    c = base;  c.mix_outer = mixes(k);
    txt = evalc('Rk = invzp_exec_s1_failure_census(Bx, '''', c);'); %#ok<NASGU>
    hB(k) = Rk.meta.hstar;  stB{k} = Rk.meta.hmf_status;  nfB(k) = Rk.meta.n_failed;
end
okB = strcmp(stB, 'ok');
spreadB = NaN;
if nnz(okB) >= 2
    spreadB = (max(hB(okB)) - min(hB(okB)))/mean(hB(okB));
end

E = struct('Bx', Bx, 'axisA', struct('h_factored', hF, 'h_defactored', hD, ...
        'rel_diff', relA, 'cmp', cmp), ...
    'axisB', struct('mix', mixes, 'hstar', hB, 'status', {stB}, 'n_failed', nfB, ...
        'rel_spread_over_closing_configs', spreadB));

fprintf('=== E3 root sensitivity at Bx = %.3f T ===\n', Bx);
fprintf('AXIS A (arithmetic coordinate, mix_outer = %.2f):\n', base.mix_outer);
fprintf('  h*(factored)   = %.12g\n  h*(defactored) = %.12g\n  relative difference = %.4g\n', ...
    hF, hD, relA);
fprintf('  first node whose |rel dh| exceeds 1e-9: ');
j = find(cmp.rel_dh > 1e-9, 1, 'first');
if isempty(j), fprintf('none\n'); else, fprintf('%d (h = %.6g)\n', j, cmp.h_fac(j)); end
fprintf('  per-node divergence (nodes where any rel diff > 1e-12):\n');
sel = cmp.rel_dh > 1e-12 | cmp.rel_dr > 1e-12 | cmp.rel_dm > 1e-12 | cmp.rel_dK0 > 1e-12;
disp(cmp(sel, :));
fprintf('AXIS B (damping window):\n');
for k = 1:numel(mixes)
    fprintf('  mix = %.2f -> %-12s failed = %2d  h* = %.12g\n', ...
        mixes(k), stB{k}, nfB(k), hB(k));
end
fprintf('  relative spread of h* over the %d closing configurations: %.4g\n', ...
    nnz(okB), spreadB);

if ~isempty(save_path), save(save_path, 'E', 'RF', 'RD'); end
end
