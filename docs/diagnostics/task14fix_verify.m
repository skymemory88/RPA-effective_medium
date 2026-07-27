function task14fix_verify()
% (1) bitwise pre/post comparison of the DEFAULT resummed bare path, and of the strict OK path;
% (2) the finding-2 field-presence contract, measured across every bare/jensen return path.
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT, 'invz_common'));
addpath(genpath(fullfile(ROOT, 'invz_projected')));

pre  = load(fullfile(ROOT, '.superpowers', 'sdd', 'task14fix_pre.mat'));
post = load(fullfile(ROOT, '.superpowers', 'sdd', 'task14fix_post.mat'));
nmis = 0;
for nm = {'B', 'C', 'D'}          % B = strict OK path, C/D = default resummed path
    a = pre.(nm{1});  b = post.(nm{1});
    fa = sort(fieldnames(a));  fb = sort(fieldnames(b));
    same = isequal(fa, fb);
    for i = 1:numel(fa)
        same = same && isequaln(a.(fa{i}), b.(fa{i}));   % NO tolerance
    end
    nmis = nmis + ~same;
    fprintf('G9 %-2s nfields=%d bitwise_equal=%d\n', nm{1}, numel(fa), same);
end
fprintf('G9COMPARE mismatches=%d\n', nmis);

% ---- finding 2: the field-presence contract -----------------------------------------------
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o   = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
ob  = o;  ob.static_medium = 'strict_1z_dyson_ref';
oh  = ob; oh.ref_margin = 2;
oj  = o;  oj.ordered_mode = 'jensen';

cases = {
    'bare_full',      invz_solve_point_ordered(ion, 0.31, [2.85 0 0], Jnu, o)
    'bare_strict_ok', invz_solve_point_ordered(ion, 0.31, [2.85 0 0], Jnu, ob)
    'bare_halt',      invz_solve_point_ordered(ion, 0.31, [2.85 0 0], Jnu, oh)
    'bare_early',     invz_solve_point_ordered(ion, 2.0,  [2.85 0 0], Jnu, o)
    'jensen_full',    invz_solve_point_ordered(ion, 0.31, [2.85 0 0], Jnu, oj)
    'jensen_early',   invz_solve_point_ordered(ion, 2.0,  [2.85 0 0], Jnu, oj)
};
jonly = {'stable_1z','crit_1z','Dq_min','D_uni','omit_mu3','omit_cubic','omit_max','path_omit_max'};
both  = {'static_medium','medium_status','medium_denom','medium_margin','Jmom'};
for i = 1:size(cases, 1)
    p = cases{i, 2};
    fprintf('F2 %-15s jensen_only_present=%d/%d  common_present=%d/%d\n', ...
        cases{i,1}, nnz(cellfun(@(f) isfield(p, f), jonly)), numel(jonly), ...
        nnz(cellfun(@(f) isfield(p, f), both)),  numel(both));
end
fprintf('VERIFY DONE\n');
end
