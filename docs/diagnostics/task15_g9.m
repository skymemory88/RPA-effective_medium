% G9 regression capture: run the SAME resummed configurations under the current worktree and
% under HEAD's invz_spectra_map/invz_solve_auto (set G9_HEAD=1 to shadow with the HEAD copies).
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
SP = ['/private/tmp/claude-503/-Users-yikaiyang-Library-CloudStorage-OneDrive-Nexus365-' ...
      'Programming-scripts-Matlab-Simulation-invZ-expansion/' ...
      'e72c63b7-dc11-4663-9d4b-867feace2d26/scratchpad'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));
isHead = ~isempty(getenv('G9_HEAD'));
if isHead
    addpath(fullfile(SP, 'g9head'), '-begin');   % HEAD copies shadow the worktree versions
end
tag = 'WORKTREE';  if isHead, tag = 'HEAD'; end
fprintf('G9 capture: %s (invz_spectra_map from %s)\n', tag, which('invz_spectra_map'));

ion = invz_ion();  T = 0.31;
w  = (0.02:0.02:0.6).';
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
th = 3*pi/180;

cfg = {};
cfg{end+1} = struct('name', 'synth_jensen',  'fields', [2.85 3.30 5.50], ...
                    'opts', struct('Jnu', Jnu, 'info', info, 'verbose', false));
cfg{end+1} = struct('name', 'synth_bare',    'fields', [2.85 3.30 5.50], ...
                    'opts', struct('Jnu', Jnu, 'info', info, 'verbose', false, ...
                                   'ordered_1z', 'bare'));
cfg{end+1} = struct('name', 'synth_tilt',    'fields', [2.85 5.50], ...
                    'opts', struct('Jnu', Jnu, 'info', info, 'verbose', false, ...
                                   'field_dir', [cos(th) 0 sin(th)]));
cfg{end+1} = struct('name', 'real_zero3_8',  'fields', [0 3 8], ...
                    'opts', struct('grid', [4 4 4], 'dpRng', 12, 'cache', false, ...
                                   'verbose', false));

out = struct();
for c = 1:numel(cfg)
    nm = cfg{c}.name;
    try
        S = invz_spectra_map(ion, T, cfg{c}.fields, w, cfg{c}.opts);
        keep = {'chiz','chirpa','Sigma0','phase','phase_1z','crit_pm','m_1z','D_ord', ...
                'suspect','Bc_auto','Bc_1z','Epeak','Epeak_rpa','fields'};
        r = struct();
        for f = keep, r.(f{1}) = S.(f{1}); end
        out.(nm) = r;
        fprintf('%-14s OK  phase_1z=%s Bc_1z=%g\n', nm, mat2str(S.phase_1z), S.Bc_1z);
    catch err
        out.(nm) = struct('errid', err.identifier, 'errmsg', err.message);
        fprintf('%-14s THREW %s\n', nm, err.identifier);
    end
end

% --- error-policy probe: a NON-whitelisted invz:* identifier reachable through the solver
% catches. A NaN in the coupling multiset makes invz_coupling_moments throw
% invz:couplingMoments inside BOTH point solvers (task 14 made it unconditional).
JnuBad = Jnu;  JnuBad(3) = NaN;
try
    invz_spectra_map(ion, T, 5.5, w, struct('Jnu', JnuBad, 'info', info, 'verbose', false));
    out.badJnu = struct('errid', '(none -- absorbed)', 'errmsg', '');
    fprintf('badJnu         ABSORBED (no error)\n');
catch err
    out.badJnu = struct('errid', err.identifier, 'errmsg', err.message);
    fprintf('badJnu         THREW %s\n', err.identifier);
end
% and the same probe one level down, at invz_solve_auto itself
try
    [~, ph, di] = invz_solve_auto(ion, T, 5.5, JnuBad, struct('J0eff', 6.4e-3));
    out.badJnuAuto = struct('errid', sprintf('(none) phase=%d ord=%s para=%s', ph, ...
                            di.ordered_err, di.para_err), 'errmsg', '');
    fprintf('badJnuAuto     ABSORBED phase=%d ordered_err=%s para_err=%s\n', ...
            ph, di.ordered_err, di.para_err);
catch err
    out.badJnuAuto = struct('errid', err.identifier, 'errmsg', err.message);
    fprintf('badJnuAuto     THREW %s\n', err.identifier);
end
% whitelisted recoverable signal must STILL be absorbed (invz:degenerateDoublet)
try
    [pt, ph, di] = invz_solve_auto(ion, 1.9, 0.04, Jnu, struct('J0eff', 6.4e-3));
    out.degen = struct('errid', sprintf('(none) phase=%d para=%s empty=%d', ph, ...
                       di.para_err, isempty(pt)), 'errmsg', '');
    fprintf('degenDoublet   ABSORBED phase=%d para_err=%s\n', ph, di.para_err);
catch err
    out.degen = struct('errid', err.identifier, 'errmsg', err.message);
    fprintf('degenDoublet   THREW %s\n', err.identifier);
end

if isHead, save(fullfile(SP, 'g9_head.mat'), 'out');
else,      save(fullfile(SP, 'g9_new.mat'),  'out');  end
