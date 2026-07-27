% G9 harness: run the same configuration matrix under EITHER the worktree invz_spectra_map.m or
% a shadowed historical copy, and dump the results. Mode + shadow dir + output file come from
% env vars so the two sides run in SEPARATE MATLAB processes (no rehash/shadowing games).
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));
shadow = getenv('G9_SHADOW');
if ~isempty(shadow), addpath(shadow, '-begin'); end
fprintf('G9 using invz_spectra_map from: %s\n', which('invz_spectra_map'));

ion = invz_ion();
w   = (0.02:0.02:0.6).';
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
syn  = struct('Jnu', Jnu, 'info', info, 'verbose', false);

cfg = {};
cfg{end+1} = {'synth_jensen', [2.85 3.30 5.50], setfield(syn, 'ordered_1z', 'jensen')};
cfg{end+1} = {'synth_bare',   [2.85 3.30 5.50], setfield(syn, 'ordered_1z', 'bare')};
tl = syn;  tl.field_dir = [cosd(3) 0 sind(3)];
cfg{end+1} = {'synth_tilt3',  [2.85 5.50], tl};
cfg{end+1} = {'real_zero3_8', [0 3 8], ...
              struct('grid',[4 4 4],'dpRng',12,'cache',false,'verbose',false)};
% DEGENERATE INPUTS -- the class the original four-configuration spot check missed (review F1)
cfg{end+1} = {'deg_empty',      [],            syn};
cfg{end+1} = {'deg_single',     [5.50],        syn};
cfg{end+1} = {'deg_unsorted',   [5.50 2.85],   syn};
cfg{end+1} = {'deg_duplicate',  [3.30 3.30],   syn};
cfg{end+1} = {'deg_zero_only',  [0],           ...
              struct('grid',[4 4 4],'dpRng',12,'cache',false,'verbose',false)};

out = struct();
for k = 1:numel(cfg)
    name = cfg{k}{1};
    try
        S = invz_spectra_map(ion, 0.31, cfg{k}{2}, w, cfg{k}{3});
        out.(name) = struct('threw', false, 'S', S);
        fprintf('G9RUN %-14s OK\n', name);
    catch e
        out.(name) = struct('threw', true, 'id', e.identifier, 'msg', e.message);
        fprintf('G9RUN %-14s THREW %s\n', name, e.identifier);
    end
end
save(getenv('G9_OUT'), 'out');
fprintf('G9RUNDONE\n');
