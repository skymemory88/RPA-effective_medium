% Task-15 final fix pass -- I2 hunt round 3: an UNCOUNTED masked reason, ideally with a
% FINITE Bc_1z (which would make the pre-fix sweep emit no warning at all).
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));

ion = invz_ion();
w   = (0.02:0.04:0.42).';
syn = @(s) struct('Jnu', s*linspace(-2e-3, 6.0e-3, 24).', 'info', struct('Jcc0', s*6.4e-3), ...
                  'verbose', false);
j   = @(c) strjoin(c, ' | ');
UNCOUNTED = {'unstable_endpoint', 'solver_failed', 'not_attempted_longitudinal', ...
             'bare_not_ordered', 'response_failed'};

function report(nm, S, UNCOUNTED, j)
    nmsk = nnz(S.phase_1z == 0);
    unc  = nnz(ismember(S.phase_1z_reason, UNCOUNTED));
    fprintf('%-26s Bc_1z=%-8.4g %-12s n_msk=%d  UNCOUNTED=%d  reason = %s\n', ...
            nm, S.Bc_1z, S.Bc_1z_status, nmsk, unc, j(S.phase_1z_reason));
    if unc > 0 && isfinite(S.Bc_1z)
        fprintf('   *** JACKPOT: uncounted masked reason WITH a finite Bc_1z ***\n');
    end
end

ws = warning('off', 'invz:spectraMapMasked');

%% ---- (a) tilt + bz_tol split: low fields transverse, tail longitudinal -------------------
fprintf('########## (a) split tilt ##########\n');
for th = [3 20 45 80]
  for Bhi = [9 20 60]
    fd = [cosd(th) 0 sind(th)];
    F  = [2.85 3.30 5.50 Bhi];
    tol = 0.5*(5.50*sind(th) + Bhi*sind(th));   % between the 5.5 T and Bhi z-components
    o = syn(1);  o.static_medium = 'strict_1z_dyson_ref';  o.field_dir = fd;  o.bz_tol = tol;
    try
        S = invz_spectra_map(ion, 0.31, F, w, o);
        report(sprintf('tilt%02d Bhi=%g', th, Bhi), S, UNCOUNTED, j);
    catch ME
        fprintf('tilt%02d Bhi=%-4g THREW %s\n', th, Bhi, ME.identifier);
    end
  end
end

%% ---- (b) weak / strong couplings with a B = 0 column -------------------------------------
fprintf('\n########## (b) coupling strength x B = 0 ##########\n');
for s = [0.02 0.1 0.3 0.6 1.0 2.0]
    o = syn(s);  o.static_medium = 'strict_1z_dyson_ref';
    try
        S = invz_spectra_map(ion, 0.31, [0 2.85 3.30 5.50], w, o);
        report(sprintf('scale=%g', s), S, UNCOUNTED, j);
    catch ME
        fprintf('scale=%-6g THREW %s\n', s, ME.identifier);
    end
end

%% ---- (c) hyp = false ----------------------------------------------------------------------
fprintf('\n########## (c) hyp = false ##########\n');
for s = [0.3 1.0]
    o = syn(s);  o.static_medium = 'strict_1z_dyson_ref';  o.hyp = false;
    try
        S = invz_spectra_map(ion, 0.31, [0 2.85 3.30 5.50], w, o);
        report(sprintf('hypfalse scale=%g', s), S, UNCOUNTED, j);
    catch ME
        fprintf('hypfalse scale=%-4g THREW %s\n', s, ME.identifier);
    end
end
warning(ws);
fprintf('DONE_SCAN3\n');
