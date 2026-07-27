ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));

ion = invz_ion();  T = 0.31;  w = (0.02:0.04:0.42).';
o = struct('grid', [4 4 4], 'dpRng', 12, 'cache', false, 'verbose', false);

fprintf('--- map [0 3 8] default opts ---\n');
try
    S = invz_spectra_map(ion, T, [0 3 8], w, o);
    fprintf('phase    = %s\n', mat2str(S.phase));
    fprintf('phase_1z = %s\n', mat2str(S.phase_1z));
    fprintf('crit_pm  = %s\n', mat2str(S.crit_pm, 6));
    fprintf('Sigma0   = %s\n', mat2str(S.Sigma0, 6));
    fprintf('Bc_1z    = %g   Bc_auto = %g\n', S.Bc_1z, S.Bc_auto);
catch err
    fprintf('MAP ERROR %s : %s\n', err.identifier, err.message);
end

% Direct probe of the B = 0 column's internals, on the SAME couplings.
fprintf('\n--- B = 0 leg-by-leg ---\n');
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 12, 'cache', false));
sopts = struct('hyp', true, 'J0eff', info.Jcc0, 'Jxx0', ion.Jxx0, 'bz_tol', 1e-9);
B0 = [0 0 0];
[pt, phase, di] = invz_solve_auto(ion, T, B0, Jnu, sopts);
fprintf('auto: phase=%d  ordered_err=%s  para_err=%s  empty(pt)=%d\n', ...
        phase, di.ordered_err, di.para_err, isempty(pt));
if ~isempty(pt)
    fprintf('auto pt: is_ordered=%d converged=%d Sigma0=%g m0=%g\n', ...
        getf(pt,'is_ordered',false), pt.converged, pt.Sigma0, getf(pt,'m0',NaN));
end
try
    ptp = invz_solve_point(ion, T, B0, Jnu, sopts);
    fprintf('PM probe: converged=%d crit=%g Sigma0=%g\n', ptp.converged, ptp.crit, ptp.Sigma0);
catch err
    fprintf('PM probe THREW %s\n', err.identifier);
end
so2 = sopts;  so2.ordered_mode = 'jensen';
try
    ptj = invz_solve_point_ordered(ion, T, B0, Jnu, so2);
    fprintf('jensen: is_ordered=%d converged=%d Sigma0=%g stable_1z=%d medium_status=%s hmf_status=%s\n', ...
        ptj.is_ordered, ptj.converged, ptj.Sigma0, getf(ptj,'stable_1z',false), ...
        getf(ptj,'medium_status','?'), getf(ptj,'hmf_status','(n/a)'));
catch err
    fprintf('jensen THREW %s : %s\n', err.identifier, err.message);
end

% and at 3 T / 8 T for the strict dispatcher's own anchors
for B = [3 8]
    fprintf('\n--- B = %g leg-by-leg ---\n', B);
    Bv = [B 0 0];
    [pt2, ph2, di2] = invz_solve_auto(ion, T, Bv, Jnu, sopts);
    fprintf('auto: phase=%d ordered_err=%s para_err=%s\n', ph2, di2.ordered_err, di2.para_err);
    try
        ptp2 = invz_solve_point(ion, T, Bv, Jnu, sopts);
        fprintf('PM: converged=%d crit=%g\n', ptp2.converged, ptp2.crit);
    catch err
        fprintf('PM THREW %s\n', err.identifier);
    end
end
