function diag11()
% Independent check of the QCP claim now consolidated in invzp_convg_diagnosis_Claude.md
% SS9.5 that the RESPONSE evaluator is
% healthy once it is handed a valid state, i.e. the failure is state construction, not chi.
%   (a) production copts carry no 'Jfull'  -> invz_chi_realaxis is single-shot (static check
%       done by grep; here we confirm the npass loop is a no-op by comparing npass 1 vs 3)
%   (b) bare ordered point at 1 T   -> count finite chi'' samples
%   (c) jensen ordered point at 4 T -> count finite chi'' samples
%   (d) same at 3.8 T (window edge) and a PM point at 4.6 T as controls
%   (e) exhibit the Gamma-gap: min_q Dq - D_uni = (J0eff - max_q Jnu)*|Gstat| > 0, so a state
%       can have D_uni < 0 (uniform mode unstable) while every GRID mode is still evaluable.
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root); addpath(fullfile(root,'invz_projected')); addpath(fullfile(root,'invz_common'));
warning('off','MATLAB:nearlySingularMatrix');
SD = '/private/tmp/claude-501/-Users-yikaiyang-Library-CloudStorage-OneDrive-Nexus365-Programming-scripts-Matlab-Simulation-invZ-expansion/fe404c1e-1dd6-4ec7-868e-75af95bf66a2/scratchpad';
LOG = fullfile(SD,'diag11.log'); if exist(LOG,'file'), delete(LOG); end

ion = invz_ion();
[Jnu, info] = invz_bz_couplings(ion, struct('grid',[16 16 16],'dpRng',30,'cache',true));
Jf = Jnu(:); Jaa0 = ion.Jxx0; if isfield(info,'Jaa0'), Jaa0 = info.Jaa0; end
T = 0.31;
w = linspace(0.0, 0.6, 60).';            % 60-point real-axis grid, as in the QCP doc
copts0 = struct('Jsel', info.Jcc0, 'eta', 5e-3, 'Jxx0', Jaa0, 'Jshape', 0, ...
                'hyp', true, 'transverse_mf', 'legacy_x');

say(LOG,'--- (b,c,d) invz_chi_realaxis given a VALID state, T = 0.31 K, 60-point w grid ---');
say(LOG,'  B   mode      solved  Sigma0        finite chi''''   finite Sigma(w)   max|chi''''|');

cases = { 1.0,'bare'; 3.8,'jensen'; 4.0,'jensen'; 4.6,'pm' };
for k = 1:size(cases,1)
    B = cases{k,1}; mode = cases{k,2};
    so = struct('J0eff', info.Jcc0, 'Jxx0', Jaa0);
    switch mode
        case 'pm'
            pt = invz_solve_point(ion, T, [B 0 0], Jf, so);
            solved = pt.converged;
        otherwise
            so.ordered_mode = mode;
            pt = invz_solve_point_ordered(ion, T, [B 0 0], Jf, so);
            solved = pt.converged && pt.is_ordered;
    end
    if ~solved
        say(LOG,sprintf(' %4.1f  %-8s  NO      (state not available; hmf_status=%s)', ...
            B, mode, getf(pt,'hmf_status','-')));
        continue
    end
    o = invz_chi_realaxis(ion, T, [B 0 0], pt, w, copts0);
    ch = imag(o.chi_cc_q(1,:)).';
    say(LOG,sprintf(' %4.1f  %-8s  yes     %-12.6g  %2d / %2d      %2d / %2d        %.6g', ...
        B, mode, pt.Sigma0, sum(isfinite(ch)), numel(ch), ...
        sum(isfinite(o.Sigma_w)), numel(o.Sigma_w), max(abs(ch))));
    % (a) npass has no effect without Jfull
    c3 = copts0; c3.npass = 12;
    o3 = invz_chi_realaxis(ion, T, [B 0 0], pt, w, c3);
    say(LOG,sprintf('        npass 3 -> 12 : max|d chi''''| = %.3e  (0 confirms no lattice loop runs)', ...
        max(abs(imag(o3.chi_cc_q(1,:)).' - ch))));
end

% ---- (e) the Gamma gap, exactly ------------------------------------------------------
say(LOG,'');
say(LOG,'--- (e) min_q Dq - D_uni = (J0eff - max_q Jnu)*|Gstat| ---');
say(LOG,sprintf('  J0eff        = %.9g meV', info.Jcc0));
say(LOG,sprintf('  max_q Jnu    = %.9g meV', max(Jf)));
say(LOG,sprintf('  J0eff-max_q  = %.9g meV   (ratio J0eff/max_q = %.6f)', ...
    info.Jcc0-max(Jf), info.Jcc0/max(Jf)));
say(LOG,sprintf('  1/J0eff      = %.5f     1/max_q Jnu = %.5f   (the record''s |G0| ladder)', ...
    1/info.Jcc0, 1/max(Jf)));
soj = struct('J0eff', info.Jcc0, 'Jxx0', Jaa0, 'ordered_mode','jensen');
pj = invz_solve_point_ordered(ion, T, [4.0 0 0], Jf, soj);
if pj.converged && pj.is_ordered
    Gstat = pj.G(1); K0 = pj.K(1);
    Dq = 1 + (Jf - K0).*Gstat;  Duni = 1 + (info.Jcc0 - K0).*Gstat;
    say(LOG,sprintf('  at the B=4.0 T jensen endpoint: Gstat = %.6g', Gstat));
    say(LOG,sprintf('    D_uni = %+.6g    min_q Dq = %+.6g    gap = %.6g  (predicted %.6g)', ...
        Duni, min(Dq), min(Dq)-Duni, (info.Jcc0-max(Jf))*abs(Gstat)));
end
say(LOG,'done');
end
function v = getf(s,f,d), if isstruct(s)&&isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end
function say(LOG,s), fid=fopen(LOG,'a'); fprintf(fid,'%s\n',s); fclose(fid); end
