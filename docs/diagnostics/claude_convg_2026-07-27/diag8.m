function diag8()
% Where is the 1/z PM boundary, really? Cross-check part (B) with the PRODUCTION
% boundary finder, and locate the temperature at which Bc = 3.03 T (the value the
% Stage-2 record and the Gate-0 control placement assume).
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root); addpath(fullfile(root,'invz_projected')); addpath(fullfile(root,'invz_common'));
warning('off','MATLAB:nearlySingularMatrix');
SD = '/private/tmp/claude-501/-Users-yikaiyang-Library-CloudStorage-OneDrive-Nexus365-Programming-scripts-Matlab-Simulation-invZ-expansion/fe404c1e-1dd6-4ec7-868e-75af95bf66a2/scratchpad';
LOG = fullfile(SD,'diag8.log'); if exist(LOG,'file'), delete(LOG); end
ion = invz_ion();
[Jnu, info] = invz_bz_couplings(ion, struct('grid',[16 16 16],'dpRng',30,'cache',true));
Jf = Jnu(:); Jaa0 = ion.Jxx0; if isfield(info,'Jaa0'), Jaa0 = info.Jaa0; end

% (1) production boundary finder at the Gate-0 temperature and at mine
say(LOG,'--- (1) invz_critical (production PM-side Bc finder), resummed default ---');
for T = [0.10 0.31]
    try
        bx = invz_critical(ion, T, Jf, struct('J0eff',info.Jcc0,'Jxx0',Jaa0,'window',[7.5 2.0]));
        say(LOG, sprintf('  T=%.2f K :  Bc = %.4f T', T, bx));
    catch ME
        say(LOG, sprintf('  T=%.2f K :  THREW %s : %s', T, ME.identifier, ME.message));
    end
end

% (2) PM crit(B) at several T -- which T puts the crit-zero at ~3.03 T?
say(LOG,'');
say(LOG,'--- (2) PM crit(B) vs T (invz_solve_point, resummed) : locate Bc(T) ---');
say(LOG,'   T(K)     Bc_interp(T)   [samples: B -> crit (conv)]');
for T = [0.10 0.31 0.60 0.90 1.20 1.50]
    Bs = 2.0:0.5:7.0;  cr = nan(size(Bs)); cv = false(size(Bs));
    for k = 1:numel(Bs)
        try
            p = invz_solve_point(ion, T, [Bs(k) 0 0], Jf, struct('J0eff',info.Jcc0,'Jxx0',Jaa0));
            if p.converged, cr(k) = p.crit; cv(k) = true; end
        catch, end
    end
    Bc = NaN;
    ok = find(cv & isfinite(cr));
    for j = 1:numel(ok)-1
        a = ok(j); b = ok(j+1);
        if cr(a) < 0 && cr(b) >= 0
            Bc = Bs(a) + (Bs(b)-Bs(a))*(-cr(a))/(cr(b)-cr(a));
        end
    end
    s = '';
    for k = 1:numel(Bs)
        if cv(k), s = [s sprintf(' %.1f:%+.3f', Bs(k), cr(k))]; else, s = [s sprintf(' %.1f:x', Bs(k))]; end %#ok<AGROW>
    end
    say(LOG, sprintf('  %5.2f    %8.4f      %s', T, Bc, s));
end
say(LOG,'done');
end
function say(LOG,s), fid=fopen(LOG,'a'); fprintf(fid,'%s\n',s); fclose(fid); end
