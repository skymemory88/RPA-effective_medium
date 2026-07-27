ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT); addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));
ion = invz_ion();
prov = struct('grid',[16 16 16],'dpRng',30,'dipole','bruteforce','cache',false);
[Jnu, info] = invz_bz_couplings(ion, prov);
mom = invz_coupling_moments(Jnu(:));
for scheme = {'strict_1z_dyson_ref','strict_1z_bare_ref','resummed'}
    s = scheme{1};
    o = struct('J0eff', info.Jcc0, 'Jxx0', ion.Jxx0, 'hyp', true, ...
               'static_medium', s, 'Jmom', mom);
    if strcmp(s,'resummed'), o = rmfield(o,'Jmom'); end
    fprintf('\n########## static_medium = %s\n', s);
    try
        [hstar, p] = invz_hmf_ordered(ion, 0.1, [1 0 0], Jnu(:), o);
        fprintf('  status      = %s\n', p.status);
        fprintf('  hstar       = %.17g   isfinite=%d  >0=%d\n', hstar, isfinite(hstar), hstar>0);
        fprintf('  crit_star   = %.17g\n', p.crit_star);
        fprintf('  slope0      = %.17g\n', p.slope0);
        fprintf('  m_star      = %.17g\n', p.m_star);
        u = unique(p.medium_status);
        for k = 1:numel(u)
            fprintf('  medium_status %-24s : %d\n', u{k}, nnz(strcmp(p.medium_status,u{k})));
        end
        if isfield(p,'omit_max') && ~isempty(p.omit_max)
            om = p.omit_max(isfinite(p.omit_max));
            if ~isempty(om)
                fprintf('  omit_max    = [%.6g , %.6g]   (omit_promote = 0.10)\n', min(om), max(om));
            end
        end
        if isfield(p,'int_Sigma0')
            fprintf('  int_Sigma0=%.6g  int_r_minus_1=%.6g  ratio=%.6g\n', ...
                p.int_Sigma0, p.int_r_minus_1, p.int_r_minus_1/p.int_Sigma0);
        end
    catch e
        fprintf('  THREW id=%s msg=%s\n', e.identifier, e.message);
    end
end
fprintf('\n### done ###\n');
