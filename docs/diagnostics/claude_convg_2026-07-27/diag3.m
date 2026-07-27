function diag3()
% Same ordered sweep under the STRICT one-shot static medium, plus a finer field scan.
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root); addpath(fullfile(root,'invz_projected')); addpath(fullfile(root,'invz_common'));
warning('off','MATLAB:nearlySingularMatrix');
SD = '/private/tmp/claude-501/-Users-yikaiyang-Library-CloudStorage-OneDrive-Nexus365-Programming-scripts-Matlab-Simulation-invZ-expansion/fe404c1e-1dd6-4ec7-868e-75af95bf66a2/scratchpad';
LOG = fullfile(SD,'diag3.log'); if exist(LOG,'file'), delete(LOG); end

ion = invz_ion();
[Jnu, info] = invz_bz_couplings(ion, struct('grid',[16 16 16],'dpRng',30,'cache',true));
Jf = Jnu(:);
Jaa0 = ion.Jxx0; if isfield(info,'Jaa0'), Jaa0 = info.Jaa0; end
T = 0.31;
Bs = [1.0 2.0 3.0 3.6 3.8 3.9 4.0 4.1];

for scheme = {'resummed','strict_1z_dyson_ref'}
    say(LOG, sprintf('===== static_medium = %s =====', scheme{1}));
    for B = Bs
        o = struct('J0eff', info.Jcc0, 'Jxx0', Jaa0, 'ordered_mode','jensen', 'static_medium', scheme{1});
        try
            pt = invz_solve_point_ordered(ion, T, [B 0 0], Jf, o);
            st = 'n/a';
            if isfield(pt,'hmf_prof'), st = pt.hmf_prof.status; end
            if isfield(pt,'hmf_status'), st = pt.hmf_status; end
            extra = '';
            if isfield(pt,'hmf_prof') && ~isempty(pt.hmf_prof.node_term_reason)
                tr = pt.hmf_prof.node_term_reason; u = unique(tr);
                extra = strjoin(arrayfun(@(k) sprintf('%s=%d',u{k},sum(strcmp(tr,u{k}))), 1:numel(u),'uni',0), ',');
            end
            say(LOG, sprintf('B=%4.2f ord=%d conv=%d status=%-20s m0=%9.4g D_uni=%10.4g stable=%g omit=%8.3g  [%s]', ...
                B, pt.is_ordered, pt.converged, st, pt.m0, gf(pt,'D_uni'), gf(pt,'stable_1z'), gf(pt,'path_omit_max'), extra));
        catch ME
            say(LOG, sprintf('B=%4.2f THREW %s : %s', B, ME.identifier, ME.message));
        end
    end
end
say(LOG,'done');
end
function say(LOG,s), fid=fopen(LOG,'a'); fprintf(fid,'%s\n',s); fclose(fid); end
function v = gf(s,f)
if isfield(s,f), v=s.(f); if islogical(v), v=double(v); end; if isempty(v), v=NaN; end, else, v=NaN; end
end
