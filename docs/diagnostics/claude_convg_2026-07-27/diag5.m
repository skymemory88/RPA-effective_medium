function diag5()
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root); addpath(fullfile(root,'invz_projected')); addpath(fullfile(root,'invz_common'));
warning('off','MATLAB:nearlySingularMatrix');
SD = '/private/tmp/claude-501/-Users-yikaiyang-Library-CloudStorage-OneDrive-Nexus365-Programming-scripts-Matlab-Simulation-invZ-expansion/fe404c1e-1dd6-4ec7-868e-75af95bf66a2/scratchpad';
LOG = fullfile(SD,'diag5.log'); if exist(LOG,'file'), delete(LOG); end
ion = invz_ion();  T = 0.31;

% ---------- (a) bare ordered mode + strict PM leg, on the default 16^3 grid ----------------
[Jnu, info] = invz_bz_couplings(ion, struct('grid',[16 16 16],'dpRng',30,'cache',true));
Jf = Jnu(:); Jaa0 = ion.Jxx0; if isfield(info,'Jaa0'), Jaa0 = info.Jaa0; end
Bs = [1.0 2.0 3.0 3.6 4.0];
say(LOG,'--- (a1) ordered_mode = BARE (endpoint-only, no H_MF path), static_medium = resummed ---');
for B = Bs
    o = struct('J0eff',info.Jcc0,'Jxx0',Jaa0,'ordered_mode','bare');
    try
        pt = invz_solve_point_ordered(ion,T,[B 0 0],Jf,o);
        say(LOG, sprintf('B=%4.2f ord=%d conv=%d m0=%9.4g Sigma0=%10.4g crit=%10.4g iters=%d medium=%s', ...
            B, pt.is_ordered, pt.converged, pt.m0, pt.Sigma0, pt.crit, gf(pt,'outer_iters'), pt.medium_status));
    catch ME, say(LOG,sprintf('B=%4.2f THREW %s: %s',B,ME.identifier,ME.message)); end
end
say(LOG,'--- (a2) strict PM leg under strict_1z_dyson_ref ---');
for B = Bs
    o = struct('J0eff',info.Jcc0,'Jxx0',Jaa0,'static_medium','strict_1z_dyson_ref');
    try
        p = invz_solve_point(ion,T,[B 0 0],Jf,o);
        say(LOG, sprintf('B=%4.2f conv=%d Sigma0=%10.4g crit=%10.4g iters=%3d medium=%s', ...
            B,p.converged,p.Sigma0,p.crit,gf(p,'outer_iters'),p.medium_status));
    catch ME, say(LOG,sprintf('B=%4.2f THREW %s: %s',B,ME.identifier,ME.message)); end
end
% ---------- (b) why strict fails at 1 T -------------------------------------------------
say(LOG,'--- (b) strict at B = 1 T: which domain event? ---');
o = struct('J0eff',info.Jcc0,'Jxx0',Jaa0,'ordered_mode','jensen','static_medium','strict_1z_dyson_ref');
try
    pt = invz_solve_point_ordered(ion,T,[1 0 0],Jf,o);
    say(LOG,sprintf('  status=%s  medium_status=%s  ref_denom=%.6g  margin=%.6g', ...
        gs(pt,'hmf_status'), pt.medium_status, pt.medium_denom, pt.medium_margin));
    if isfield(pt,'hmf_prof')
        pr = pt.hmf_prof;
        say(LOG,sprintf('  profile medium_status tally: %s', tally(pr.medium_status)));
        say(LOG,sprintf('  ref_denom over nodes: min=%.5g max=%.5g ; Sigma0 min=%.5g max=%.5g', ...
            min(pr.ref_denom), max(pr.ref_denom), min(pr.Sigma0), max(pr.Sigma0)));
    end
catch ME, say(LOG,sprintf('  THREW %s: %s',ME.identifier,ME.message)); end

% ---------- (c) grid-refinement test of the "working sliver" -----------------------------
say(LOG,'--- (c) grid refinement: max(Jq) -> J0eff, so the resummed sliver must shrink ---');
for N = [8 12 16 20]
    try
        [Jn, inf2] = invz_bz_couplings(ion, struct('grid',[N N N],'dpRng',30,'cache',true));
        Jv = Jn(:);
        say(LOG, sprintf('  N=%2d  nq=%6d  J0eff=%.6g  max(Jq)=%.6g  gap=J0-maxJ=%.4g', ...
            N, numel(Jv), inf2.Jcc0, max(Jv), inf2.Jcc0-max(Jv)));
        Ja = ion.Jxx0; if isfield(inf2,'Jaa0'), Ja = inf2.Jaa0; end
        for B = [3.4 3.6 3.8 4.0]
            o = struct('J0eff',inf2.Jcc0,'Jxx0',Ja,'ordered_mode','jensen');
            try
                pt = invz_solve_point_ordered(ion,T,[B 0 0],Jv,o);
                say(LOG,sprintf('     B=%4.2f  resummed jensen: ord=%d conv=%d status=%s m0=%.4g', ...
                    B, pt.is_ordered, pt.converged, gs2(pt), pt.m0));
            catch ME, say(LOG,sprintf('     B=%4.2f THREW %s',B,ME.identifier)); end
        end
    catch ME
        say(LOG, sprintf('  N=%2d coupling build failed: %s', N, ME.message));
    end
end
say(LOG,'done');
end
function say(LOG,s), fid=fopen(LOG,'a'); fprintf(fid,'%s\n',s); fclose(fid); end
function v=gf(s,f), if isfield(s,f), v=s.(f); if islogical(v),v=double(v);end; if isempty(v),v=NaN;end, else, v=NaN; end, end
function v=gs(s,f), if isfield(s,f), v=s.(f); else, v='n/a'; end, end
function v=gs2(pt)
v='n/a'; if isfield(pt,'hmf_prof'), v=pt.hmf_prof.status; end
if isfield(pt,'hmf_status'), v=pt.hmf_status; end
end
function s=tally(c)
u=unique(c); s=strjoin(arrayfun(@(k) sprintf('%s=%d',u{k},sum(strcmp(c,u{k}))),1:numel(u),'uni',0),',');
end
