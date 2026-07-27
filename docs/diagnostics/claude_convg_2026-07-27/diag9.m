function diag9()
% The Gate-0 field set [0.05 .. 3.0] T + controls [3.1 3.5] T is the field set of a
% boundary at ~3.03 T, which diag8 locates at T ~ 1.0 K -- not at the registered
% T = 0.10 K (Bc = 4.71 T there). Re-run the SAME field set at T = 1.00 K, where it
% is the field set it was evidently designed to be, and report clause (a)/(c)/(e)
% quantities.
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root); addpath(fullfile(root,'invz_projected')); addpath(fullfile(root,'invz_common'));
warning('off','MATLAB:nearlySingularMatrix');
SD = '/private/tmp/claude-501/-Users-yikaiyang-Library-CloudStorage-OneDrive-Nexus365-Programming-scripts-Matlab-Simulation-invZ-expansion/fe404c1e-1dd6-4ec7-868e-75af95bf66a2/scratchpad';
LOG = fullfile(SD,'diag9.log'); if exist(LOG,'file'), delete(LOG); end
ion = invz_ion();  T = 1.00;
[Jnu, info] = invz_bz_couplings(ion, struct('grid',[16 16 16],'dpRng',30,'cache',true));
Jf = Jnu(:); Jaa0 = ion.Jxx0; if isfield(info,'Jaa0'), Jaa0 = info.Jaa0; end
say(LOG, sprintf('T = %.2f K  (diag8: Bc_1z(1.00 K) ~ 3.03 T -- the boundary the Gate-0 field set fits)', T));

say(LOG,'');
say(LOG,'--- Gate-0 ORDERED field set at T = 1.00 K, nH = 33 ---');
say(LOG,'   B    resummed                             strict_1z_dyson_ref');
for B = [0.05 0.25 0.5 1 2 2.5 2.9 3.0]
    r1 = one(ion,T,B,Jf,info,Jaa0,'resummed');
    r2 = one(ion,T,B,Jf,info,Jaa0,'strict_1z_dyson_ref');
    say(LOG, sprintf('%5.2f   %-36s %-36s', B, r1, r2));
end

say(LOG,'');
say(LOG,'--- Gate-0 PM CONTROLS at T = 1.00 K (prereg requires converged, finite, POSITIVE mass) ---');
for B = [3.1 3.5]
    s1 = pmone(ion,T,B,Jf,info,Jaa0,'resummed');
    s2 = pmone(ion,T,B,Jf,info,Jaa0,'strict_1z_dyson_ref');
    say(LOG, sprintf('  B=%.2f  resummed: %-46s | strict: %s', B, s1, s2));
end
say(LOG,'done');
end

function s = one(ion,T,B,Jf,info,Jaa0,scheme)
o = struct('J0eff',info.Jcc0,'Jxx0',Jaa0,'ordered_mode','jensen','static_medium',scheme);
try
    pt = invz_solve_point_ordered(ion,T,[B 0 0],Jf,o);
    st='n/a'; if isfield(pt,'hmf_prof'), st=pt.hmf_prof.status; end
    if isfield(pt,'hmf_status'), st=pt.hmf_status; end
    om = NaN; if isfield(pt,'path_omit_max'), om = pt.path_omit_max; end
    s = sprintf('%-20s m0=%7.4g omit=%7.4g', st, pt.m0, om);
catch ME
    s = sprintf('THREW %s', ME.identifier);
end
end
function s = pmone(ion,T,B,Jf,info,Jaa0,scheme)
o = struct('J0eff',info.Jcc0,'Jxx0',Jaa0,'static_medium',scheme);
try
    p = invz_solve_point(ion,T,[B 0 0],Jf,o);
    s = sprintf('conv=%d crit=%+12.6g it=%3d %s', p.converged, p.crit, p.outer_iters, p.medium_status);
catch ME
    s = sprintf('THREW %s', ME.identifier);
end
end
function say(LOG,s), fid=fopen(LOG,'a'); fprintf(fid,'%s\n',s); fclose(fid); end
