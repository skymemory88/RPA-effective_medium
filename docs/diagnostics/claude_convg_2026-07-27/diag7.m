function diag7()
% Commensurate with the frozen Gate-0 fixture: T = 0.10 K, N=16, dpRng=30, bruteforce,
% Ecut=40, nH=33, transverse_mf='legacy_x', demag=0.
%  (A) head-to-head resummed vs strict on the Gate-0 ORDERED field set  -> closes
%      "Not established: that the strict scheme extends the ordered domain at 2-3 T".
%  (B) the strict candidate's OWN PM boundary                          -> closes
%      "Not established: the strict candidate's own PM phase boundary" (their open Q3).
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root); addpath(fullfile(root,'invz_projected')); addpath(fullfile(root,'invz_common'));
warning('off','MATLAB:nearlySingularMatrix');
SD = '/private/tmp/claude-501/-Users-yikaiyang-Library-CloudStorage-OneDrive-Nexus365-Programming-scripts-Matlab-Simulation-invZ-expansion/fe404c1e-1dd6-4ec7-868e-75af95bf66a2/scratchpad';
LOG = fullfile(SD,'diag7.log'); if exist(LOG,'file'), delete(LOG); end

ion = invz_ion();  T = 0.10;
[Jnu, info] = invz_bz_couplings(ion, struct('grid',[16 16 16],'dpRng',30,'cache',true));
Jf = Jnu(:); Jaa0 = ion.Jxx0; if isfield(info,'Jaa0'), Jaa0 = info.Jaa0; end
say(LOG, sprintf('T = %.2f K (Gate-0 fixture temperature).  J0eff=%.10g  n=%d', T, info.Jcc0, numel(Jf)));

% ---------- (A) head-to-head on the Gate-0 ordered field set + PM controls ---------------
Bs = [0.25 0.5 1 2 2.5 2.9 3.0 3.1 3.5];
say(LOG,'');
say(LOG,'--- (A) ordered (jensen) leg, HEAD-TO-HEAD, nH = 33 ---');
say(LOG,'   B    resummed                             strict_1z_dyson_ref');
for B = Bs
    r1 = one(ion,T,B,Jf,info,Jaa0,'resummed');
    r2 = one(ion,T,B,Jf,info,Jaa0,'strict_1z_dyson_ref');
    say(LOG, sprintf('%5.2f   %-36s %-36s', B, r1, r2));
end

% ---------- (B) the strict candidate's own PM boundary -----------------------------------
say(LOG,'');
say(LOG,'--- (B) PM leg crit(B): where does each scheme put its own boundary at T = 0.10 K? ---');
say(LOG,'   B    resummed: conv  crit          iters | strict: conv  crit          iters  medium');
for B = [2.5 3.0 3.1 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0]
    s1 = pmone(ion,T,B,Jf,info,Jaa0,'resummed');
    s2 = pmone(ion,T,B,Jf,info,Jaa0,'strict_1z_dyson_ref');
    say(LOG, sprintf('%5.2f   %-32s | %-44s', B, s1, s2));
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
    s = sprintf('%-20s m0=%7.4g omit=%7.3g', st, pt.m0, om);
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
