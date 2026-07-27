function diag2(Bin)
% Deep trace of the jensen ordered node loop at one (T,B).
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root); addpath(fullfile(root,'invz_projected')); addpath(fullfile(root,'invz_common'));
warning('off','MATLAB:nearlySingularMatrix');
SD = '/private/tmp/claude-501/-Users-yikaiyang-Library-CloudStorage-OneDrive-Nexus365-Programming-scripts-Matlab-Simulation-invZ-expansion/fe404c1e-1dd6-4ec7-868e-75af95bf66a2/scratchpad';
if nargin<1, Bin = 1.0; end
LOG = fullfile(SD, sprintf('diag2_B%g.log', Bin));
if exist(LOG,'file'), delete(LOG); end

ion = invz_ion();
[Jnu, info] = invz_bz_couplings(ion, struct('grid',[16 16 16],'dpRng',30,'cache',true));
Jf = Jnu(:);
Jaa0 = ion.Jxx0; if isfield(info,'Jaa0'), Jaa0 = info.Jaa0; end
T = 0.31;  B = [Bin 0 0];
say(LOG, sprintf('T=%.3g K  B=%.3g T   J0eff=%.6g  max(Jq)=%.6g  meanJ=%.6g  mu2=%.4g', ...
    T, Bin, info.Jcc0, max(Jf), mean(Jf), mean((Jf-mean(Jf)).^2)));

o = struct('J0eff', info.Jcc0, 'Jxx0', Jaa0, 'trace', true);
[hstar, prof, trc] = invz_hmf_ordered(ion, T, B, Jf, o);
say(LOG, sprintf('hmf_star=%.6g   status=%s   slope0=%.6g   n_extend=%d  redens=%d  nnodes(profile)=%d', ...
    hstar, prof.status, prof.slope0, prof.n_extend, prof.redensified, numel(prof.hgrid)));
save(fullfile(SD, sprintf('diag2_B%g.mat',Bin)), 'prof','trc','hstar','info','-v7.3');

% ---- profile node table
say(LOG,'');
say(LOG,'idx      h(meV)    acc  term_reason        r          m       Sigma0        K0        D_uni      Dq_min    gstat_den   Delta      crit');
for k = 1:numel(prof.hgrid)
    say(LOG, sprintf('%3d %12.5g %4d  %-18s %9.4g %9.4g %11.4g %11.4g %11.4g %11.4g %11.4g %9.4g %10.4g', ...
        k, prof.hgrid(k), prof.node_conv(k), prof.node_term_reason{k}, prof.r(k), prof.m(k), ...
        prof.Sigma0(k), prof.K0(k), prof.D_uni(k), prof.Dq_min(k), prof.gstat_local_denom(k), ...
        prof.Delta(k), prof.crit(k)));
end

% ---- trace node ledger
say(LOG,'');
say(LOG,'--- trc.nodes ledger (id, phase, h, outer_iters, hit_max, dS_break, resid_static, ok) ---');
for k = 1:numel(trc.nodes)
    n = trc.nodes(k);
    say(LOG, sprintf('%3d %-10s h=%12.5g  outer=%3d hitmax=%d dSbreak=%d resid_static=%10.3g ok=%d term=%s  K0=%.5g D_uni=%.5g Dqmin=%.5g', ...
        n.id, n.phase, n.h, n.outer_iters, n.outer_hit_max, n.dS_break, n.resid_static, n.ok_final, n.term_reason, n.K0, n.D_uni, n.Dq_min));
end

% ---- per-iteration trace of the FIRST failing node and of node 1
fail_ids = [];
for k = 1:numel(trc.nodes), if ~trc.nodes(k).ok_final, fail_ids(end+1) = trc.nodes(k).id; end, end %#ok<AGROW>
say(LOG,'');
say(LOG, sprintf('failing node ids: %s  (of %d)', mat2str(fail_ids), numel(trc.nodes)));
ids = unique([1, fail_ids(1:min(2,numel(fail_ids)))]);
for id = ids
    say(LOG,'');
    say(LOG, sprintf('--- outer-iteration trace, node %d (h=%.5g, phase=%s) ---', id, trc.nodes(id).h, trc.nodes(id).phase));
    sel = [trc.iters.node_id] == id;
    it = trc.iters(sel);
    say(LOG,'  out   resid_map   resid_static        K0        D_uni      Dq_min     Dq_max  Dqneg  cflag');
    for k = 1:numel(it)
        say(LOG, sprintf('%5d %11.4g %14.4g %11.5g %11.4g %11.4g %11.4g %6d %5d', ...
            it(k).outer, it(k).resid_map, it(k).resid_static, it(k).K0, it(k).D_uni, ...
            it(k).Dq_min, it(k).Dq_max, it(k).Dq_neg_count, it(k).converged_flag));
    end
end
say(LOG,'done');
end

function say(LOG,s)
fid = fopen(LOG,'a'); fprintf(fid,'%s\n',s); fclose(fid);
end
