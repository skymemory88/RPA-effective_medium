function diag1()
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root); addpath(fullfile(root,'invz_projected')); addpath(fullfile(root,'invz_common'));
warning('off','MATLAB:nearlySingularMatrix');
LOG = '/private/tmp/claude-501/-Users-yikaiyang-Library-CloudStorage-OneDrive-Nexus365-Programming-scripts-Matlab-Simulation-invZ-expansion/fe404c1e-1dd6-4ec7-868e-75af95bf66a2/scratchpad/diag1.log';
if exist(LOG,'file'), delete(LOG); end
say(LOG, 'start');

t0 = tic;
ion = invz_ion();
[Jnu, info] = invz_bz_couplings(ion, struct('grid',[16 16 16],'dpRng',30,'cache',true));
say(LOG, sprintf('coupling build: %.1f s ; nq=%d ; Jcc0=%.6g meV ; Jaa0=%.6g', toc(t0), size(Jnu,1), info.Jcc0, info.Jaa0));
Jf = Jnu(:);
say(LOG, sprintf('Jnu(flat): n=%d  min=%.6g  max=%.6g  mean=%.6g  std=%.6g', numel(Jf), min(Jf), max(Jf), mean(Jf), std(Jf,1)));
save(fullfile(fileparts(LOG),'jnu.mat'),'Jnu','info','-v7.3');

T = 0.31;
Bs = [1.0 2.0 3.0 3.6 4.0 4.2 4.4];
Jaa0 = ion.Jxx0; if isfield(info,'Jaa0'), Jaa0 = info.Jaa0; end

say(LOG, '--- jensen ordered leg (default static_medium = resummed) ---');
for B = Bs
    o = struct('J0eff', info.Jcc0, 'Jxx0', Jaa0, 'ordered_mode','jensen');
    tb = tic;
    try
        pt = invz_solve_point_ordered(ion, T, [B 0 0], Jf, o);
        st = 'n/a';
        if isfield(pt,'hmf_prof'), st = pt.hmf_prof.status; end
        if isfield(pt,'hmf_status'), st = pt.hmf_status; end
        say(LOG, sprintf(['B=%4.2f  ord=%d conv=%d  status=%-20s m0=%9.4g  D_uni=%9.3g  crit=%9.3g' ...
                 '  stable=%g  final_resid=%9.3g  (%.1f s)'], ...
            B, pt.is_ordered, pt.converged, st, pt.m0, gf(pt,'D_uni'), pt.crit, ...
            gf(pt,'stable_1z'), gf(pt,'final_resid'), toc(tb)));
        if isfield(pt,'hmf_prof') && ~isempty(pt.hmf_prof.node_term_reason)
            tr = pt.hmf_prof.node_term_reason;
            u = unique(tr);
            cnt = cellfun(@(s) sum(strcmp(tr,s)), u);
            say(LOG, ['    node_term_reason: ' strjoin(arrayfun(@(k) sprintf('%s=%d',u{k},cnt(k)), 1:numel(u),'uni',0), ', ')]);
            say(LOG, sprintf('    nodes=%d accepted=%d  slope0=%.4g  crit_star=%.4g  Dq_min_star=%.4g  n_extend=%d redens=%d', ...
                numel(tr), sum(pt.hmf_prof.node_conv), pt.hmf_prof.slope0, pt.hmf_prof.crit_star, pt.hmf_prof.Dq_min_star, pt.hmf_prof.n_extend, pt.hmf_prof.redensified));
            say(LOG, sprintf('    Dq_min over nodes: min=%.4g  #(<=0)=%d/%d', min(pt.hmf_prof.Dq_min), sum(pt.hmf_prof.Dq_min<=0), numel(pt.hmf_prof.Dq_min)));
        end
    catch ME
        say(LOG, sprintf('B=%4.2f  THREW %s : %s', B, ME.identifier, ME.message));
    end
end

say(LOG, '--- strict PM leg (invz_solve_point) ---');
for B = Bs
    o = struct('J0eff', info.Jcc0, 'Jxx0', Jaa0);
    try
        p = invz_solve_point(ion, T, [B 0 0], Jf, o);
        say(LOG, sprintf('B=%4.2f  conv=%d  Sigma0=%9.4g  crit=%9.4g  iters=%d', B, p.converged, p.Sigma0, p.crit, gf(p,'outer_iters')));
    catch ME
        say(LOG, sprintf('B=%4.2f  THREW %s : %s', B, ME.identifier, ME.message));
    end
end
say(LOG,'done');
end

function say(LOG, s)
fid = fopen(LOG,'a'); fprintf(fid,'%s\n',s); fclose(fid);
end

function v = gf(s,f)
if isfield(s,f), v = s.(f); if islogical(v), v = double(v); end
    if isempty(v), v = NaN; end
else, v = NaN; end
end
