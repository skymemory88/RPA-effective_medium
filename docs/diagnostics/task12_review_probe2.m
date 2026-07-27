% Task-12 review probe 2: force the ADAPTIVE (extension + redensify) path and check that the
% export block and BOTH path integrals describe the FINAL adapted grid.
ROOT = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(genpath(fullfile(ROOT,'invz_projected'))); addpath(fullfile(ROOT,'invz_common'));
addpath(ROOT);
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
base = struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true, 'transverse_mf', 'legacy_x');
ARRS = {'crit','r_minus_1','Delta','Dq_min','ref_denom','ref_margin', ...
    'gstat_local_denom','omit_mu3','omit_cubic','omit_max'};

Bcand = [3.05 3.00 2.95];
h_ref = NaN;
for b = Bcand
    [hh, pp] = invz_hmf_ordered(ion, T, [b 0 0], Jnu, base);
    fprintf('scan B=%.3f -> h=%.6g status=%s slope0=%.4g\n', b, hh, pp.status, pp.slope0);
    if isfinite(hh) && hh > 0, h_ref = hh; p_ref = pp; Btest = b; break; end
end
if ~isfinite(h_ref), error('no below-Bc field found in the scan'); end

frac = min(0.5, 4*h_ref/max(p_ref.hgrid));
oa = base;  oa.hmin_frac = frac;  oa.nH = 17;
[ha, pa] = invz_hmf_ordered(ion, T, [Btest 0 0], Jnu, oa);
fprintf('\nADAPT B=%.3f frac=%.6g: h=%.6g status=%s n_extend=%d redensified=%d nNodes=%d\n', ...
    Btest, frac, ha, pa.status, pa.n_extend, pa.redensified, numel(pa.hgrid));
bad = 0;
for k = 1:numel(ARRS)
    if numel(pa.(ARRS{k})) ~= numel(pa.hgrid)
        fprintf('  ALIGN FAIL %s: %d vs %d\n', ARRS{k}, numel(pa.(ARRS{k})), numel(pa.hgrid)); bad = bad+1;
    end
end
if numel(pa.medium_status)    ~= numel(pa.hgrid), fprintf('  ALIGN FAIL medium_status\n');    bad=bad+1; end
if numel(pa.node_term_reason) ~= numel(pa.hgrid), fprintf('  ALIGN FAIL node_term_reason\n'); bad=bad+1; end
fprintf('align failures = %d ; numel(r)=%d numel(h0)=%d numel(F)=%d\n', bad, ...
    numel(pa.r), numel(pa.h0), numel(pa.F));
lhs = pa.int_r_minus_1;  rhs = pa.h0(end) - pa.hgrid(end);
fprintf('int_r_minus_1 = %.17g\n h0(end)-hgrid(end) = %.17g   reldiff = %.3g\n', lhs, rhs, ...
    abs(lhs-rhs)/max(abs(rhs),realmin));
fprintf('int_Sigma0 = %.17g   trapz-recheck reldiff = %.3g\n', pa.int_Sigma0, ...
    abs(pa.int_Sigma0 - trapz([0 pa.hgrid],[pa.Sigma0_pm0 pa.Sigma0]))/abs(pa.int_Sigma0));
% integrals must NOT match a superseded grid: compare against the pre-redensify hmax scale
fprintf('hgrid(1)=%.6g hgrid(end)=%.6g hmin_initial=%.6g (redensified grid differs from initial: %d)\n', ...
    pa.hgrid(1), pa.hgrid(end), pa.hmin_initial, pa.hgrid(1) ~= pa.hmin_initial);
fprintf('crit vs r+J0eff*G0bare : max reldiff = %.3g\n', ...
    max(abs(pa.crit - (pa.r + oa.J0eff*pa.G0bare))./max(abs(pa.crit),realmin)));
fprintf('crit_star = %.17g  r_star+J0eff*G0bare_star = %.17g  crit_star>0 = %d\n', ...
    pa.crit_star, pa.r_star + oa.J0eff*pa.G0bare_star, pa.crit_star > 0);
fprintf('F sign changes = %d ; r_minus_1==r-1 bitwise = %d\n', ...
    nnz(diff(sign(pa.F))~=0), isequaln(pa.r_minus_1, pa.r - 1));
ot = oa;  ot.trace = true;
[hb, pb, trcA] = invz_hmf_ordered(ion, T, [Btest 0 0], Jnu, ot);
fprintf('traced: h bitwise = %d ; prof isequaln = %d ; trcNodes=%d nFields=%d\n', ...
    isequaln(ha,hb), isequaln(pa,pb), numel(trcA.nodes), numel(fieldnames(trcA.nodes)));
fprintf('phases = %s ; ids contiguous = %d\n', strjoin(unique({trcA.nodes.phase}),','), ...
    isequal([trcA.nodes.id], 1:numel(trcA.nodes)));
nsweep = nnz(strcmp({trcA.nodes.phase},'sweep'));
next   = nnz(strcmp({trcA.nodes.phase},'extend'));
nred   = nnz(strcmp({trcA.nodes.phase},'redensify'));
nbis   = nnz(strcmp({trcA.nodes.phase},'bisect'));
nroot  = nnz(strcmp({trcA.nodes.phase},'root'));
npred  = nnz(strcmp({trcA.nodes.phase},'predictor'));
fprintf('counts: pred=%d sweep=%d extend=%d redensify=%d bisect=%d root=%d  sum=%d\n', ...
    npred, nsweep, next, nred, nbis, nroot, npred+nsweep+next+nred+nbis+nroot);
fprintf('expected extend = 3*n_extend = %d\n', 3*pa.n_extend);
fprintf('\nPROBE2 DONE\n');
