% Task-12 review probe (read-only measurement; no repo files modified).
% (1) ADAPTIVE path: does the export block + both integrals describe the FINAL grid?
% (2) fbare exit: full-schema trace record, one per eval_node call?
% (3) strict DOMAIN exit: full-schema trace record + prof arrays?
ROOT = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(genpath(fullfile(ROOT,'invz_projected'))); addpath(fullfile(ROOT,'invz_common'));
addpath(ROOT);
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
base = struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true, 'transverse_mf', 'legacy_x');

REC_FIELDS = {'r','m','Sigma0','K0','D_uni','Dq_min','G0bare','Gstat','accepted', ...
    'crit','Delta','medium_status','term_reason','ref_denom','ref_margin', ...
    'gstat_local_denom','omit_mu3','omit_cubic','omit_max'};
PRE_TRACE = {'id','h','phase','seed_kind','seed_from','outer_iters','outer_hit_max', ...
    'dS_break','ok_final','term_reason','K0','D_uni','resid_static'};
ARRS = {'crit','r_minus_1','Delta','Dq_min','ref_denom','ref_margin', ...
    'gstat_local_denom','omit_mu3','omit_cubic','omit_max'};

% ---------- (1) ADAPTIVE: extension + redensify -------------------------------------------
fprintf('\n== (1) ADAPTIVE ==\n');
Btest = 3.09;
o = base;
[h_ref, p_ref] = invz_hmf_ordered(ion, T, [Btest 0 0], Jnu, o);
fprintf('ref: h=%.6g status=%s\n', h_ref, p_ref.status);
frac = min(0.5, 4*h_ref/max(p_ref.hgrid));
oa = base;  oa.hmin_frac = frac;  oa.nH = 17;
[ha, pa] = invz_hmf_ordered(ion, T, [Btest 0 0], Jnu, oa);
fprintf('ADAPT: h=%.6g status=%s n_extend=%d redensified=%d nNodes=%d nH=%d\n', ...
    ha, pa.status, pa.n_extend, pa.redensified, numel(pa.hgrid), oa.nH);
bad = 0;
for k = 1:numel(ARRS)
    if numel(pa.(ARRS{k})) ~= numel(pa.hgrid)
        fprintf('  ALIGN FAIL %s: %d vs %d\n', ARRS{k}, numel(pa.(ARRS{k})), numel(pa.hgrid)); bad = bad+1;
    end
end
if numel(pa.medium_status) ~= numel(pa.hgrid), fprintf('  ALIGN FAIL medium_status\n'); bad=bad+1; end
if numel(pa.node_term_reason) ~= numel(pa.hgrid), fprintf('  ALIGN FAIL node_term_reason\n'); bad=bad+1; end
fprintf('ADAPT align failures = %d\n', bad);
% Analytic identity: int(r-1)dh == h0(end) - hgrid(end) IFF the integral was built on the
% FINAL grid with h0's own first-panel seeding. Independent of re-running trapz.
lhs = pa.int_r_minus_1;  rhs = pa.h0(end) - pa.hgrid(end);
fprintf('int_r_minus_1 = %.17g\n  h0(end)-hgrid(end) = %.17g\n  reldiff = %.3g\n', ...
    lhs, rhs, abs(lhs-rhs)/max(abs(rhs),realmin));
fprintf('int_Sigma0 = %.17g   trapz-recheck reldiff = %.3g\n', pa.int_Sigma0, ...
    abs(pa.int_Sigma0 - trapz([0 pa.hgrid],[pa.Sigma0_pm0 pa.Sigma0]))/abs(pa.int_Sigma0));
fprintf('crit vs r+J0eff*G0bare : max reldiff = %.3g\n', ...
    max(abs(pa.crit - (pa.r + oa.J0eff*pa.G0bare))./max(abs(pa.crit),realmin)));
fprintf('crit_star = %.17g  r_star+J0eff*G0bare_star = %.17g\n', ...
    pa.crit_star, pa.r_star + oa.J0eff*pa.G0bare_star);
fprintf('slope0 = %.17g  r_pm0+J0eff*G0bare_pm0 = %.17g  bitwise=%d\n', ...
    pa.slope0, pa.r_pm0 + oa.J0eff*pa.G0bare_pm0, ...
    isequaln(pa.slope0, pa.r_pm0 + oa.J0eff*pa.G0bare_pm0));
% traced adaptive run: one trc.nodes entry per eval_node call, full schema on every entry
ot = oa;  ot.trace = true;
[hb, pb, trcA] = invz_hmf_ordered(ion, T, [Btest 0 0], Jnu, ot);
fprintf('traced ADAPT: h bitwise-equal untraced = %d ; prof isequaln = %d\n', ...
    isequaln(ha,hb), isequaln(pa,pb));
fn = fieldnames(trcA.nodes);
fprintf('trc.nodes: n=%d nFields=%d  missing REC=%s  missing PRE=%s\n', ...
    numel(trcA.nodes), numel(fn), strjoin(setdiff(REC_FIELDS, fn), ','), ...
    strjoin(setdiff(PRE_TRACE, fn), ','));
fprintf('phases = %s\n', strjoin(unique({trcA.nodes.phase}), ','));
fprintf('ids contiguous 1..n = %d\n', isequal([trcA.nodes.id], 1:numel(trcA.nodes)));
fprintf('term_reason set = %s\n', strjoin(unique({trcA.nodes.term_reason}), ','));
fprintf('any not_evaluated/not_applicable term_reason = %d\n', ...
    any(ismember({trcA.nodes.term_reason}, {'not_evaluated'})));

% ---------- (2) fbare exit ------------------------------------------------------------------
fprintf('\n== (2) FBARE ==\n');
ob = base;  ob.force_bare = true;  ob.trace = true;
[hf, pf, trcF] = invz_hmf_ordered(ion, T, [2.85 0 0], Jnu, ob);
fprintf('fbare: h=%.6g status=%s nNodes=%d trcNodes=%d\n', hf, pf.status, ...
    numel(pf.hgrid), numel(trcF.nodes));
fnf = fieldnames(trcF.nodes);
fprintf('missing REC=%s  missing PRE=%s\n', strjoin(setdiff(REC_FIELDS,fnf),','), ...
    strjoin(setdiff(PRE_TRACE,fnf),','));
fprintf('term_reason set = %s ; K0 set = %s ; D_uni allNaN = %d\n', ...
    strjoin(unique({trcF.nodes.term_reason}),','), mat2str(unique([trcF.nodes.K0])), ...
    all(isnan([trcF.nodes.D_uni])));
fprintf('prof.crit allNaN = %d ; prof.medium_status set = %s\n', all(isnan(pf.crit)), ...
    strjoin(unique(pf.medium_status),','));
badf = 0;
for k = 1:numel(ARRS)
    if numel(pf.(ARRS{k})) ~= numel(pf.hgrid), badf = badf+1; end
end
fprintf('fbare align failures = %d ; Delta finite = %d/%d\n', badf, ...
    nnz(isfinite(pf.Delta)), numel(pf.Delta));

% ---------- (3) strict DOMAIN exit ----------------------------------------------------------
fprintf('\n== (3) STRICT DOMAIN (ref_margin = 1e9) ==\n');
od = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
            'static_medium', 'strict_1z_dyson_ref', 'ref_margin', 1e9, 'trace', true);
[hd, pd, trcD] = invz_hmf_ordered(ion, T, [2.85 0 0], Jnu, od);
fprintf('domain: h=%.6g status=%s static_medium=%s nNodes=%d trcNodes=%d\n', ...
    hd, pd.status, pd.static_medium, numel(pd.hgrid), numel(trcD.nodes));
fnd = fieldnames(trcD.nodes);
fprintf('missing REC=%s  missing PRE=%s  nFields=%d\n', ...
    strjoin(setdiff(REC_FIELDS,fnd),','), strjoin(setdiff(PRE_TRACE,fnd),','), numel(fnd));
fprintf('medium_status set = %s\n', strjoin(unique(pd.medium_status),','));
fprintf('node_term_reason set = %s\n', strjoin(unique(pd.node_term_reason),','));
fprintf('trc term_reason set = %s\n', strjoin(unique({trcD.nodes.term_reason}),','));
fprintf('ref_denom(1)=%.17g  ref_margin(1)=%.17g  gstat_local_denom allNaN=%d\n', ...
    pd.ref_denom(1), pd.ref_margin(1), all(isnan(pd.gstat_local_denom)));
fprintf('Dq_min allNaN=%d  crit allNaN=%d  Delta finite=%d/%d\n', ...
    all(isnan(pd.Dq_min)), all(isnan(pd.crit)), nnz(isfinite(pd.Delta)), numel(pd.Delta));
badd = 0;
for k = 1:numel(ARRS)
    if numel(pd.(ARRS{k})) ~= numel(pd.hgrid), badd = badd+1; end
end
fprintf('domain align failures = %d ; int_Sigma0=%.6g int_r_minus_1=%.6g\n', badd, ...
    pd.int_Sigma0, pd.int_r_minus_1);
fprintf('domain nodes trace field-order == adaptive field-order : %d\n', isequal(fnd, fn));
fprintf('\nPROBE DONE\n');
