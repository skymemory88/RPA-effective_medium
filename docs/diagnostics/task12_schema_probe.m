root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root); addpath(fullfile(root,'invz_projected')); addpath(fullfile(root,'invz_common'));
ion = invz_ion(); T = 0.31; Bx = [2.85 0 0]; Jnu = linspace(-2e-3,6.0e-3,24).';
J0eff = 6.42e-3;
recF = {'r','m','Sigma0','K0','D_uni','Dq_min','G0bare','Gstat','accepted','crit','Delta', ...
        'medium_status','term_reason','ref_denom','ref_margin','gstat_local_denom', ...
        'omit_mu3','omit_cubic','omit_max'};
oldTrcF = {'id','h','phase','seed_kind','seed_from','outer_iters','outer_hit_max', ...
           'dS_break','ok_final','term_reason','K0','D_uni','resid_static'};

% --- (a) BARE-SHORTCUT record, traced: append_trace_node must emit the FULL schema ----------
ob = struct('J0eff',J0eff,'Jxx0',ion.Jxx0,'hyp',true,'force_bare',true,'nH',5,'trace',true);
[hb, pb, tb] = invz_hmf_ordered(ion, T, Bx, Jnu, ob);
fn = fieldnames(tb.nodes);
fprintf('PROBE(a) fbare+trace: nNodes=%d  missing_record_fields=%d  missing_old_fields=%d\n', ...
    numel(tb.nodes), nnz(~ismember(recF, fn)), nnz(~ismember(oldTrcF, fn)));
fprintf('PROBE(a)   term_reason set = {%s}  ok_final all true = %d\n', ...
    strjoin(unique({tb.nodes.term_reason},'stable'),','), all([tb.nodes.ok_final]));
fprintf('PROBE(a)   prof.node_term_reason = {%s}  medium_status={%s}\n', ...
    strjoin(unique(pb.node_term_reason,'stable'),','), strjoin(unique(pb.medium_status,'stable'),','));
fprintf('PROBE(a)   prof arrays aligned with hgrid(%d): ', numel(pb.hgrid));
alnF = {'crit','r_minus_1','Delta','Dq_min','ref_denom','ref_margin','gstat_local_denom', ...
        'omit_mu3','omit_cubic','omit_max','medium_status','node_term_reason'};
ok = true; for k=1:numel(alnF), ok = ok && numel(pb.(alnF{k}))==numel(pb.hgrid); end
fprintf('%d\n', ok);
fprintf('PROBE(a)   fbare rec: r=%g Sigma0=%g K0=%g crit=%g Delta=%g (G0bare NaN by design)\n', ...
    pb.r(1), pb.Sigma0(1), pb.K0(1), pb.crit(1), pb.Delta(1));

% --- (b) STRICT-mode measured diagnostics (the diagnostics fixture) ------------------------
os = struct('J0eff',J0eff,'Jxx0',ion.Jxx0,'hyp',true,'static_medium','strict_1z_dyson_ref');
[hs, ps] = invz_hmf_ordered(ion, T, Bx, Jnu, os);
fprintf('\nPROBE(b) STRICT status=%s hstar=%.17g scheme=%s nNodes=%d n_extend=%d redens=%d\n', ...
    ps.status, hs, ps.static_medium, numel(ps.hgrid), ps.n_extend, ps.redensified);
fprintf('PROBE(b) slope0=%.17g  r_pm0=%.17g  G0bare_pm0=%.17g\n', ps.slope0, ps.r_pm0, ps.G0bare_pm0);
fprintf('PROBE(b) crit: min=%.6g max=%.6g nNeg=%d nNaN=%d\n', min(ps.crit), max(ps.crit), ...
    nnz(ps.crit<0), nnz(isnan(ps.crit)));
fprintf('PROBE(b) crit_star=%.17g  r_star=%.17g  G0bare_star=%.17g  Dq_min_star=%.6g\n', ...
    ps.crit_star, ps.r_star, ps.G0bare_star, ps.Dq_min_star);
k = find(abs(ps.m)>1e-3 & isfinite(ps.Sigma0) & isfinite(ps.r), 1, 'last');
fprintf('PROBE(b) finite-m node k=%d h=%.6g m=%.6g: r-1=%.17g Sigma0=%.17g |diff|=%.6g\n', ...
    k, ps.hgrid(k), ps.m(k), ps.r_minus_1(k), ps.Sigma0(k), abs(ps.r_minus_1(k)-ps.Sigma0(k)));
fprintf('PROBE(b) int_Sigma0=%.17g  int_r_minus_1=%.17g  ratio=%.6g\n', ...
    ps.int_Sigma0, ps.int_r_minus_1, ps.int_r_minus_1/ps.int_Sigma0);
fprintf('PROBE(b) medium_status set = {%s}\n', strjoin(unique(ps.medium_status,'stable'),','));
fprintf('PROBE(b) node_term_reason set = {%s}\n', strjoin(unique(ps.node_term_reason,'stable'),','));
fprintf('PROBE(b) ref_denom: min=%.6g max=%.6g | ref_margin: min=%.6g max=%.6g\n', ...
    min(ps.ref_denom), max(ps.ref_denom), min(ps.ref_margin), max(ps.ref_margin));
fprintf('PROBE(b) ref_margin == ref_denom - 1e-6 everywhere: %d\n', ...
    all(abs(ps.ref_margin - (ps.ref_denom - 1e-6)) < 1e-18));
fprintf('PROBE(b) omit_mu3 max=%.6g  omit_cubic max=%.6g  omit_max max=%.6g\n', ...
    max(ps.omit_mu3), max(ps.omit_cubic), max(ps.omit_max));
fprintf('PROBE(b) gstat_local_denom: min=%.6g max=%.6g | Dq_min: min=%.6g\n', ...
    min(ps.gstat_local_denom), max(ps.gstat_local_denom), min(ps.Dq_min));

% --- (c) RESUMMED reference values (notes SS3 comparison) -----------------------------------
orr = struct('J0eff',J0eff,'Jxx0',ion.Jxx0,'hyp',true);
[hr, pr] = invz_hmf_ordered(ion, T, Bx, Jnu, orr);
fprintf('\nPROBE(c) RESUMMED scheme=%s hstar=%.17g crit_star=%.17g\n', pr.static_medium, hr, pr.crit_star);
kr = find(abs(pr.m)>1e-3 & isfinite(pr.Sigma0) & isfinite(pr.r), 1, 'last');
fprintf('PROBE(c) crit min=%.6g max=%.6g nNeg=%d | finite-m |r-1 - Sigma0|=%.6g\n', ...
    min(pr.crit), max(pr.crit), nnz(pr.crit<0), abs(pr.r_minus_1(kr)-pr.Sigma0(kr)));
fprintf('PROBE(c) int_Sigma0=%.17g int_r_minus_1=%.17g\n', pr.int_Sigma0, pr.int_r_minus_1);
fprintf('PROBE(c) medium_status set = {%s} ref_denom all NaN=%d ref_margin all NaN=%d\n', ...
    strjoin(unique(pr.medium_status,'stable'),','), all(isnan(pr.ref_denom)), all(isnan(pr.ref_margin)));

% --- (d) no_bare_order + node_failed exits: arrays still aligned ---------------------------
[hz, pz] = invz_hmf_ordered(ion, 0.31, [10 0 0], Jnu, orr);       % far PM: bare does not order
fprintf('\nPROBE(d) no_bare_order: status=%s hgrid=%d crit=%d medium_status=%d term=%d cell=%d sm=%s\n', ...
    pz.status, numel(pz.hgrid), numel(pz.crit), numel(pz.medium_status), ...
    numel(pz.node_term_reason), iscell(pz.medium_status), pz.static_medium);
of = orr; of.max_outer = 1; of.nH = 5;
[hf, pf] = invz_hmf_ordered(ion, T, Bx, Jnu, of);
alignedF = all(cellfun(@(f) numel(pf.(f))==numel(pf.hgrid), alnF));
fprintf('PROBE(d) node_failed: status=%s h=%g nHgrid=%d aligned=%d term={%s}\n', ...
    pf.status, hf, numel(pf.hgrid), alignedF, strjoin(unique(pf.node_term_reason,'stable'),','));
fprintf('PROBE(d) node_failed int_Sigma0=%g int_r_minus_1=%g (NaN when ok0 false)\n', ...
    pf.int_Sigma0, pf.int_r_minus_1);
fprintf('\nPROBE DONE\n');
