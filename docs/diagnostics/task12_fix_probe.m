% Fix-round probe for F1 (domain + bare-shortcut schema tests) and F11 (K0 investigation).
% Read-only: no committed file touched. Reuses the EXACT fixture() the diagnostics test uses.
ROOT = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(ROOT); addpath(fullfile(ROOT,'invz_projected')); addpath(fullfile(ROOT,'invz_common'));

ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
           'static_medium', 'strict_1z_dyson_ref');

ALIGN = {'crit', 'r_minus_1', 'Delta', 'Dq_min', 'ref_denom', 'ref_margin', ...
         'gstat_local_denom', 'omit_mu3', 'omit_cubic', 'omit_max'};

fprintf('=== DOMAIN: fixture() + ref_margin=1e9 ===\n');
od = o;  od.ref_margin = 1e9;
[hd, pd] = invz_hmf_ordered(ion, T, Bx, Jnu, od);
fprintf('status=%s hstar=%s n=%d\n', pd.status, mat2str(hd), numel(pd.hgrid));
bad = 0;
for k = 1:numel(ALIGN)
    nk = numel(pd.(ALIGN{k}));
    if nk ~= numel(pd.hgrid), fprintf('  ALIGN FAIL %s: %d vs %d\n', ALIGN{k}, nk, numel(pd.hgrid)); bad=bad+1; end
end
fprintf('align failures=%d (of %d fields) ; medium_status iscell=%d numel=%d ; node_term_reason numel=%d\n', ...
    bad, numel(ALIGN), iscell(pd.medium_status), numel(pd.medium_status), numel(pd.node_term_reason));
fprintf('medium_status set = {%s}\n', strjoin(unique(pd.medium_status,'stable'), ','));
fprintf('node_term_reason set = {%s}\n', strjoin(unique(pd.node_term_reason,'stable'), ','));
fprintf('ref_denom(1)=%.17g  ref_margin(1)=%.17g  ref_margin(1)==ref_denom(1)-1e9: %d\n', ...
    pd.ref_denom(1), pd.ref_margin(1), abs(pd.ref_margin(1)-(pd.ref_denom(1)-1e9))<1e-6);
fprintf('all ref_margin == ref_denom - 1e9: %d\n', all(abs(pd.ref_margin - (pd.ref_denom-1e9)) < 1e-6));
fprintf('crit allNaN=%d  Delta finite=%d/%d  r allNaN=%d\n', ...
    all(isnan(pd.crit)), nnz(isfinite(pd.Delta)), numel(pd.Delta), all(isnan(pd.r)));
fprintf('node_conv (accepted) any true=%d\n', any(pd.node_conv));

fprintf('\n=== BARE-SHORTCUT: fixture() + force_bare=true ===\n');
ob = o;  ob.force_bare = true;
[hb, pb] = invz_hmf_ordered(ion, T, Bx, Jnu, ob);
fprintf('status=%s hstar=%s n=%d\n', pb.status, mat2str(hb), numel(pb.hgrid));
bad = 0;
for k = 1:numel(ALIGN)
    nk = numel(pb.(ALIGN{k}));
    if nk ~= numel(pb.hgrid), fprintf('  ALIGN FAIL %s: %d vs %d\n', ALIGN{k}, nk, numel(pb.hgrid)); bad=bad+1; end
end
fprintf('align failures=%d (of %d fields) ; medium_status iscell=%d numel=%d ; node_term_reason numel=%d\n', ...
    bad, numel(ALIGN), iscell(pb.medium_status), numel(pb.medium_status), numel(pb.node_term_reason));
fprintf('medium_status set = {%s}\n', strjoin(unique(pb.medium_status,'stable'), ','));
fprintf('node_term_reason set = {%s}\n', strjoin(unique(pb.node_term_reason,'stable'), ','));
fprintf('crit allNaN=%d  r all==1: %d  Sigma0 all==0: %d  K0 all==0: %d  Delta finite=%d/%d\n', ...
    all(isnan(pb.crit)), all(pb.r==1), all(pb.Sigma0==0), all(pb.K0==0), ...
    nnz(isfinite(pb.Delta)), numel(pb.Delta));

fprintf('\n=== F11: K0 provenance on a TRACED force_bare run ===\n');
obt = ob;  obt.trace = true;
[hbt, pbt, tbt] = invz_hmf_ordered(ion, T, Bx, Jnu, obt);
fprintf('trc.nodes numel=%d ; K0 set = %s ; term_reason set = {%s}\n', numel(tbt.nodes), ...
    mat2str(unique([tbt.nodes.K0])), strjoin(unique({tbt.nodes.term_reason},'stable'), ','));
fprintf('prof.K0 set on this path = %s\n', mat2str(unique(pbt.K0)));
fprintf('trc.schema_version = %d\n', tbt.schema_version);

fprintf('\nPROBE DONE\n');
