function out = invzt_a3d_convergence(opts)
%INVZT_A3D_CONVERGENCE  7C Step-3b VERTEX-BASIS CONVERGENCE study (required before
%   production): compact cc;cc dominant Vcc at Nv = 16 vs 24 (and 32 if the budget guard
%   admits it), on a REDUCED frequency grid at >= 1 ordered anchor.
%
%   WHY (2026-07-20 refinement). The rank-16 field-adapted basis captures ~0.98 of the
%   static cc response but only 0.665-0.838 of the equal-time Jz variance (measured, 7B):
%   high-energy states enter the four-point vertex as virtual INTERMEDIATES with no static
%   denominator suppression, so static capture is NOT vertex convergence. The actual check
%   is Vcc(16) vs Vcc(24): agreement at the few-percent scale supports production at rank
%   16; a large spread is a STOP (sanctioned fallbacks: (1) enlarge the field-adapted
%   cluster through the next robust energy gap, or (2) a response-augmented/block-Krylov
%   basis seeded by the low-energy cluster + the omitted Q*Jz*P directions -- NEVER the
%   zero-field projector, a moving Esplit cut, or an untracked eigenstate ranking).
%
%   For each anchor and each Nv this builds the fixed-rank field-adapted vertex basis
%   (INVZT_ORDERED_VERTEX_BASIS), passes it as opts.dom_basis to INVZT_SIGMA_TENSOR
%   (dress = 'dominant'), and records the converged cc self-energy Vcc(n) = squeeze(
%   Vmat(3,3,:)). The state cap max_vertex_states is raised to admit Nv > 16; the
%   max_vertex_work / max_vertex_bytes guards are LEFT AT DEFAULT so the budget guard
%   itself decides whether Nv = 32 fits (a refusal is caught and reported, not fatal).
%
%   DEV-TIME driver (never run by a test suite): returns a data struct and prints the
%   comparison for the controller to log to the ODD-LOG. No file writes. Runtime is
%   dominated by the one-time vertex build per (anchor, Nv) ~ 6*Nv^4*nwn*nl kernel evals;
%   reduce Ecut to speed it up (the convergence comparison is grid-shape agnostic).
%
%   opts (all optional):
%     .anchors   {[0.1 3.0], [1.0 3.0]}  {[T_K Bx_T], ...} ordered anchors
%     .Nv_list   [16 24 32]  vertex ranks to compare (32 gated by the budget guard)
%     .Ecut      10          reduced Matsubara cutoff (meV)
%     .odd       true        false -> invzt_odd_mask(lat.Jt) (drop c<->a,b mediation)
%     .mix_outer/.tol_outer/.max_outer   forwarded to invzt_sigma_tensor
%     .rel_floor 1e-3        relative-difference denominator floor (fraction of max|Vcc|)
%     .grid_label '16^3 halfopen, dpRng 30'
%
%   ACCEPTANCE (7E row, reported here, judged by the controller): |Vcc(16)-Vcc(24)| at the
%   few-percent scale supports rank 16; a large spread STOPs. This script REPORTS, it does
%   not gate.
%
%   See also INVZT_SIGMA_TENSOR, INVZT_ORDERED_VERTEX_BASIS, INVZT_A3D_BENCHMARK.
if nargin < 1 || isempty(opts), opts = struct(); end
gf = @(f, d) def(opts, f, d);
anchors    = gf('anchors',   {[0.1 3.0], [1.0 3.0]});
Nv_list    = gf('Nv_list',   [16 24 32]);
Ecut       = gf('Ecut',      10);
odd        = gf('odd',       true);
rel_floor  = gf('rel_floor', 1e-3);
grid_label = gf('grid_label','16^3 halfopen, dpRng 30');

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));                 % invz_tensor
addpath(fullfile(here, '..', '..', '..'));           % repo root
addpath(fullfile(here, '..', '..', '..', 'invz_common'));

ion = invz_ion();
g   = invzt_qgrid(16, 'halfopen');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 30, 'cache', true));
J0z = lat.info.Jcc0;  Jxx0 = lat.info.Jaa0;
lat_eff = lat;
if ~odd, lat_eff.Jt = invzt_odd_mask(lat.Jt); end

st_opts_base = struct('dress', 'dominant', 'max_vertex_states', max(Nv_list) + 8);
for f = {'mix_outer', 'tol_outer', 'max_outer'}
    if isfield(opts, f{1}), st_opts_base.(f{1}) = opts.(f{1}); end
end

fprintf('\n=== A3d VERTEX-BASIS CONVERGENCE STUDY (7C Step 3b) ===\n');
fprintf('lattice %s ; Ecut = %g meV ; odd = %d ; Nv = %s\n', ...
    grid_label, Ecut, odd, mat2str(Nv_list));

na = numel(anchors);
res = struct('anchor', {}, 'T', {}, 'Bx', {}, 'm0', {}, 'nwn', {}, 'nl', {}, ...
    'Nv', {}, 'ran', {}, 'converged', {}, 'chi_share', {}, 'var_share', {}, ...
    'p_mass', {}, 'gap_16_17', {}, 'Vcc', {}, 'cmp', {});
for ia = 1:na
    a = anchors{ia};  T = a(1);  Bx = a(2);
    siopts  = struct('hyp', true, 'order', true, 'J0z', J0z, 'Jxx0', Jxx0, 'transverse_mf', 'legacy_x');
    si_full = invz_single_ion(ion, T, [Bx 0 0], siopts);
    m0 = si_full.Jexp(3);
    [wn, ~, beta] = invz_matsubara(T, Ecut);
    nwn = numel(wn);  nl = 2*(nwn - 1) + 1;
    fprintf('\n-- anchor %d: T = %.3g K, Bx = %.3g T | m0 = %.4f | nwn = %d, nl = %d\n', ...
        ia, T, Bx, m0, nwn, nl);

    Vcc = cell(1, numel(Nv_list));
    ran = false(1, numel(Nv_list));  conv = false(1, numel(Nv_list));
    fprintf('  %4s  %5s  %10s  %10s  %8s  %11s  %8s\n', ...
        'Nv', 'ran', 'chi_share', 'var_share', 'p_mass', 'gap16_17', 'converged');
    for iv = 1:numel(Nv_list)
        Nv = Nv_list(iv);
        vb = invzt_ordered_vertex_basis(ion, T, si_full, struct('n_vertex', Nv));
        db = struct('E', vb.E, 'p', vb.p, 'Mz', vb.Mz);
        st_opts = st_opts_base;  st_opts.dom_basis = db;
        try
            stt = invzt_sigma_tensor(si_full, T, lat_eff, wn, beta, st_opts);
            Vcc{iv} = squeeze(stt.Vmat(3, 3, :));
            ran(iv) = true;  conv(iv) = stt.converged;
        catch ME
            if strcmp(ME.identifier, 'invzt:orderedA3Budget')
                ran(iv) = false;
                fprintf('  %4d  %5s  (budget-refused: %s)\n', Nv, 'no', ME.message);
                idx = numel(res) + 1;
                res(idx) = pack(ia, T, Bx, m0, nwn, nl, Nv, false, false, vb, [], struct());
                continue;
            else
                rethrow(ME);
            end
        end
        fprintf('  %4d  %5s  %10.5f  %10.5f  %8.5f  %11.4g  %8d\n', ...
            Nv, 'yes', vb.chi_share, vb.var_share, vb.p_mass, vb.gap_16_17, conv(iv));
        res(numel(res)+1) = pack(ia, T, Bx, m0, nwn, nl, Nv, true, conv(iv), vb, Vcc{iv}, struct());
    end

    % --- pairwise convergence vs the Nv=16 reference ---------------------------------
    i16 = find(Nv_list == 16, 1);
    if ~isempty(i16) && ran(i16)
        Vref = Vcc{i16};
        peak = max(abs(Vref));
        fprintf('  convergence vs Nv=16 (static |dVcc(0)|/|Vcc(0)| ; max-over-freq |dVcc|/max|Vcc16|):\n');
        for iv = 1:numel(Nv_list)
            if iv == i16 || ~ran(iv), continue; end
            dv = Vcc{iv} - Vref;
            static_rel  = abs(dv(1)) / max(abs(Vref(1)), rel_floor*peak);
            maxfreq_rel = max(abs(dv)) / max(peak, realmin);
            fprintf('    Nv=%2d vs 16:  static_rel = %.3e   max-freq_rel = %.3e\n', ...
                Nv_list(iv), static_rel, maxfreq_rel);
            % attach to the matching res row
            for k = 1:numel(res)
                if res(k).anchor == ia && res(k).Nv == Nv_list(iv)
                    res(k).cmp = struct('static_rel', static_rel, 'maxfreq_rel', maxfreq_rel);
                end
            end
        end
    else
        fprintf('  (Nv=16 reference did not run -- cannot form the convergence comparison)\n');
    end
end

fprintf(['\nREPORT (controller judges 7E): Vcc(16)-vs-Vcc(24) at the FEW-PERCENT scale ' ...
    'supports production at rank 16; a large spread STOPs (enlarge the field-adapted ' ...
    'cluster through the next gap, or a response-augmented/block-Krylov basis).\n']);

out = struct();
out.res = res;
out.params = struct('anchors', {anchors}, 'Nv_list', Nv_list, 'Ecut', Ecut, 'odd', odd, ...
    'rel_floor', rel_floor, 'grid_label', grid_label);
end

% ---------------------------------------------------------------------- %
function r = pack(anchor, T, Bx, m0, nwn, nl, Nv, ran, converged, vb, Vcc, cmp)
r = struct('anchor', anchor, 'T', T, 'Bx', Bx, 'm0', m0, 'nwn', nwn, 'nl', nl, ...
    'Nv', Nv, 'ran', ran, 'converged', converged, 'chi_share', vb.chi_share, ...
    'var_share', vb.var_share, 'p_mass', vb.p_mass, 'gap_16_17', vb.gap_16_17, ...
    'Vcc', Vcc, 'cmp', cmp);
end

function v = def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
