function vb = invzt_ordered_vertex_basis(ion, T, si_full, opts)
%INVZT_ORDERED_VERTEX_BASIS  Fixed-RANK, FIELD-ADAPTED spectral subspace for the ordered a3d.
%   vb = INVZT_ORDERED_VERTEX_BASIS(ion, T, si_full, opts) builds the 16-state dominant vertex
%   basis for the ordered-phase full-response solve as the LOWEST n_vertex eigenstates of the
%   ORDERED MEAN-FIELD Hamiltonian at (T, B) -- i.e. a DIRECT TRUNCATION of the converged full
%   electronuclear single ion si_full. Fixed RANK (16 by default), field-ADAPTED CONTENT: the
%   basis tracks the field-rotated low-energy manifold instead of a frozen zero-field subspace.
%
%   FIXED-RANK FRAMING. The object is a SPECTRAL SUBSPACE (the span of the nv lowest eigenvectors),
%   not a labelled state list: it is INVARIANT under any internal rotation/relabelling of those nv
%   states (only their span enters the truncated response). The truncation is therefore
%   well-defined and continuous EXCEPT where the 16/17 gap closes -- the ONLY event that can
%   reshuffle which eigenvectors span the subspace. vb.gap_16_17 / vb.gap_ratio / vb.gap_kBT
%   monitor exactly that edge; vb.V16 lets a caller measure adjacent-field subspace overlap
%   (principal angles) to confirm smooth tracking.
%
%   WHY FIELD-ADAPTED (executed falsification, 2026-07-20). Round-2 review P0-3 recommended a
%   FIELD-INDEPENDENT content-defined basis: the zero-field electronic ground doublet tensored
%   with I8 (INVZT_RUNG_BASIS 'e2xI8', dim 16). Task 7B MEASURED that this captures ~NONE of the
%   ordered static cc response and its subspace drifts from the true low manifold under the
%   transverse field:
%
%       B (T)   overlap_deficit   chi_share (cc, elastic)
%       3.00        1.56264        2.8e-13
%       3.50        1.67335        ~0
%       4.00        1.76804        ~0
%       4.40        1.83487        1e-5
%       4.65        1.87331        1.3e-4
%
%   Root cause: the zero-field ground doublet is an Ising doublet, <+|Jz|-> = 0, and the ordering
%   field hz = J0z<Jz> fully polarizes it (reduced Jz-fluctuation -> 0). The entire longitudinal
%   (cc) fluctuation lives in VIRTUAL transitions to CF excited states, which a doublet-only basis
%   cannot represent. But 16 states ARE enough IF field-adapted: the 16 LOWEST eigenstates of the
%   ordered Hamiltonian capture chi_share = 0.977 at B = 3.0 T (measured), rising to >0.99 by
%   n = 32. So 7B fixes the RANK and adapts the CONTENT.
%
%   COUNT-DRIFT (P0-3) RESOLUTION. The objection to the instantaneous E < Esplit ENERGY cut was
%   that it selects a drifting 13/13/11/10/9/8 -state set across the field sweep -- a moving
%   DIMENSION that injects cutoff jumps. This basis fixes the dimension at n_vertex (16, constant
%   by construction) and keeps only the CONTENT field-adapted; continuity is then a measurable
%   diagnostic (chi_share smoothness + subspace overlap of vb.V16), not a dimensional discontinuity.
%
%   CONVERGENCE CAVEAT (verbatim-in-spirit, do NOT over-read chi_share). A high STATIC
%   chi_share (~0.98) does NOT prove vertex convergence: high-energy states are DENOMINATOR-
%   SUPPRESSED in the static (z = 0) susceptibility but enter the four-point vertex as virtual
%   INTERMEDIATES with no such suppression. The equal-time connected Jz variance captured by the
%   truncation (vb.var_share) is only ~0.665-0.838 across 3.0-4.65 T -- materially below the
%   static 0.98. Task 7C's Vcc(16)-vs-Vcc(24) study is the actual vertex-convergence check;
%   var_share here is the early-warning diagnostic, chi_share alone is NOT sufficient.
%
%   CONSUMER CENTERING RULE. The four-point vertex is built from the CONNECTED (mean-subtracted)
%   Jz. Consumers must center vb.Mz before use:
%       Mz_centered = vb.Mz - real(sum(vb.p .* diag(vb.Mz)))*eye(vb.n_vertex);
%   (subtract the population-weighted <Jz> of the truncated basis, not the full-ion <Jz>).
%
%   SIGNATURE DEVIATION (documented): the original Task-7B brief wrote (ion, si_full, opts);
%   vb.p needs T (see populations below), so T is threaded in as the 2nd argument. si_full is the
%   CONVERGED ordered mean-field single ion, exactly the struct INVZT_SOLVE_POINT_ORDERED builds
%   (invz_single_ion with hyp+order+J0z/Jxx0).
%
%   OPTIONS (getf): n_vertex (default 16) -- the fixed truncation rank; errors invzt:vertexCount
%   if > n_full. vb_prev (default []) -- a previous-field vb struct; when supplied the constructor
%   fills vb.min_subspace_overlap / vb.projector_distance from svd(vb_prev.V16' * vb.V16).
%
%   OUTPUT struct vb:
%     vb.E          [nv x 1]  si_full.E(1:nv), already ground-shifted (E(1)=0).
%     vb.Mx,My,Mz   [nv x nv] si_full.M{x,y,z}(1:nv,1:nv): the full 3x3 truncated response tensor
%                   (a3d's chi_dom hybrid needs all three, not just Jz).
%     vb.p          [nv x 1]  populations of the CONVERGED state restricted to the subspace,
%                   si_full.P(1:nv)/vb.p_mass (same-state provenance, not a recomputed Boltzmann).
%     vb.p_mass     scalar    sum(si_full.P(1:nv)): thermal weight the subspace carries (~1 when
%                   the 16/17 gap >> kB*T). Soft coverage gate (> 0.99).
%     vb.V16        [n_full x nv] si_full.V(:,1:nv): the vertex eigenvectors in the full space.
%     vb.n_full     numel(si_full.E) (136 for the full electronuclear ion).
%     vb.n_vertex   nv (= 16 by default).
%     vb.chi_share  static cc chi0 of the TRUNCATED system / full (elastic ON both). GATE (> 0.9).
%     vb.var_share  truncated / full equal-time CONNECTED Jz variance (the vertex-relevant weight;
%                   < chi_share by construction -- see CONVERGENCE CAVEAT).
%     vb.gap_16_17  si_full.E(nv+1) - si_full.E(nv): truncation-edge gap; NaN if nv == n_full.
%     vb.cluster_width vb.E(end) - vb.E(1): energy span of the retained subspace.
%     vb.gap_ratio  vb.gap_16_17 / max(vb.cluster_width, eps): edge isolation vs internal spread.
%     vb.gap_kBT    vb.gap_16_17 / (kB*T): edge isolation vs temperature (population leak monitor).
%     vb.min_subspace_overlap  min svd(vb_prev.V16' * vb.V16) if opts.vb_prev given, else NaN.
%     vb.projector_distance    sqrt(1 - min_subspace_overlap^2) (largest principal-angle sine),
%                   else NaN.
%   Diagnostics (not asserted): vb.Jexp, vb.chi0_cc_vertex, vb.chi0_cc_full, vb.si_vertex, vb.T.
%
%   See also INVZT_SOLVE_POINT_ORDERED, INVZ_SINGLE_ION, INVZ_CHI0Z, INVZT_SOLVE_POINT
%   (local_rung_si), INVZT_RUNG_BASIS.
if nargin < 4 || isempty(opts), opts = struct(); end

C = invz_const();

% --- guard: si_full must be a full converged single-ion struct -------------------------
for f = {'E', 'V', 'P', 'Mx', 'My', 'Mz'}
    if ~isfield(si_full, f{1})
        error('invzt:vertexSi', ['si_full must be a full INVZ_SINGLE_ION struct (missing ' ...
            'field ''%s''): the ordered vertex basis truncates the converged mean-field ' ...
            'single ion.'], f{1});
    end
end
n_full = numel(si_full.E);

nv = getf(opts, 'n_vertex', 16);
if ~(isscalar(nv) && nv == round(nv) && nv >= 1)
    error('invzt:vertexCount', 'n_vertex must be a positive integer; got %s.', mat2str(nv));
end
if nv > n_full
    error('invzt:vertexCount', ['n_vertex = %d exceeds n_full = %d: the vertex basis is a ' ...
        'TRUNCATION of si_full, it cannot exceed the full space.'], nv, n_full);
end

% --- direct truncation of the ordered mean-field single ion (fixed rank, adapted content) --
idx = 1:nv;
E   = si_full.E(idx);   E = E(:);                     % already ground-shifted
Mx  = si_full.Mx(idx, idx);
My  = si_full.My(idx, idx);
Mz  = si_full.Mz(idx, idx);
V16 = si_full.V(:, idx);
% Populations from the CONVERGED state restricted to the subspace (same-state provenance;
% mathematically identical to a fresh Boltzmann over E for the Boltzmann si_full.P, but
% guarantees consistency with the state si_full actually evaluated and yields p_mass).
p_mass = sum(si_full.P(idx));
p      = si_full.P(idx) / p_mass;   p = p(:);

% --- reduced si-like struct (INVZ_CHI0Z field surface) fed to the truncated-leg chi0 ---
si_red.E  = E;
si_red.P  = p;
si_red.Mx = Mx;  si_red.My = My;  si_red.Mz = Mz;
si_red.Jexp = [real(diag(Mx)).'*p; real(diag(My)).'*p; real(diag(Mz)).'*p];

% --- static cc susceptibility ratio (elastic ON both, matched convention) --------------
chi_v = invz_chi0z(si_red,  T, 0, struct('elastic', true));
chi_f = invz_chi0z(si_full, T, 0, struct('elastic', true));
chi0_cc_vertex = real(chi_v(3,3,1));
chi0_cc_full   = real(chi_f(3,3,1));
chi_share = chi0_cc_vertex / chi0_cc_full;

% --- equal-time connected Jz variance captured (the vertex-relevant weight) -------------
var_trunc = real(diag(Mz*Mz)).'*p          - (real(diag(Mz)).'*p)^2;
var_full  = real(diag(si_full.Mz*si_full.Mz)).'*si_full.P ...
                                           - (real(diag(si_full.Mz)).'*si_full.P)^2;
var_share = var_trunc / var_full;

% --- isolation / coverage diagnostics --------------------------------------------------
if nv < n_full
    gap_16_17 = si_full.E(nv+1) - si_full.E(nv);
else
    gap_16_17 = NaN;
end
cluster_width = E(end) - E(1);
gap_ratio = gap_16_17 / max(cluster_width, eps);
gap_kBT   = gap_16_17 / (C.kB*T);

% --- optional adjacent-field continuity (exercises the caller's sweep path) ------------
min_subspace_overlap = NaN;
projector_distance   = NaN;
if isfield(opts, 'vb_prev') && ~isempty(opts.vb_prev) && isfield(opts.vb_prev, 'V16')
    s = svd(opts.vb_prev.V16' * V16);
    min_subspace_overlap = min(s);
    projector_distance   = sqrt(max(0, 1 - min_subspace_overlap^2));
end

% --- assemble (contract fields first, then diagnostics) --------------------------------
vb.E          = E;
vb.Mx         = Mx;
vb.My         = My;
vb.Mz         = Mz;
vb.p          = p;
vb.p_mass     = p_mass;
vb.V16        = V16;
vb.n_full     = n_full;
vb.n_vertex   = nv;
vb.chi_share  = chi_share;
vb.var_share  = var_share;
vb.gap_16_17  = gap_16_17;
vb.cluster_width = cluster_width;
vb.gap_ratio  = gap_ratio;
vb.gap_kBT    = gap_kBT;
vb.min_subspace_overlap = min_subspace_overlap;
vb.projector_distance   = projector_distance;
% diagnostics (not asserted by the gates)
vb.Jexp           = si_red.Jexp;
vb.chi0_cc_vertex = chi0_cc_vertex;
vb.chi0_cc_full   = chi0_cc_full;
vb.si_vertex      = si_red;
vb.T              = T;
end
