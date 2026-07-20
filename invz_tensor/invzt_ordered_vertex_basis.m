function vb = invzt_ordered_vertex_basis(ion, T, si_full, opts)
%INVZT_ORDERED_VERTEX_BASIS  Fixed-count, FIELD-ADAPTED ordered dominant-vertex basis (a3d).
%   vb = INVZT_ORDERED_VERTEX_BASIS(ion, T, si_full, opts) builds the 16-state dominant
%   vertex basis for the ordered-phase full-response solve as the LOWEST n_vertex eigenstates
%   of the ORDERED MEAN-FIELD Hamiltonian at (T, B) -- i.e. a DIRECT TRUNCATION of the
%   converged full electronuclear single ion si_full. Fixed COUNT (16 by default),
%   field-ADAPTED CONTENT: the basis tracks the field-rotated low-energy manifold instead of
%   a frozen zero-field subspace.
%
%   WHY FIELD-ADAPTED (executed falsification, 2026-07-20). Round-2 review P0-3 recommended a
%   FIELD-INDEPENDENT content-defined basis: the zero-field electronic ground doublet tensored
%   with I8 (INVZT_RUNG_BASIS 'e2xI8', dim 16). Task 7B MEASURED that this captures ~NONE of
%   the ordered static cc response and its subspace drifts from the true low manifold under
%   the transverse field:
%
%       B (T)   overlap_deficit   chi_share (cc, elastic)
%       3.00        1.56264        2.8e-13
%       3.50        1.67335        ~0
%       4.00        1.76804        ~0
%       4.40        1.83487        1e-5
%       4.65        1.87331        1.3e-4
%
%   Root cause: the zero-field ground doublet is an Ising doublet, <+|Jz|-> = 0, and the
%   ordering field hz = J0z<Jz> fully polarizes it (reduced Jz-fluctuation -> 0). The entire
%   longitudinal (cc) fluctuation lives in VIRTUAL transitions to CF excited states, which a
%   doublet-only basis cannot represent. But 16 states ARE enough IF field-adapted: the 16
%   LOWEST eigenstates of the ordered Hamiltonian capture chi_share = 0.977 at B = 3.0 T
%   (measured), rising to >0.99 by n = 32. So 7B fixes the COUNT and adapts the CONTENT.
%
%   COUNT-DRIFT (P0-3) RESOLUTION. The objection to the instantaneous E < Esplit ENERGY cut
%   was that it selects a drifting 13/13/11/10/9/8 -state set across the field sweep -- a
%   moving DIMENSION that injects cutoff jumps. This basis fixes the dimension at n_vertex
%   (16, constant by construction) and keeps only the CONTENT field-adapted; continuity is
%   then a measurable diagnostic (vb.chi_share smoothness + adjacent-field subspace overlap of
%   vb.V16), not a dimensional discontinuity. vb.gap_16_17 monitors the truncation edge for
%   level crossings that would make the fixed-count truncation ambiguous.
%
%   SIGNATURE DEVIATION (documented): the original Task-7B brief wrote (ion, si_full, opts);
%   vb.p (Boltzmann populations of vb.E) needs beta = 1/(kB T), so T is threaded in as the 2nd
%   argument. si_full is the CONVERGED ordered mean-field single ion, exactly the struct
%   INVZT_SOLVE_POINT_ORDERED builds (invz_single_ion with hyp+order+J0z/Jxx0).
%
%   OPTIONS (getf): n_vertex (default 16) -- the fixed truncation count; errors
%   invzt:vertexCount if > n_full.
%
%   OUTPUT struct vb:
%     vb.E          [nv x 1]  si_full.E(1:nv), already ground-shifted (E(1)=0).
%     vb.Mz         [nv x nv] si_full.Mz(1:nv,1:nv): Jz matrix elements in the truncated basis.
%     vb.p          [nv x 1]  Boltzmann populations of vb.E at T, renormalized to sum 1.
%     vb.V16        [n_full x nv] si_full.V(:,1:nv): the vertex eigenvectors in the full
%                   space (for cross-field subspace-overlap diagnostics).
%     vb.n_full     numel(si_full.E) (136 for the full electronuclear ion).
%     vb.n_vertex   nv (= 16 by default).
%     vb.chi_share  static cc chi0 of the TRUNCATED system / full (elastic ON both): the
%                   fraction of chi0_cc(0) captured by the truncation. DIRECTION GATE (> 0.9).
%     vb.gap_16_17  si_full.E(nv+1) - si_full.E(nv): truncation-edge gap (level-crossing
%                   monitor); NaN if nv == n_full.
%   Diagnostics (not asserted): vb.Mx, vb.My, vb.Jexp, vb.chi0_cc_vertex, vb.chi0_cc_full,
%   vb.si_vertex (the reduced si-like struct fed to chi0), vb.T.
%
%   See also INVZT_SOLVE_POINT_ORDERED, INVZ_SINGLE_ION, INVZ_CHI0Z, INVZT_SOLVE_POINT
%   (local_rung_si), INVZT_RUNG_BASIS.
if nargin < 4 || isempty(opts), opts = struct(); end

C    = invz_const();
beta = 1/(C.kB*T);

% --- guard: si_full must be a full converged single-ion struct -------------------------
for f = {'E', 'V', 'Mx', 'My', 'Mz'}
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

% --- direct truncation of the ordered mean-field single ion (fixed count, adapted content)
idx = 1:nv;
E   = si_full.E(idx);   E = E(:);                     % already ground-shifted
Mz  = si_full.Mz(idx, idx);
Mx  = si_full.Mx(idx, idx);
My  = si_full.My(idx, idx);
V16 = si_full.V(:, idx);
p   = exp(-beta*E);  p = p/sum(p);                    % Boltzmann over vb.E, renormalized

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

% --- truncation-edge gap (level-crossing monitor) --------------------------------------
if nv < n_full
    gap_16_17 = si_full.E(nv+1) - si_full.E(nv);
else
    gap_16_17 = NaN;
end

% --- assemble (contract fields first, then diagnostics) --------------------------------
vb.E          = E;
vb.Mz         = Mz;
vb.p          = p;
vb.V16        = V16;
vb.n_full     = n_full;
vb.n_vertex   = nv;
vb.chi_share  = chi_share;
vb.gap_16_17  = gap_16_17;
% diagnostics (not asserted by the gates)
vb.Mx             = Mx;
vb.My             = My;
vb.Jexp           = si_red.Jexp;
vb.chi0_cc_vertex = chi0_cc_vertex;
vb.chi0_cc_full   = chi0_cc_full;
vb.si_vertex      = si_red;
vb.T              = T;
end
