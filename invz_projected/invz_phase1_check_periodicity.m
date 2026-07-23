function res = invz_phase1_check_periodicity(ion, dpRng, qsample, Gs)
%INVZ_PHASE1_CHECK_PERIODICITY Item 2, reciprocal periodicity (docs/invzp_phase1_quadrature_
% prereg.md): "For a sample of integer reciprocal vectors G, the sorted eigenvalue-branch
% spectrum ... at q+G equals that at q within |DJ| <= AbsTol_J + RelTol_J*|J(q)|, AbsTol_J=1e-10
% meV, RelTol_J=1e-8. Test the sorted branch spectrum / coupling multiset -- NOT elementwise
% equality of raw sublattice matrices (which may differ by a reciprocal-space gauge
% transformation)." This is a property of invz_jq_modes + the real-space lattice sum at a given
% dpRng; it is NOT offset/convention-specific (no Phase-1 grid is built here) -- the driver runs
% it once per dpRng actually used in the audit.
%
% Calls invz_jq_modes DIRECTLY at q and q+G (both may lie outside [-0.5,0.5]^3 -- invz_jq_modes
% imposes no BZ-range restriction; it is a real-space exponential sum, well-defined and manifestly
% periodic in q with integer period). Reuses invz_jq_modes' own SORTED [1x4] output
% (sort(real(eig(Jcc))), invz_jq_modes.m:89) -- a gauge-invariant eigenvalue spectrum, never the
% raw 4x4 matrix. Does NOT reimplement invz_jq_modes/MF_dipole. All sample points (q and every
% q+G) are batched into ONE invz_jq_modes call so the q-independent real-space geometry is built
% only once.
%
% ion        invz_ion()-shaped struct.
% dpRng      real-space dipole cutoff.
% qsample    [ns x 3] sample of q-points (reduced coords); default: a small fixed, deterministic
%            5-point sample spanning generic (non-symmetry-special) reduced coordinates.
% Gs         [nG x 3] integer reciprocal vectors; default: a small fixed set covering single-axis,
%            negative, and combined integer shifts.
%
% OUTPUT res (struct): .AbsTol_J, .RelTol_J, .n_pairs, .pass (logical, true iff EVERY sampled
% (q,G,branch) triple agrees within tolerance), .max_violation_margin (max over all triples of
% max(|DJ|-tol); <=0 everywhere means PASS), .detail ([n_pairs x 1] struct array: q, G,
% maxAbsDiff, tol, pass).
if nargin < 3 || isempty(qsample)
    qsample = [0.10 0.20 0.30; -0.30 0.15 0.40; 0.00 0.25 -0.10; 0.45 -0.45 0.05; -0.20 -0.20 0.30];
end
if nargin < 4 || isempty(Gs)
    Gs = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 1 1 0; 1 1 1];
end
AbsTol_J = 1e-10;
RelTol_J = 1e-8;
ns = size(qsample, 1);
nG = size(Gs, 1);

qall = qsample;
for j = 1:nG
    qall = [qall; qsample + Gs(j,:)]; %#ok<AGROW> small (ns*(nG+1) rows), fine
end
[Jall, ~] = invz_jq_modes(ion, qall, struct('dpRng', dpRng, 'cache', true));
Jq = Jall(1:ns, :);   % [ns x 4], sorted branches at the unshifted samples

detail = repmat(struct('q', [0 0 0], 'G', [0 0 0], 'maxAbsDiff', 0, 'tol', 0, 'pass', true), ns*nG, 1);
k = 0;
worst_margin = -Inf;   % max(|DJ| - tol) across every (q,G,branch); <=0 everywhere means PASS
for j = 1:nG
    JqG = Jall(ns*j + 1 : ns*(j+1), :);   % [ns x 4]
    for i = 1:ns
        d      = abs(JqG(i,:) - Jq(i,:));
        tol_ij = AbsTol_J + RelTol_J * abs(Jq(i,:));
        k = k + 1;
        detail(k).q          = qsample(i,:);
        detail(k).G          = Gs(j,:);
        detail(k).maxAbsDiff = max(d);
        detail(k).tol        = max(tol_ij);
        detail(k).pass       = all(d <= tol_ij);
        worst_margin = max(worst_margin, max(d - tol_ij));
    end
end
res.AbsTol_J = AbsTol_J;
res.RelTol_J = RelTol_J;
res.n_pairs  = k;
res.max_violation_margin = worst_margin;
res.pass   = all([detail.pass]);
res.detail = detail;
end
