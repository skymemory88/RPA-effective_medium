function rb = invzt_rung_basis(ion, rung, opts)
%INVZT_RUNG_BASIS  Basis-content-defined state-space rung -> projector + basis energies.
%
%   rb = INVZT_RUNG_BASIS(ion, rung, opts) builds the reduced-basis PROJECTOR for one
%   rung of the A4 state-space ladder (docs/superpowers/plans/2026-07-17-invz-tensor-
%   full.md, Task 13). A rung is defined by BASIS CONTENT, NEVER by a lowest-N cut of
%   the 136-dim electronuclear spectrum:
%
%     'e3' | 'e6' | 'e17'          the lowest 3 / 6 / 17 CF states of the ZERO-FIELD
%                                  crystal-field Hamiltonian as an ELECTRONIC basis
%                                  (no nuclear space). 'e17' is the full J = 8 electronic
%                                  space (17 states).
%     'e3xI8' | 'e6xI8' | 'e17xI8' the same electronic subspace tensored with the COMPLETE
%                                  I = 7/2 nuclear space (8 states). 'e17xI8' = 136 = full.
%
%   MULTIPLET-COMPLETENESS (LOCKED, v3 review Other 6). The nominal count eN is NOMINAL:
%   if the lowest-N cut would SPLIT a degeneracy multiplet of the CF spectrum, the WHOLE
%   multiplet is included and the ACTUAL electronic dimension is recorded in rb.dim_el
%   (and rb.dim_actual = rb.dim_el * nuclear_dim). Completion extends UPWARD from the
%   nominal edge until the next CF gap exceeds opts.degtol: any multiplet the edge lands
%   inside is thereby completed (all states below the edge are already in). The ground
%   doublet is thus never split by a cut.
%
%   OUTPUT struct rb:
%     rb.projector          [n_full, dim_actual] complex. Columns are an orthonormal basis
%                           of the rung subspace in the FULL single-ion space
%                           (n_full = 17 electronic, or 136 for xI8). For xI8 the projector
%                           is kron(P_el, eye(8)) (Cartesian/CF-major, nuclear-minor -- the
%                           SAME ordering invz_single_ion's kron(M, eye(nI)) uses).
%     rb.dim_actual         size(rb.projector, 2) = rb.dim_el * rb.nuclear_dim.
%     rb.dim_el             electronic basis dimension (>= nominal; larger iff a multiplet
%                           completed at the edge).
%     rb.nuclear_dim        8 for a xI8 rung, else 1.
%     rb.E_basis            [dim_el, 1] the ground-shifted CF energies of the basis states
%                           (the CF structure that DEFINES the rung; nuclear multiplicity
%                           is recorded separately in nuclear_dim/dim_actual).
%     rb.multiplet_complete logical: the electronic basis has NO partial multiplet at its
%                           top edge (true by construction; a genuine post-condition check
%                           against the first omitted CF state).
%     rb.nominal            nominal count (3/6/17). rb.base_label 'eN' (nuclear-stripped).
%     rb.hyp                logical (rb.nuclear_dim > 1); rb.rung the echoed label.
%     rb.degtol             the degeneracy tolerance used.
%
%   The projector is a DIAGNOSTIC/SOLVER input, not physics on its own: INVZT_SOLVE_POINT
%   (nlevels rung) projects the field/MF single-ion Hamiltonian into this basis and
%   diagonalizes there ("build the Hamiltonian IN the reduced basis"). The rung's static
%   chi0 vs the full-136 chi0 (a driver-level diagnostic) DIAGNOSES -- does NOT bound --
%   the missing virtual intermediate states.
%
%   See also INVZT_SOLVE_POINT, INVZT_RUN_LADDER, INVZ_SINGLE_ION, STEVENS_OPS.
if nargin < 3 || isempty(opts), opts = struct(); end
degtol = getf(opts, 'degtol', 1e-6);   % 1e-6 meV: separates the ~1e-14 numerical doublet
                                       % degeneracy from the smallest real CF gap (~0.48 meV)

% --- parse the rung label: 'eN' with an optional 'xI8' nuclear-product suffix ----------
[base, nomN, nucdim] = parse_rung(rung);

% --- zero-field electronic crystal-field Hamiltonian (EXACTLY invz_single_ion's Hcf) ---
oJ = stevens_ops(ion.J);
B44s = 0;  if isfield(ion, 'B44s'), B44s = ion.B44s; end
Hcf = ion.B20*oJ.O20 + ion.B40*oJ.O40 + ion.B44*oJ.O44 + B44s*oJ.O44s ...
    + ion.B60*oJ.O60 + ion.B64c*oJ.O64c + ion.B64s*oJ.O64s;
Hcf = (Hcf + Hcf')/2;
nel = size(Hcf, 1);                                    % 17 for J = 8
if nomN > nel
    error('invzt:rungTooLarge', 'rung ''%s'' nominal %d exceeds the electronic dimension %d.', ...
        char(rung), nomN, nel);
end
[V, Dv] = eig(Hcf, 'vector');
[Eall, ix] = sort(real(Dv));  V = V(:, ix);
Eall = Eall - Eall(1);                                 % ground-shifted

% --- multiplet-complete the nominal cut (extend upward past any degeneracy) -------------
dim_el = nomN;
while dim_el < nel && (Eall(dim_el+1) - Eall(dim_el)) < degtol
    dim_el = dim_el + 1;
end
% post-condition: the basis top edge is NOT inside a multiplet (checked, not assumed).
multiplet_complete = (dim_el == nel) || ((Eall(dim_el+1) - Eall(dim_el)) >= degtol);

Pel = V(:, 1:dim_el);                                  % [nel, dim_el] electronic projector
E_basis = Eall(1:dim_el);

% --- nuclear product (xI8): kron(P_el, eye(nI)) into the full 136-dim space -------------
if nucdim > 1
    oI = stevens_ops(ion.I);
    nI = size(oI.Jz, 1);                               % 8 for I = 7/2
    if nI ~= nucdim
        error('invzt:rungNuclear', 'xI8 expects I=7/2 (dim 8); got nuclear dim %d.', nI);
    end
    projector = kron(Pel, eye(nI));                    % [nel*nI, dim_el*nI]
else
    projector = Pel;
end

rb = struct();
rb.rung               = char(rung);
rb.base_label         = base;
rb.nominal            = nomN;
rb.projector          = projector;
rb.dim_el             = dim_el;
rb.nuclear_dim        = nucdim;
rb.dim_actual         = size(projector, 2);            % = dim_el * nucdim
rb.E_basis            = E_basis;
rb.multiplet_complete = multiplet_complete;
rb.hyp                = nucdim > 1;
rb.degtol             = degtol;
end

% ------------------------------------------------------------------------------------- %
function [base, nomN, nucdim] = parse_rung(rung)
% 'e3' -> ('e3', 3, 1);  'e17xI8' -> ('e17', 17, 8). Errors on any other shape.
s = char(rung);
nucdim = 1;
tok = regexp(s, '^e(\d+)(xI8)?$', 'tokens', 'once');
if isempty(tok)
    error('invzt:rungLabel', ...
        ['unknown rung ''%s''. Basis-content rungs are ''eN'' (lowest-N CF states) with ' ...
         'an optional ''xI8'' nuclear product, e.g. e3, e6, e17, e3xI8, e17xI8. ' ...
         '(The ''three'' toy is handled by invzt_threestate, not invzt_rung_basis.)'], s);
end
nomN = str2double(tok{1});
if ~isempty(tok{2}), nucdim = 8; end
base = sprintf('e%d', nomN);
end
