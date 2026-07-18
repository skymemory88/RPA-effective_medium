function Jm = invzt_odd_mask(J)
%INVZT_ODD_MASK Zero the Cartesian-off-diagonal entries of every 3x3
%sublattice block of a [12,12] or [12,12,nq] coupling tensor (the odd=false
%rule).
%   Jm = INVZT_ODD_MASK(J) returns a copy of J with every Cartesian
%   cross-term (mu ~= nu, mu/nu in {a,b,c}) zeroed, for ALL sublattice pairs
%   (s,s' = 1..4) and ALL pages: the cc sector then no longer sees the
%   transverse (ab/ac/bc) blocks. J indexes i = 3*(s-1)+mu (mu = 1(a),2(b),
%   3(c), s = 1..4 sublattice; see INVZT_JQ_TENSOR), so keeping same-
%   Cartesian entries only (aa/bb/cc, any s,s' pair) is equivalent to zeroing
%   the off-diagonal of every 3x3 (mu,nu) sublattice block -- both
%   descriptions are the same [12,12] logical mask since it depends only on
%   mu = mod(i-1,3), not on s.
%
%   J : [12,12] or [12,12,nq] coupling tensor. Jm has the same size; the
%       [12,12] mask broadcasts over pages.
%
%   PROVENANCE: extracted from INVZT_SOLVE_POINT's inline odd=false rule
%   (R1, 2026-07-18 review) into a shared helper so INVZT_CHI_REALAXIS's
%   explicit-q continuation can apply the IDENTICAL rule (the continuation
%   previously consumed the unmasked tensor, giving an odd=false point
%   ODD-on couplings at finite q).
%
%   See also INVZT_SOLVE_POINT, INVZT_CHI_REALAXIS, INVZT_JQ_TENSOR.
cart = mod((0:11).', 3);
keepmask = (cart == cart.');                             % [12,12] logical
Jm = J .* keepmask;                                       % broadcasts over pages
end
