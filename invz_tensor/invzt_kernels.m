function k = invzt_kernels()
%INVZT_KERNELS  KMS-scaled, exprel-stable divided-difference kernels (A3 tensor vertex).
%
%   k = INVZT_KERNELS() returns a struct of function handles implementing the
%   imaginary-time divided-difference kernels phi / I2 / I3 that appear in the
%   component-labelled four-/three-point path sums (invzt_vertex4 / invzt_vertex3).
%   Every handle takes beta EXPLICITLY (no captured global), so the same struct
%   evaluates rows at any temperature.  All handles are ELEMENTWISE on their last
%   frequency/energy argument (accept a scalar OR a vector), which lets the vertex
%   engine vectorise the innermost path-sum index.
%
%   PLAIN kernels (population-free divided differences; reference of record is the
%   mpmath oracle verify_tensor_vertex.py -> vertex_oracle.json):
%       k.phi(a,beta)          phi(a)   = (exp(beta a) - 1)/a           (exprel)
%       k.dphi(a,beta)         phi'(a),  k.d2phi(a,beta)  phi''(a)
%       k.I2(x,y,beta)         first divided difference of phi at {x, x+y}
%       k.I3(x,y,z,beta)       second divided difference at {x, x+y, x+y+z}
%
%   KMS-SCALED kernels (hard constraint 4).  The path sum weights each kernel by
%   the INITIAL-STATE population p_r.  At large beta*Delta the raw kernel carries
%   an exp(beta*(E_r-E_j)) endpoint layer ~ 1e86 while p_r ~ 1e-87; forming
%   p_r * (huge) risks overflow AND loses the sign/phase of the divided-difference
%   cancellations.  The *_s handles take the initial log-weight  logpr = log(p_r)
%   and return  p_r * (plain kernel)  computed in log-REGROUPED form: every bare
%   exp(beta a) becomes exp(logpr + beta a) and every bare 1 becomes exp(logpr).
%   Because  logpr + beta*Re(a) = log(p_r) + beta*(E_r - E_j) = log(p_j) <= 0, each
%   folded exponential has magnitude <= 1 -- no overflow, no catastrophic
%   cancellation, and the complex Matsubara phase (Im part of a) is carried
%   exactly.  (A bare real exp(log p_r + beta*x_acc) would drop that phase -- the
%   error the constraint forbids.)
%       k.phis(logpr,a,beta)   = p_r*phi(a),  k.dphis, k.d2phis
%       k.I2s(logpr,x,y,beta)  = p_r*I2(x,y)
%       k.I3s(logpr,x,y,z,beta)= p_r*I3(x,y,z)
%   The plain handles are exactly the scaled ones at logpr = 0 (p_r = 1).
%
%   Repeated-node (Hermite) branches (constraint 5): I2(x,0), I3(x,y,-y), I3(x,y,0),
%   I3(x,0,0) and complex arguments are all regular by construction via the exprel
%   Taylor limits below -- so the Bx=0 doublet (exact degeneracy) is NaN-free.
%
%   Branch thresholds match the oracle generator VERBATIM (phi 1e-12, dphi 1e-10,
%   d2phi 1e-8, divided-difference denominators 1e-10) so repeated-node rows select
%   the identical series limb.
%
%   See also INVZT_VERTEX4, INVZT_VERTEX3.

k.phi    = @(a,beta)     phis(0.0, a, beta);
k.dphi   = @(a,beta)     dphis(0.0, a, beta);
k.d2phi  = @(a,beta)     d2phis(0.0, a, beta);
k.I2     = @(x,y,beta)   I2s(0.0, x, y, beta);
k.I3     = @(x,y,z,beta) I3s(0.0, x, y, z, beta);

k.phis   = @phis;
k.dphis  = @dphis;
k.d2phis = @d2phis;
k.I2s    = @I2s;
k.I3s    = @I3s;
end

% ========================================================================= %
%  KMS-scaled primitives.  logpr is a real scalar (= log p_r); the frequency/
%  energy argument a (or the last of x/y/z) may be a scalar or an array.  Each
%  returns  p_r * (plain value), overflow-safe and phase-correct.
% ========================================================================= %

function out = phis(logpr, a, beta)
% p_r * phi(a) = (exp(logpr + beta a) - exp(logpr)) / a  ,  Hermite limit for a~0.
a  = complex(a);
Ea = exp(logpr + beta.*a);
Ep = exp(logpr);
out = (Ea - Ep) ./ a;
ser = Ep .* (beta + beta.^2 .* a ./ 2 + beta.^3 .* a.^2 ./ 6);
small = abs(a) <= 1e-12;
out(small) = ser(small);
end

function out = dphis(logpr, a, beta)
% p_r * phi'(a).
a  = complex(a);
Ea = exp(logpr + beta.*a);
Ep = exp(logpr);
out = (beta .* a .* Ea - (Ea - Ep)) ./ a.^2;
ser = Ep .* (beta.^2 ./ 2 + beta.^3 .* a ./ 3 + beta.^4 .* a.^2 ./ 8);
small = abs(a) <= 1e-10;
out(small) = ser(small);
end

function out = d2phis(logpr, a, beta)
% p_r * phi''(a).
a  = complex(a);
Ea = exp(logpr + beta.*a);
Ep = exp(logpr);
out = (beta.^2 .* a.^2 .* Ea - 2 .* (beta .* a .* Ea - (Ea - Ep))) ./ a.^3;
ser = Ep .* (beta.^3 ./ 3 + beta.^4 .* a ./ 4);
small = abs(a) <= 1e-8;
out(small) = ser(small);
end

function out = I2hs(logpr, x, Y, beta)
% First divided difference of (scaled) phi at {x, x+Y}.  x scalar; Y scalar/array.
out = (phis(logpr, x + Y, beta) - phis(logpr, x, beta)) ./ Y;
ser = dphis(logpr, x, beta) + Y .* d2phis(logpr, x, beta) ./ 2;
small = abs(Y) <= 1e-10;
out(small) = ser(small);
end

function out = dI2ys(logpr, x, Y, beta)
% d/dY of I2h(x,Y)  (the y-repeated-node limit of I3).  x scalar; Y scalar/array.
out = (dphis(logpr, x + Y, beta) ...
       - (phis(logpr, x + Y, beta) - phis(logpr, x, beta)) ./ Y) ./ Y;
ser = d2phis(logpr, x, beta) ./ 2;                 % scalar
small = abs(Y) <= 1e-10;
out(small) = ser;                                  % broadcast scalar into masked slots
end

function out = I2s(logpr, x, y, beta)
% p_r * I2(x,y).
out = I2hs(logpr, x, y, beta);
end

function out = I3s(logpr, x, y, z, beta)
% p_r * I3(x,y,z) = second divided difference of (scaled) phi at {x, x+y, x+y+z}.
% x,y scalar; z scalar OR array (elementwise, for the vectorised innermost index).
out   = (I2hs(logpr, x, y + z, beta) - I2hs(logpr, x, y, beta)) ./ z;
serz  = dI2ys(logpr, x, y, beta);                  % scalar (z -> 0 limit)
small = abs(z) <= 1e-10;
out(small) = serz;                                 % broadcast scalar into masked slots
end
