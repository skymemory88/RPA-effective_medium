function out = invzf_scalar_functional_local(loc, lattice, m, h, H, ring)
%INVZF_SCALAR_FUNCTIONAL_LOCAL Assemble the common functional from local data.
%
%   OUT = INVZF_SCALAR_FUNCTIONAL_LOCAL(LOC,LATTICE,M,H,HPHYS)
%   evaluates the strict scalar common functional using a precomputed local
%   generator LOC.  LOC may come from the analytic two-level helper or the
%   full electronuclear oracle.  An already evaluated RING may be supplied
%   as the sixth argument to avoid repeating an expensive local/root trial.

if nargin < 6 || isempty(ring)
    if isfield(lattice,'mode_weights'), qw = lattice.mode_weights; else, qw = []; end
    ring = invzf_ring_scalar(loc,lattice.Jmodes,qw);
end
validateattributes(m, {'numeric'}, {'real','scalar','finite'});
validateattributes(h, {'numeric'}, {'real','scalar','finite'});
validateattributes(H, {'numeric'}, {'real','scalar','finite'});
validateattributes(lattice.J0, {'numeric'}, {'real','scalar','finite'});
if ~strcmp(ring.status,'ok')
    out = struct('status',ring.status,'f',NaN,'u',NaN, ...
        'grad',nan(2,1),'hessian',nan(2),'loc',loc,'ring',ring);
    return
end

J0 = lattice.J0;
f = loc.f0+(h-H)*m-0.5*J0*m^2+ring.f;
u = loc.u0+(h-H)*m-0.5*J0*m^2+ring.dBetaf_dBeta;
grad = [h-H-J0*m; -loc.m+m+ring.dfdh];
hessian = [-J0,1;1,-loc.chi+ring.d2fdh2];
out = struct('status','ok','f',f,'u',u,'grad',grad, ...
    'hessian',hessian,'loc',loc,'ring',ring);
end
