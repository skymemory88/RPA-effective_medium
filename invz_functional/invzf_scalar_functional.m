function out = invzf_scalar_functional(model, lattice, m, h, H, wn)
%INVZF_SCALAR_FUNCTIONAL Independently varied strict scalar (m,h) functional.
%
%   OUT = INVZF_SCALAR_FUNCTIONAL(MODEL, LATTICE, M, HREF, HPHYS, WN)
%   evaluates per site
%
%     f = f0(href) + (href-Hphys)*m - J0*m^2/2 + f_ring(href).
%
%   MODEL has Delta, M, beta.  LATTICE has J0, Jmodes and optional
%   mode_weights.  OUT.grad is ordered as [df/dm; df/dh], and OUT.hessian
%   uses the same variable order.

required_model = {'Delta','M','beta'};
required_lattice = {'J0','Jmodes'};
for k = 1:numel(required_model)
    if ~isfield(model, required_model{k})
        error('invzf:modelField', 'Missing model.%s.', required_model{k});
    end
end
for k = 1:numel(required_lattice)
    if ~isfield(lattice, required_lattice{k})
        error('invzf:latticeField', 'Missing lattice.%s.', required_lattice{k});
    end
end
validateattributes(m, {'numeric'}, {'real','scalar','finite'});
validateattributes(h, {'numeric'}, {'real','scalar','finite'});
validateattributes(H, {'numeric'}, {'real','scalar','finite'});
validateattributes(lattice.J0, {'numeric'}, {'real','scalar','finite'});
if isfield(lattice, 'mode_weights')
    qweights = lattice.mode_weights;
else
    qweights = [];
end

loc = invzf_twolevel_local(model.Delta, model.M, h, model.beta, wn);
ring = invzf_ring_scalar(loc,lattice.Jmodes,qweights);
out = invzf_scalar_functional_local(loc,lattice,m,h,H,ring);
end
