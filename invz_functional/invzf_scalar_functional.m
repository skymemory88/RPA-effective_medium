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
ring = invzf_ring_scalar(loc, lattice.Jmodes, qweights);
if ~strcmp(ring.status, 'ok')
    out = struct('status', ring.status, 'f', NaN, 'u', NaN, ...
        'grad', nan(2,1), 'hessian', nan(2), 'loc', loc, 'ring', ring);
    return
end

J0 = lattice.J0;
f = loc.f0 + (h-H)*m - 0.5*J0*m^2 + ring.f;
u = loc.u0 + (h-H)*m - 0.5*J0*m^2 + ring.dBetaf_dBeta;
grad = [h-H-J0*m; -loc.m+m+ring.dfdh];
hessian = [-J0, 1; 1, -loc.chi+ring.d2fdh2];

out = struct('status', 'ok', 'f', f, 'u', u, 'grad', grad, ...
    'hessian', hessian, 'loc', loc, 'ring', ring);
end
