function Jq = build_jq(qvec, Rvec, J_of_R, options)
% BUILD_JQ Build interaction matrix J(q) from real-space couplings J(R).
%
% Inputs:
%   qvec   : n_q x 3 reduced reciprocal coordinates
%   Rvec   : n_R x 3 real-space lattice vectors in matching reduced units
%   J_of_R : 3 x 3 x n_R coupling tensor for each R
%   options (optional struct):
%       .phase_factor      (default 2*pi)
%       .symmetrize        (default true)
%       .demag_tensor      (default [])
%       .demag_prefactor   (default 1)
%
% Output:
%   Jq : 3 x 3 x n_q complex interaction matrices

if nargin < 4
    options = struct();
end
if ~isfield(options, 'phase_factor')
    options.phase_factor = 2 * pi;
end
if ~isfield(options, 'symmetrize')
    options.symmetrize = true;
end
if ~isfield(options, 'demag_tensor')
    options.demag_tensor = [];
end
if ~isfield(options, 'demag_prefactor')
    options.demag_prefactor = 1;
end

n_q = size(qvec, 1);
n_R = size(Rvec, 1);
if size(J_of_R, 3) ~= n_R
    error('build_jq:badShape', ...
        'J_of_R third dimension (%d) must match size(Rvec,1) (%d).', ...
        size(J_of_R, 3), n_R);
end

Jq = zeros(3, 3, n_q, 'like', J_of_R);
for iq = 1:n_q
    accum = zeros(3, 3, 'like', J_of_R);
    phase = options.phase_factor * (Rvec * qvec(iq, :).');
    ephase = exp(1i * phase);

    for iR = 1:n_R
        accum = accum + J_of_R(:,:,iR) * ephase(iR);
    end

    if ~isempty(options.demag_tensor)
        accum = accum - options.demag_prefactor * options.demag_tensor;
    end

    if options.symmetrize
        accum = 0.5 * (accum + accum');
    end

    Jq(:,:,iq) = accum;
end

end
