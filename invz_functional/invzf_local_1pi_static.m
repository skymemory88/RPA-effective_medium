function v = invzf_local_1pi_static(loc)
%INVZF_LOCAL_1PI_STATIC Static scalar 1PI vertices from one local generator.
%
%   V = INVZF_LOCAL_1PI_STATIC(LOC) amputates the zero-frequency connected
%   response and its first two source derivatives:
%
%     gamma2 = 1/C
%     gamma3 = -C'/C^3
%     gamma4 = (3*C'^2 - C*C'')/C^5.
%
%   These are derivatives of the local Legendre transform with respect to
%   its static moment.  LOC may come from invzf_twolevel_local or
%   invzf_electronuclear_local, but all entries must come from the same
%   source-biased Hamiltonian.

required = {'status','wn','C2','dC2dh','d2C2dh2'};
for k = 1:numel(required)
    if ~isfield(loc,required{k})
        error('invzf:local1piField', 'Missing loc.%s.',required{k});
    end
end
iz = find(loc.wn(:) == 0);
if numel(iz) ~= 1
    error('invzf:local1piZero', 'loc.wn must contain exactly one zero-frequency entry.');
end
C = loc.C2(iz);
Cp = loc.dC2dh(iz);
Cpp = loc.d2C2dh2(iz);
if ~strcmp(loc.status,'ok') || ~all(isfinite([C Cp Cpp])) || C <= 0
    v = struct('status','invalid_local','C2',C,'C3static',Cp, ...
        'C4static',Cpp,'gamma2',NaN,'gamma3',NaN,'gamma4',NaN);
    return
end
v = struct('status','ok','C2',C,'C3static',Cp,'C4static',Cpp, ...
    'gamma2',1/C,'gamma3',-Cp/C^3, ...
    'gamma4',(3*Cp^2-C*Cpp)/C^5);
end
