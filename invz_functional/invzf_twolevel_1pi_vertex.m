function out = invzf_twolevel_1pi_vertex(Delta, M, h, beta, labels)
%INVZF_TWOLEVEL_1PI_VERTEX Exact scalar dynamic gamma2/gamma3/gamma4 gate.
%
%   OUT = INVZF_TWOLEVEL_1PI_VERTEX(...,LABELS) amputates the exact
%   connected two-level cumulant for two, three, or four conserved
%   Matsubara labels.  For rank four, the three gamma3*C2*gamma3 channels
%   required by local Legendre inversion are included explicitly.

validateattributes(labels, {'numeric'}, {'real','vector','finite','integer'});
labels = labels(:).';
rank = numel(labels);
if rank < 2 || rank > 4 || sum(labels) ~= 0
    error('invzf:twolevel1piLabels', ...
        'labels must have rank two through four and sum to zero.');
end

cumulant = invzf_twolevel_cumulant(Delta,M,h,beta,labels);
Cleg = zeros(1,rank);
for k = 1:rank
    c2 = invzf_twolevel_cumulant(Delta,M,h,beta,[labels(k),-labels(k)]);
    Cleg(k) = real(c2.connected);
end
if ~strcmp(cumulant.status,'ok') || any(~isfinite(Cleg)) || any(Cleg <= 0)
    out = struct('status','invalid_local','rank',rank,'labels',labels, ...
        'connected',cumulant.connected,'C2_legs',Cleg,'gamma',NaN, ...
        'conversion_channels',[]);
    return
end

conversion = [];
switch rank
    case 2
        gamma = 1/Cleg(1);
    case 3
        gamma = -cumulant.connected/prod(Cleg);
    case 4
        pairing = [1 2 3 4;1 3 2 4;1 4 2 3];
        conversion = complex(zeros(3,1));
        for k = 1:3
            a = pairing(k,1);
            b = pairing(k,2);
            c = pairing(k,3);
            d = pairing(k,4);
            internal = -(labels(a)+labels(b));
            left = invzf_twolevel_1pi_vertex( ...
                Delta,M,h,beta,[labels(a),labels(b),internal]);
            right = invzf_twolevel_1pi_vertex( ...
                Delta,M,h,beta,[-internal,labels(c),labels(d)]);
            line = invzf_twolevel_cumulant( ...
                Delta,M,h,beta,[internal,-internal]);
            if ~strcmp(left.status,'ok') || ~strcmp(right.status,'ok') ...
                    || ~strcmp(line.status,'ok')
                out = struct('status','invalid_conversion','rank',rank, ...
                    'labels',labels,'connected',cumulant.connected, ...
                    'C2_legs',Cleg,'gamma',NaN, ...
                    'conversion_channels',conversion);
                return
            end
            conversion(k) = left.gamma*line.connected*right.gamma;
        end
        gamma = -cumulant.connected/prod(Cleg)+sum(conversion);
end

scale = max(1,abs(gamma));
if abs(imag(gamma)) <= 1024*eps(scale), gamma = real(gamma); end
out = struct('status','ok','rank',rank,'labels',labels, ...
    'connected',cumulant.connected,'C2_legs',Cleg,'gamma',gamma, ...
    'conversion_channels',conversion,'cumulant',cumulant);
end
