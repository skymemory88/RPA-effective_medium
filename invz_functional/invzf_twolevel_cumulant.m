function out = invzf_twolevel_cumulant(Delta, M, h, beta, labels)
%INVZF_TWOLEVEL_CUMULANT Exact connected C2/C3/C4 Matsubara vertex.
%
%   OUT = INVZF_TWOLEVEL_CUMULANT(DELTA,M,H,BETA,LABELS) evaluates the
%   connected local cumulant of X=M*sigma_x for two, three, or four integer
%   bosonic Matsubara labels whose sum is zero.  The convention is
%
%     X(tau)=beta^(-1)*sum_n X_n*exp(-i*omega_n*tau).
%
%   One operator time is fixed at zero.  Every ordering of the remaining
%   times is integrated exactly.  An ordered k-simplex is evaluated as the
%   (1,k+1) entry of a small upper-bidiagonal matrix exponential; repeated
%   energies and zero frequencies therefore use the analytic Hermite limit
%   without explicit divided denominators.
%
%   At fourth order, all three disconnected C2*C2 pairings, including their
%   beta*delta frequency factors, are subtracted explicitly.

validateattributes(Delta, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(M, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(h, {'numeric'}, {'real','scalar','finite'});
validateattributes(beta, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(labels, {'numeric'}, {'real','vector','finite','integer'});
labels = labels(:).';
rank = numel(labels);
if rank < 2 || rank > 4
    error('invzf:twolevelCumulantRank', ...
        'labels must contain two, three, or four Matsubara indices.');
end
if sum(labels) ~= 0
    error('invzf:twolevelCumulantFrequency', ...
        'The Matsubara labels must sum exactly to zero.');
end

X = M*[0 1;1 0];
Hloc = [0 0;0 Delta]-h*X;
[V,D] = eig((Hloc+Hloc')/2,'vector');
[E,ord] = sort(real(D));
V = V(:,ord);
boltz_log = -beta*(E-min(E));
logZ = max(boltz_log)+log(sum(exp(boltz_log-max(boltz_log))));
logp = boltz_log-logZ;
Xe = V'*X*V;
p = exp(logp);
m = real(sum(p.*diag(Xe)));
A = Xe-m*eye(2);

full = ordered_full(E,logp,A,beta,labels);
connected = full;
disconnected = complex(0);
if rank == 4
    pairing = [1 2 3 4;1 3 2 4;1 4 2 3];
    for k = 1:3
        a = pairing(k,1);
        b = pairing(k,2);
        c = pairing(k,3);
        d = pairing(k,4);
        if labels(a)+labels(b) == 0 && labels(c)+labels(d) == 0
            Cleft = ordered_full(E,logp,A,beta,labels([a b]));
            Cright = ordered_full(E,logp,A,beta,labels([c d]));
            disconnected = disconnected+beta*Cleft*Cright;
        end
    end
    connected = full-disconnected;
end

scale = max(1,abs(connected));
if abs(imag(connected)) <= 512*eps(scale)
    connected = real(connected);
end
out = struct('status','ok','rank',rank,'labels',labels, ...
    'full',full,'disconnected',disconnected,'connected',connected, ...
    'm',m,'E',E,'logp',logp,'operator_centered',A, ...
    'convention','X(tau)=beta^-1 sum_n X_n exp(-i*omega_n*tau)');
if any(~isfinite([real(full),imag(full),real(connected),imag(connected)]))
    out.status = 'nonfinite';
end
end

function value = ordered_full(E,logp,A,beta,labels)
rank = numel(labels);
timed = rank-1;
orderings = perms(1:timed);
omega = 2*pi*labels/beta;
nstates = numel(E);
sequences = dec2base(0:nstates^rank-1,nstates)-'0'+1;
value = complex(0);
for ip = 1:size(orderings,1)
    order = orderings(ip,:);
    for is = 1:size(sequences,1)
        state = sequences(is,:);
        amplitude = complex(1);
        for k = 1:rank-1
            amplitude = amplitude*A(state(k),state(k+1));
        end
        amplitude = amplitude*A(state(rank),state(1));
        if amplitude == 0, continue; end
        increment = zeros(1,timed);
        for k = 1:timed
            increment(k) = E(state(k))-E(state(k+1)) ...
                +1i*omega(order(k));
        end
        nodes = [0,cumsum(increment)];
        value = value+amplitude*simplex_exp(nodes,beta,logp(state(1)));
    end
end
end

function value = simplex_exp(nodes,beta,log_weight)
% log_weight is folded into every divided-difference node.  This evaluates
% p_a times the simplex integral without separately forming a tiny p_a or a
% potentially large exp(beta*(E_a-E_b)).
n = numel(nodes);
T = diag(nodes+log_weight/beta);
T((1:n-1)*(n+1)) = 1;
R = expm(beta*T);
value = R(1,n);
end
