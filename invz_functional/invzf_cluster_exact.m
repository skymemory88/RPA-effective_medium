function out = invzf_cluster_exact(Delta, M, h, J, beta, wn, sites)
%INVZF_CLUSTER_EXACT Dense exact scalar two-level finite-cluster oracle.
%
%   OUT = INVZF_CLUSTER_EXACT(DELTA,M,H,J,BETA,WN,SITES) diagonalizes
%
%     H = sum_i [DELTA*|1><1|_i - H_i*X_i]
%         - sum_{i<j} J(i,j)*X_i*X_j,  X=M*sigma_x.
%
%   H may be scalar or one source per site.  J must be finite, real
%   symmetric, and have zero diagonal.  The initial dense oracle is capped
%   at eight sites.  OUT contains F, U, site moments, the total uniform
%   static susceptibility, and connected local Matsubara susceptibilities
%   for SITES (default all sites).  KMS/Hermite divided differences are
%   evaluated without subtracting nearly equal populations.

validateattributes(Delta, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(M, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(h, {'numeric'}, {'real','vector','finite'});
validateattributes(beta, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(J, {'numeric'}, {'real','finite','2d','square'});
validateattributes(wn, {'numeric'}, {'real','vector','finite','integer'});
N = size(J,1);
if N < 1 || N > 8
    error('invzf:clusterSize', 'Dense exact clusters require 1 <= N <= 8.');
end
if isscalar(h)
    hsite = repmat(h,N,1);
elseif numel(h) == N
    hsite = h(:);
else
    error('invzf:clusterSource', ...
        'h must be scalar or contain one source for each of the %d sites.',N);
end
scale = max(1,max(abs(J),[],'all'));
if norm(J-J.',inf) > 64*eps(scale) || any(abs(diag(J)) > 64*eps(scale))
    error('invzf:clusterCoupling', 'J must be symmetric with zero diagonal.');
end
J = (J+J.')/2;
J(1:N+1:end) = 0;
if nargin < 7 || isempty(sites), sites = 1:N; end
validateattributes(sites, {'numeric'}, ...
    {'real','vector','finite','integer','>=',1,'<=',N});
sites = sites(:).';
if numel(unique(sites)) ~= numel(sites)
    error('invzf:clusterSites', 'sites must not contain duplicates.');
end

d = 2^N;
I2 = eye(2);
X1 = M*[0 1;1 0];
Xops = cell(N,1);
H = zeros(d);
for i = 1:N
    Hloc1 = [0 0;0 Delta]-hsite(i)*X1;
    Xops{i} = embed_one(X1,I2,i,N);
    H = H+embed_one(Hloc1,I2,i,N);
end
for i = 1:N
    for j = i+1:N
        if J(i,j) ~= 0
            H = H-J(i,j)*(Xops{i}*Xops{j});
        end
    end
end
H = (H+H')/2;
[V,D] = eig(H,'vector');
[E,ord] = sort(real(D));
V = V(:,ord);
E0 = E(1);
boltz = exp(-beta*(E-E0));
Zshift = sum(boltz);
p = boltz/Zshift;
F = E0-log(Zshift)/beta;
U = sum(p.*E);

m_site = zeros(N,1);
for i = 1:N
    Ae = V'*Xops{i}*V;
    m_site(i) = real(sum(p.*diag(Ae)));
end
Xtot = zeros(d);
for i = 1:N, Xtot = Xtot+Xops{i}; end
Xte = V'*Xtot*V;
m_total = sum(m_site);
Xtc = Xte-m_total*eye(d);
chi_uniform = real(lehmann(E,p,Xtc,beta,0));

wn = wn(:);
C_local = zeros(numel(wn),numel(sites));
for a = 1:numel(sites)
    i = sites(a);
    Ae = V'*Xops{i}*V;
    Ac = Ae-m_site(i)*eye(d);
    for k = 1:numel(wn)
        C_local(k,a) = real(lehmann(E,p,Ac,beta,wn(k)));
    end
end

out = struct('F',F,'U',U,'m_site',m_site,'m_total',m_total, ...
    'chi_uniform',chi_uniform,'C_local',C_local,'sites',sites, ...
    'wn',wn,'E',E,'p',p,'J',J,'h_site',hsite);
end

function A = embed_one(op,I2,site,N)
A = 1;
for k = 1:N
    if k == site, A = kron(A,op); else, A = kron(A,I2); end
end
end

function C = lehmann(E,p,A,beta,n)
omega = 2*pi*n/beta;
C = complex(0);
N = numel(E);
for r = 1:N
    for s = 1:N
        amp = A(r,s)*A(s,r);
        if amp == 0, continue; end
        d = E(r)-E(s);
        if n == 0
            if d == 0
                ratio = beta*p(r);
            elseif d < 0
                z = beta*d;
                if z < -50, ratio = p(r)/(-d);
                else,       ratio = p(r)*beta*exprel_local(z); end
            else
                z = -beta*d;
                if z < -50, ratio = p(s)/d;
                else,       ratio = p(s)*beta*exprel_local(z); end
            end
            C = C+ratio*amp;
        elseif d ~= 0
            if d < 0, dp = p(r)*expm1(beta*d);
            else,     dp = -p(s)*expm1(-beta*d); end
            C = C+dp*amp/(1i*omega+d);
        end
    end
end
end

function y = exprel_local(x)
if abs(x) < 1e-5
    y = 1+x/2+x^2/6+x^3/24+x^4/120+x^5/720;
else
    y = expm1(x)/x;
end
end
