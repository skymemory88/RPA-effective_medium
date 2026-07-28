function fac = invz_ordered_node_jacobian_factors(node, u)
%INVZ_ORDERED_NODE_JACOBIAN_FACTORS Exact bordered low-rank Jacobian.
%
%   For u=[Sigma(:);K0], the fixed-h Jacobian is represented as
%
%     J = [diag(d)+U*V'   b
%          c'             e].
%
%   U and V have two columns. No approximation or rank truncation is used.

nw = numel(node.wn);
if ~isnumeric(u) || ~isreal(u) || ~isvector(u) || numel(u) ~= nw+1 || ...
        any(~isfinite(u),'all')
    error('invzp:OrderedFactors:InvalidInput', ...
        'u must be a finite real vector with numel(node.wn)+1 entries.');
end
u = u(:);
Sigma = u(1:nw);
K0 = u(end);
G0 = node.G0(:);
Jf = node.Jnu_flat(:);
D = 1+Sigma;

den = D.'+Jf*G0.';
A = mean(Jf./den,1).';
Ap = -mean(Jf./(den.^2),1).';
Hk = 1-A.*G0;
Nk = A.*D;
kp = ((Ap.*D+A).*Hk+Nk.*Ap.*G0)./(Hk.^2);

nvar = nw+1;
dK = zeros(nw,nvar);
dK(2:nw,2:nw) = diag(kp(2:nw));
dK(1,end) = 1;
dK0 = zeros(1,nvar);
dK0(end) = 1;

g = node.g(:);
wts = node.wts(:);
L = zeros(3,nvar);
for p = 1:3
    L(p,:) = ((wts.*g.^p).'/node.beta)*dK;
end

tl = node.tl;
pref = tl.M2/tl.n01^2;
c0 = 1-tl.n01^2;
c1 = 0.5*(tl.g0+node.beta*c0);
amPref = tl.m^2/tl.n01^2;
amK = c0*(1+0.5*node.beta*tl.Delta*tl.n01)*tl.g0;
mfac = 2*tl.m^2/tl.M2;

rowA = (pref-amPref)*L(2,:)+ ...
    (-pref*c1+amPref*tl.g0)*L(1,:)- ...
    amPref*(4/tl.g0)*L(3,:)+amPref*amK*dK0;
rowG = pref*((1-mfac)*L(1,:)+mfac*c0*dK0);

diagonal = -ones(nw,1);
diagonal(2:nw) = diagonal(2:nw)-pref*c0*g(2:nw).*kp(2:nw);
left = [ones(nw,1),g];
right = [rowA(1:nw).',rowG(1:nw).'];
K0Column = rowA(end)*ones(nw,1)+rowG(end)*g;
K0Column(1) = K0Column(1)-pref*c0*g(1);

Kdyn = A.*D./Hk;
lam = invz_lambdas([K0;Kdyn(2:end)],g,wts,node.beta,[1 2 3]);
a = node.G0inel0;
b = node.G0el0;
d0 = 1+Sigma(1)+K0*a;
dd0 = zeros(1,nvar);
dd0(1) = 1;
dd0(end) = a;
t = tl.m^2*tl.n01^2*node.beta*K0-tl.M2*node.beta*lam(1);
numXi = 1+tanh(t);
denXi = 1+(4*tl.n01^2*K0*tl.g0+2*lam(2)+ ...
    tl.g0*lam(1))*tl.M2/tl.n01^2;
dt = tl.m^2*tl.n01^2*node.beta*dK0-tl.M2*node.beta*L(1,:);
dnumXi = (1-tanh(t)^2)*dt;
ddenXi = (tl.M2/tl.n01^2)*(4*tl.n01^2*tl.g0*dK0+ ...
    2*L(2,:)+tl.g0*L(1,:));
xi = numXi/denXi;
dxi = (dnumXi*denXi-numXi*ddenXi)/(denXi^2);
Hz = a+xi*b*d0;
dHz = b*(d0*dxi+xi*dd0);
z = d0/Hz;
dz = (dd0*Hz-d0*dHz)/(Hz^2);

Jscale = max(abs(Jf));
if isinf(z) || abs(z) > Jscale/sqrt(eps)
    rho = Hz/d0;
    drho = (dHz*d0-Hz*dd0)/(d0^2);
    couplingOffset = Jf-K0;
    F = 1./(1+rho*couplingOffset);
    Qrho = mean(F);
    Nrho = mean(Jf.*F);
    dF = -(F.^2).*(couplingOffset*drho-rho*dK0);
    dQ = mean(dF,1);
    dN = mean(Jf.*dF,1);
    dJloc = (Qrho*dN-Nrho*dQ)/(Qrho^2);
else
    E = z+Jf-K0;
    Gbar = mean(1./E);
    N = mean(Jf./E);
    dE = dz-dK0;
    dGbar = -mean(1./(E.^2))*dE;
    dN = -mean(Jf./(E.^2))*dE;
    dJloc = (Gbar*dN-N*dGbar)/(Gbar^2);
end
staticRow = (dK0-dJloc)/Jscale;

fac = struct( ...
    'schema','invzp_ordered_node_jacobian_factors/v1', ...
    'sigma_diagonal',diagonal, ...
    'sigma_left',left, ...
    'sigma_right',right, ...
    'K0_column',K0Column, ...
    'static_sigma_row',staticRow(1:nw), ...
    'static_K0',staticRow(end), ...
    'z',z, ...
    'z_gradient',dz, ...
    'w',z-K0, ...
    'w_gradient',dz-dK0);
end
