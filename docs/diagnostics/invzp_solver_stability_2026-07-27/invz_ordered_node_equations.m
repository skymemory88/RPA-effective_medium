function [R, diagout, state, J] = invz_ordered_node_equations(node, u)
%INVZ_ORDERED_NODE_EQUATIONS Defactored residual and analytic Jacobian.
%
%   [R,DIAG,STATE,J] = INVZ_ORDERED_NODE_EQUATIONS(NODE,U) exposes the
%   exact diagnostic equations used by INVZ_ORDERED_NODE_NEWTON so a
%   branch driver can form a bordered pseudo-arclength system. Unknowns
%   are U=[Sigma(:);K0]. J is dR/dU and is nargout-gated.

if ~isstruct(node) || ~isscalar(node)
    error('invzp:OrderedEquations:InvalidInput','node must be a scalar struct.');
end
required = {'wn','G0','Jnu_flat','eopts','g','wts','beta','tl', ...
    'G0inel0','G0el0','eso'};
for k = 1:numel(required)
    if ~isfield(node,required{k})
        error('invzp:OrderedEquations:InvalidInput', ...
            'node.%s is required.',required{k});
    end
end
nw = numel(node.wn);
if ~isnumeric(u) || ~isreal(u) || ~isvector(u) || numel(u) ~= nw+1 || ...
        any(~isfinite(u),'all')
    error('invzp:OrderedEquations:InvalidInput', ...
        'u must be a finite real vector with numel(node.wn)+1 entries.');
end
u = u(:);
Sigma = u(1:nw);
K0 = u(end);
med = invz_emt_scalar(node.G0,Sigma,node.Jnu_flat,node.eopts);
K = [K0;med.K(2:end)];
lam = invz_lambdas(K,node.g,node.wts,node.beta,[1 2 3]);
sg = invz_sigma_ordered(node.tl,lam,K,node.g,node.beta);
[Gstat,go] = invz_gstat_ordered(node.tl,lam,K0,Sigma(1),node.beta, ...
    node.G0inel0,node.G0el0,struct('stable_form',true));

d0 = go.gstat_local_denom;
Hz = node.G0inel0+go.xi*node.G0el0*d0;
z = d0/Hz;
Gtil0 = 1/(z-K0);
r = go.G0bare*(z-K0);

Jf = node.Jnu_flat(:);
Jscale = max(abs(Jf));
if ~(isfinite(Jscale) && Jscale > 0)
    error('invz:orderedNewtonCoupling', ...
        'The static coupling scale must be positive and finite.');
end
if isinf(z)
    Gbar = 0;
    Jloc = mean(Jf);
    poleMargin = 1;
    meanMargin = 0;
    qCancel = 1;
    Q = 1;
else
    scale = max(abs(z),Jscale);
    E = z+Jf-K0;
    weights = scale./E;
    meanWeights = mean(weights);
    Gbar = meanWeights/scale;
    Jloc = mean(Jf.*weights)/meanWeights;
    poleMargin = min(abs(E))/scale;
    meanMargin = abs(Gbar)*Jscale;
    qCancel = abs(meanWeights)/mean(abs(weights));
    Q = z*Gbar;
end

R = [sg.Sigma-Sigma;(K0-Jloc)/Jscale];
state = struct('Sigma',Sigma,'K',K,'lam',lam,'K0s',K0);
diagout = struct('z',z,'Gstat',Gstat,'Gtil0',Gtil0,'r',r, ...
    'Gbar',Gbar,'pole_margin',poleMargin,'mean_margin',meanMargin, ...
    'gstat_local_denom',d0,'xi',go.xi,'h0',go.h0, ...
    'G0bare',go.G0bare,'q_cancel',qCancel,'Q',Q, ...
    'Hz',Hz,'rho',1/z,'Jloc',Jloc,'Jscale',Jscale);

if nargout >= 4
    J = analyticJacobian(u,node);
end
end

function J = analyticJacobian(u,node)
nw = numel(node.wn);
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

dalpha = pref*(L(2,:)-c1*L(1,:));
dalphaM = amPref*(L(2,:)-tl.g0*L(1,:)+ ...
    (4/tl.g0)*L(3,:)-amK*dK0);
dgamma = pref*(repmat(L(1,:),nw,1)-c0*dK);
dgamma0 = pref*(L(1,:)-c0*dK0);
dSigmaMap = repmat(dalpha-dalphaM,nw,1)+ ...
    g.*(dgamma-mfac*repmat(dgamma0,nw,1));
Jsigma = dSigmaMap-[eye(nw),zeros(nw,1)];

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
    c = Jf-K0;
    F = 1./(1+rho*c);
    Qrho = mean(F);
    Nrho = mean(Jf.*F);
    dF = -(F.^2).*(c*drho-rho*dK0);
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
J = [Jsigma;(dK0-dJloc)/Jscale];
end
