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
fac = invz_ordered_node_jacobian_factors(node,u);
Jsigma = diag(fac.sigma_diagonal)+fac.sigma_left*fac.sigma_right.';
J = [Jsigma,fac.K0_column;fac.static_sigma_row,fac.static_K0];
end
