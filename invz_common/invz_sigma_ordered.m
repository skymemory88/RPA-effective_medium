function sig = invz_sigma_ordered(tl, lam, K, g, beta)
%INVZ_SIGMA_ORDERED Ordered-phase (m~=0) Jensen self-energy, framework SS9.2, J 2.26-2.27.
%   Sigma(iwn) = (alpha - alpha_m) + [ gamma(iwn) - (2 m^2 / M^2) gamma(0) ] g(iwn)   (eq 37)
%   alpha_m = (m^2/n01^2){ lambda_2 - g0*lambda_1 + (4/g0) lambda_3
%                          - (1-n01^2)(1 + 1/2 beta Delta n01) K(0) g0 }               (eq 38)
% Reduces to paramagnetic invz_sigma at m = 0 (alpha/gamma, eqs 23/27, are the SAME formula,
% so they are computed by calling invz_sigma directly rather than re-derived here). Requires
% lam = [lambda_1; lambda_2; lambda_3] (invz_lambdas(...,[1 2 3])). K(:) is on the bosonic
% n>=0 grid so K(1) is the omega_n = 0 slot.
s0 = invz_sigma(tl, lam, K, g, beta);
K = K(:);  g = g(:);
K0 = K(1);  g0 = tl.g0;
gamma0 = s0.gamma(1);                                                  % gamma at omega_n = 0

alpha_m = (tl.m^2 / tl.n01^2) * ( lam(2) - g0*lam(1) + (4/g0)*lam(3) ...
            - (1 - tl.n01^2)*(1 + 0.5*beta*tl.Delta*tl.n01)*K0*g0 );   % eq 38

sig.alpha   = s0.alpha;
sig.alpha_m = alpha_m;
sig.gamma   = s0.gamma;
sig.Sigma   = (s0.alpha - alpha_m) + (s0.gamma - (2*tl.m^2/tl.M2)*gamma0) .* g;
end
