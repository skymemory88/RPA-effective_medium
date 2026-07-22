function [Gstat, out] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0)
%INVZ_GSTAT_ORDERED Jensen ordered elastic static single-site function (framework SS9.2, J 2.28-2.29),
% in the BOUNDARY-PRESERVING HYBRID parametrization (stage-2 plan SS3, round-2 P0-A):
%   Gstat = G0inel0/(1 + Sigma0 + K0*G0inel0)  +  xi*G0el0
%   xi    = (1 + tanh(m^2*n01^2*beta*K0 - M2*beta*lam(1)))
%           / (1 + (4*n01^2*K0*g0 + 2*lam(2) + g0*lam(1))*M2/n01^2)         (J 2.29)
% G-convention (G = -chi, meV^-1). The CALLER supplies the static weights: the production
% loops pass the FULL-electronuclear split (invz_chi0z elastic:false static, and
% elastic:true minus elastic:false), so the m -> 0 closure fixed point is the PM solver's
% own (round-2 P0-A); the transcription tests pass the two-level weights
% G0inel0 = -M2*g0, G0el0 = -m^2*h0, under which the first term's denominator becomes
% 1 + Sigma0 - M2*K0*g0 and Gstat is J 2.28 VERBATIM. Only xi is two-level (SS6).
% At w = 0 the K-cancellation does NOT go through -- Gstat inserts directly into the
% effective-medium form G(q,0) = Gstat/[1 + (J(q)-K0)*Gstat], and K0 itself requires the
% static closure of INVZ_EMT_STATIC_ORDERED (invz_emt_scalar's direct solve is
% ordinary-Dyson only).
% Identities pinned by test_invz_gstat_ordered (NEVER adjust signs to pass them, SS9):
%   bare (Sigma0=K0=lam=0):  xi = 1,  Gstat = G0bare = G0inel0 + G0el0,
%                            -G0bare = d<Jz>/d(hmf)                          (J 2.31)
%   m = 0 (G0el0 = 0):  Gstat = G0inel0/(1+Sigma0+K0*G0inel0);
%                       Gtil0 = G0inel0/(1+Sigma0);  r = 1+Sigma0   (any G0inel0)
% out: xi, h0 = beta*(1-n01^2), G0bare, Gtil0 = Gstat/(1-K0*Gstat), r = G0bare/Gtil0
% (the H_MF integrand of J 2.33).
m = tl.m;  M2 = tl.M2;  n01 = tl.n01;  g0 = tl.g0;
h0 = beta*(1 - n01^2);
xi = (1 + tanh(m^2*n01^2*beta*K0 - M2*beta*lam(1))) / ...
     (1 + (4*n01^2*K0*g0 + 2*lam(2) + g0*lam(1))*M2/n01^2);
Gstat  = G0inel0/(1 + Sigma0 + K0*G0inel0) + xi*G0el0;
G0bare = G0inel0 + G0el0;
Gtil0  = Gstat/(1 - K0*Gstat);
out = struct('xi', xi, 'h0', h0, 'G0bare', G0bare, 'Gtil0', Gtil0, 'r', G0bare/Gtil0);
end
