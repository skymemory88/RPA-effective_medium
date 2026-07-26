function [Gstat, out] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0, opts)
%INVZ_GSTAT_ORDERED Jensen ordered elastic static single-site function (framework SS9.2, J 2.28-2.29),
% in the BOUNDARY-PRESERVING HYBRID parametrization (stage-2 plan SS3, round-2 P0-A):
%   Gstat = G0inel0/(1 + Sigma0 + K0*G0inel0)  +  xi*G0el0
%   xi    = (1 + tanh(m^2*n01^2*beta*K0 - M2*beta*lam(1)))
%           / (1 + (4*n01^2*K0*g0 + 2*lam(2) + g0*lam(1))*M2/n01^2)         (J 2.29)
% G-convention (G = -chi, meV^-1; ferromagnetic positive J). The CALLER supplies the static weights: the production
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
% opts.stable_form (default false) selects the arithmetic used for out.Gtil0/out.r:
%   false (default -- EVERY existing/resummed call site, preserved bit-identically, G9):
%     out.Gtil0 = Gstat/(1-K0*Gstat);  out.r = G0bare/Gtil0 -- the historical arithmetic.
%   true (strict mode only; wired by invz_emt_static_ordered under strict static_medium):
%     the EXACT reassociation 1/Gtil0 = 1/Gstat - K0, i.e.
%     out.Gtil0 = 1/(1/Gstat - K0);  out.r = G0bare*(1/Gstat - K0).
%     Algebraically identical to the false branch everywhere away from a pole (measured 0 ulp
%     over a representative sweep); NOT a regulariser -- no broadening, no added tolerance, no
%     sign change. Where Gstat itself diverges (its own local denominator
%     gstat_local_denom -> 0), the divergence CANCELS in this arrangement --
%     Gtil0 -> -1/K0, r -> -G0bare*K0 -- while the false-branch form evaluates
%     Inf/(-Inf) = NaN there. See test_invz_gstat_removable_pole (spec G17).
% out: xi, h0 = beta*(1-n01^2), G0bare, Gtil0, r (per opts.stable_form above -- the H_MF
%   integrand of J 2.33), and gstat_local_denom = 1 + Sigma0 + K0*G0inel0, the signed local
%   denominator of the first Gstat term (identical on both branches; exposed for G17/Gate 0).
m = tl.m;  M2 = tl.M2;  n01 = tl.n01;  g0 = tl.g0;
h0 = beta*(1 - n01^2);
xi = (1 + tanh(m^2*n01^2*beta*K0 - M2*beta*lam(1))) / ...
     (1 + (4*n01^2*K0*g0 + 2*lam(2) + g0*lam(1))*M2/n01^2);
gstat_local_denom = 1 + Sigma0 + K0*G0inel0;
Gstat  = G0inel0/gstat_local_denom + xi*G0el0;
G0bare = G0inel0 + G0el0;
if nargin >= 8 && getf(opts, 'stable_form', false)
% EXACT REASSOCIATION (not a regulariser -- identical value, no broadening, no tolerance):
%   1/Gtil0 = 1/Gstat - K0   <=>   Gtil0 = Gstat/(1 - K0*Gstat)
% Written this way because Gstat has its own pole where 1 + Sigma0 + K0*G0inel0 -> 0, and the
% divergence CANCELS in the quantities Jensen's condition actually integrates:
%   Gtil0 -> -1/K0,   r -> -G0bare*K0,   crit = r + J0eff*G0bare -> G0bare*(J0eff - K0)
% all finite, with the same limit from both sides. The former arrangement
% (Gtil0 = Gstat/(1-K0*Gstat); r = G0bare/Gtil0) evaluates Inf/(-Inf) = NaN at the crossing and
% loses precision approaching it, which would turn a removable singularity into a node failure.
% Pinned by test_invz_gstat_removable_pole.m (spec G17).
    invGtil0 = 1/Gstat - K0;
    Gtil0    = 1/invGtil0;
    r        = G0bare*invGtil0;
else
    % G9 compatibility: preserve the historical arithmetic on every existing/resummed call.
    Gtil0 = Gstat/(1 - K0*Gstat);
    r     = G0bare/Gtil0;
end
out = struct('xi', xi, 'h0', h0, 'G0bare', G0bare, 'Gtil0', Gtil0, 'r', r, ...
             'gstat_local_denom', gstat_local_denom);
end
