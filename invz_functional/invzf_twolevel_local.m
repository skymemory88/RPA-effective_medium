function loc = invzf_twolevel_local(Delta, M, h, beta, wn)
%INVZF_TWOLEVEL_LOCAL Exact source-biased scalar two-level local data.
%
%   LOC = INVZF_TWOLEVEL_LOCAL(DELTA, M, H, BETA, WN) evaluates
%
%       Hloc = DELTA*|1><1| - H*M*(|0><1| + |1><0|)
%
%   on integer bosonic Matsubara indices WN.  C2 is the positive connected
%   susceptibility.  dC2dh, d2C2dh2, and dC2dbeta are analytic derivatives
%   of that same C2, including its elastic zero-frequency contribution.
%
%   This is an isolated strict-functional helper.  It does not use the
%   production convention G=-C2 and does not apply an ordered static split.

validateattributes(Delta, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(M, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(h, {'numeric'}, {'real','scalar','finite'});
validateattributes(beta, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(wn, {'numeric'}, {'real','vector','finite','integer'});

wn = wn(:);
delta = Delta/2;
y = h*M;
E = hypot(delta, y);
x = beta*E;
t = tanh(x);
s2 = sech(x)^2;

% f0 includes the irrelevant-to-response DELTA/2 shift so that it equals
% -(1/beta)*log(1 + exp(-beta*DELTA)) at h=0.
log2cosh = abs(x) + log1p(exp(-2*abs(x)));
f0 = delta - log2cosh/beta;
u0 = delta-E*t;
m = h*M^2*t/E;
chi = M^2*delta^2*t/E^3 + beta*h^2*M^4*s2/E^2;

omega = 2*pi*wn/beta;
D = 4*E^2 + omega.^2;
b = delta^2/E^2;
db = -2*b/E;
d2b = 6*b/E^2;

tp = beta*s2;
tpp = -2*beta^2*s2*t;
u = 4*E*t;
up = 4*(t + E*tp);
upp = 4*(2*tp + E*tpp);
Dp = 8*E;
Dpp = 8;
L = u./D;
Lp = up./D - u.*Dp./D.^2;
Lpp = upp./D - 2*up.*Dp./D.^2 - u.*Dpp./D.^2 ...
       + 2*u.*Dp.^2./D.^3;

sp = -2*beta*s2*t;
spp = 4*beta^2*s2*t^2 - 2*beta^2*s2^2;
iszero = (wn == 0);
R = beta*(1-b)*s2;
Rp = beta*(-db*s2 + (1-b)*sp);
Rpp = beta*(-d2b*s2 - 2*db*sp + (1-b)*spp);

q = b*L + iszero.*R;
qp = db*L + b*Lp + iszero.*Rp;
qpp = d2b*L + 2*db*Lp + b*Lpp + iszero.*Rpp;

% Explicit beta derivative at fixed source and integer Matsubara index.
% omega_n=2*pi*n/beta therefore contributes through the denominator.
Lbeta = 4*E*(E*s2./D + 2*t*omega.^2./(beta*D.^2));
Rbeta = (1-b)*s2*(1-2*beta*E*t);
qbeta = b*Lbeta + iszero.*Rbeta;

E_h = h*M^2/E;
E_hh = M^2*delta^2/E^3;
C2 = M^2*q;
dC2dh = M^2*qp*E_h;
d2C2dh2 = M^2*(qpp*E_h^2 + qp*E_hh);
dC2dbeta = M^2*qbeta;

loc = struct();
loc.Delta = Delta;
loc.M = M;
loc.h = h;
loc.beta = beta;
loc.wn = wn;
loc.omega = omega;
loc.E = E;
loc.f0 = f0;
loc.u0 = u0;
loc.m = m;
loc.chi = chi;
loc.C2 = real(C2);
loc.dC2dh = real(dC2dh);
loc.d2C2dh2 = real(d2C2dh2);
loc.dC2dbeta = real(dC2dbeta);
loc.zero_mode_identity = loc.C2(wn == 0) - chi;
derived = [E;f0;u0;m;chi;loc.C2(:);loc.dC2dh(:); ...
    loc.d2C2dh2(:);loc.dC2dbeta(:)];
if all(isfinite(derived(:)))
    loc.status = 'ok';
else
    loc.status = 'nonfinite';
end
end
