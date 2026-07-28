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

% At the exact Z2 point, the general chain-rule form for d2C2/dh2 contains
% separately divergent terms when beta*Delta is tiny.  Evaluate the same
% Hermite limit in the dimensionless x=beta*Delta/2 variable.  For n=0,
% q2=3*(x*sech(x)^2-tanh(x))/x^3 -> -2.  For n~=0 the leading term is
% O(x^2); the series avoids subtracting two O(1) terms.
if h == 0
    x0 = beta*delta;
    q2 = zeros(size(wn));
    iz = (wn == 0);
    if abs(x0) <= 0.1
        x2 = x0^2;
        q2(iz) = -2+x2*(8/5+x2*(-34/35+x2*(496/945 ...
            +x2*(-2764/10395+x2*87376/675675))));
    else
        q2(iz) = 3*(x0*sech(x0)^2-tanh(x0))/x0^3;
    end
    inz = ~iz;
    if any(inz)
        k = 2*pi*wn(inz);
        if abs(x0) <= 0.1
            x2 = x0^2;
            a2 = -8*(12+k.^2)./(3*k.^4);
            a4 = 32*(120+10*k.^2+k.^4)./(15*k.^6);
            a6 = -8*(20160+1680*k.^2+168*k.^4+17*k.^6)./(105*k.^8);
            a8 = 64*(362880+30240*k.^2+3024*k.^4+306*k.^6 ...
                +31*k.^8)./(2835*k.^10);
            q2(inz) = x2.*(a2+x2.*(a4+x2.*(a6+x2.*a8)));
        else
            den0 = 4*x0^2+k.^2;
            f0x = 4*x0*tanh(x0)./den0;
            fp0x = 4*(tanh(x0)+x0*sech(x0)^2)./den0 ...
                -32*x0^2*tanh(x0)./den0.^2;
            q2(inz) = fp0x/x0-2*f0x/x0^2;
        end
    end
    dC2dh(:) = 0;
    d2C2dh2 = M^4*beta^3*q2;
end

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
