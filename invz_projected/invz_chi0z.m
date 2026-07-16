function chi = invz_chi0z(si, T, z, opts)
%INVZ_CHI0Z Single-ion susceptibility tensor at arbitrary complex frequencies.
% chi(mu,nu,iz) = sum_{a,b inelastic} (p_a-p_b) M_mu(a,b) M_nu(b,a) / (E_b-E_a - z(iz))
%              [+ elastic beta-term on entries with |z|<ztol].
if nargin < 4, opts = struct(); end
degtol = getf(opts, 'degtol', 1e-8);
ztol   = getf(opts, 'ztol', 1e-12);
elast  = getf(opts, 'elastic', true);
C = invz_const();  beta = 1/(C.kB*T);
E = si.E;  p = si.P;  n = numel(E);
dE = E.' - E;                    % dE(a,b) = E(b)-E(a)
dp = p - p.';                    % dp(a,b) = p(a)-p(b)
inel = abs(dE) > degtol;
M = {si.Mx, si.My, si.Mz};
z = z(:); nz = numel(z);
chi = zeros(3,3,nz);
for mu = 1:3
    for nu = 1:3
        Nmn = M{mu} .* (M{nu}.');            % M_mu(a,b)*M_nu(b,a)
        w   = Nmn .* dp;  w(~inel) = 0;
        dEi = dE(inel);  wi = w(inel);
        for iz = 1:nz
            chi(mu,nu,iz) = sum(wi ./ (dEi - z(iz)));
        end
        if elast
            P2 = repmat(p, 1, n);                % P2(a,b) = p(a)
            el = beta*(sum(Nmn(~inel).*P2(~inel)) - si.Jexp(mu)*si.Jexp(nu));
            idx0 = abs(z) < ztol;
            chi(mu,nu,idx0) = chi(mu,nu,idx0) + el;
        end
    end
end
end
