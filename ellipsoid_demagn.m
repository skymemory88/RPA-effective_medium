function N = ellipsoid_demagn(alpha)
% Demagnetization tensor for a rotational ellipsoid (spheroid) with a = b.
%
% Definition:
%   alpha = a/c  (equatorial-to-polar axis ratio)
%     alpha = 0   → needle (c >> a)
%     alpha = 1   → sphere
%     alpha = Inf → disk   (a >> c)
%
% References:
%   A. Hubert and R. Schäfer, Magnetic Domains, p. 121
%   (standard prolate/oblate spheroid demagnetization factors)
if alpha == 1 % sphere
    Nz = 1/3; Nx = (1-Nz)/2; Ny = Nx;
elseif alpha == 0 % needle (prolate limit)
    Nz = 0.0;  Nx = 0.5;     Ny = 0.5;
elseif isinf(alpha) % disk (oblate limit)
    Nz = 1.0;  Nx = 0.0;     Ny = 0.0;
elseif alpha < 1 % prolate ("cigar")
    Nz = (alpha^2/(1-alpha^2)) * ((1/sqrt(1-alpha^2))*asinh(sqrt(1-alpha^2)/alpha) - 1);
    Nx = 0.5*(1-Nz); Ny = Nx;
else % oblate ("disk/smarties")
    Nz = (alpha^2/(alpha^2-1)) * (1 - (1/sqrt(alpha^2-1))*asin(sqrt(alpha^2-1)/alpha));
    Nx = 0.5*(1-Nz); Ny = Nx;
end

N = [Nx  0   0
      0  Ny  0
      0  0  Nz];
end
