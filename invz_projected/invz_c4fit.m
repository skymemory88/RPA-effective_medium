function [A, phi0, resid] = invz_c4fit(phi_deg, y)
%INVZ_C4FIT Least-squares C4 harmonic fit y = A1 + A2*cos4phi + A3*sin4phi + A4*cos8phi.
% phi0 (deg) = principal-axis angle of the leading harmonic, atan2d(A3,A2)/4.
% resid = max |fit - y| (worst-case pointwise residual, not RMS).
phi_deg = phi_deg(:);  y = y(:);
M = [ones(size(phi_deg)), cosd(4*phi_deg), sind(4*phi_deg), cosd(8*phi_deg)];
A = M \ y;
phi0 = atan2d(A(3), A(2))/4;
resid = max(abs(M*A - y));
end
