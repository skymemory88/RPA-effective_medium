function C = invz_const()
%INVZ_CONST Physical constants for the invz module (meV / K / T units).
C.kB     = 0.08617333;    % meV/K
C.muB    = 0.05788382;    % meV/T
C.Gh2mV  = 4.135667e-3;   % GHz -> meV
gL = 5/4;
C.gfac   = 0.08388;       % mu0/4pi*(gL*muB)^2 [meV*Ang^3]; check: gfac*4/Vc = 1.1654e-3 meV = J_D
end
