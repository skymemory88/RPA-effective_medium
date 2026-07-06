function ion = invz_ion()
%INVZ_ION LiHoF4 parameters (Rønnow et al., PRB 75, 054426 (2007), Table I last row).
ion.J   = 8;    ion.gL = 5/4;   ion.I = 3.5;
ion.A   = 3.36e-3;                                 % meV
ion.B20 = -0.06;      ion.B40 = 0.35e-3;  ion.B44 = 3.6e-3;    % meV
ion.B60 = 4.0e-7;     ion.B64c = 7.0e-5;  ion.B64s = 0.98e-5;  % meV
ion.a   = [5.175 0 0; 0 5.175 0; 0 0 10.75];       % Å (rows = lattice vectors)
ion.tau = [0 0 0; 0 0.5 0.25; 0.5 0.5 0.5; 0.5 0 0.75];  % Ho sites, fractional
ion.Vc  = abs(det(ion.a));                          % 287.9 Å^3
ion.J12 = -0.1e-3;                                  % meV, nn exchange (4 neighbours)
% Uniform-mode couplings in meV (validated in Task 5): J0eff = J_D*D_cc(0) + 4*J12 = 6.821e-3 + 4*(-0.1e-3) = 6.421e-3 meV; Jxx0 = J_D*D_aa(0) + 4*J12 = 3.912e-3 + 4*(-0.1e-3) = 3.512e-3 meV.
ion.J0eff = 6.421e-3;   % meV: J_D*D_cc(0) + 4*J12  (R 2007, after eq 11)
ion.Jxx0  = 3.512e-3;   % meV: J_D*D_aa(0) + 4*J12, transverse MF channel
end
