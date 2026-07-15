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
% Uniform-mode couplings in meV -- REFERENCE/FALLBACK values only. The AUTHORITATIVE live
% values are the derived info.Jcc0 / info.Jaa0 returned by invz_jq_modes (dpRng-dependent:
% e.g. Jcc0 ~ 6.424e-3 at dpRng 30); every driver and solver passes those via opts.J0eff /
% opts.Jxx0. These pasted constants (from R2007's after-eq-11 D_cc/D_aa(0), validated in
% Task 5) are used only when a caller supplies no derived value -- keep them as documented
% fallbacks, do not treat them as the source of truth. J0eff = J_D*D_cc(0) + 4*J12 =
% 6.821e-3 + 4*(-0.1e-3); Jxx0 = J_D*D_aa(0) + 4*J12 = 3.912e-3 + 4*(-0.1e-3).
ion.J0eff = 6.421e-3;   % meV: fallback J_D*D_cc(0) + 4*J12  (R 2007, after eq 11)
ion.Jxx0  = 3.512e-3;   % meV: fallback J_D*D_aa(0) + 4*J12, transverse MF channel
% Sample-shape (demagnetization) knob. demag = 0 (default): intrinsic couplings, the R2007
% benchmark. demag ~= 0: the ellipsoid shape (ellipsoid_demagn(alpha)) enters ONLY as
%   (a) info.Jshape_cc -- strict-uniform observable correction applied in invz_chi_realaxis
%       (chi_meas = chi/(1 + Jshape_cc*chi)), and
%   (b) demag-aware info.Jaa0 -- the transverse mean-field channel.
% Consequences: info.Jcc0/Jnu and the ordering-channel criticality are demag-INVARIANT
% (R2007: the demagnetizing field cancels from the critical condition; ordering at q -> 0+);
% Tc(B=0) is exactly demag-invariant (<Jx> = 0 there); Bc(T) vs APPLIED field can still
% shift through (b) -- the internal-vs-applied transverse field relation.
ion.demag = 0;          % demag factor (0 = off/intrinsic; 1 = full ellipsoid)
ion.alpha = 1;          % spheroid aspect ratio a/c for ellipsoid_demagn (1 sphere, 0 c-needle, Inf disk)
end
