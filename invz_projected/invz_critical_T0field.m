function Tc = invz_critical_T0field(ion, Sc, J0eff)
%INVZ_CRITICAL_T0FIELD Zero-field Tc from 1 + Sigma_c = J(0)*chi0_cc(0;T).
% Sc is the closed-form critical self-energy (invz_sigma_crit); Sc=0 gives the MF-RPA Tc.
f = @(T) J0eff*static_chi_cc(ion, T) - (1 + Sc);
Tlo = 0.8;  Thi = 3.0;
assert(f(Tlo) > 0 && f(Thi) < 0, 'invz:bracket', 'Tc not bracketed in [0.8, 3.0] K');
tol = 1e-9;                                     % bracket width (K); ~31 bisections vs the old fixed 60
for it = 1:60
    if Thi - Tlo < tol, break; end             % each step is a 136-state single-ion solve; stop early
    Tm = 0.5*(Tlo + Thi);
    if f(Tm) > 0, Tlo = Tm; else, Thi = Tm; end
end
Tc = 0.5*(Tlo + Thi);
end

function c = static_chi_cc(ion, T)
si = invz_single_ion(ion, T, [0 0 0], struct('hyp', true));
cc = real(invz_chi0z(si, T, 0, struct('elastic', true)));
c  = cc(3,3);
end
