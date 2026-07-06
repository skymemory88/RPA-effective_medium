function tl = invz_twolevel(ion, T, Bx)
%INVZ_TWOLEVEL Electronic two-level (split doublet) parameters for the Jensen self-energy.
C  = invz_const();
si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', false));
tl.Delta = si.E(2) - si.E(1);
if tl.Delta < 1e-4
    error('invz:degenerateDoublet', ...
        'Doublet splitting %.2e meV too small: Bx=0 limit needs the closed-form Sigma_c (invz_sigma_crit).', tl.Delta);
end
tl.M2  = abs(si.Mz(1,2))^2;
tl.m   = real(si.Mz(1,1));
if abs(tl.m) > 1e-3
    error('invz:orderedPhase', 'Nonzero diagonal moment m=%.3g: outside paramagnetic-phase scope.', tl.m);
end
tl.n01 = tanh(tl.Delta/(2*C.kB*T));
tl.g0  = 2*tl.n01/tl.Delta;
end
