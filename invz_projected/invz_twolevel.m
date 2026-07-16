function tl = invz_twolevel(ion, T, Bx, opts)
%INVZ_TWOLEVEL Electronic two-level (split doublet) parameters for the Jensen self-energy.
% Bx: scalar (transverse, historical) or [Bx By Bz] vector (T).
% opts.Jxx0 (optional): transverse MF coupling forwarded to invz_single_ion (default ion.Jxx0).
if nargin < 4, opts = struct(); end
Jxx0 = ion.Jxx0;  if isfield(opts,'Jxx0'), Jxx0 = opts.Jxx0; end
C  = invz_const();
si = invz_single_ion(ion, T, invz_field_vec(Bx), struct('hyp', false, 'Jxx0', Jxx0));
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
