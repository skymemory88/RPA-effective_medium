function tl = invz_twolevel_ordered(ion, T, Bx, hz, opts)
%INVZ_TWOLEVEL_ORDERED Electronic two-level (split doublet) params in the ordered phase.
% Like invz_twolevel, but solves the doublet with a FIXED longitudinal mean field hz (meV) --
% the same hz produced by invz_single_ion(...,'order',true).hz. tl.m is nonzero here (that is
% the point of the ordered phase), so unlike invz_twolevel this does NOT error on m~=0.
% Returns, in addition to the paramagnetic fields:
%   tl.h0 = beta*(1 - n01^2)  the static elastic (Curie) weight h(0); exponentially small at
%                             low T, where the elastic (m^2 h) sector of Sigma is negligible.
% opts.Jxx0 (optional): transverse MF coupling forwarded to invz_single_ion (default ion.Jxx0).
if nargin < 5, opts = struct(); end
Jxx0 = ion.Jxx0;  if isfield(opts,'Jxx0'), Jxx0 = opts.Jxx0; end
C  = invz_const();
si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', false, 'hz_fixed', hz, 'Jxx0', Jxx0));
tl.Delta = si.E(2) - si.E(1);
if tl.Delta < 1e-4
    error('invz:degenerateDoublet', ...
        'Doublet splitting %.2e meV too small for the two-level Sigma.', tl.Delta);
end
tl.M2  = abs(si.Mz(1,2))^2;
tl.m   = real(si.Mz(1,1));
tl.n01 = tanh(tl.Delta/(2*C.kB*T));
tl.g0  = 2*tl.n01/tl.Delta;
tl.h0  = (1 - tl.n01^2)/(C.kB*T);         % beta*(1 - n01^2), the elastic weight h(0)
tl.hz  = hz;
end
