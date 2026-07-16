function R = invz_chi_tensor_ref(ion, T, Bvec, w, opts)
%INVZ_CHI_TENSOR_REF Sigma = 0 scalar-vs-tensor RPA cross-check at one (T, B).
% Quantifies the cross-channel (xz/yz/xy) error the scalar cc pipeline omits under
% a c-axis tilt (spec 2026-07-16 SR3). Both sides share ONE single-ion state
% (order mode, sign-aware seed) and ONE bare chi0 tensor from invz_chi0z; they
% differ ONLY in propagation:
%   scalar:  chi_cc = chi0_cc / (1 - Jsel*chi0_cc)              (cc channel only)
%   tensor:  chi    = chi0 * inv(I3 - J*chi0),  J = diag(Jaa0, Jaa0, Jsel)
% NOTE: in the SPONTANEOUSLY ordered phase (m0 ~= 0 at theta = 0) the cross
% channels are nonzero even without tilt -- a pre-existing property of the scalar
% pipeline that this reference makes visible; expect eps_spec(theta=0) > 0 there.
% Intrinsic response only: asserts ion.demag == 0.
%
% R fields: chi_sc, chi_ten [nw x 1] (chi''_cc), eps_spec (floored spectral L2
% relative error), Epeak_sc/Epeak_ten/dE_peak (invz_peak_energy, opts.peak_wmin
% default 0.05 meV to skip the hyperfine line), amp_sc/amp_ten/eps_amp (peak
% amplitude error over the SAME wmin-masked w range as the peak search, so the
% gate reads the electronic mode, not a sub-wmin hyperfine feature -- the
% GATED metric; eps_spec is a reported diagnostic only, since it is dominated
% by the zero-tilt yz peak-offset artifact for sharp lines, not by lineshape
% fidelity), w, B.
if nargin < 5, opts = struct(); end
eta   = getf(opts, 'eta', 5e-3);
Jsel  = getf(opts, 'Jsel', ion.J0eff);
Jaa0  = getf(opts, 'Jaa0', ion.Jxx0);
hyp   = getf(opts, 'hyp', true);
wmin  = getf(opts, 'peak_wmin', 0.05);
if ion.demag ~= 0
    error('invz:tensorRef', 'reference defined for the intrinsic response (ion.demag = 0).');
end
B  = invz_field_vec(Bvec);
si = invz_single_ion(ion, T, B, struct('hyp', hyp, 'order', true, 'J0z', Jsel, 'Jxx0', Jaa0));
w  = w(:);
z  = w + 1i*eta;
c0 = invz_chi0z(si, T, z, struct('elastic', false));
J  = diag([Jaa0, Jaa0, Jsel]);
nw = numel(z);
chi_sc = zeros(nw, 1);  chi_ten = zeros(nw, 1);
for k = 1:nw
    X = c0(:, :, k);
    chi_sc(k)  = imag(X(3, 3) / (1 - Jsel*X(3, 3)));
    Xt = X / (eye(3) - J*X);                     % chi0 * inv(I - J*chi0), full 3x3
    chi_ten(k) = imag(Xt(3, 3));
end
R.chi_sc = chi_sc;  R.chi_ten = chi_ten;  R.w = w;  R.B = B;
floorv = 1e-12 * max(abs(chi_ten)) * sqrt(nw);   % guards the metric at spectral zeros
R.eps_spec  = norm(chi_sc - chi_ten, 2) / max(norm(chi_ten, 2), floorv);
R.Epeak_sc  = invz_peak_energy(chi_sc,  w, wmin);
R.Epeak_ten = invz_peak_energy(chi_ten, w, wmin);
R.dE_peak   = abs(R.Epeak_sc - R.Epeak_ten);
% Peak-observable amplitude error (the GATED intensity metric; L2 lineshape
% metrics are positional-artifact-dominated for sharp lines -- see spec sec. 7):
% Same wmin mask as the peak search: the gate must measure the ELECTRONIC mode,
% not a sub-wmin hyperfine feature.
msk = w >= wmin;
R.amp_sc  = max(chi_sc(msk));
R.amp_ten = max(chi_ten(msk));
R.eps_amp = abs(R.amp_sc - R.amp_ten) / max(R.amp_ten, floorv);
end
