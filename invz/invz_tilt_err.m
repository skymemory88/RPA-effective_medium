function eps = invz_tilt_err(R, R0)
%INVZ_TILT_ERR Tilt-referenced spectral error of the scalar cc chain.
% Measures how well the scalar pipeline captures the TILT-INDUCED change (the
% spec's accuracy statement), by differencing against the theta = 0 reference
% R0 at the same field -- this removes the theta-independent baseline
% discrepancy from the symmetry-allowed yz cross channel (B64s), which exists
% at zero tilt and is not a tilt error.
Dsc  = R.chi_sc  - R0.chi_sc;
Dten = R.chi_ten - R0.chi_ten;
floorv = 1e-12 * max(abs(R.chi_ten)) * sqrt(numel(R.chi_ten));
eps = norm(Dsc - Dten, 2) / max(norm(Dten, 2), floorv);
end
