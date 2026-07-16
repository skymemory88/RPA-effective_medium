%INVZ_RUN_TENSOR_REF Sigma=0 scalar-vs-tensor cross-channel error vs tilt angle.
% Two-layer metric (spec section 7, amended after the Task-9 measurement):
%   eps_spec (invz_chi_tensor_ref): RAW discrepancy incl. the theta-independent
%     baseline from the symmetry-allowed yz cross channel (B64s; measured
%     yz/zz = 0.183 at 6 T even at theta = 0). REPORTED, never gated.
%   eps_tilt (invz_tilt_err): error in the TILT-INDUCED change, differenced
%     against the theta = 0 reference at the same field. Diagnostic only.
% Comparison eta = 0.02 meV (4 pts/HWHM on the 0.005 grid): at the production
% eta = 5e-3 the L2 norm is dominated by sub-linewidth peak misalignment (a
% metric instability, not physics).
% GATE (final, user-approved): peak observables only -- every L2 lineshape
% metric is dominated by the zero-tilt yz peak-offset artifact (delta0/eta;
% verified 0.11 vs 2.2/20 at 6T, ~0.28 vs 7.7/20 at 2T). eps_spec/eps_tilt are
% reported diagnostics; lineshape fidelity is explicitly not claimed.
%   supported(theta > 0) <=> eps_amp <= 0.10 AND dE_peak <= max(0.02*Epeak_ten, eta)
%   at EVERY tested field. Copy the printed table into
%   docs/SESSION-2026-07-16-field-angle.md and the constants of
%   invz/tests/test_invz_tensor_ref.m (reproducibility assertion).
addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
ion = invz_ion();
T       = 0.1;                       % K -- the spectra-driver default temperature
angles  = [0 0.5 1 2 5];             % deg (spec grid; 0 = baseline row, not gated)
fieldsB = [2 4.95 6];                % T: ordered / near-crossover / paramagnetic
w       = (0:0.005:0.6).';           % meV
eta     = 0.02;                      % meV comparison broadening (metric stability)

% live couplings, as in the production drivers (cached lattice sum)
[~, info, Jaa0] = invz_bz_couplings(ion);   % shared BZ-grid coupling branches (Jaa0-aware)
ropts = struct('Jsel', info.Jcc0, 'Jaa0', Jaa0, 'eta', eta);

fprintf('%8s %8s %12s %12s %12s %12s %10s %10s\n', 'theta', '|B| (T)', 'eps_spec', 'eps_tilt', 'eps_amp', 'dE_peak', 'Ep_sc', 'Ep_ten');
supported = true(size(angles));
for ib = 1:numel(fieldsB)
    B = fieldsB(ib);
    R0 = invz_chi_tensor_ref(ion, T, [B 0 0], w, ropts);   % theta = 0 reference
    for ia = 1:numel(angles)
        a = angles(ia);
        if a == 0
            R = R0;  et = 0;
        else
            R = invz_chi_tensor_ref(ion, T, B*[cosd(a) 0 sind(a)], w, ropts);
            et = invz_tilt_err(R, R0);
        end
        ok = a == 0 || (R.eps_amp <= 0.10 && ...
             ( (isnan(R.dE_peak) && isnan(R.Epeak_sc) == isnan(R.Epeak_ten)) || ...
               R.dE_peak <= max(0.02*R.Epeak_ten, eta) ));
        if a > 0, supported(ia) = supported(ia) && ok; end
        if a == 0, vs = 'base'; elseif ok, vs = 'ok'; else, vs = 'FAIL'; end
        fprintf('%8.2f %8.2f %12.4g %12.4g %12.4g %12.4g %10.4g %10.4g   %s\n', ...
                a, B, R.eps_spec, et, R.eps_amp, R.dE_peak, R.Epeak_sc, R.Epeak_ten, vs);
    end
end
fprintf('Supported tilt range (peak-observable criterion): theta_c <= %.2g deg\n', ...
        max([0, angles(supported & angles > 0)]));
