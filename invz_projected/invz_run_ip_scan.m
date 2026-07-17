%INVZ_RUN_IP_SCAN In-plane rotation diagnostics (IP0) + Sigma=0 tensor reference (IP3).
% Section A -- single-ion angular scans at T = 0.31 K, hyp = false, B = [2 4 6] T,
%   phi_ab = [0:5:90 union 11 union 79] deg, for the three transverse-MF models
%   ('none' = bare CF+Zeeman, 'legacy_x', 'vector_ab'):
%   Delta(phi) = E2-E1 and its C4 harmonic fit (invz_c4fit), span (max-min)/mean,
%   principal-axis angle phi0. The legacy_x rows exist to display the C4 violation,
%   never for production.
% Section B -- Sigma=0 scalar-vs-tensor cc comparison (invz_chi_tensor_ref,
%   transverse_mf = 'vector_ab') at T = 0.1 K, fields [2 4.95 6] T,
%   phi_ab = [0 5 11 15 30 45 60 75 79 90] deg, w = (0:0.005:0.6), eta = 0.02:
%   per row dE_peak, eps_amp, eps_W, Epeak_sc/ten; gate per the tilt criterion
%   (eps_amp <= 0.10 AND dE_peak <= max(0.02*Epeak_ten, eta)).
% Copy both printed tables into docs/SESSION-2026-07-16-inplane-rotation.md.
addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')), '..'));  addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_common'));
ion = invz_ion();

% ---- knobs -----------------------------------------------------------------
% Section A: single-ion angular scans
T_A   = 0.31;  fieldsA = [2 4 6];  phis = unique([0:5:90, 11, 79]);
modes = {'none', 'legacy_x', 'vector_ab'};
% Section B: Sigma=0 scalar-vs-tensor comparison
T_B = 0.1;  fieldsB = [2 4.95 6];  phisB = [0 5 11 15 30 45 60 75 79 90];
w = (0:0.005:0.6).';  eta = 0.02;
% -----------------------------------------------------------------------------

% ---- Section A: single-ion angular scans ----
fprintf('%10s %6s %10s %10s %10s %8s\n', 'model', 'B(T)', 'span(%)', 'A4/A0(%)', 'A8/A0(%)', 'phi0');
for im = 1:numel(modes)
    for ib = 1:numel(fieldsA)
        B = fieldsA(ib);  d = zeros(numel(phis), 1);
        for k = 1:numel(phis)
            s = invz_single_ion(ion, T_A, B*[cosd(phis(k)) sind(phis(k)) 0], ...
                                struct('hyp', false, 'transverse_mf', modes{im}));
            d(k) = s.E(2) - s.E(1);
        end
        [A, phi0, ~] = invz_c4fit(phis, d);
        fprintf('%10s %6.1f %10.3f %10.3f %10.3f %8.2f\n', modes{im}, B, ...
            100*(max(d)-min(d))/mean(d), 100*hypot(A(2),A(3))/A(1), 100*abs(A(4))/A(1), phi0);
    end
end

% ---- Section B: Sigma=0 scalar-vs-tensor over angle (production couplings) ----
[~, info, Jaa0] = invz_bz_couplings(ion);   % shared BZ-grid coupling branches (Jaa0-aware)
ropts = struct('Jsel', info.Jcc0, 'Jaa0', Jaa0, 'eta', eta, 'transverse_mf', 'vector_ab');
fprintf('\n%8s %8s %12s %12s %12s %10s %10s %5s\n', 'phi', '|B| (T)', 'dE_peak', 'eps_amp', 'eps_W', 'Ep_sc', 'Ep_ten', 'ok');
supported = true(size(phisB));
for ib = 1:numel(fieldsB)
    for ia = 1:numel(phisB)
        ph = phisB(ia);
        R = invz_chi_tensor_ref(ion, T_B, fieldsB(ib)*[cosd(ph) sind(ph) 0], w, ropts);
        ok = R.eps_amp <= 0.10 && ...
             ( (isnan(R.dE_peak) && isnan(R.Epeak_sc) == isnan(R.Epeak_ten)) || ...
               R.dE_peak <= max(0.02*R.Epeak_ten, eta) );
        supported(ia) = supported(ia) && ok;
        fprintf('%8.1f %8.2f %12.4g %12.4g %12.4g %10.4f %10.4f %5d\n', ...
            ph, fieldsB(ib), R.dE_peak, R.eps_amp, R.eps_W, R.Epeak_sc, R.Epeak_ten, ok);
    end
end
fprintf('\nsupported in-plane angles (all fields): %s deg\n', mat2str(phisB(supported)));
