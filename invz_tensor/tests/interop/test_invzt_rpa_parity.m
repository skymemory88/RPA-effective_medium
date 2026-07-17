function tests = test_invzt_rpa_parity
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));                                 % invz_tensor
addpath(fullfile(here, '..', '..', '..'));                           % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', '..', 'invz_common'));            % shared single-ion engine
addpath(fullfile(here, '..', '..', '..', 'invz_projected'));         % parity targets
addpath(fullfile(here, '..', '..', '..', 'invz_projected', 'tests'));% invz_odd_anchors fixture
end

function test_uniform_mode_rpa_parity_with_tensor_ref(testCase)
% Sigma = 0 uniform-mode parity vs invz_chi_tensor_ref (invz_projected). The
% reference's "tensor" side is a SINGLE-SITE 3x3 RPA closure Xt = chi0*inv(I3 -
% J*chi0) with J = diag(Jaa0, Jaa0, Jsel). The 12x12 lattice analogue: Jd puts
% the SAME per-Cartesian-diagonal value on EVERY one of the 16 sublattice pairs,
% scaled by 1/4 so the S4-uniform sublattice mode (eigenvalue 1 of ones(4)/4)
% picks up exactly Jaa0/Jsel -- matching the lat.info.Jaa0/Jcc0 uniform-mode
% convention (Global Constraints). u = kron(ones(4,1)/2, eye(3)) is the
% isometric embedding of that normalized uniform mode (u'*u = eye(3)); Xu =
% u'*X*u is the tensor RPA response projected onto it. Algebraically (resolvent
% identity (I-AB)^-1 A = A(I-BA)^-1 plus the rank-1 sublattice projector ones(4)/4
% = v*v', v = ones(4,1)/2), Xu reduces EXACTLY to Xt for ANY chi0 -- so this is a
% linear-algebra identity, not a physics claim; RelTol 1e-8 is a generous margin
% over the ~1e-8..1e-12 numerical floor of two independently-coded 12x12 vs 3x3
% linear solves.
%
% READ invz_chi_tensor_ref.m FIRST (done): its si is built with
% struct('hyp',hyp,'order',true,'J0z',Jsel,'Jxx0',Jaa0,'transverse_mf',tmf),
% Jsel = getf(opts,'Jsel',ion.J0eff), Jaa0 = getf(opts,'Jaa0',ion.Jxx0), eta =
% getf(opts,'eta',5e-3), z = w + 1i*eta, c0 = invz_chi0z(si,T,z,struct('elastic',
% false)). The local wrapper below mirrors that EXACTLY (same getf defaults,
% same field/transverse-MF normalization calls) and self-checks against
% R.chi_sc before the tensor-side comparison.
ion = invz_ion();
T = 1.6;  Bvec = [0.1 0 0];
w = (0:0.05:0.4).';
opts = struct('hyp', true);                     % Jsel/Jaa0 left at defaults: ion.J0eff/ion.Jxx0

R = invz_chi_tensor_ref(ion, T, Bvec, w, opts);

% --- local wrapper mirroring invz_chi_tensor_ref's si/c0 construction exactly ---
eta  = getf(opts, 'eta', 5e-3);
Jsel = getf(opts, 'Jsel', ion.J0eff);
Jaa0 = getf(opts, 'Jaa0', ion.Jxx0);
hyp  = getf(opts, 'hyp', true);
B    = invz_field_vec(Bvec);
tmf  = invz_check_transverse_mf(opts, B(2));
si = invz_single_ion(ion, T, B, struct('hyp', hyp, 'order', true, 'J0z', Jsel, ...
    'Jxx0', Jaa0, 'transverse_mf', tmf));
wc = w(:);  z = wc + 1i*eta;
c0 = invz_chi0z(si, T, z, struct('elastic', false));      % [3,3,nw]
nw = numel(z);

% Self-check: our local si/c0 reproduces R.chi_sc (same scalar-branch formula)
% BEFORE trusting the tensor-side comparison below -- confirms the local wrapper
% has not drifted from invz_chi_tensor_ref's internal si construction.
chi_sc_local = zeros(nw, 1);
for k = 1:nw
    Xk = c0(:,:,k);
    chi_sc_local(k) = imag(Xk(3,3) / (1 - Jsel*Xk(3,3)));
end
verifyEqual(testCase, chi_sc_local, R.chi_sc, 'RelTol', 1e-12);

% --- 12x12 uniform-mode tensor RPA (invzt_chi_rpa), projected with u ---------
Jd = kron(ones(4)/4, diag([ion.Jxx0, ion.Jxx0, ion.J0eff]));
u  = kron(ones(4,1)/2, eye(3));
chi_ten_lat = zeros(nw, 1);
for k = 1:nw
    X  = invzt_chi_rpa(c0(:,:,k), Jd);
    Xu = u' * X * u;
    chi_ten_lat(k) = imag(Xu(3,3));
end
verifyEqual(testCase, chi_ten_lat, R.chi_ten, 'RelTol', 1e-8);
end
