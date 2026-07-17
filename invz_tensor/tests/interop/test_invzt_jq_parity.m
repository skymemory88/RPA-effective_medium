function tests = test_invzt_jq_parity
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

function test_block_parity_with_projected_odd_blocks(testCase)
% ca/cb/cc block equality with invz_odd_blocks — inherits its DS2023 geometry
% guards transitively. Exact same assembly path expected.
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09; 0 0 0];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
[Vca, Vcb, Vcc, infoS] = invz_odd_blocks(ion, q, struct('dpRng', 15, 'cache', false));
for iq = 1:3
    verifyEqual(testCase, lat.Jt(3:3:12, 1:3:12, iq), Vca(:,:,iq), 'AbsTol', 1e-15);
    verifyEqual(testCase, lat.Jt(3:3:12, 2:3:12, iq), Vcb(:,:,iq), 'AbsTol', 1e-15);
    verifyEqual(testCase, lat.Jt(3:3:12, 3:3:12, iq), Vcc(:,:,iq), 'AbsTol', 1e-14);
end
verifyEqual(testCase, lat.info.Jcc0, infoS.Jcc0, 'RelTol', 1e-12);
verifyEqual(testCase, lat.info.Jaa0, infoS.Jaa0, 'RelTol', 1e-12);
end

function test_cc_eigen_parity_with_jq_modes(testCase)
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
[Jnu, ~] = invz_jq_modes(ion, q, struct('dpRng', 15, 'cache', false));
for iq = 1:2
    Jcc = (lat.Jt(3:3:12, 3:3:12, iq) + lat.Jt(3:3:12, 3:3:12, iq)')/2;
    verifyEqual(testCase, sort(real(eig(Jcc))).', Jnu(iq,:), 'AbsTol', 1e-12);
end
end
