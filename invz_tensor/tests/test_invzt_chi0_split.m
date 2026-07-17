function tests = test_invzt_chi0_split
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion engine
end

function test_split_exact_and_groupsize(testCase)
ion = invz_ion();  T = 1.6;
si = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', true));
z = [0, 1i*0.05, 1i*0.5, 0.1 + 1i*5e-3];
full = invz_chi0z(si, T, z, struct('elastic', true));
[cdom, crest, mspec] = invzt_chi0_split(si, T, z, struct());
verifyEqual(testCase, cdom + crest, full, 'AbsTol', 1e-14*max(abs(full(:))));
verifyEqual(testCase, mspec.ndom, 16);
si0 = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', false));
[~, ~, m0] = invzt_chi0_split(si0, T, 0, struct());
verifyEqual(testCase, m0.ndom, 2);
end

function test_dominant_sector_is_longitudinal(testCase)
ion = invz_ion();  T = 1.6;
si = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', true));
[~, crest, mspec] = invzt_chi0_split(si, T, 0, struct());
verifyGreaterThan(testCase, mspec.fdom_cc0, 0.90);
verifyLessThan(testCase, mspec.fdom_perp0, 0.10);
verifyGreaterThan(testCase, real(crest(1,1)), 0);
verifyLessThan(testCase, mspec.elastic_conv_share, 0.02);   % convention residual small
end

function test_esplit_insensitivity_band(testCase)
% Insensitivity band is (0.13, 0.93) meV: above the full ground-doublet
% hyperfine manifold (0 -> 0.130 meV) and below the first excited CF level
% (10.84 K = 0.934 meV). Both probes below sit inside that band (0.1 does
% not -- it splits the ground hyperfine manifold itself, measured
% ndom(0.1)=12 vs ndom(0.7)=16; corrected per controller 2026-07-17).
ion = invz_ion();  T = 1.6;
si = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', true));
c1 = invzt_chi0_split(si, T, 1i*0.05, struct('Esplit', 0.2));
c2 = invzt_chi0_split(si, T, 1i*0.05, struct('Esplit', 0.7));
verifyEqual(testCase, c1, c2, 'AbsTol', 1e-15*max(abs(c1(:))));
end
