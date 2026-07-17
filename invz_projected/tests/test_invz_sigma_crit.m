function tests = test_invz_sigma_crit
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here,'..','..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_fcc_watson(testCase)
% R 2007 below eq (10): fcc nearest-neighbour value 0.3447 (= Watson integral 1.34466 - 1).
S40 = invz_sigma_crit(12, fcc_jq(40));
S80 = invz_sigma_crit(12, fcc_jq(80));
Sc  = 2*S80 - S40;   % Richardson: raw grid mean converges O(1/n) due to ~1/q^2 integrable singularity at Gamma
verifyEqual(testCase, Sc, 0.3447, 'AbsTol', 0.001);
end

function Jq = fcc_jq(n)
% Nearest-neighbour fcc J(q) on an n^3 grid over the BZ (a=1), Gamma point dropped.
f = (0:n-1)/n;
[F1,F2,F3] = ndgrid(f,f,f);
b = 2*pi*[-1 1 1; 1 -1 1; 1 1 -1];          % fcc primitive reciprocal vectors (a=1)
qx = F1*b(1,1)+F2*b(2,1)+F3*b(3,1);
qy = F1*b(1,2)+F2*b(2,2)+F3*b(3,2);
qz = F1*b(1,3)+F2*b(2,3)+F3*b(3,3);
Jq = 4*(cos(qx/2).*cos(qy/2) + cos(qy/2).*cos(qz/2) + cos(qz/2).*cos(qx/2));
Jq = Jq(:);
Jq = Jq(2:end);                             % drop Gamma (f=0 origin is element 1) so no exclusion warning fires
end

function test_lihof4_sigma_crit(testCase)
% R 2007 eq (10): Sigma_c(0; H=0) = 0.3004 with J12 = -0.1 ueV.  SLOW (dipole sums on 3D grids).
% Raw grid means converge sublinearly (integrable Gamma singularity): 0.2504/0.2639/0.2721/0.2809
% at 8/12/16/24^3, so we Richardson-extrapolate the (12^3, 24^3) pair; controller study gives
% 2*S24 - S12 = 0.2980 vs published 0.3004 (0.8%, different lattice-sum methods).
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
S = zeros(1,2); ns = [12 24];
for k = 1:2
    [T, qvec] = evalc("qVec_generator(ion.a, 'mode', 'grid', 'grid', [ns(k) ns(k) ns(k)], 'range', [-0.5 0.5])"); %#ok<ASGLU> capture fprintf noise
    qvec = qvec(any(abs(qvec) > 1e-12, 2), :);
    [Jnu, info] = invz_jq_modes(ion, qvec, struct('dpRng', 30, 'cache', true));
    S(k) = invz_sigma_crit(info.Jcc0, Jnu(:));
end
Sc = 2*S(2) - S(1);   % Richardson over (12^3, 24^3)
verifyEqual(testCase, Sc, 0.3004, 'AbsTol', 0.006);
verifyEqual(testCase, info.Jcc0, 6.421e-3, 'RelTol', 0.03);
end
