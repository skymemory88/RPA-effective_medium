function tests = test_invzt_emt_matrix
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion engine
end

function test_scalar_limit_parity_with_emt_scalar(testCase)
% Nonsingular spectator channels (v2: avoids the null-channel inverse), zero
% a/b couplings -> K_aa = K_bb = 0 exactly, K_cc == med.K (POSITIVE sign).
% S4-symmetric synthetic cc blocks: circulant eigenvectors (commute with the
% sublattice cyclic shift), so the four sites stay equivalent. v3 (review Other
% 3): KEEP the COMPLEX Hermitian DFT result -- `F` is complex, so `real(F*diag*F')`
% would change the eigenvalues away from `Jnu`. `F*diag(Jnu)*F'` is complex
% Hermitian, circulant, and has eigenvalues EXACTLY `Jnu` (F unitary). The full
% 12x12 page is then complex Hermitian, which the matrix EMT handles natively.
%
% invz_emt_scalar (the 1-channel closed form / parity target) lives in
% invz_projected, which is deliberately OFF the CORE path. Reach the standalone
% (dependency-free) function via a guarded, self-cleaning local addpath so the
% CORE suite stays isolated (setupOnce never adds invz_projected) while this one
% parity test still calls the reference verbatim.
projdir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'invz_projected');
if isempty(which('invz_emt_scalar'))
    addpath(projdir);
    cleaner = onCleanup(@() rmpath(projdir));  %#ok<NASGU> removed on test exit
end
rng(7);
nq = 40;  nz = 5;
F = dftmtx(4)/2;                                    % unitary circulant eigenbasis
Jnu = 6e-3 * (rand(nq, 4) - 0.3);
Jt = zeros(12, 12, nq);
for iq = 1:nq
    Jt(3:3:12, 3:3:12, iq) = F*diag(Jnu(iq,:))*F';  % complex Hermitian, eig == Jnu
end
lat = struct('Jt', Jt, 'qvec', zeros(nq,3), 'w', ones(nq,1)/nq, 'conv', 'explicit', ...
    'JtGamma', zeros(12), 'info', struct());
g0 = -[0.8; 1.1; 1.6; 2.2; 3.0];
Sigma = 0.1*ones(nz, 1);
ctil = zeros(3, 3, nz);
ctil(1,1,:) = 1e-3;  ctil(2,2,:) = 1e-3;            % nonsingular spectators
ctil(3,3,:) = -g0 ./ (1 + Sigma);
[K, chi_bar] = invzt_emt_matrix(ctil, lat, struct());
med = invz_emt_scalar(g0, Sigma, sort(real(eig_all(Jt))), struct());   % v3: eig_all reads ONLY the 4x4 cc blocks (not the full 12x12 -- that would inject 8 spurious zero transverse branches into invz_emt_scalar)
verifyEqual(testCase, squeeze(K(3,3,:)), med.K, 'RelTol', 1e-8);       % POSITIVE sign (v2)
verifyEqual(testCase, squeeze(-chi_bar(3,3,:)), med.G, 'RelTol', 1e-8);
verifyEqual(testCase, max(abs(squeeze(K(1,1,:)))), 0, 'AbsTol', 1e-12);% decoupled spectator
end

function test_direct_closure_and_transpose_symmetry(testCase)
% Defining identities on physical input: chi_lat is the bare RPA of ctil (K
% cancels), chi_imp(K) == chi_bar on the active subspace.
% v3 amend (Task-9 review, verified empirically): chi0(iwn) is NON-Hermitian at
% wn~=0 -- the gyrotropic (~B) part is a PHYSICAL imaginary-antisymmetric term
% (measured rel anti-Hermitian part ~6-17% at Bx=0.5 T). So chi_bar and K are
% non-Hermitian OFF the static slot; they obey the LOCKED transpose relation
% (constraint 9) X(-iwn) = X(iwn).', and are Hermitian ONLY at wn=0. The v2
% asserts (chi_bar == (bz+bz')/2; K == K') contradicted constraint 9 and are
% corrected below. The gyrotropic part is NOT symmetrized away. (Tolerances here
% are magnitude-scaled sketches -- the implementer adapts to the measured norms.)
ion = invz_ion();  T = 1.6;
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
si = invz_single_ion(ion, T, [0.5 0 0], struct('hyp', true));
% One call over +wn AND -wn so the frequency-consistent active-subspace projector
% is shared (the transpose relation is exact only under a shared projector):
zset = 1i*[0.01 0.1 -0.01 -0.1];
ctil = invz_chi0z(si, T, zset, struct('elastic', true));
[K, chi_bar, info] = invzt_emt_matrix(ctil, lat, struct());
for k = 1:2                                              % positive-wn slots
    X = invzt_chi_rpa(ctil(:,:,k), lat.Jt);
    acc = zeros(3, 3, size(X, 3));
    for s = 1:4, acc = acc + X(3*(s-1)+(1:3), 3*(s-1)+(1:3), :) / 4; end
    bz = sum(acc .* reshape(lat.w, 1, 1, []), 3);
    verifyEqual(testCase, chi_bar(:,:,k), bz, 'RelTol', 1e-10);   % RAW average (NOT Hermitized)
    P = info.projector;                                 % [3 x r] orthonormal columns
    R = ctil(:,:,k) \ chi_bar(:,:,k) - K(:,:,k)*chi_bar(:,:,k);
    verifyEqual(testCase, P'*R*P, eye(size(P, 2)), 'AbsTol', 1e-8);
    % transpose relation (constraint 9): slot 2+k is -wn of slot k
    verifyEqual(testCase, K(:,:,2+k),       K(:,:,k).',       'AbsTol', 1e-9);
    verifyEqual(testCase, chi_bar(:,:,2+k), chi_bar(:,:,k).', 'AbsTol', 1e-9);
end
% off the static slot K is genuinely non-Hermitian (asserting Hermiticity here is
% the v2 bug): confirm the anti-Hermitian part is O(1), not O(1e-12):
verifyGreaterThan(testCase, norm(K(:,:,2) - K(:,:,2)'), 1e-3*norm(K(:,:,2)));
% static slot IS Hermitian (real symmetric):
c0 = invz_chi0z(si, T, 0, struct('elastic', true));
[K0, cb0] = invzt_emt_matrix(c0, lat, struct());
verifyLessThan(testCase, norm(K0(:,:,1)  - K0(:,:,1)'),  1e-9*max(norm(K0(:,:,1)), 1e-6));
verifyLessThan(testCase, norm(cb0(:,:,1) - cb0(:,:,1)'), 1e-9*max(norm(cb0(:,:,1)), 1e-6));
end

function test_a2_map_equals_a1_on_diagonal(testCase)
% v3 amend (Task-9 review): the exact identity "matrix EMT == scalar EMT on
% DIAGONAL input" is a property of the K MAP, not the converged Sigma0. At this
% chi0_diag point the T6 self-consistent solve is MULTI-ROOT (both roots are
% valid a1 fixed points; Anderson amplifies a ~1e-14 map difference into a
% Delta~0.09 root split), so a fixed-point Sigma0 equality is fragile. Instead
% verify the a2-converged root is ALSO an a1 fixed point (shared fixed-point set
% == identical map on diagonal input), which is the real content of the identity.
ion = invz_ion();  T = 1.6;  B = [0.5 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
p2 = invzt_solve_point(ion, T, B, lat, struct('mode', 'a2', 'chi0_diag', true, 'odd', false));
% seed a1 at the a2 root: a1 must STAY there (the a2 root is an a1 fixed point):
pchk = invzt_solve_point(ion, T, B, lat, ...
    struct('mode', 'a1', 'chi0_diag', true, 'odd', false, 'Sigma_seed', p2.Sigma));
verifyEqual(testCase, pchk.Sigma0, p2.Sigma0, 'AbsTol', 1e-9);
% physical (non-diagonal) a2 - a1 difference is a REPORT (the matrix-medium content):
p3 = invzt_solve_point(ion, T, B, lat, struct('mode', 'a2', 'odd', true));
p4 = invzt_solve_point(ion, T, B, lat, struct('mode', 'a1', 'odd', true));
fprintf('A2 vs A1 (physical, odd on): dSigma0 = %+.3e, dcrit = %+.3e\n', ...
    p3.Sigma0 - p4.Sigma0, p3.crit - p4.crit);
verifyTrue(testCase, p3.converged && p4.converged);
end

% ------------------------------- local helpers ----------------------------------

function e = eig_all(Jt)
% Flatten the eigenvalues of the 4x4 cc BLOCKS Jt(3:3:12,3:3:12,:) over all pages
% (the 4 cc branches per q). NOT the full 12x12 (which would inject 8 spurious
% zero transverse branches into invz_emt_scalar).
nq = size(Jt, 3);
e = zeros(4*nq, 1);
for iq = 1:nq
    Jcc = Jt(3:3:12, 3:3:12, iq);
    e((iq-1)*4 + (1:4)) = eig((Jcc + Jcc')/2);
end
end
