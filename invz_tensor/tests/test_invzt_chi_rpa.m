function tests = test_invzt_chi_rpa
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion engine
end

function test_decoupled_cc_closes_to_scalar_branches(testCase)
% Cartesian-diagonal chi0 + Cartesian-diagonal-only Jt: the cc sector must equal
% scalar RPA over the cc-block eigenvalues COMPUTED LOCALLY from Jt (no
% projected dependency).
ion = invz_ion();  T = 1.6;
q = [0.25 0 0; 0.31 0.17 0.09];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
si = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', true));
c0 = invz_chi0z(si, T, 1i*0.05, struct('elastic', true));
c0d = diag(diag(c0(:,:,1)));
Jz = zero_offdiag_blocks(lat.Jt);
X = invzt_chi_rpa(c0d, Jz);
x0cc = real(c0d(3,3));
for iq = 1:2
    Jcc = (lat.Jt(3:3:12, 3:3:12, iq) + lat.Jt(3:3:12, 3:3:12, iq)')/2;
    Jnu = sort(real(eig(Jcc)));
    Xcc = X(3:3:12, 3:3:12, iq);
    got = sort(real(eig((Xcc + Xcc')/2)));
    verifyEqual(testCase, got, x0cc ./ (1 - Jnu * x0cc), 'RelTol', 1e-10);
end
end

function test_gcc_lattice_consistency(testCase)
% Brute-force equality of the weighted average (exact identity). S4 sublattice
% EQUALITY is NOT asserted here — 3 arbitrary q points are not a symmetry-
% complete set; the S4 check lives in the solver tests on full grids (report).
ion = invz_ion();  T = 1.6;
q = [0.25 0 0; 0.31 0.17 0.09; 0.1 0.2 0.3];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 10, 'cache', false));
si = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', true));
zn = 1i*[0.01 0.05 0.2];
c0 = invz_chi0z(si, T, zn, struct('elastic', true));
[Gcc, diag4] = invzt_gcc_lattice(c0, lat);
for k = 1:3
    X = invzt_chi_rpa(c0(:,:,k), lat.Jt);
    brute = zeros(4,1);
    for s = 1:4
        brute(s) = sum(lat.w(:).' .* squeeze(real(X(3*(s-1)+3, 3*(s-1)+3, :))).');
    end
    verifyEqual(testCase, diag4(:,k), brute, 'RelTol', 1e-12);
    verifyEqual(testCase, Gcc(k), mean(brute), 'RelTol', 1e-12);
end
end

function test_schur_complement_equals_E1_direct(testCase)
% THE convention gate: transverse-transverse blocks zeroed + synthetic diagonal
% chi0 -> Schur elimination of (a,b) equals scalar RPA with Jcc + dJpre, where
% dJpre is built DIRECTLY (xperp*(Vca*Vca' + Vcb*Vcb')). No invz_odd_deltaJ
% call: its caller contract requires a full uniform BZ mesh, which one generic
% q is not (v2 review finding 6a).
ion = invz_ion();
q = [0.31 0.17 0.09];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
xperp = 11.05;  xcc = 40.0;
Jz = lat.Jt;
for s = 1:4, for sp = 1:4
    Jz(3*(s-1)+(1:2), 3*(sp-1)+(1:2), 1) = 0;
end, end
X0 = kron(eye(4), diag([xperp, xperp, xcc]));
X = (eye(12) - X0*Jz(:,:,1)) \ X0;
Xcc = X(3:3:12, 3:3:12);
Vca = Jz(3:3:12, 1:3:12, 1);  Vcb = Jz(3:3:12, 2:3:12, 1);
Jcc = Jz(3:3:12, 3:3:12, 1);
dJpre = xperp*(Vca*Vca') + xperp*(Vcb*Vcb');
Xcc_e1 = xcc * ((eye(4) - (Jcc + dJpre)*xcc) \ eye(4));
verifyEqual(testCase, Xcc, Xcc_e1, 'RelTol', 1e-10);
end

function test_full_schur_enhancement_reported(testCase)
% Exact Schur with real transverse-transverse blocks kept; E1 is the
% Jperp,perp -> 0 limit. REPORT the ratio (finite, < 1).
ion = invz_ion();
q = [0.31 0.17 0.09];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
Xp = invz_chiperp(ion, 1.53, [0 0 0], struct());        % invz_common function
P = zeros(8);  Vc = zeros(4, 8);
pidx = @(s) 2*(s-1)+(1:2);
for s = 1:4, for sp = 1:4
    P(pidx(s), pidx(sp)) = lat.Jt(3*(s-1)+(1:2), 3*(sp-1)+(1:2), 1);
    Vc(s, pidx(sp)) = lat.Jt(3*(s-1)+3, 3*(sp-1)+(1:2), 1);
end, end
Jcc = lat.Jt(3:3:12, 3:3:12, 1);
S = Jcc + Vc * ((kron(eye(4), inv(Xp)) - P) \ Vc');
dJpre = Vc * kron(eye(4), Xp) * Vc';
r = norm(S - (Jcc + dJpre), 'fro') / max(norm(dJpre, 'fro'), 1e-30);
fprintf('full-Schur vs E1 correction ratio: %.3g (expect ~ Xp*Jaa0 ~ 0.05 scale)\n', r);
verifyTrue(testCase, isfinite(r) && r < 1);
end

% ------- local helper (not a test: no testCase arg, name has no 'test' prefix) -------
function Jz = zero_offdiag_blocks(Jt)
%ZERO_OFFDIAG_BLOCKS Zero every Cartesian-off-diagonal (mu ~= nu) block of a
% [12,12,nq] tensor, for ALL sublattice pairs, keeping only the mu==nu (aa/bb/cc)
% blocks. Leaves a Cartesian-diagonal-only coupling: the cc channel then
% decouples exactly from (a,b), so RPA over Jz reduces to scalar RPA over the
% Jcc eigenvalues alone (the identity test_decoupled_cc_closes_to_scalar_branches
% exercises).
Jz = Jt;
for mu = 1:3
    for nu = 1:3
        if mu ~= nu
            Jz(mu:3:12, nu:3:12, :) = 0;
        end
    end
end
end
