function tests = test_invzt_ordered_vertex_basis
%TEST_INVZT_ORDERED_VERTEX_BASIS  Task 7B: FIXED-COUNT, FIELD-ADAPTED ordered vertex basis.
%   The a3d vertex basis is the 16 LOWEST eigenstates of the ORDERED mean-field Hamiltonian
%   at each point -- a direct truncation of si_full (fixed COUNT, adapted CONTENT). This
%   replaces the falsified field-INDEPENDENT zero-field e2xI8 projection (review P0-3), which
%   captured ~0 of the ordered static cc response (see test_falsification_record below).
%
%   (1) B = 3.0 T gates: n_vertex == 16, chi_share > 0.9, gap_16_17 finite positive.
%   (2) field-continuity over B = [3.0 3.5 4.0 4.4 4.65]: adjacent |Δchi_share| < 0.05,
%       adjacent-field subspace overlap min(svd(V16_k' V16_{k+1})) > 0.9, all E finite.
%   (3) falsification record (executable): the OLD zero-field e2xI8 projection chi_share
%       < 0.01 at B = 3.0 T -- the recorded reason field-adaptation is REQUIRED.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
ion = invz_ion();
g   = invzt_qgrid(16, 'halfopen');
% CACHED lattice (loads from disk in seconds) -- only Jcc0/Jaa0 are needed for the
% ordered single-ion mean field; no per-point lattice solve in these tests.
lat = invzt_jq_tensor(ion, g, struct('dpRng', 30, 'cache', true));
tc.TestData.ion  = ion;
tc.TestData.J0z  = lat.info.Jcc0;
tc.TestData.Jxx0 = lat.info.Jaa0;
tc.TestData.T    = 0.1;
end

function si = ordered_si(tc, Bx)
% Ordered (spontaneous-moment) FULL electronuclear single ion at (T, Bx) -- the SAME
% si_full invzt_solve_point_ordered builds (lines 108-113): hyp + order + J0z/Jxx0 threaded.
siopts = struct('hyp', true, 'order', true, 'J0z', tc.TestData.J0z, ...
    'Jxx0', tc.TestData.Jxx0, 'transverse_mf', 'legacy_x');
si = invz_single_ion(tc.TestData.ion, tc.TestData.T, [Bx 0 0], siopts);
end

function test_vertex_gates_at_3T(tc)
% (1) DIRECTION GATES at the deep ordered anchor B = 3.0 T (evidence, not knobs).
si = ordered_si(tc, 3.0);
vb = invzt_ordered_vertex_basis(tc.TestData.ion, tc.TestData.T, si, struct());
verifyEqual(tc, vb.n_full, numel(si.E));
verifyEqual(tc, vb.n_vertex, 16);
verifyEqual(tc, size(vb.Mz), [16 16]);
verifyEqual(tc, numel(vb.E), 16);
verifyEqual(tc, numel(vb.p), 16);
verifyEqual(tc, size(vb.V16), [numel(si.E) 16]);
verifyTrue(tc, all(isfinite(vb.E)));
verifyLessThan(tc, abs(sum(vb.p) - 1), 1e-12);
verifyTrue(tc, isfinite(vb.gap_16_17) && vb.gap_16_17 > 0);
fprintf('7B @3T: chi_share = %.5f (> 0.9), gap_16_17 = %.6f meV\n', ...
    vb.chi_share, vb.gap_16_17);
verifyGreaterThan(tc, vb.chi_share, 0.9);          % DIRECTION GATE (do NOT loosen)
end

function test_field_continuity(tc)
% (2) field-continuity: no cutoff-induced jumps + smooth subspace tracking to the QPT.
Bs = [3.0 3.5 4.0 4.4 4.65];
cs = zeros(size(Bs));
Vs = cell(size(Bs));
for k = 1:numel(Bs)
    si = ordered_si(tc, Bs(k));
    vb = invzt_ordered_vertex_basis(tc.TestData.ion, tc.TestData.T, si, struct());
    verifyEqual(tc, vb.n_vertex, 16);
    verifyTrue(tc, all(isfinite(vb.E)));
    cs(k) = vb.chi_share;
    Vs{k} = vb.V16;
    fprintf('B=%.2f T: chi_share=%.5f gap_16_17=%.6f\n', Bs(k), vb.chi_share, vb.gap_16_17);
end
ov = zeros(1, numel(Bs)-1);
for k = 1:numel(Bs)-1
    ov(k) = min(svd(Vs{k}' * Vs{k+1}));            % smallest principal-angle cosine
    fprintf('subspace overlap B=%.2f->%.2f : min svd = %.5f\n', Bs(k), Bs(k+1), ov(k));
end
fprintf('max adjacent |Δchi_share| = %.5f (< 0.05); min subspace overlap = %.5f (> 0.9)\n', ...
    max(abs(diff(cs))), min(ov));
verifyLessThan(tc, max(abs(diff(cs))), 0.05);
verifyGreaterThan(tc, min(ov), 0.9);
end

function test_falsification_record_zerofield_e2xI8(tc)
% (3) EXECUTABLE falsification of review P0-3's field-INDEPENDENT recommendation: the
% zero-field e2xI8 (ground doublet (x) I8) projection captures ~nothing of the ordered
% static cc response, because the Ising doublet has <+|Jz|-> = 0 and the ordering field hz
% fully polarizes it (reduced Jz-fluctuation -> 0). This is WHY 7B is field-adapted.
si = ordered_si(tc, 3.0);
s_old = old_zerofield_e2xI8_chishare(tc.TestData.ion, tc.TestData.T, si);
fprintf('falsification record: zero-field e2xI8 chi_share @3T = %.3e (< 0.01)\n', s_old);
verifyLessThan(tc, s_old, 0.01);
end

% ------------------------------------------------------------------------------------- %
function s = old_zerofield_e2xI8_chishare(ion, T, si_full)
%OLD_ZEROFIELD_E2XI8_CHISHARE  Recorded (dead) projection: chi0_cc of the fixed zero-field
% ground-doublet (x) I8 subspace / full. Kept ONLY as an executable falsification record;
% the shipping constructor no longer uses this projector machinery.
rb = invzt_rung_basis(ion, 'e2xI8');  P = rb.projector;
C  = invz_const();  beta = 1/(C.kB*T);
oJ = stevens_ops(ion.J);  oI = stevens_ops(ion.I);  nI = size(oI.Jz, 1);
Jx = kron(oJ.Jx, eye(nI));  Jy = kron(oJ.Jy, eye(nI));  Jz = kron(oJ.Jz, eye(nI));
Hfull = si_full.V * diag(si_full.E + si_full.E0) * si_full.V';  Hfull = (Hfull + Hfull')/2;
Hred  = P'*Hfull*P;  Hred = (Hred + Hred')/2;
[Vr, Dr] = eig(Hred, 'vector');  [E, ix] = sort(real(Dr));  Vr = Vr(:, ix);  E = E - E(1);
p = exp(-beta*E);  p = p/sum(p);
sr.E = E;  sr.P = p;
sr.Mx = Vr'*(P'*Jx*P)*Vr;  sr.My = Vr'*(P'*Jy*P)*Vr;  sr.Mz = Vr'*(P'*Jz*P)*Vr;
sr.Jexp = [real(diag(sr.Mx)).'*p; real(diag(sr.My)).'*p; real(diag(sr.Mz)).'*p];
cv = invz_chi0z(sr,      T, 0, struct('elastic', true));
cf = invz_chi0z(si_full, T, 0, struct('elastic', true));
s  = real(cv(3,3,1)) / real(cf(3,3,1));
end
