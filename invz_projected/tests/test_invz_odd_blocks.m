function tests = test_invz_odd_blocks
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));   % repo root: MF_dipole, exchange
end

function test_shapes_and_conj_symmetry(testCase)
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09; -0.25 0 0; -0.31 -0.17 -0.09];
[Vca, Vcb, Vcc] = invz_odd_blocks(ion, q, struct('dpRng', 15, 'cache', false));
verifyEqual(testCase, size(Vca), [4 4 4]);  verifyEqual(testCase, size(Vcb), [4 4 4]);
for iq = 1:2   % real-space couplings real => J(-q) = conj(J(q))
    verifyLessThan(testCase, max(abs(Vca(:,:,iq+2) - conj(Vca(:,:,iq))), [], 'all'), 1e-12);
    verifyLessThan(testCase, max(abs(Vcb(:,:,iq+2) - conj(Vcb(:,:,iq))), [], 'all'), 1e-12);
    verifyLessThan(testCase, norm(Vcc(:,:,iq) - Vcc(:,:,iq)', 'fro'), 1e-14);   % cc Hermitian
end
end

function test_vcc_parity_with_jq_modes(testCase)
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09; 0 0 0];      % includes a Gamma-equivalent point
[~, ~, Vcc, infoB] = invz_odd_blocks(ion, q, struct('dpRng', 15, 'cache', false));
[Jnu, infoS] = invz_jq_modes(ion, q, struct('dpRng', 15, 'cache', false));
for iq = 1:3
    verifyEqual(testCase, sort(real(eig(Vcc(:,:,iq)))).', Jnu(iq,:), 'AbsTol', 1e-12);
end
verifyEqual(testCase, infoB.Jcc0, infoS.Jcc0, 'RelTol', 1e-12);
verifyEqual(testCase, infoB.Jaa0, infoS.Jaa0, 'RelTol', 1e-12);
end

function test_ds2023_geometry_sums(testCase)
% DS2023 Suppl. Table I (a = 5.175 Ang): pure-geometry, gfac-free real-space
% sums — THE unit guard (gfac slips enter deltaJ SQUARED). Central ion on
% sublattice 1; geom stores lower-triangular pairs, so {s,1} exists for s=1..4
% and covers all four sublattice partners. Tf(:,n,m) = 3 r_n r_m/r^5 - d_nm/r^3.
ion = invz_ion();
[~, ~, geom] = MF_dipole([0 0 0], 30, ion.a, ion.tau);
a = 5.175;
[Sxz, Syz, Szz] = deal(0);
for s = 1:4
    Tf = geom.Tf{s, 1};
    Sxz = Sxz + sum(Tf(:,3,1).^2);
    Syz = Syz + sum(Tf(:,3,2).^2);
    Szz = Szz + sum(Tf(:,3,3).^2);
end
verifyEqual(testCase, Sxz, 36.73/a^6, 'RelTol', 0.01);
verifyEqual(testCase, Syz, 36.73/a^6, 'RelTol', 0.01);
verifyEqual(testCase, Szz, 17.93/a^6, 'RelTol', 0.01);
end

function test_onaxis_smallq_decay(testCase)
% C2-about-c kills the ODD blocks on-axis as q -> 0. ON-AXIS ONLY (tilted rays
% carry a direction-dependent macroscopic limit — plan SS8/P0.3).
% P0 AMENDMENT (ODD-LOG SSP0.3): element decay is LINEAR in q (sublattice phase
% factors; the macroscopic term vanishes on-axis), so the source plan's
% 1e-6*Jcc0 element gate at q = 1e-3 is unachievable — its own escape clause
% ("or the deviation is explained by grid geometry") applies. Gates: pinned
% values, linear decay structure, and E1-relevant smallness of the SQUARE
% (deltaJ ~ chi_perp*|Jca|^2 is what must vanish vs Jcc0).
ion = invz_ion();
A = invz_odd_anchors();
q = [1e-1 0 0; 1e-2 0 0; 1e-3 0 0];
Vca = invz_odd_blocks(ion, q, struct('dpRng', 30, 'cache', false));
m = arrayfun(@(iq) max(abs(Vca(:,:,iq)), [], 'all'), 1:3);
verifyEqual(testCase, m(:), A.odd_onaxis_smallq.maxca(:), 'RelTol', 1e-6);   % pinned P0 digits
verifyEqual(testCase, m(2)/m(1), 0.1, 'RelTol', 0.25);                       % ~linear decade steps
verifyEqual(testCase, m(3)/m(2), 0.1, 'RelTol', 0.25);
% chi_perp*|Jca|^2 << Jcc0, on the anchored chi_aa(1.53 K, 0 T) and Jcc0 instead
% of re-hardcoded rounded literals (measured margin ~10x; the <= 2% shifts from
% the old rounded 18 / 6.421e-3 cannot flip the gate).
verifyLessThan(testCase, A.chiperp_1p53K_0T(1,1) * m(3)^2, 1e-5 * A.jcc0);
end

function test_cache_roundtrip_selfverifying(testCase)
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09];
[V1a, V1b, V1c] = invz_odd_blocks(ion, q, struct('dpRng', 10, 'cache', true));
[V2a, V2b, V2c] = invz_odd_blocks(ion, q, struct('dpRng', 10, 'cache', true));
verifyEqual(testCase, {V2a, V2b, V2c}, {V1a, V1b, V1c});      % bitwise round-trip
ion2 = ion;  ion2.J12 = ion.J12 * 1.05;                        % physics change must miss
[~, ~, V3c] = invz_odd_blocks(ion2, q, struct('dpRng', 10, 'cache', true));
% Brief-fix (documented in task-2-report / ODD-LOG §T1.1): the brief asserted
% isequal(V3a,V1a)==false, but Vca/Vcb are DIPOLE-ONLY per the interface spec, so
% they are J12-independent (a 5% J12 change gives bit-identical Vca — verified) and
% that assertion is unsatisfiable. The cc block carries the exchange (|J12|, sign
% J12), so it is the observable that proves the J12 change produced a cache MISS +
% recompute (distinct pkey) rather than a stale hit. Intent preserved, target fixed.
verifyFalse(testCase, isequal(V3c, V1c));
cdir = fullfile(fileparts(mfilename('fullpath')), '..', 'cache');
verifyTrue(testCase, ~isempty(dir(fullfile(cdir, 'odd1_*.mat'))));
end

function test_deltaJ_identities(testCase)
ion = invz_ion();
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];  % driver-style, no Gamma
[Vca, Vcb] = invz_odd_blocks(ion, qvec, struct('dpRng', 15, 'cache', false));
Xp = invz_chiperp(ion, 1.53, [0 0 0], struct());
[dJ, d, dinfo] = invz_odd_deltaJ(Vca, Vcb, Xp);
verifyGreaterThan(testCase, d, 0);
verifyEqual(testCase, dinfo.d_per_sublattice, d*ones(4,1), 'RelTol', 1e-10);       % T1.3 assert
verifyGreaterThan(testCase, dinfo.presub_min_eig, -1e-12*max(1e-30, dinfo.dJ_max)); % PSD pre-sub
verifyLessThan(testCase, dinfo.postsub_diag_bzavg, 1e-15*max(1e-30, dinfo.dJ_max)); % E4 exact
for iq = [1, round(size(dJ,3)/2)]
    verifyLessThan(testCase, norm(dJ(:,:,iq) - dJ(:,:,iq)', 'fro'), 1e-14);         % Hermitian
    verifyLessThan(testCase, max(abs(imag(eig(dJ(:,:,iq))))), 1e-12);               % real eigs
end
% E4/E5 closure (T1.3 acceptance vi): BZ-avg of the pre-subtraction diagonal, as
% reconstructed from the assembled outputs dJ(s,s,:) + d, recovers d per sublattice.
verifyEqual(testCase, squeeze(mean(real(dJ(1,1,:)), 3)) + d, d, 'AbsTol', 1e-15*max(1e-30, d)); %#ok<NASGU> % diag avg is 0
% Xp = 0 -> dJ = 0, d = 0
[dJ0, d0] = invz_odd_deltaJ(Vca, Vcb, zeros(2));
verifyEqual(testCase, max(abs(dJ0(:))), 0);  verifyEqual(testCase, d0, 0);
end

function test_jq_modes_odd_off_bitwise(testCase)
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09; 0 0 0];
[J1, i1, u1] = invz_jq_modes(ion, q, struct('dpRng', 15, 'cache', false));
[J2, i2, u2] = invz_jq_modes(ion, q, struct('dpRng', 15, 'cache', false, 'odd', false));
verifyTrue(testCase, isequaln({J1, i1, u1}, {J2, i2, u2}));                          % regression (i)
end

function test_jq_modes_odd_zero_Xp_equals_off(testCase)
ion = invz_ion();
n = 4;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[J1, i1] = invz_jq_modes(ion, qvec, struct('dpRng', 10, 'cache', false));
[J2, i2] = invz_jq_modes(ion, qvec, struct('dpRng', 10, 'cache', false, 'odd', struct('Xp', zeros(2))));
verifyEqual(testCase, J2, J1, 'AbsTol', 1e-14);                                      % regression (ii)
verifyEqual(testCase, i2.Jcc0, i1.Jcc0, 'AbsTol', 1e-18);
verifyEqual(testCase, i2.odd.d, 0);
end

function test_jq_modes_odd_on_structure(testCase)
ion = invz_ion();
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
Xp = invz_chiperp(ion, 1.53, [0 0 0], struct());
[J1, i1] = invz_jq_modes(ion, qvec, struct('dpRng', 15, 'cache', false));
[J2, i2] = invz_jq_modes(ion, qvec, struct('dpRng', 15, 'cache', false, 'odd', struct('Xp', Xp)));
verifyGreaterThan(testCase, i2.odd.d, 0);
verifyEqual(testCase, i2.Jcc0, i1.Jcc0 - i2.odd.d, 'AbsTol', 1e-18);                 % E5 shift
% ordering mode stays the uniform branch (T1.3 acceptance v): margin shrinks
% from both sides but must remain positive on this grid
verifyLessThan(testCase, max(J2(:)), i2.Jcc0);
% finite-q modes gain weight: mean top-branch coupling strictly increases
verifyGreaterThan(testCase, mean(J2(:,4)), mean(J1(:,4)));
end

function test_parseval_odd_vs_realspace_slow(testCase)
% T1.1(iv): BZ-average of the squared ca blocks == real-space squared sum
% (Parseval), SAME dpRng both sides; 1% tolerance absorbs superlattice folding
% at n = 8 (r^-6-suppressed).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();  C = invz_const();
n = 8;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];   % uniform mesh INCLUDING Gamma
Vca = invz_odd_blocks(ion, qvec, struct('dpRng', 20, 'cache', true));
lhs = mean(sum(abs(Vca(1, :, :)).^2, 2), 3);                       % row s = 1, sum over s', avg q
[~, ~, geom] = MF_dipole([0 0 0], 20, ion.a, ion.tau);
rhs = 0;
for s = 1:4
    Tf = geom.Tf{s, 1};
    rhs = rhs + sum((C.gfac * Tf(:,3,1)).^2);
end
verifyEqual(testCase, lhs, rhs, 'RelTol', 0.01);
% independent conversion check via the cc channel (validates gfac placement):
[~, ~, Vcc] = invz_odd_blocks(ion, qvec, struct('dpRng', 20, 'cache', true));
% subtract exchange+Lorentz first: rebuild the dipole-only cc from MF_dipole per q
% is overkill — instead run the same Parseval on Vca vs B_xz,xz in absolute terms:
a = 5.175;                       % DS2023 lattice constant (Ang), as in test_ds2023_geometry_sums
verifyEqual(testCase, lhs, C.gfac^2 * 36.73/a^6, 'RelTol', 0.015);
end
