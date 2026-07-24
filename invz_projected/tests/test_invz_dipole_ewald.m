function tests = test_invz_dipole_ewald
% Gate-A Task-1 (strict contract / geometry / exact candidate union / preflight)
% and Task-2 (real+reciprocal+self assembly + extended-zone gauge restore) tests
% for the opt-in Ewald dipolar primitive invz_dipole_ewald.m.
%
% Structural/exact-identity checks use M_id (max component |A-B| <= 1e-12*T_scale,
% T_scale = max component of max(|A|,|B|)) per docs/invzp_ewald_prereg.md sec 3.
% This file inlines the helpers it needs (Task-1/2 must not depend on the later
% invz_ewald_fixtures / invz_ewald_metrics deliverables). Every reconstruction
% builds all ordered (n,m) blocks independently -- no conjugate block is ever
% filled by assignment -- so Hermiticity/reconstruction stay non-tautological.
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));            % invz_projected
addpath(fullfile(here, '..', '..'));      % repo root: invz_dipole_ewald, MF_dipole
addpath(fullfile(here, '..', '..', 'invz_common'));  % invz_ion, invz_const
end

% =====================================================================
% shared helpers (names deliberately free of "test" so functiontests skips them)
% =====================================================================
function a0 = alpha0_of(a)
a0 = sqrt(pi)/abs(det(a))^(1/3);
end

function eo = mk_eopts(alpha, r_cut, g_cut, boundary)
eo = struct('alpha', alpha, 'r_cut', r_cut, 'g_cut', g_cut, 'boundary', boundary);
end

function eo = default_eopts(a)
a0 = alpha0_of(a);
eo = mk_eopts(a0, 5.5/a0, 11*a0, 'conducting_k0_omitted');   % frozen production defaults
end

function m = mid_margin(A, B)
% M_id residual: <= 0 means PASS. Symmetric in A,B.
Tscale = max(max(abs(A(:))), max(abs(B(:))));
m = max(abs(A(:) - B(:))) - 1e-12*Tscale;
end

function x = indep_real_cells(a, tau, n, m, r_cut)
% Independent ndgrid enumeration of retained real cells for ordered pair (n,m),
% in the SAME deterministic order the primitive uses (self cell R=0 dropped by
% cell index only on the diagonal pair).
taucart = tau*a;
d = taucart(m,:) - taucart(n,:);
sa = min(svd(a));
nmax = ceil((r_cut + norm(d))/sa);
rng = -nmax:nmax;
[H, Kk, L] = ndgrid(rng, rng, rng);
hkl = [H(:) Kk(:) L(:)];
xx = hkl*a + d;
r = vecnorm(xx, 2, 2);
keep = (r <= r_cut);
if n == m
    keep = keep & ~all(hkl == 0, 2);
end
x = xx(keep, :);
end

function [dR, dG, dS, dG0] = reconstruct_parts(geom, qrow, alpha, g_cut)
% Independent (test-only) reconstruction from PUBLIC geom data of the real,
% full-reciprocal, self, and isolated-G=0 contributions. Every ordered (n,m)
% block is formed on its own; no conjugate block is copied.
ntau = geom.ntau; B = geom.B; Vc = geom.Vc; taucart = geom.taucart; tau = geom.tau;
Gcart = geom.Gcart; Ghkl = geom.Ghkl;
K = floor(qrow + 0.5); qbar = qrow - K; qcart = qbar*B;
k = Gcart + qcart; kk = sum(k.^2, 2); keep = (kk <= g_cut^2) & (kk > 0);
isG0 = all(Ghkl == 0, 2); g0k = isG0(keep);
selfval = 4*alpha^3/(3*sqrt(pi));
z = complex(zeros(3,3,ntau,ntau)); dR = z; dG = z; dS = z; dG0 = z;
ksel = k(keep,:); kk2 = kk(keep); kernel = (4*pi/Vc)*exp(-kk2/(4*alpha^2))./kk2;
Gk = Gcart(keep,:);
for n = 1:ntau
    for m = 1:ntau
        d = taucart(m,:) - taucart(n,:);
        gph = exp(-1i*2*pi*(K*(tau(m,:) - tau(n,:)).'));
        X = geom.real{n,m}.x; gab = geom.real{n,m}.gab; ph = exp(-1i*(X*qcart.'));
        bR = zeros(3,3);
        for aa = 1:3, for bb = 1:3, bR(aa,bb) = -sum(ph.*gab(:,aa,bb)); end, end %#ok<ALIGN>
        dR(:,:,n,m) = gph*bR;
        phG = exp(1i*(Gk*d.')); bG = zeros(3,3);
        for aa = 1:3, for bb = 1:3, bG(aa,bb) = sum(kernel.*ksel(:,aa).*ksel(:,bb).*phG); end, end %#ok<ALIGN>
        dG(:,:,n,m) = gph*bG;
        if any(g0k)
            ker0 = kernel(g0k); k0 = ksel(g0k,:); ph0 = exp(1i*(Gk(g0k,:)*d.')); b0 = zeros(3,3);
            for aa = 1:3, for bb = 1:3, b0(aa,bb) = sum(ker0.*k0(:,aa).*k0(:,bb).*ph0); end, end %#ok<ALIGN>
            dG0(:,:,n,m) = gph*b0;
        end
        if n == m, dS(:,:,n,m) = gph*(-selfval*eye(3)); end
    end
end
end

function u = recip_used_indep(geom, qrow, g_cut)
% Independent per-q retained reciprocal count.
K = floor(qrow + 0.5); qbar = qrow - K; qcart = qbar*geom.B;
k = geom.Gcart + qcart; kk = sum(k.^2, 2);
u = sum((kk <= g_cut^2) & (kk > 0));
end

% =====================================================================
% Task 1 -- contract, geometry, candidate union, preflight
% =====================================================================
function test_output_geometry_counts_schema(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
[dip, counts, geom] = invz_dipole_ewald([0.137 0.291 0.453], a, tau, default_eopts(a));
% output rank for nq==1
verifyEqual(testCase, size(dip), [3 3 ntau ntau]);
verifyTrue(testCase, isa(dip,'double') && ~isreal(dip));   % complex-typed double output
% counts schema
verifyEqual(testCase, size(counts.real_pair), [ntau ntau]);
verifyTrue(testCase, isscalar(counts.recip_candidates));
verifyEqual(testCase, size(counts.recip_used), [1 1]);     % [nq,1] with nq==1
verifyEqual(testCase, size(counts.preflight.real_cube_bound), [ntau ntau]);
verifyTrue(testCase, isscalar(counts.preflight.recip_cube_bound));
verifyTrue(testCase, isscalar(counts.preflight.estimated_peak_bytes));
verifyTrue(testCase, isstruct(counts.preflight.array_manifest));
% geom schema (fields the later reconstructions rely on)
for f = {'B','Vc','ntau','real','Ghkl','Gcart','recip','real_pair_count', ...
        'real_cube_bound','recip_cube_bound','fingerprint','tau','taucart', ...
        'alpha','r_cut','g_cut','boundary'}
    verifyTrue(testCase, isfield(geom, f{1}), sprintf('geom missing field %s', f{1}));
end
verifyEqual(testCase, geom.ntau, ntau);
end

function test_alpha0_rule_not_rounded(testCase)
ion = invz_ion(); a = ion.a;
a0 = alpha0_of(a);
% exact live rule, NOT the displayed rounded value 0.268431
verifyEqual(testCase, a0, sqrt(pi)/ion.Vc^(1/3), 'RelTol', 1e-14);
verifyGreaterThan(testCase, abs(a0 - 0.268431), 0);          % distinct from rounded
verifyLessThan(testCase, abs(a0 - 0.268431), 5e-5);          % but consistent with it
% the primitive accepts the live alpha0
verifyWarningFree(testCase, @() invz_dipole_ewald([0.1 0.2 0.3], a, ion.tau, default_eopts(a)));
end

function test_real_cutoff_and_self_cell_omission(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
eo = default_eopts(a);
[~, ~, geom] = invz_dipole_ewald([0 0 0], a, tau, eo);
for n = 1:ntau
    for m = 1:ntau
        x = geom.real{n,m}.x;
        r = vecnorm(x, 2, 2);
        verifyLessThanOrEqual(testCase, max(r), eo.r_cut);         % all within cutoff
        if n == m
            verifyGreaterThan(testCase, min(r), 0);                % self cell x=0 omitted
        end
        % exact set/order equality with independent enumeration
        xind = indep_real_cells(a, tau, n, m, eo.r_cut);
        verifyEqual(testCase, x, xind, 'AbsTol', 0);
    end
end
end

function test_deterministic_ordering_real_and_recip(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
eo = default_eopts(a);
[~, ~, gA] = invz_dipole_ewald([0.2 0.1 0.4], a, tau, eo);
[~, ~, gB] = invz_dipole_ewald([-0.3 0.05 0.2], a, tau, eo);   % different q, same geom
for n = 1:ntau
    for m = 1:ntau
        verifyEqual(testCase, gA.real{n,m}.x, gB.real{n,m}.x, 'AbsTol', 0);
        verifyEqual(testCase, gA.real{n,m}.gab, gB.real{n,m}.gab, 'AbsTol', 0);
    end
end
verifyEqual(testCase, gA.Gcart, gB.Gcart, 'AbsTol', 0);
verifyEqual(testCase, gA.Ghkl, gB.Ghkl, 'AbsTol', 0);
end

function test_G0_in_candidate_union(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau;
[~, ~, geom] = invz_dipole_ewald([0.1 0.2 0.3], a, tau, default_eopts(a));
verifyTrue(testCase, any(all(geom.Ghkl == 0, 2)), 'G=[0,0,0] must be in the cached candidate union');
end

function test_retained_counts_within_bounds(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
[~, counts] = invz_dipole_ewald([0.137 0.291 0.453], a, tau, default_eopts(a));
for n = 1:ntau
    for m = 1:ntau
        verifyLessThanOrEqual(testCase, counts.real_pair(n,m), counts.preflight.real_cube_bound(n,m));
    end
end
verifyLessThanOrEqual(testCase, counts.recip_candidates, counts.preflight.recip_cube_bound);
end

function test_manifest_row_fields(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau;
[~, counts] = invz_dipole_ewald([0.1 0.2 0.3], a, tau, default_eopts(a));
man = counts.preflight.array_manifest;
verifyEqual(testCase, sort(fieldnames(man)), sort({'name';'class';'is_complex';'size';'numel';'bytes'}));
for i = 1:numel(man)
    row = man(i);
    verifyTrue(testCase, ischar(row.name) && ~isempty(row.name));
    verifyTrue(testCase, ischar(row.class));
    verifyTrue(testCase, islogical(row.is_complex));
    verifyEqual(testCase, row.numel, prod(row.size));
    switch row.class
        case 'double'
            w = 8; if row.is_complex, w = 16; end
        case 'logical'
            w = 1;
        otherwise
            w = 8;
    end
    verifyEqual(testCase, row.bytes, row.numel*w);
end
end

function test_estimate_is_125x_sum_and_conservative(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
[~, counts] = invz_dipole_ewald([0.137 0.291 0.453], a, tau, default_eopts(a));
man = counts.preflight.array_manifest;
% exact 1.25 * sum(bytes)
verifyEqual(testCase, counts.preflight.estimated_peak_bytes, 1.25*sum([man.bytes]), 'RelTol', 0);
% the estimate is built from CONSERVATIVE cube bounds, not actual retained counts
names = {man.name};
outrow = man(strcmp(names, 'dip_output'));
verifyEqual(testCase, outrow.numel, 9*ntau^2*1);              % [3,3,ntau,ntau,nq=1]
verifyTrue(testCase, outrow.is_complex);
xrow = man(strcmp(names, 'real_x'));
verifyEqual(testCase, xrow.numel, sum(counts.preflight.real_cube_bound(:))*3);  % conservative, not actual
end

function test_nq_updates_qwork_output_estimate(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
eo = default_eopts(a);
[~, c1, geom] = invz_dipole_ewald([0.137 0.291 0.453], a, tau, eo);
q3 = [0.137 0.291 0.453; 0.25 0 0.1; -0.311 0.173 -0.227];
[dip3, c3] = invz_dipole_ewald(q3, a, tau, eo, geom);   % reuse geom, different nq
verifyEqual(testCase, size(dip3), [3 3 ntau ntau 3]);
verifyEqual(testCase, size(c3.recip_used), [3 1]);
% the q-work/output part of the estimate grew with nq
verifyGreaterThan(testCase, c3.preflight.estimated_peak_bytes, c1.preflight.estimated_peak_bytes);
o1 = local_named(c1.preflight.array_manifest, 'dip_output');
o3 = local_named(c3.preflight.array_manifest, 'dip_output');
verifyEqual(testCase, o3.numel, 3*o1.numel);                 % output scales exactly with nq
end

function row = local_named(man, name)
row = man(strcmp({man.name}, name));
end

function test_geom_reuse_identical_and_fingerprint_guard(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau;
eo = default_eopts(a);
q = [0.137 0.291 0.453];
[dref, ~, geom] = invz_dipole_ewald(q, a, tau, eo);
dreuse = invz_dipole_ewald(q, a, tau, eo, geom);
verifyEqual(testCase, dreuse, dref, 'AbsTol', 0);            % bit-identical reuse

% changing any single fingerprint field (via the CALL) must reject the stale geom
a2 = a; a2(1,1) = a2(1,1) + 1e-6;
verifyError(testCase, @() invz_dipole_ewald(q, a2, tau, eo, geom), 'invz:ewaldGeomReuse');
tau2 = tau; tau2(2,1) = tau2(2,1) + 1e-6;
verifyError(testCase, @() invz_dipole_ewald(q, a, tau2, eo, geom), 'invz:ewaldGeomReuse');
verifyError(testCase, @() invz_dipole_ewald(q, a, tau, mk_eopts(eo.alpha*1.01, eo.r_cut, eo.g_cut, eo.boundary), geom), 'invz:ewaldGeomReuse');
verifyError(testCase, @() invz_dipole_ewald(q, a, tau, mk_eopts(eo.alpha, eo.r_cut*1.01, eo.g_cut, eo.boundary), geom), 'invz:ewaldGeomReuse');
verifyError(testCase, @() invz_dipole_ewald(q, a, tau, mk_eopts(eo.alpha, eo.r_cut, eo.g_cut*1.01, eo.boundary), geom), 'invz:ewaldGeomReuse');

% tampering the stored convention / schema fields also rejects
gbad = geom; gbad.fingerprint.qconv = 'tampered';
verifyError(testCase, @() invz_dipole_ewald(q, a, tau, eo, gbad), 'invz:ewaldGeomReuse');
gbad2 = geom; gbad2.fingerprint.schema = 'other/v9';
verifyError(testCase, @() invz_dipole_ewald(q, a, tau, eo, gbad2), 'invz:ewaldGeomReuse');
end

function test_input_validation_errors(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau; eo = default_eopts(a);
% q shape / finiteness / realness
verifyError(testCase, @() invz_dipole_ewald([1 2], a, tau, eo), 'invz:ewaldArgs');
verifyError(testCase, @() invz_dipole_ewald([Inf 0 0], a, tau, eo), 'invz:ewaldArgs');
verifyError(testCase, @() invz_dipole_ewald([1+1i 0 0], a, tau, eo), 'invz:ewaldArgs');
% a shape / singular / finiteness
verifyError(testCase, @() invz_dipole_ewald([0 0 0], eye(2), tau, eo), 'invz:ewaldArgs');
verifyError(testCase, @() invz_dipole_ewald([0 0 0], [1 0 0;2 0 0;0 0 1], tau, eo), 'invz:ewaldArgs');
verifyError(testCase, @() invz_dipole_ewald([0 0 0], [Inf 0 0;0 1 0;0 0 1], tau, eo), 'invz:ewaldArgs');
% tau shape / finiteness
verifyError(testCase, @() invz_dipole_ewald([0 0 0], a, [1 2], eo), 'invz:ewaldArgs');
verifyError(testCase, @() invz_dipole_ewald([0 0 0], a, [NaN 0 0], eo), 'invz:ewaldArgs');
% eopts field set
verifyError(testCase, @() invz_dipole_ewald([0 0 0], a, tau, rmfield(eo,'g_cut')), 'invz:ewaldArgs');
eox = eo; eox.extra = 1;
verifyError(testCase, @() invz_dipole_ewald([0 0 0], a, tau, eox), 'invz:ewaldArgs');
% eopts control positivity / finiteness
verifyError(testCase, @() invz_dipole_ewald([0 0 0], a, tau, mk_eopts(0, eo.r_cut, eo.g_cut, eo.boundary)), 'invz:ewaldArgs');
verifyError(testCase, @() invz_dipole_ewald([0 0 0], a, tau, mk_eopts(eo.alpha, -1, eo.g_cut, eo.boundary)), 'invz:ewaldArgs');
verifyError(testCase, @() invz_dipole_ewald([0 0 0], a, tau, mk_eopts(eo.alpha, eo.r_cut, Inf, eo.boundary)), 'invz:ewaldArgs');
% boundary type / value
verifyError(testCase, @() invz_dipole_ewald([0 0 0], a, tau, mk_eopts(eo.alpha, eo.r_cut, eo.g_cut, 5)), 'invz:ewaldArgs');
verifyError(testCase, @() invz_dipole_ewald([0 0 0], a, tau, mk_eopts(eo.alpha, eo.r_cut, eo.g_cut, 'lorentz')), 'invz:ewaldBoundary');
end

function test_boundary_and_accuracy_guards(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau;
% alpha*r_cut just below 4.5
verifyError(testCase, @() invz_dipole_ewald([0 0 0], a, tau, mk_eopts(1, 4.4, 9, 'conducting_k0_omitted')), 'invz:ewaldGuard');
% g_cut/(2 alpha) just below 4.5
verifyError(testCase, @() invz_dipole_ewald([0 0 0], a, tau, mk_eopts(1, 5, 8.9, 'conducting_k0_omitted')), 'invz:ewaldGuard');
% exact boundary (alpha*r_cut = 4.5, g_cut/(2 alpha) = 4.5) is ACCEPTED
[dip, ~] = invz_dipole_ewald([0.1 0.2 0.3], a, tau, mk_eopts(1, 4.5, 9, 'conducting_k0_omitted'));
verifyEqual(testCase, size(dip), [3 3 4 4]);
end

% =====================================================================
% Task 2 -- assembly, shape, reconstruction, Hermiticity, recip_used, G0
% =====================================================================
function test_shape_one_q_and_multi_q(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
eo = default_eopts(a);
[d1, c1] = invz_dipole_ewald([0.137 0.291 0.453], a, tau, eo);
verifyEqual(testCase, size(d1), [3 3 ntau ntau]);
verifyEqual(testCase, size(c1.recip_used), [1 1]);
q = [0.137 0.291 0.453; 0.25 0 0.1; -0.311 0.173 -0.227];
[d3, c3] = invz_dipole_ewald(q, a, tau, eo);
verifyEqual(testCase, size(d3), [3 3 ntau ntau 3]);
verifyEqual(testCase, size(c3.recip_used), [3 1]);
end

function test_reconstruction_sum_matches(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau;
eo = default_eopts(a);
qset = [0 0 0
        0.137 0.291 0.453
        0.25 0 0.1
       -0.311 0.173 -0.227
        1.137 0.291 -1.547];       % last has nonzero K -> exercises gauge restore
[dip, ~, geom] = invz_dipole_ewald(qset, a, tau, eo);
for i = 1:size(qset,1)
    [dR, dG, dS] = reconstruct_parts(geom, qset(i,:), eo.alpha, eo.g_cut);
    recon = dR + dG + dS;                     % real + full reciprocal + self
    verifyLessThanOrEqual(testCase, mid_margin(recon, dip(:,:,:,:,i)), 0);
end
end

function test_returned_hermiticity(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
eo = default_eopts(a);
qset = [0.137 0.291 0.453; 0.25 0 0.1; -0.311 0.173 -0.227];
dip = invz_dipole_ewald(qset, a, tau, eo);
for i = 1:size(qset,1)
    A = dip(:,:,:,:,i);
    herm = complex(zeros(3,3,ntau,ntau));
    for n = 1:ntau
        for m = 1:ntau
            herm(:,:,n,m) = conj(A(:,:,m,n));   % emergent identity, not how it was built
        end
    end
    verifyLessThanOrEqual(testCase, mid_margin(A, herm), 0);
end
end

function test_recip_used_matches_enumeration(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau;
eo = default_eopts(a);
qset = [0 0 0; 0.137 0.291 0.453; 0.25 0 0.1; 1.137 0.291 -1.547];
[~, counts, geom] = invz_dipole_ewald(qset, a, tau, eo);
for i = 1:size(qset,1)
    verifyEqual(testCase, counts.recip_used(i), recip_used_indep(geom, qset(i,:), eo.g_cut));
end
end

function test_G0_used_at_nonzero_q_and_omitted_at_gamma(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau;
eo = default_eopts(a);
[~, ~, geom] = invz_dipole_ewald([0 0 0], a, tau, eo);
isG0 = all(geom.Ghkl == 0, 2);
gmag = vecnorm(geom.Gcart, 2, 2);

% nonzero reduced q: the G=0 candidate participates and its isolated term is nonzero
qnz = [0.137 0.291 0.453];
qcart = (qnz - floor(qnz+0.5))*geom.B;
knz = geom.Gcart + qcart; kknz = sum(knz.^2,2);
keepnz = (kknz <= eo.g_cut^2) & (kknz > 0);
verifyTrue(testCase, keepnz(isG0), 'G=0 candidate must be used at nonzero reduced q');
[~, ~, ~, dG0] = reconstruct_parts(geom, qnz, eo.alpha, eo.g_cut);
verifyGreaterThan(testCase, max(abs(dG0(:))), 0);

% exact Gamma: G=0 is the ONE omitted k among |Gcart|<=g_cut
k0 = geom.Gcart; kk0 = sum(k0.^2,2); keep0 = (kk0 <= eo.g_cut^2) & (kk0 > 0);
verifyFalse(testCase, keep0(isG0));                              % G=0 omitted at Gamma
verifyEqual(testCase, sum(~keep0 & (gmag <= eo.g_cut)), 1);      % exactly one omission
end

function test_gauge_restore_sign(testCase)
% dip_nm(qbar+K) = exp(-i*2*pi*K.(tau_m-tau_n)) * dip_nm(qbar). Catches a wrong
% gauge sign (which would leave an O(1) residual on pairs with K.(tau_m-tau_n)
% non-integer, e.g. tau4-tau1 = [0.5 0 0.75] with K = [0 0 1]).
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
eo = default_eopts(a);
qbar = [0.137 0.291 0.453];
K = [0 0 1];
[~, ~, geom] = invz_dipole_ewald(qbar, a, tau, eo);
dbar = invz_dipole_ewald(qbar, a, tau, eo, geom);
dK   = invz_dipole_ewald(qbar + K, a, tau, eo, geom);
pred = complex(zeros(3,3,ntau,ntau));
for n = 1:ntau
    for m = 1:ntau
        gph = exp(-1i*2*pi*(K*(tau(m,:) - tau(n,:)).'));
        pred(:,:,n,m) = gph*dbar(:,:,n,m);
    end
end
verifyLessThanOrEqual(testCase, mid_margin(dK, pred), 0);
end

function test_gamma_tensor_is_real(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau;
eo = default_eopts(a);
dip0 = invz_dipole_ewald([0 0 0], a, tau, eo);
verifyLessThanOrEqual(testCase, max(abs(imag(dip0(:)))), 1e-12*max(abs(dip0(:))));
end

function test_mfdipole_secondary_crosscheck(testCase)
% SECONDARY diagnostic only (not a Gate-A acceptance): the conducting_k0_omitted
% Ewald tensor at a generic q away from Gamma is cross-checked against the
% brute-force MF_dipole spherical truncation (|R+d|<=N*min|a|), N=6..20. We record
% that the finite-N brute force converges toward the Ewald value; we do NOT
% promote its truncation residual to a frozen tolerance.
ion = invz_ion(); a = ion.a; tau = ion.tau;
eo = default_eopts(a);
q = [0.137 0.291 0.453];
dew = invz_dipole_ewald(q, a, tau, eo);
Ns = [6 10 16 20];
d = zeros(size(Ns));
for j = 1:numel(Ns)
    dmf = MF_dipole(q, Ns(j), a, tau);
    d(j) = max(abs(dmf(:) - dew(:)));
end
% overall convergence toward Ewald + finest brute force is close (loose bound)
verifyLessThan(testCase, d(end), d(1));
verifyLessThan(testCase, d(end), 1e-3);
end

% =====================================================================
% Task 1 -- resource caps are a hard pre-ALLOCATION gate (distinct ids)
% =====================================================================
function test_real_cap_fires_before_allocation(testCase)
% eye lattice + huge r_cut drives the per-pair real cube bound past 3.0e6 while
% the guards still pass; must raise invz:ewaldRealCap without allocating.
eo = mk_eopts(0.05, 200, 0.45, 'conducting_k0_omitted');   % a*r=10>=4.5, g/(2a)=4.5
verifyError(testCase, @() invz_dipole_ewald([0 0 0], eye(3), [0 0 0; 0.2 0.2 0.2], eo), ...
    'invz:ewaldRealCap');
end

function test_recip_cap_fires_before_allocation(testCase)
% large lattice shrinks svd_min(B) so the reciprocal cube bound exceeds 3.0e6
% while the real bound stays tiny; must raise invz:ewaldRecipCap.
eo = mk_eopts(1, 5, 10, 'conducting_k0_omitted');          % a*r=5>=4.5, g/(2a)=5
verifyError(testCase, @() invz_dipole_ewald([0 0 0], 100*eye(3), [0 0 0; 0.2 0.2 0.2], eo), ...
    'invz:ewaldRecipCap');
end

function test_memory_cap_fires_before_allocation(testCase)
% dimension-only lower bound (complex output 9*ntau^2*nq) exceeds 4 GiB; must
% raise invz:ewaldMemoryCap before any geometry is built.
ntau = 40; nq = 20000;
tau = 0.1*ones(ntau,3) + (1:ntau).'*1e-3;   % valid finite basis
q = zeros(nq,3);
eo = mk_eopts(1, 5, 9, 'conducting_k0_omitted');
verifyError(testCase, @() invz_dipole_ewald(q, 3*eye(3), tau, eo), 'invz:ewaldMemoryCap');
end
