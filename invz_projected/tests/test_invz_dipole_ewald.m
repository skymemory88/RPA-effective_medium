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

% =====================================================================
% Task 3 -- Gate-A test 1: screened-Hessian finite-difference validation
% =====================================================================
function test_gateA1_screened_hessian(testCase)
% Gate-A #1 (docs/invzp_ewald_prereg.md sec 3, item 1): the primitive's
% precomputed real-space screened tensor g_ab must equal the analytic
% Hessian of f(y) = erfc(alpha*|y|)/|y|, validated three ways:
%   (i)   Richardson-extrapolated central-difference Hessians (steps 4e-3,
%         2e-3, 1e-3 Angstrom -> R12, R23) agree with each other and with an
%         INDEPENDENTLY-derived closed g_ab (indep_closed_gab below), over
%         all 200 frozen hess_x samples at each of the 5 frozen alpha values;
%   (ii)  the closed g_ab reduces EXACTLY to the bare
%         (3 x_a x_b - r^2 delta_ab)/r^5 tensor at alpha = 0 (P->1, Q->0),
%         not merely to within a loose 1e-8 floor;
%   (iii) BRIDGE (anti-tautology): the primitive's own geom.real{n,m}.gab,
%         for every retained real vector of every ordered sublattice pair,
%         matches this SEPARATELY-written closed g_ab at M_id -- this is
%         what would catch a wrong PRODUCTION formula (a test that only
%         re-checks its own duplicated formula cannot).
%
% indep_closed_gab (below) is derived fresh from calculus --
%   h(r) = erfc(alpha r)/r,
%   d_a d_b h = (h''(r) - h'(r)/r) x_a x_b/r^2 + (h'(r)/r) delta_ab  (the
%   general radial-function Hessian identity) --
% NOT copied from invz_dipole_ewald's private P/Q-form local_gab: it is a
% different algebraic decomposition (delta/xaxb coefficients A(r), B(r)
% derived directly from h, h', h'' rather than P(r), Q(r) applied to the
% bare tensor), independently typed, so a coefficient/sign/exponent typo
% specific to the production code would generally NOT be replicated here.
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();

% ---- fixture cardinalities this test exercises (fail fast, non-tautological)
verifyEqual(testCase, size(fx.q_int), [3 3]);
verifyEqual(testCase, size(fx.QA), [30 3]);
verifyEqual(testCase, size(fx.Kset), [8 3]);
verifyEqual(testCase, size(fx.hess_x), [200 3]);
verifyEqual(testCase, size(unique(fx.hess_x, 'rows'), 1), 200);   % no duplicate sample
verifyGreaterThan(testCase, fx.alpha0, 0);
verifyEqual(testCase, fx.alpha0, sqrt(pi)/abs(det(invz_ion().a))^(1/3), 'RelTol', 1e-14);
verifyGreaterThan(testCase, abs(fx.alpha0 - 0.268431), 0);   % never the rounded literal

h1 = 4e-3; h2 = 2e-3; h3 = 1e-3;
mult = [0.6 0.8 1.0 1.2 1.5];
nX = size(fx.hess_x, 1);

worst_analytic_ratio   = 0;
worst_richardson_ratio = 0;

for ai = 1:numel(mult)
    alpha = mult(ai)*fx.alpha0;
    for k = 1:nX
        x = fx.hess_x(k,:);
        closed = squeeze(indep_closed_gab(x, alpha));

        Hh1 = fd_hessian_gab(x, alpha, h1);
        Hh2 = fd_hessian_gab(x, alpha, h2);
        Hh3 = fd_hessian_gab(x, alpha, h3);
        R12 = (4*Hh2 - Hh1)/3;
        R23 = (4*Hh3 - Hh2)/3;

        mrich = M.mhfd(R12, R23, x);
        verifyLessThanOrEqual(testCase, mrich.worst_margin, 0, sprintf( ...
            'Gate-A1: R12/R23 Richardson estimates disagree by > 1e-7*H_scale at alpha=%.6g (x%.2g), sample %d.', ...
            alpha, mult(ai), k));

        mclosed = M.mhfd(closed, R23, x);
        verifyLessThanOrEqual(testCase, mclosed.worst_margin, 0, sprintf( ...
            'Gate-A1: R23 fails M_HFD vs the closed g_ab at alpha=%.6g (x%.2g), sample %d.', ...
            alpha, mult(ai), k));
        verifyTrue(testCase, mclosed.sign_ok, sprintf( ...
            'Gate-A1: sign mismatch on a gated component (closed vs R23) at alpha=%.6g (x%.2g), sample %d.', ...
            alpha, mult(ai), k));

        worst_analytic_ratio   = max(worst_analytic_ratio, mclosed.worst_ratio);
        worst_richardson_ratio = max(worst_richardson_ratio, mrich.worst_ratio);
    end
end
fprintf(['Gate-A1 screened-Hessian: worst normalized analytic ratio = %.3e, ' ...
    'worst adjacent-Richardson ratio = %.3e (tol 1e-7 each).\n'], ...
    worst_analytic_ratio, worst_richardson_ratio);

% ---- alpha = 0: closed g_ab must reduce EXACTLY to the bare tensor ----------
for k = 1:nX
    x = fx.hess_x(k,:);
    closed0 = squeeze(indep_closed_gab(x, 0));
    r = norm(x);
    bare = (3*(x.'*x) - r^2*eye(3))/r^5;
    m0 = M.mid(closed0, bare);
    verifyLessThanOrEqual(testCase, m0.worst_margin, 0, sprintf( ...
        'Gate-A1: alpha=0 closed g_ab is not EXACTLY the bare tensor at sample %d.', k));
end

% ---- BRIDGE: production geom.real{n,m}.gab vs the independent closed g_ab --
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
eo = default_eopts(a);
[~, ~, geom] = invz_dipole_ewald(fx.q_int(1,:), a, tau, eo);
worst_bridge_margin = -inf;
for n = 1:ntau
    for m = 1:ntau
        X       = geom.real{n,m}.x;
        prodgab = geom.real{n,m}.gab;
        verifyGreaterThan(testCase, size(X,1), 0, sprintf( ...
            'Gate-A1 bridge: pair (%d,%d) retains zero real vectors -- cannot bridge.', n, m));
        indepgab = indep_closed_gab(X, eo.alpha);
        mbridge = M.mid(prodgab, indepgab);
        verifyLessThanOrEqual(testCase, mbridge.worst_margin, 0, sprintf( ...
            'Gate-A1 bridge: geom.real{%d,%d}.gab disagrees with the independent closed g_ab (M_id).', n, m));
        worst_bridge_margin = max(worst_bridge_margin, mbridge.worst_margin);
    end
end
fprintf(['Gate-A1 bridge (production geom.real.gab vs independent closed g_ab): ' ...
    'worst M_id margin = %.3e.\n'], worst_bridge_margin);
end

% =====================================================================
% Task 3 helpers -- INDEPENDENT closed-form screened Hessian + FD stencils
% (do not reuse invz_dipole_ewald's private local_gab; see test docstring)
% =====================================================================
function gab = indep_closed_gab(x, alpha)
% Independent closed-form screened dipolar Hessian g_ab(x) = d_a d_b h(r),
% h(r) = erfc(alpha*r)/r, r = |x|, derived via the general radial-function
% Hessian identity
%   d_a d_b h = (h''(r) - h'(r)/r) x_a x_b/r^2 + (h'(r)/r) delta_ab.
% x: [K,3] Cartesian vectors (K>=1, all r>0). Returns gab: [K,3,3].
x = reshape(x, [], 3);
K = size(x,1);
r = vecnorm(x, 2, 2);
if any(r <= 0)
    error('invz:ewaldMetricShape', 'indep_closed_gab: all sample points must have |x|>0.');
end
r2 = r.^2; r3 = r2.*r; r4 = r2.^2; r5 = r4.*r;
z  = alpha*r;
E  = exp(-z.^2);
Ec = erfc(z);
sp = sqrt(pi);
% delta_ab coefficient: h'(r)/r
Acoef = -( Ec./r3 + (2*alpha/sp)*E./r2 );
% x_a x_b coefficient: (h''(r) - h'(r)/r)/r^2
Bcoef = (4*alpha^3/sp)*E./r2 + (6*alpha/sp)*E./r4 + 3*Ec./r5;
gab = zeros(K,3,3);
for aa = 1:3
    for bb = 1:3
        gab(:,aa,bb) = Bcoef.*(x(:,aa).*x(:,bb)) + Acoef.*double(aa==bb);
    end
end
end

function H = fd_hessian_gab(x, alpha, h)
% Central-difference Hessian of f(y)=erfc(alpha*|y|)/|y| at y=x, step h.
% Diagonal (a=b): 3-point [f(x+h e_a)-2f(x)+f(x-h e_a)]/h^2.
% Off-diagonal (a~=b): 4-point mixed
%   [f(x+h e_a+h e_b)-f(x+h e_a-h e_b)-f(x-h e_a+h e_b)+f(x-h e_a-h e_b)]/(4h^2).
% The mixed formula is NEVER used for a=b -- that would silently turn the
% labelled step h into a 2h diagonal stencil.
x = reshape(x, 1, 3);
H = zeros(3,3);
f0 = screened_f(x, alpha);
for aa = 1:3
    ea = zeros(1,3); ea(aa) = 1;
    for bb = 1:3
        if aa == bb
            H(aa,bb) = (screened_f(x+h*ea, alpha) - 2*f0 + screened_f(x-h*ea, alpha))/h^2;
        else
            eb = zeros(1,3); eb(bb) = 1;
            H(aa,bb) = (screened_f(x+h*ea+h*eb, alpha) - screened_f(x+h*ea-h*eb, alpha) ...
                      - screened_f(x-h*ea+h*eb, alpha) + screened_f(x-h*ea-h*eb, alpha))/(4*h^2);
        end
    end
end
end

function f = screened_f(y, alpha)
r = norm(y);
f = erfc(alpha*r)/r;
end

% =====================================================================
% Task 4 -- Gate-A tests 2, 6, 7, 8 (small-shell real signs / self term /
% boundary + no-surface-term / counts + caps)
% =====================================================================
function test_gateA2_small_shell_real_signs(testCase)
% Gate-A #2 (docs/invzp_ewald_prereg.md sec 3, item 2): explicit small-shell
% real-space signs and phase. Uses the frozen two-site fx.shell_fixture.
% (i)  independently enumerate every integer cell with |R+d|<=r_cut and build
%      the complete [3,3,2,2] screened REAL tensor from the defining formula
%      (leading minus, total-displacement phase e^{-i q_cart.(R+d)});
% (ii) independently reconstruct the primitive's real contribution from its
%      retained geom.real{n,m}.{x,gab} (reconstruct_parts' dR, ignoring
%      dG/dS/dG0 -- the reciprocal part is deliberately never added here).
% Compare the two COMPLETE tensors at M_T (do NOT re-add a reciprocal copy;
% A2 is specifically the finite screened real shell).
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();
sf = fx.shell_fixture;
a = sf.a; tau = sf.tau; alpha = sf.alpha; r_cut = sf.r_cut; q = sf.q;
ntau = size(tau,1);
g_cut = 11*alpha;   % valid synthetic reciprocal cutoff (guard g_cut/(2a)=5.5>=4.5);
                    % A2 never uses the reciprocal part -- only needed to form a valid eopts.
eo = mk_eopts(alpha, r_cut, g_cut, 'conducting_k0_omitted');
[~, ~, geom] = invz_dipole_ewald(q, a, tau, eo);

K = floor(q + 0.5); qbar = q - K; qcart = qbar*geom.B;

% ---- (i) independent reference: fresh cell enumeration + defining formula --
ref = complex(zeros(3,3,ntau,ntau));
for n = 1:ntau
    for m = 1:ntau
        x = indep_real_cells(a, tau, n, m, r_cut);
        verifyGreaterThan(testCase, size(x,1), 0, sprintf( ...
            'Gate-A2: pair (%d,%d) retains zero real vectors in the small-shell fixture.', n, m));
        gab = indep_closed_gab(x, alpha);
        ph  = exp(-1i*(x*qcart.'));                      % total-displacement phase e^{-i q.(R+d)}
        gph = exp(-1i*2*pi*(K*(tau(m,:) - tau(n,:)).'));
        blk = zeros(3,3);
        for aa = 1:3, for bb = 1:3, blk(aa,bb) = -sum(ph.*gab(:,aa,bb)); end, end %#ok<ALIGN>
        ref(:,:,n,m) = gph*blk;                           % leading minus sign
    end
end

% ---- (ii) independent reconstruction of the primitive's own retained geometry
[dR, ~, ~, ~] = reconstruct_parts(geom, q, alpha, g_cut);

mres = M.mt(ref, dR);
verifyTrue(testCase, mres.pass, sprintf( ...
    'Gate-A2: small-shell real-space tensor mismatch (signs/phase), worst_margin=%.3e.', mres.worst_margin));
fprintf('Gate-A2 small-shell real signs: worst M_T margin = %.3e (ratio %.3e).\n', ...
    mres.worst_margin, mres.worst_ratio);
end

function test_gateA6_self_term(testCase)
% Gate-A #6 (prereg sec 3, item 6): isolate the self term by subtracting an
% independent geom-reconstruction of the real + reciprocal NON-SELF parts from
% the returned primitive; compare the isolated actual self tensor against one
% complete expected self tensor (-4*alpha^3/(3*sqrt(pi))*eye(3) on every
% same-site block, zero off-site) at M_T, at every frozen positive alpha. Uses
% a generic nonzero frozen q (fx.q_int(1,:)) throughout -- NOT Gamma-only --
% so this is a genuine reconstruction/subtraction test, not an
% imaginary-trace/Gamma-realness placeholder.
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();
q  = fx.q_int(1,:);
mult = [0.6 0.8 1.0 1.2 1.5];

sameMask = false(3,3,ntau,ntau);
for n = 1:ntau
    sameMask(:,:,n,n) = true;
end

worst_margin = -inf;
for ai = 1:numel(mult)
    alpha = mult(ai)*fx.alpha0;
    eo = mk_eopts(alpha, 5.5/alpha, 11*alpha, 'conducting_k0_omitted');
    [dip, ~, geom] = invz_dipole_ewald(q, a, tau, eo);

    [dR, dG, ~, ~] = reconstruct_parts(geom, q, alpha, eo.g_cut);   % real + reciprocal NON-SELF
    dS_actual = dip - dR - dG;                                      % isolated actual self term

    selfval = 4*alpha^3/(3*sqrt(pi));
    dS_expected = complex(zeros(3,3,ntau,ntau));
    for n = 1:ntau
        dS_expected(:,:,n,n) = -selfval*eye(3);
    end

    mres = M.mt(dS_actual, dS_expected);
    verifyTrue(testCase, mres.pass, sprintf( ...
        'Gate-A6: isolated self tensor mismatch at alpha=%.6g (x%.2g), worst_margin=%.3e.', ...
        alpha, mult(ai), mres.worst_margin));

    sameMax = max(abs(dS_actual(sameMask) - dS_expected(sameMask)));
    offMax  = max(abs(dS_actual(~sameMask) - dS_expected(~sameMask)));
    fprintf(['Gate-A6 self term (alpha=%.6g, x%.2g): same-site max|err|=%.3e, ' ...
        'off-site max|err|=%.3e (M_T worst_margin=%.3e, complete-tensor scale unchanged).\n'], ...
        alpha, mult(ai), sameMax, offMax, mres.worst_margin);

    worst_margin = max(worst_margin, mres.worst_margin);
end
fprintf('Gate-A6 self term: worst M_T margin over all frozen alpha = %.3e.\n', worst_margin);
end

function test_gateA7_boundary_and_no_surface_term(testCase)
% Gate-A #7 (prereg sec 3, item 7): (1) every non-frozen boundary string
% raises invz:ewaldBoundary; (2) an unknown surface/demag-flavored eopts
% control raises invz:ewaldArgs (an unknown primitive control, exactly as any
% other unknown field would); (3) a geom-reconstruction (real+recip+self)
% equals the primitive at M_id across every fx.QA point (Gamma + q_int + all
% face/edge/corner candidate-boundary probes), proving the output contains
% EXACTLY these three terms and no fourth macroscopic (surface/shape/demag)
% term.
ion = invz_ion(); a = ion.a; tau = ion.tau;
eo = default_eopts(a);

% ---- (1) every non-frozen boundary string errors with the SPECIFIC id -----
bad_boundaries = {'lorentz','vacuum','open','tinfoil','shape_dependent','surface_dipole'};
for i = 1:numel(bad_boundaries)
    eob = mk_eopts(eo.alpha, eo.r_cut, eo.g_cut, bad_boundaries{i});
    verifyError(testCase, @() invz_dipole_ewald([0.1 0.2 0.3], a, tau, eob), ...
        'invz:ewaldBoundary', sprintf('Gate-A7: boundary ''%s'' must raise invz:ewaldBoundary.', bad_boundaries{i}));
end

% ---- (2) unknown surface/demag-flavored controls error as unknown controls -
extra_fields = {'surface_correction','demag_shape','lorentz_radius'};
for i = 1:numel(extra_fields)
    eox = eo; eox.(extra_fields{i}) = 1;
    verifyError(testCase, @() invz_dipole_ewald([0.1 0.2 0.3], a, tau, eox), ...
        'invz:ewaldArgs', sprintf('Gate-A7: unknown control ''%s'' must raise invz:ewaldArgs.', extra_fields{i}));
end

% ---- (3) no fourth macroscopic term: geom-reconstruction == primitive, ----
% ---- at M_id, over every fx.QA candidate-boundary probe -------------------
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();
[dip, ~, geom] = invz_dipole_ewald(fx.QA, a, tau, eo);
worst_margin = -inf;
for i = 1:size(fx.QA,1)
    [dR, dG, dS, ~] = reconstruct_parts(geom, fx.QA(i,:), eo.alpha, eo.g_cut);
    recon = dR + dG + dS;
    mres = M.mid(recon, dip(:,:,:,:,i));
    verifyTrue(testCase, mres.pass, sprintf( ...
        ['Gate-A7: geom-reconstruction (real+recip+self) fails M_id at fx.QA row %d -- ' ...
         'implies a fourth macroscopic term, worst_margin=%.3e.'], i, mres.worst_margin));
    worst_margin = max(worst_margin, mres.worst_margin);
end
fprintf(['Gate-A7 no-fourth-term reconstruction: worst M_id margin over all %d fx.QA points ' ...
    '= %.3e.\n'], size(fx.QA,1), worst_margin);
end

function test_gateA8_counts_and_caps(testCase)
% Gate-A #8 (prereg sec 3, item 8): (1) counts.real_pair / recip_candidates /
% recip_used each match an independent enumeration; (2) all three hard caps
% fire their SPECIFIC error id, with synthetic small-input configurations,
% and without the prohibited large allocation (timing kept diagnostic only).
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
eo = default_eopts(a);
fx = invz_ewald_fixtures();
qset = fx.q_int;
[~, counts, geom] = invz_dipole_ewald(qset, a, tau, eo);

% ---- counts.real_pair vs independent integer-cell enumeration -------------
for n = 1:ntau
    for m = 1:ntau
        xind = indep_real_cells(a, tau, n, m, eo.r_cut);
        verifyEqual(testCase, counts.real_pair(n,m), size(xind,1), sprintf( ...
            'Gate-A8: counts.real_pair(%d,%d) disagrees with independent cell enumeration.', n, m));
    end
end

% ---- counts.recip_candidates vs independent closed-box+slack enumeration ---
cnt_ind = indep_recip_candidates(a, eo.g_cut);
verifyEqual(testCase, counts.recip_candidates, cnt_ind, ...
    'Gate-A8: counts.recip_candidates disagrees with the independent closed-box(+1e-12 slack) enumeration.');

% ---- each counts.recip_used(iq) vs independent per-q filtering ------------
for i = 1:size(qset,1)
    verifyEqual(testCase, counts.recip_used(i), recip_used_indep(geom, qset(i,:), eo.g_cut), sprintf( ...
        'Gate-A8: counts.recip_used(%d) disagrees with independent per-q filtering.', i));
end
fprintf('Gate-A8 counts: real_pair / recip_candidates (=%d) / recip_used all match independent enumeration.\n', ...
    counts.recip_candidates);

% ---- three hard caps: SPECIFIC id, no prohibited large allocation ----------
% (timing is a diagnostic only -- never turned into a pass/fail assertion)
eo_real = mk_eopts(0.05, 400, 0.45, 'conducting_k0_omitted');    % a*r=20>=4.5, g/(2a)=4.5
tR0 = tic;
verifyError(testCase, @() invz_dipole_ewald([0 0 0], eye(3), [0 0 0; 0.2 0.2 0.2], eo_real), ...
    'invz:ewaldRealCap');
tR = toc(tR0);

eo_recip = mk_eopts(1, 5, 20, 'conducting_k0_omitted');          % a*r=5>=4.5, g/(2a)=10
tG0 = tic;
verifyError(testCase, @() invz_dipole_ewald([0 0 0], 400*eye(3), [0 0 0; 0.2 0.2 0.2], eo_recip), ...
    'invz:ewaldRecipCap');
tG = toc(tG0);

ntauM = 50; nqM = 15000;
tauM = 0.1*ones(ntauM,3) + (1:ntauM).'*1e-3;
qM = zeros(nqM,3);
eo_mem = mk_eopts(1, 5, 9, 'conducting_k0_omitted');
tM0 = tic;
verifyError(testCase, @() invz_dipole_ewald(qM, 3*eye(3), tauM, eo_mem), 'invz:ewaldMemoryCap');
tM = toc(tM0);

fprintf(['Gate-A8 caps (timing diagnostic only, not asserted): RealCap %.3fs, RecipCap %.3fs, ' ...
    'MemoryCap %.3fs.\n'], tR, tG, tM);
end

% =====================================================================
% Task 4 helper -- independent reciprocal candidate-union enumeration
% =====================================================================
function cnt = indep_recip_candidates(a, g_cut)
% Independent count of the reciprocal candidate union, replicating the EXACT
% frozen convention (docs/invzp_ewald_prereg.md): closed box qbar in
% [-0.5,0.5]^3, keep every integer G with
%   dmin(G) = min_{qbar in box} |(qbar+G)*B| <= g_cut*(1+1e-12).
% Uses a CLOSED-FORM per-axis minimization -- NOT the primitive's 27
% free/lower/upper active-set enumeration -- valid because `a` (hence
% B=2*pi*inv(a)') is EXACTLY DIAGONAL for every lattice this helper is called
% with (asserted below): the quadratic |(qbar+G)*B|^2 then separates
% axis-by-axis, and each 1-D box minimization has the closed form
% max(|G_j|-0.5,0)*B_jj (0 exactly when G_j==0, since 0 is then inside the
% box; otherwise the nearest box point to the origin on that axis).
assert(isequal(a, diag(diag(a))), ...
    'indep_recip_candidates: this closed-form independent check requires a diagonal lattice matrix.');
B  = 2*pi*inv(a).';  %#ok<MINV>
Bd = diag(B).';                                    % [1,3], B is exactly diagonal
corners01 = [-0.5 0.5];
[c1,c2,c3] = ndgrid(corners01, corners01, corners01);
corners = [c1(:) c2(:) c3(:)];
qmax = max(vecnorm(corners*B, 2, 2));
nmax_G = ceil((g_cut + qmax)/min(abs(Bd))) + 1;    % generous conservative enumeration range
rng = -nmax_G:nmax_G;
[H,Kk,L] = ndgrid(rng, rng, rng);
Ghkl_all = [H(:) Kk(:) L(:)];
gapv = max(abs(Ghkl_all) - 0.5, 0);                % [N,3] per-axis closed-form box-min |v_j|
dmin = sqrt(sum((gapv.*Bd).^2, 2));
keep = dmin <= g_cut*(1 + 1e-12);
cnt = nnz(keep);
end

% =====================================================================
% Task 5 -- Gate-A tests 3, 4 (reciprocal gauge covariance + periodic
% coupling spectrum; half-open reduction + candidate completeness) and the
% shared test-only coupling mapper (reused by a later task's Gate-A#9)
% =====================================================================
function test_gateA3_reciprocal_gauge_and_periodic_spectrum(testCase)
% Gate-A #3 (docs/invzp_ewald_prereg.md sec 3, item 3): reciprocal gauge
% covariance + periodic coupling spectrum. For EVERY q in fx.QA (30 points:
% Gamma + 3 generic + 6 face + 12 edge + 8 corner) and EVERY frozen K in
% fx.Kset (8 translations, 240 (q,K) pairs total):
%   (i)  the full gauge-predicted tensor
%        pred_nm(q,K) = exp(-i*2*pi*K.(tau_m-tau_n)) * dip_nm(q)
%        must equal dip_nm(q+K) EXACTLY at M_id. Per the Task-3 reviewer's
%        LOAD-BEARING forward-note, M_id is called ONCE PER (q,K) on that
%        single q's COMPLETE [3,3,ntau,ntau] tensor -- never on a stacked
%        multi-q array (which would silently loosen T_scale to the worst
%        component over the whole stack);
%   (ii) via the shared ewald_coupling_mapper (mirrors the non-ODD Ewald
%        algebra through `exchange`, NO Ewald Lorentz term added -- the
%        Ewald dip already carries the isotropic G=0-omitted term), the four
%        sorted cc branches Jnu(q) and Jnu(q+K) agree under a DIRECT
%        AbsTol=1e-10 meV / RelTol=1e-8 tolerance (the frozen Phase-1 item-2
%        coupling tolerance -- NOT invz_ewald_metrics' M_J, whose internal
%        AbsTol_J=1e-8*J_ref literal is a DIFFERENT frozen number). This meV
%        tolerance is applied ONLY to the coupling branches -- never directly
%        to the raw inverse-cubic-Angstrom dip eigenvalues.
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
eo = default_eopts(a);
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();

[~, geomX] = exchange([0 0 0], abs(ion.J12), a, tau);      % prime exchange geom ONCE
[dipQ, ~, geom] = invz_dipole_ewald(fx.QA, a, tau, eo);    % ONE geom, ONE batched dip(Q_A) call

nQ = size(fx.QA,1); nK = size(fx.Kset,1);
outQ = cell(nQ,1);
for qi = 1:nQ
    outQ{qi} = ewald_coupling_mapper(fx.QA(qi,:), dipQ(:,:,:,:,qi), geomX);
end

AbsTol_J3 = 1e-10; RelTol_J3 = 1e-8;   % frozen Phase-1 item-2 tolerance (prereg sec 3 item 3)

worst_tensor_margin   = -inf;
worst_coupling_margin = -inf;
for ki = 1:nK
    K = fx.Kset(ki,:);
    dipQK = invz_dipole_ewald(fx.QA + K, a, tau, eo, geom);   % batched over all 30 q, geom reused

    for qi = 1:nQ
        q   = fx.QA(qi,:);
        dq  = dipQ(:,:,:,:,qi);
        dqK = dipQK(:,:,:,:,qi);

        % ---- (i) gauge-predicted tensor vs dip(q+K), M_id per (q,K) --------
        pred = gauge_predict(dq, K, tau, ntau);
        mres = M.mid(pred, dqK);
        verifyTrue(testCase, mres.pass, sprintf( ...
            'Gate-A3: gauge-predicted tensor fails M_id at fx.QA row %d, K=[%g %g %g] (worst_margin=%.3e).', ...
            qi, K(1), K(2), K(3), mres.worst_margin));
        worst_tensor_margin = max(worst_tensor_margin, mres.worst_margin);

        % ---- (ii) sorted cc coupling branches: DIRECT meV tolerance --------
        outqK = ewald_coupling_mapper(q + K, dqK, geomX);
        diffv = abs(outQ{qi}.Jnu - outqK.Jnu);
        allowed = AbsTol_J3 + RelTol_J3*max(abs(outQ{qi}.Jnu), abs(outqK.Jnu));
        margin = max(diffv - allowed);
        verifyLessThanOrEqual(testCase, margin, 0, sprintf( ...
            ['Gate-A3: sorted cc Jnu branches fail the direct meV tolerance at fx.QA row %d, ' ...
             'K=[%g %g %g] (margin=%.3e meV).'], qi, K(1), K(2), K(3), margin));
        worst_coupling_margin = max(worst_coupling_margin, margin);
    end
end
fprintf(['Gate-A3 reciprocal gauge + periodic spectrum: worst tensor M_id margin = %.3e, ' ...
    'worst coupling (direct meV) margin = %.3e over %d q x %d K = %d pairs.\n'], ...
    worst_tensor_margin, worst_coupling_margin, nQ, nK, nQ*nK);
end

function test_gateA4_reduction_and_completeness(testCase)
% Gate-A #4 (docs/invzp_ewald_prereg.md sec 3, item 4): canonical half-open
% reduction + extended-zone gauge restoration + reciprocal-candidate
% completeness.
%   (a) exact q=+0.5 (each axis alone, and all three at once) reduces via the
%       frozen convention K=floor(q+0.5), qbar=q-K EXACTLY to qbar=-0.5 (both
%       the raw arithmetic AND the physical gauge-restored dip, at M_id);
%   (b) all EIGHT frozen near-upper-face corners (fx corners at
%       0.5-fx.delta_q) restore correctly under EVERY one of the 8 frozen K
%       (64 combos total, M_id) -- an earlier draft's three-x-axis-shift-only
%       loop was insufficient;
%   (c) for a representative Gamma/face/edge/corner/generic q, INDEPENDENTLY
%       enumerate (fresh generous bounding box, never touching geom.Ghkl in
%       the enumeration) every integer reciprocal G with |qcart+Gcart|<=g_cut
%       (k~=0) and assert every such G is present in the cached geom.Ghkl
%       candidate union -- explicitly covering G=0 at nonzero q.
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
eo = default_eopts(a);
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();
verifyEqual(testCase, size(fx.QA), [30 3]);   % documented Gamma+q_int(3)+face(6)+edge(12)+corner(8)

[~, ~, geom] = invz_dipole_ewald([0 0 0], a, tau, eo);   % prime ONE geom, reused throughout

% ---- (a) exact +0.5 reduces EXACTLY to -0.5 --------------------------------
axis_hi  = {[0.5 0 0], [0 0.5 0], [0 0 0.5], [0.5 0.5 0.5]};
axis_lo  = {[-0.5 0 0], [0 -0.5 0], [0 0 -0.5], [-0.5 -0.5 -0.5]};
axis_lbl = {'x','y','z','xyz'};
for ci = 1:numel(axis_hi)
    qhi = axis_hi{ci}; qlo = axis_lo{ci}; Kexp = qhi - qlo;
    verifyEqual(testCase, floor(qhi + 0.5), Kexp, 'AbsTol', 0, sprintf( ...
        'Gate-A4a: K=floor(q+0.5) reduction wrong for exact +0.5 (%s).', axis_lbl{ci}));
    verifyEqual(testCase, qhi - floor(qhi + 0.5), qlo, 'AbsTol', 0, sprintf( ...
        'Gate-A4a: qbar=q-K does not reduce exact +0.5 to -0.5 (%s).', axis_lbl{ci}));

    dhi  = invz_dipole_ewald(qhi, a, tau, eo, geom);
    dlo  = invz_dipole_ewald(qlo, a, tau, eo, geom);
    pred = gauge_predict(dlo, Kexp, tau, ntau);
    mres = M.mid(pred, dhi);
    verifyTrue(testCase, mres.pass, sprintf( ...
        'Gate-A4a: exact +0.5->-0.5 (%s) gauge-restored tensor fails M_id (worst_margin=%.3e).', ...
        axis_lbl{ci}, mres.worst_margin));
end

% ---- (b) all EIGHT near-upper-face corners x every frozen K ----------------
corners = fx.QA(23:30,:);
verifyEqual(testCase, size(corners), [8 3]);
dipC = invz_dipole_ewald(corners, a, tau, eo, geom);
worst_corner_margin = -inf;
for ki = 1:size(fx.Kset,1)
    K = fx.Kset(ki,:);
    dipCK = invz_dipole_ewald(corners + K, a, tau, eo, geom);
    for ci = 1:size(corners,1)
        pred = gauge_predict(dipC(:,:,:,:,ci), K, tau, ntau);
        mres = M.mid(pred, dipCK(:,:,:,:,ci));
        verifyTrue(testCase, mres.pass, sprintf( ...
            'Gate-A4b: corner %d fails M_id gauge restoration at K=[%g %g %g] (worst_margin=%.3e).', ...
            ci, K(1), K(2), K(3), mres.worst_margin));
        worst_corner_margin = max(worst_corner_margin, mres.worst_margin);
    end
end
fprintf('Gate-A4b corners x K gauge restoration: worst M_id margin over 8x8=64 combos = %.3e.\n', ...
    worst_corner_margin);

% ---- (c) independent reciprocal candidate completeness ---------------------
reps = struct( ...
    'label', {'Gamma','face','edge','corner','generic'}, ...
    'q',     {fx.QA(1,:), fx.QA(5,:), fx.QA(11,:), fx.QA(23,:), fx.q_int(1,:)});
for ri = 1:numel(reps)
    Gneed = indep_needed_G(reps(ri).q, geom.B, eo.g_cut);
    verifyGreaterThan(testCase, size(Gneed,1), 0, sprintf( ...
        'Gate-A4c: independent enumeration found zero needed G at %s -- check the bounding box.', ...
        reps(ri).label));
    present = ismember(Gneed, geom.Ghkl, 'rows');
    verifyTrue(testCase, all(present), sprintf( ...
        'Gate-A4c: %d of %d independently-enumerated needed G rows are MISSING from geom.Ghkl at %s.', ...
        nnz(~present), numel(present), reps(ri).label));
    isG0needed = any(all(Gneed == 0, 2));
    if any(reps(ri).q ~= 0)
        verifyTrue(testCase, isG0needed, sprintf( ...
            'Gate-A4c: G=0 unexpectedly not needed at nonzero q (%s).', reps(ri).label));
    end
    verifyTrue(testCase, ismember([0 0 0], geom.Ghkl, 'rows'), ...
        'Gate-A4c: G=[0,0,0] must be present in the cached candidate union.');
    fprintf('Gate-A4c completeness (%s): %d needed G, all present in geom.Ghkl (G=0 needed: %d).\n', ...
        reps(ri).label, size(Gneed,1), isG0needed);
end
end

% =====================================================================
% Task 5 helpers -- shared gauge-prediction loop, the test-only coupling
% mapper (reused by a later task's Gate-A#9), and independent per-q
% reciprocal-candidate enumeration
% =====================================================================
function pred = gauge_predict(dq, K, tau, ntau)
% Full extended-zone gauge-predicted tensor:
%   pred_nm = exp(-i*2*pi*K.(tau_m-tau_n)) * dq_nm,   dq = dip(q) [3,3,ntau,ntau].
pred = complex(zeros(3,3,ntau,ntau));
for n = 1:ntau
    for m = 1:ntau
        gph = exp(-1i*2*pi*(K*(tau(m,:) - tau(n,:)).'));
        pred(:,:,n,m) = gph*dq(:,:,n,m);
    end
end
end

function out = ewald_coupling_mapper(q, dipq, geomX)
% ewald_coupling_mapper -- TEST-ONLY shared helper (Step-4 Task 5; the plan's
% "For primitive-level coupling convergence..." paragraph). Mirrors the
% non-ODD Ewald coupling algebra EXACTLY via the production `exchange`
% primitive, WITHOUT calling or modifying invz_jq_modes. Adds NO Lorentz term
% on the Ewald branch: invz_dipole_ewald's returned dip already carries the
% isotropic G=0-omitted term (frozen prereg sec 5 Gate-C decision), so an
% extra Lorentz add-on here would double-count it.
%
% INPUTS
%   q      [1,3] reduced reciprocal (Miller) coordinates -- the SAME raw q
%          used to produce dipq (extended-zone; not pre-reduced to qbar).
%   dipq   [3,3,ntau,ntau] complex Ewald dipolar tensor at this q (a
%          caller-computed single-q slice of invz_dipole_ewald's output --
%          this mapper never calls invz_dipole_ewald itself).
%   geomX  primed q-independent `exchange` geometry, from a ONE-TIME
%          [~, geomX] = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau);
%          reused across the whole q-sweep for speed.
%
% OUTPUT out (struct):
%   Jnu        [4,1] sorted real eigenvalues of the Hermitized cc coupling
%              Jcc(q), in meV.
%   Juni       scalar real uniform-mode projection v'*Jcc*v, v=[1 1 1 1]/2
%              (meV) -- the same convention as this repo's other P.Juni.
%   Jcc, Jaa   [4,4] complex Hermitized cc/aa coupling matrices at q (meV),
%              for a caller that needs more than the sorted branches/uniform
%              mode.
%   Jcc0, Jaa0 the uniform-mode projections of Jcc/Jaa (Jcc0 == Juni),
%              populated ONLY at exact Gamma (q==[0 0 0]); [] otherwise.
ion = invz_ion(); C = invz_const();
ex  = exchange(q, abs(ion.J12), ion.a, ion.tau, geomX);
Jcc = -C.gfac*squeeze(dipq(3,3,:,:)) + sign(ion.J12)*squeeze(ex(3,3,:,:));
Jaa = -C.gfac*squeeze(dipq(1,1,:,:)) + sign(ion.J12)*squeeze(ex(1,1,:,:));
Jcc = (Jcc+Jcc')/2; Jaa = (Jaa+Jaa')/2;      % Hermitize exactly as the caller will
v = [1 1 1 1]/2;
out = struct();
out.Jnu  = sort(real(eig(Jcc)));
out.Juni = real(v*Jcc*v');
out.Jcc  = Jcc; out.Jaa = Jaa;
if all(q == 0)
    out.Jcc0 = out.Juni;
    out.Jaa0 = real(v*Jaa*v');
else
    out.Jcc0 = []; out.Jaa0 = [];
end
end

function Gneed = indep_needed_G(qrow, B, g_cut)
% Independent per-q enumeration of every integer reciprocal G=[h,k,l] with
% |qcart+Gcart|<=g_cut, k~=0 (Gcart=G*B), from a FRESH generous bounding box
% -- geom.Ghkl/geom.Gcart are never read here, so this cannot be tautological
% against the primitive's own cached candidate union.
K = floor(qrow + 0.5); qbar = qrow - K;      % same frozen extended-zone convention
qcart = qbar*B;
sb = min(svd(B));
nmax = ceil((g_cut + norm(qcart))/sb) + 1;   % generous conservative box, +1 slack
rng = -nmax:nmax;
[H, Kk, L] = ndgrid(rng, rng, rng);
Ghkl_all  = [H(:) Kk(:) L(:)];
Gcart_all = Ghkl_all*B;
k  = Gcart_all + qcart;
kk = sum(k.^2, 2);
keep = (kk <= g_cut^2) & (kk > 0);
Gneed = Ghkl_all(keep,:);
end

% =====================================================================
% Task 6 -- Gate-A test 5 (structural identities: conjugation/Hermiticity,
% Gamma-realness, common origin-shift raw invariance, per-representative
% Bravais-shift raw invariance, and cell-phase-gauge covariance)
% =====================================================================
function test_gateA5_hermiticity_and_gamma_realness(testCase)
% Gate-A #5 part 1/2 (docs/invzp_ewald_prereg.md sec 3, item 5): for EVERY
% q in fx.QA (30 points: Gamma + 3 generic + 6 face + 12 edge + 8 corner),
% the returned complete tensor is pair-conjugate/Hermitian,
% dip(:,:,m,n)==conj(dip(:,:,n,m)), at M_id -- called PER q on that q's own
% complete [3,3,ntau,ntau] slice (never on the stacked 5-D array, which
% would silently loosen T_scale to the worst component over the whole
% stack). At EXACT Gamma the complete tensor, INCLUDING off-site blocks,
% must additionally be real at M_id (max|imag|<=1e-12*T_scale): tested via
% M.mid(dip0, real(dip0)), whose |A-B| is exactly |imag(dip0)| componentwise
% and whose T_scale is exactly max(|dip0|) (since |real(dip0)|<=|dip0|
% pointwise), i.e. precisely the frozen realness criterion.
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
eo = default_eopts(a);
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();

dipQ = invz_dipole_ewald(fx.QA, a, tau, eo);   % ONE geom (internal), ONE batched call
nQ = size(fx.QA,1);

worst_herm_margin = -inf;
for qi = 1:nQ
    A = dipQ(:,:,:,:,qi);
    herm = conj_transpose_pred(A, ntau);
    mres = M.mid(A, herm);
    verifyTrue(testCase, mres.pass, sprintf( ...
        'Gate-A5: pair conjugation/Hermiticity fails M_id at fx.QA row %d (worst_margin=%.3e).', ...
        qi, mres.worst_margin));
    worst_herm_margin = max(worst_herm_margin, mres.worst_margin);
end
fprintf('Gate-A5 Hermiticity: worst M_id margin over %d fx.QA points = %.3e.\n', nQ, worst_herm_margin);

gidx = find(all(fx.QA == 0, 2));
verifyEqual(testCase, numel(gidx), 1, 'Gate-A5: fx.QA must contain exactly one exact-Gamma row.');
dip0 = dipQ(:,:,:,:,gidx);
mreal = M.mid(dip0, real(dip0));
verifyTrue(testCase, mreal.pass, sprintf( ...
    ['Gate-A5: exact-Gamma complete tensor (incl. off-site blocks) fails realness at M_id ' ...
     '(worst_margin=%.3e).'], mreal.worst_margin));
fprintf('Gate-A5 Gamma-realness: M_id margin = %.3e.\n', mreal.worst_margin);
end

function test_gateA5_origin_shift_invariance(testCase)
% Gate-A #5 part 3/5: applying the frozen COMMON fractional-origin shift
% fx.origin_shift (added to EVERY representative alike) leaves the complete
% RAW (total-displacement-gauge) tensor invariant at M_id, for every q in
% fx.QA. Shifting every tau row by the SAME vector leaves every pairwise
% difference tau_m-tau_n -- hence every real-space displacement, reciprocal
% phase, and extended-zone gauge-restore phase the primitive computes --
% UNCHANGED up to floating round-off; this is a genuine structural
% invariance of the primitive, not a numerical accident. tau's fingerprint
% differs from the baseline, so a FRESH geometry is built (geom reuse across
% a tau change would raise invz:ewaldGeomReuse).
ion = invz_ion(); a = ion.a; tau = ion.tau;
eo = default_eopts(a);
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();

dipQ  = invz_dipole_ewald(fx.QA, a, tau, eo);                 % baseline (own geom)
tau_os = tau + fx.origin_shift;                               % common shift, ALL reps
dipOS  = invz_dipole_ewald(fx.QA, a, tau_os, eo);              % fresh geom (fingerprint differs)

nQ = size(fx.QA,1);
worst_margin = -inf;
for qi = 1:nQ
    mres = M.mid(dipOS(:,:,:,:,qi), dipQ(:,:,:,:,qi));
    verifyTrue(testCase, mres.pass, sprintf( ...
        'Gate-A5: common origin-shift breaks raw invariance at fx.QA row %d (worst_margin=%.3e).', ...
        qi, mres.worst_margin));
    worst_margin = max(worst_margin, mres.worst_margin);
end
fprintf('Gate-A5 common origin-shift raw invariance: worst M_id margin over %d fx.QA points = %.3e.\n', ...
    nQ, worst_margin);
end

function test_gateA5_bravais_shift_invariance_and_cellphase_covariance(testCase)
% Gate-A #5 part 4/5 (the load-bearing sign check): for EVERY basis
% representative rep and BOTH frozen Bravais shifts L (fx.bravais_shifts),
% rebuild geometry with tau2 = tau; tau2(rep,:) = tau2(rep,:) + L (an
% INTEGER lattice vector -- this re-indexes which integer real cell counts
% as "R=0" for that one representative while leaving every other
% representative untouched). Two requirements, for every q in fx.QA:
%   (a) RAW (total-displacement-gauge) invariance at M_id vs the unshifted
%       dip -- the shift only relabels which integer cell each real-space
%       term is enumerated under (the retained {x=R+d} SET is unchanged) and
%       multiplies every reciprocal phase by exp(+/-i*2*pi*(integer
%       Ghkl).(integer L)) == 1 EXACTLY (reciprocal-lattice duality), so the
%       returned tensor is unchanged;
%   (b) in the CELL-PHASE gauge Dcell_nm=exp(+i*2*pi*q.(tau_m-tau_n))*dip_nm
%       (each result converted using its OWN tau -- unshifted tau for the
%       baseline, tau2 for the shifted result), the shifted result obeys
%       Dcell_shifted = U' * Dcell_unshifted * U at M_id, with U diagonal in
%       sublattice, U(rep,rep)=exp(+i*2*pi*q.L), U(k,k)=1 for k~=rep.
%       *** U' is on the LEFT: U'*Dcell*U, NOT U*Dcell*U' -- this exact sign
%       was corrected and verified numerically against MF_dipole at freeze;
%       the reverse fails by O(1) (8.2e-2 vs 1.0e-6). See
%       docs/invzp_ewald_prereg.md sec 3 item 5. ***
% Every M_id call is made PER (q,rep,L) on that combination's own complete
% [3,3,ntau,ntau] tensor, built in full before the call (never a stacked
% array, never an ad hoc block-relative tolerance).
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
eo = default_eopts(a);
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();

dipQ = invz_dipole_ewald(fx.QA, a, tau, eo);   % baseline (own geom)
nQ = size(fx.QA,1);

worst_raw_margin = -inf;
worst_cov_margin = -inf;
for rep = 1:ntau
    for li = 1:numel(fx.bravais_shifts)
        L = fx.bravais_shifts{li};
        tau2 = tau;
        tau2(rep,:) = tau2(rep,:) + L;
        dip2 = invz_dipole_ewald(fx.QA, a, tau2, eo);   % fresh geom (re-indexed real cells)

        for qi = 1:nQ
            q  = fx.QA(qi,:);
            d0 = dipQ(:,:,:,:,qi);
            d2 = dip2(:,:,:,:,qi);

            % ---- (a) raw total-displacement-gauge invariance ---------------
            mraw = M.mid(d2, d0);
            verifyTrue(testCase, mraw.pass, sprintf( ...
                ['Gate-A5: per-representative Bravais shift (rep=%d, L=[%g %g %g]) breaks raw ' ...
                 'invariance at fx.QA row %d (worst_margin=%.3e).'], ...
                rep, L(1), L(2), L(3), qi, mraw.worst_margin));
            worst_raw_margin = max(worst_raw_margin, mraw.worst_margin);

            % ---- (b) cell-phase covariance, U' * Dcell * U (U' on the LEFT) --
            Dcell0 = to_cellphase(d0, q, tau,  ntau);    % unshifted result, its OWN tau
            Dcell2 = to_cellphase(d2, q, tau2, ntau);    % shifted result, its OWN tau2
            pred   = covariance_pred(Dcell0, q, L, rep, ntau);
            mcov = M.mid(Dcell2, pred);
            verifyTrue(testCase, mcov.pass, sprintf( ...
                ['Gate-A5: cell-phase covariance U''*Dcell*U fails (rep=%d, L=[%g %g %g]) at ' ...
                 'fx.QA row %d (worst_margin=%.3e).'], ...
                rep, L(1), L(2), L(3), qi, mcov.worst_margin));
            worst_cov_margin = max(worst_cov_margin, mcov.worst_margin);
        end
    end
end
fprintf(['Gate-A5 Bravais-shift raw invariance + cell-phase covariance: worst raw M_id margin = ' ...
    '%.3e, worst covariance M_id margin = %.3e, over %d reps x %d shifts x %d fx.QA points.\n'], ...
    worst_raw_margin, worst_cov_margin, ntau, numel(fx.bravais_shifts), nQ);
end

% =====================================================================
% Task 6 helpers -- pair-conjugation prediction, cell-phase gauge conversion,
% and the U'*Dcell*U covariance prediction (U' on the LEFT -- see the
% docstring of test_gateA5_bravais_shift_invariance_and_cellphase_covariance)
% =====================================================================
function herm = conj_transpose_pred(A, ntau)
% Predicted pair-conjugate tensor: herm(:,:,n,m) = conj(A(:,:,m,n)). Every
% block is built independently from A -- no conjugate block is ever filled
% by copying the block it will be compared against.
herm = complex(zeros(3,3,ntau,ntau));
for n = 1:ntau
    for m = 1:ntau
        herm(:,:,n,m) = conj(A(:,:,m,n));
    end
end
end

function Dcell = to_cellphase(dipslice, q, tau, ntau)
% Dcell_nm = exp(+i*2*pi*q.(tau_m-tau_n)) * dip_nm, using the CALLER-supplied
% tau (each result's OWN tau -- unshifted or per-representative-shifted).
Dcell = complex(zeros(3,3,ntau,ntau));
for n = 1:ntau
    for m = 1:ntau
        ph = exp(1i*2*pi*(q*(tau(m,:) - tau(n,:)).'));
        Dcell(:,:,n,m) = ph*dipslice(:,:,n,m);
    end
end
end

function pred = covariance_pred(Dcell0, q, L, rep, ntau)
% pred = U' * Dcell0 * U, U diagonal in sublattice with
% U(rep,rep)=exp(+i*2*pi*q.L), U(k,k)=1 for k~=rep. Block (n,m) scales by
% conj(U(n,n))*U(m,m): U' (conjugate-transpose) acts on the LEFT/row index
% n, U acts on the RIGHT/column index m. THE SIGN IS U'*Dcell*U, NOT
% U*Dcell*U' -- reversing it fails by O(1) (see the caller's docstring).
Udiag = ones(1,ntau);
Udiag(rep) = exp(1i*2*pi*(q*L.'));
pred = complex(zeros(3,3,ntau,ntau));
for n = 1:ntau
    for m = 1:ntau
        pred(:,:,n,m) = conj(Udiag(n))*Udiag(m)*Dcell0(:,:,n,m);
    end
end
end

% =====================================================================
% Task 7 -- Gate-A test 9 (alpha-bracket independence + separate-axis
% cutoff-ladder convergence of the raw tensor at M_T and the derived
% couplings at M_J, plus the default-vs-joint-refinement >=3-orders-margin
% check). This is the LAST Gate-A test. The frozen config set (prereg sec 2/
% sec 3 item 9) has 13 members in four groups, ALL evaluated at EVERY q in
% fx.QA (30 points), reusing ONE geom per config:
%   (1) alpha bracket:      {0.6,0.8,1.0,1.2,1.5}*alpha0, alpha-matched
%                            cutoffs r_cut=6.5/alpha, g_cut=13*alpha    [5]
%   (2) real-axis ladder:   C_r in {4.5,5.0,5.5} at fixed C_g=13        [3]
%   (3) recip-axis ladder:  C_g in {9,10,11}    at fixed C_r=6.5        [3]
%   (4) default (5.5,11) vs joint refinement (6.0,12)                  [2]
% =====================================================================
function test_gateA9_alpha_bracket_independence(testCase)
% Gate-A #9 part 1/3 (docs/invzp_ewald_prereg.md sec 3 item 9 / sec 2):
% ALPHA-BRACKET independence. Multipliers {0.6,0.8,1.0,1.2,1.5}*alpha0, each
% with the alpha-MATCHED generous cutoffs r_cut=6.5/alpha, g_cut=13*alpha
% (sec 2 "For alpha-independence..."), so the DIMENSIONLESS truncation
% bounds (alpha*r_cut=6.5, g_cut/(2*alpha)=6.5) stay FIXED across the whole
% bracket -- only the physical splitting point alpha itself changes. ALL 10
% pairwise combinations of the 5 members are compared, at EVERY q in fx.QA
% (not just q_int(1,:)), via M_T (raw complete tensor, per q) and M_J (the
% shared Task-5 ewald_coupling_mapper's Jnu/Juni/Jcc0/Jaa0, per q) -- full
% pairwise rather than each-vs-one-reference, so a cancellation specific to
% any single pair is still caught and localized (retains the pairwise
% diagnostics the frozen acceptance wording calls for).
ion = invz_ion(); a = ion.a; tau = ion.tau;
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();
verifyEqual(testCase, size(fx.QA), [30 3]);              % Gamma+q_int(3)+face(6)+edge(12)+corner(8)
[~, geomX] = exchange([0 0 0], abs(ion.J12), a, tau);     % prime exchange geom ONCE, reused everywhere

mult = [0.6 0.8 1.0 1.2 1.5];
nA = numel(mult);
cfgs = cell(nA,1);
tags = cell(nA,1);
for ai = 1:nA
    alpha = mult(ai)*fx.alpha0;
    cfgs{ai} = gateA9_eval_config(a, tau, eopts_alpha_matched(alpha), fx, geomX);
    tags{ai} = sprintf('alpha=x%.2g', mult(ai));
    fprintf('Gate-A9 alpha-bracket config %s (alpha=%.6g): elapsed %.3fs (diagnostic only).\n', ...
        tags{ai}, alpha, cfgs{ai}.elapsed);
end

worst_T_margin = -inf; worst_T_ratio = 0; worst_J_margin = -inf; worst_J_ratio = 0;
npairs = 0;
for i = 1:nA
    for j = i+1:nA
        res = gateA9_pair_check(testCase, M, fx, cfgs{i}, cfgs{j}, tags{i}, tags{j}, 'alpha-bracket');
        worst_T_margin = max(worst_T_margin, res.worst_T_margin);
        worst_T_ratio  = max(worst_T_ratio,  res.worst_T_ratio);
        worst_J_margin = max(worst_J_margin, res.worst_J_margin);
        worst_J_ratio  = max(worst_J_ratio,  res.worst_J_ratio);
        npairs = npairs + 1;
    end
end
fprintf(['Gate-A9 alpha-bracket independence: %d pairwise combos over %d members x %d fx.QA points -- ' ...
    'worst M_T margin=%.3e (ratio %.3e), worst M_J margin=%.3e (ratio %.3e).\n'], ...
    npairs, nA, size(fx.QA,1), worst_T_margin, worst_T_ratio, worst_J_margin, worst_J_ratio);
end

function test_gateA9_cutoff_ladders_separate_axes(testCase)
% Gate-A #9 part 2/3 (docs/invzp_ewald_prereg.md sec 3 item 9 / sec 2):
% SEPARATE-AXIS cutoff-ladder convergence, alpha=alpha0 throughout. TWO
% independent one-axis ladders, evaluated and pairwise-compared entirely
% separately -- the two cutoff axes are NEVER combined into one sweep here
% (that is part 3/3 below, the default-vs-joint-refinement pair, the ONE
% place both axes move together by design):
%   real axis:       C_r in {4.5,5.0,5.5} at FIXED generous C_g=13;
%   reciprocal axis: C_g in {9,10,11}    at FIXED generous C_r=6.5.
% Within each ladder, ALL 3 pairwise combinations of its 3 rungs (coarsest-
% mid, mid-finest, coarsest-finest) are compared at every q in fx.QA via M_T
% and M_J -- so the coarsest-to-finest convergence bracket is asserted
% directly (coarsest-vs-finest is itself one of the pairs), not merely
% inferred from adjacent-rung agreement alone.
ion = invz_ion(); a = ion.a; tau = ion.tau;
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();
[~, geomX] = exchange([0 0 0], abs(ion.J12), a, tau);
alpha0 = fx.alpha0;

Cr_list = [4.5 5.0 5.5];
eo_real  = arrayfun(@(Cr) eopts_Cr_Cg(alpha0, Cr, 13), Cr_list, 'UniformOutput', false);
tag_real = arrayfun(@(Cr) sprintf('C_r=%.1f,C_g=13', Cr), Cr_list, 'UniformOutput', false);
gateA9_run_ladder(testCase, M, fx, a, tau, geomX, eo_real, tag_real, 'real-axis ladder');

Cg_list = [9 10 11];
eo_recip  = arrayfun(@(Cg) eopts_Cr_Cg(alpha0, 6.5, Cg), Cg_list, 'UniformOutput', false);
tag_recip = arrayfun(@(Cg) sprintf('C_r=6.5,C_g=%.0f', Cg), Cg_list, 'UniformOutput', false);
gateA9_run_ladder(testCase, M, fx, a, tau, geomX, eo_recip, tag_recip, 'reciprocal-axis ladder');
end

function test_gateA9_default_vs_joint_refinement(testCase)
% Gate-A #9 part 3/3 (docs/invzp_ewald_prereg.md sec 3 item 9 / sec 2):
% production default (C_r,C_g)=(5.5,11) versus joint refinement (6.0,12),
% alpha=alpha0. This is the ONE Gate-A9 pair where BOTH cutoff axes move
% together (by design -- the "joint" refinement), in contrast to the two
% one-axis ladders above which never combine axes. Beyond the ordinary M_T/
% M_J PASS requirement (asserted inside gateA9_pair_check, per q, over every
% fx.QA point), the default-vs-joint pair must additionally retain the
% frozen >=3-orders-of-magnitude margin: WORST_RATIO (the fraction of the
% AbsTol+RelTol*max(|A|,|B|) budget actually used) for BOTH the raw M_T
% comparison and the coupling M_J comparison must be <=1e-3 (prereg sec 3
% item 9: "the default must retain the calibrated >=3-orders margin
% relative to the controlling tolerance"; sec 1 records the calibration
% default self-convergence residual as raw 3.6e-13 / coupling 2.1e-13 --
% several more orders below even this 1e-3 utilization requirement).
ion = invz_ion(); a = ion.a; tau = ion.tau;
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();
[~, geomX] = exchange([0 0 0], abs(ion.J12), a, tau);
alpha0 = fx.alpha0;

eo_default = eopts_Cr_Cg(alpha0, 5.5, 11);   % production default (matches default_eopts(a))
eo_joint   = eopts_Cr_Cg(alpha0, 6.0, 12);   % joint refinement -- BOTH axes move together

cfgD = gateA9_eval_config(a, tau, eo_default, fx, geomX);
cfgJ = gateA9_eval_config(a, tau, eo_joint,   fx, geomX);
fprintf(['Gate-A9 default-vs-joint config elapsed: default(5.5,11) %.3fs, joint(6.0,12) %.3fs ' ...
    '(diagnostic only).\n'], cfgD.elapsed, cfgJ.elapsed);

res = gateA9_pair_check(testCase, M, fx, cfgD, cfgJ, 'default(5.5,11)', 'joint(6.0,12)', ...
    'default-vs-joint');

verifyLessThanOrEqual(testCase, res.worst_T_ratio, 1e-3, sprintf( ...
    ['Gate-A9: default-vs-joint RAW (M_T) controlling-tolerance utilization exceeds the frozen ' ...
     '3-orders margin (worst_ratio=%.3e, required <=1e-3).'], res.worst_T_ratio));
verifyLessThanOrEqual(testCase, res.worst_J_ratio, 1e-3, sprintf( ...
    ['Gate-A9: default-vs-joint COUPLING (M_J) controlling-tolerance utilization exceeds the frozen ' ...
     '3-orders margin (worst_ratio=%.3e, required <=1e-3).'], res.worst_J_ratio));

fprintf(['Gate-A9 default-vs-joint refinement: worst M_T margin=%.3e (ratio %.3e, <=1e-3 required), ' ...
    'worst M_J margin=%.3e (ratio %.3e, <=1e-3 required).\n'], ...
    res.worst_T_margin, res.worst_T_ratio, res.worst_J_margin, res.worst_J_ratio);
end

% =====================================================================
% Task 7 helpers -- alpha/cutoff config builders, the per-config batched
% evaluator (ONE geom, all 30 fx.QA points, per-q coupling via the shared
% Task-5 ewald_coupling_mapper), and the shared pairwise M_T/M_J comparator.
% LOAD-BEARING: gateA9_pair_check calls M.mt/M.mj PER Q on that single q's
% own complete tensor / coupling row -- NEVER on a stacked multi-q array
% (which would silently loosen M_T's T_scale to the worst component over
% the whole stack).
% =====================================================================
function eo = eopts_alpha_matched(alpha)
% alpha-matched generous cutoffs r_cut=6.5/alpha, g_cut=13*alpha (prereg
% sec 2), so the dimensionless truncation bounds alpha*r_cut=6.5 and
% g_cut/(2*alpha)=6.5 stay fixed -- both comfortably above the GUARD=4.5
% floor for every frozen multiplier.
eo = mk_eopts(alpha, 6.5/alpha, 13*alpha, 'conducting_k0_omitted');
end

function eo = eopts_Cr_Cg(alpha, Cr, Cg)
% General (C_r,C_g) -> eopts map, r_cut=C_r/alpha, g_cut=C_g*alpha (prereg
% sec 2 "r_cut=C_r/alpha0, g_cut=C_g*alpha0"), the same pattern as this
% file's own default_eopts (5.5/a0, 11*a0 for (C_r,C_g)=(5.5,11)).
eo = mk_eopts(alpha, Cr/alpha, Cg*alpha, 'conducting_k0_omitted');
end

function out = gateA9_eval_config(a, tau, eo, fx, geomX)
% Evaluate ONE Ewald configuration over EVERY q in fx.QA (30 points): build
% geom ONCE via a single batched invz_dipole_ewald call over all 30 q, then
% run the shared Task-5 ewald_coupling_mapper per q, reusing the
% q-independent exchange geometry geomX. Returns:
%   dipQ  [3,3,ntau,ntau,30] complex raw Ewald tensor, this config's OWN geom
%   Jnu   [30,4] sorted real cc eigenvalues (meV), one row per fx.QA point
%   Juni  [30,1] real uniform-mode projection v'*Jcc*v (meV)
%   Jcc0, Jaa0   scalar real uniform-mode couplings AT EXACT GAMMA (meV) --
%                fx.QA(1,:) is always [0 0 0] so these are always populated
%   elapsed      wall time for this whole config (DIAGNOSTIC ONLY -- prereg
%                sec 2's "operational timing target" is explicitly not a
%                numerical pass/fail condition; never asserted on below)
t0 = tic;
dipQ = invz_dipole_ewald(fx.QA, a, tau, eo);
nQ = size(fx.QA,1);
Jnu = zeros(nQ,4); Juni = zeros(nQ,1); Jcc0 = []; Jaa0 = [];
for qi = 1:nQ
    r = ewald_coupling_mapper(fx.QA(qi,:), dipQ(:,:,:,:,qi), geomX);
    Jnu(qi,:) = r.Jnu(:).';
    Juni(qi)  = r.Juni;
    if ~isempty(r.Jcc0)
        Jcc0 = r.Jcc0; Jaa0 = r.Jaa0;
    end
end
out = struct('dipQ', dipQ, 'Jnu', Jnu, 'Juni', Juni, 'Jcc0', Jcc0, 'Jaa0', Jaa0, 'elapsed', toc(t0));
end

function res = gateA9_pair_check(testCase, M, fx, cfgA, cfgB, tagA, tagB, ctx)
% Pairwise Gate-A9 comparison between two EVALUATED configurations (from
% gateA9_eval_config): M_T on the RAW complete tensor and M_J on
% Jnu/Juni/Jcc0/Jaa0 (the shared coupling mapper's four required outputs).
% Every pair gets its OWN per-q assertions (never collapsed into one
% aggregate check), so a cancellation between any two members is still
% caught and localized to the (pair, q) that produced it.
nQ = size(fx.QA,1);
verifyTrue(testCase, ~isempty(cfgA.Jcc0) && ~isempty(cfgB.Jcc0), sprintf( ...
    'Gate-A9 (%s): Jcc0 unexpectedly empty for %s or %s (fx.QA(1,:) must be exact Gamma).', ...
    ctx, tagA, tagB));

worst_T_margin = -inf; worst_T_ratio = 0;
worst_J_margin = -inf; worst_J_ratio = 0;
for qi = 1:nQ
    A = cfgA.dipQ(:,:,:,:,qi); B = cfgB.dipQ(:,:,:,:,qi);
    mT = M.mt(A, B);
    verifyTrue(testCase, mT.pass, sprintf( ...
        'Gate-A9 (%s): M_T raw-tensor comparison fails between %s and %s at fx.QA row %d (worst_margin=%.3e).', ...
        ctx, tagA, tagB, qi, mT.worst_margin));
    worst_T_margin = max(worst_T_margin, mT.worst_margin);
    worst_T_ratio  = max(worst_T_ratio,  mT.worst_ratio);

    mJnu  = M.mj(cfgA.Jnu(qi,:), cfgB.Jnu(qi,:));
    mJuni = M.mj(cfgA.Juni(qi),  cfgB.Juni(qi));
    verifyTrue(testCase, mJnu.pass, sprintf( ...
        'Gate-A9 (%s): M_J Jnu comparison fails between %s and %s at fx.QA row %d (worst_margin=%.3e).', ...
        ctx, tagA, tagB, qi, mJnu.worst_margin));
    verifyTrue(testCase, mJuni.pass, sprintf( ...
        'Gate-A9 (%s): M_J Juni comparison fails between %s and %s at fx.QA row %d (worst_margin=%.3e).', ...
        ctx, tagA, tagB, qi, mJuni.worst_margin));
    worst_J_margin = max([worst_J_margin, mJnu.worst_margin, mJuni.worst_margin]);
    worst_J_ratio  = max([worst_J_ratio,  mJnu.worst_ratio,  mJuni.worst_ratio]);
end

mJcc0 = M.mj(cfgA.Jcc0, cfgB.Jcc0);
mJaa0 = M.mj(cfgA.Jaa0, cfgB.Jaa0);
verifyTrue(testCase, mJcc0.pass, sprintf( ...
    'Gate-A9 (%s): M_J Jcc0 (Gamma) comparison fails between %s and %s (worst_margin=%.3e).', ...
    ctx, tagA, tagB, mJcc0.worst_margin));
verifyTrue(testCase, mJaa0.pass, sprintf( ...
    'Gate-A9 (%s): M_J Jaa0 (Gamma) comparison fails between %s and %s (worst_margin=%.3e).', ...
    ctx, tagA, tagB, mJaa0.worst_margin));
worst_J_margin = max([worst_J_margin, mJcc0.worst_margin, mJaa0.worst_margin]);
worst_J_ratio  = max([worst_J_ratio,  mJcc0.worst_ratio,  mJaa0.worst_ratio]);

fprintf(['Gate-A9 (%s) pair [%s] vs [%s]: worst M_T margin=%.3e (ratio %.3e), ' ...
    'worst M_J margin=%.3e (ratio %.3e) over %d fx.QA points.\n'], ...
    ctx, tagA, tagB, worst_T_margin, worst_T_ratio, worst_J_margin, worst_J_ratio, nQ);

res = struct('worst_T_margin', worst_T_margin, 'worst_T_ratio', worst_T_ratio, ...
             'worst_J_margin', worst_J_margin, 'worst_J_ratio', worst_J_ratio);
end

function gateA9_run_ladder(testCase, M, fx, a, tau, geomX, eo_list, tag_list, ctx)
% Shared driver for ONE one-axis cutoff ladder: evaluate every rung, then
% pairwise-compare ALL combinations of that ladder's own rungs. The caller
% passes exactly one axis' eo_list per call, so this never crosses into the
% other axis' rungs (the "WITHOUT combining the two cutoff axes" wording).
n = numel(eo_list);
cfgs = cell(n,1);
for i = 1:n
    cfgs{i} = gateA9_eval_config(a, tau, eo_list{i}, fx, geomX);
    fprintf('Gate-A9 %s config %s: elapsed %.3fs (diagnostic only).\n', ctx, tag_list{i}, cfgs{i}.elapsed);
end
worst_T_margin = -inf; worst_T_ratio = 0; worst_J_margin = -inf; worst_J_ratio = 0;
npairs = 0;
for i = 1:n
    for j = i+1:n
        res = gateA9_pair_check(testCase, M, fx, cfgs{i}, cfgs{j}, tag_list{i}, tag_list{j}, ctx);
        worst_T_margin = max(worst_T_margin, res.worst_T_margin);
        worst_T_ratio  = max(worst_T_ratio,  res.worst_T_ratio);
        worst_J_margin = max(worst_J_margin, res.worst_J_margin);
        worst_J_ratio  = max(worst_J_ratio,  res.worst_J_ratio);
        npairs = npairs + 1;
    end
end
fprintf(['Gate-A9 %s coarsest-to-finest bracket: %d pairwise combos over %d rungs x %d fx.QA points -- ' ...
    'worst M_T margin=%.3e (ratio %.3e), worst M_J margin=%.3e (ratio %.3e).\n'], ...
    ctx, npairs, n, size(fx.QA,1), worst_T_margin, worst_T_ratio, worst_J_margin, worst_J_ratio);
end
