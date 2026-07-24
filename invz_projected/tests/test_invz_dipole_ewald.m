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
