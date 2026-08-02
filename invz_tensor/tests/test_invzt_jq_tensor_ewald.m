function tests = test_invzt_jq_tensor_ewald
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
tc.TestData.ion = invz_ion();
tc.TestData.eopts = invzt_ewald_defaults(tc.TestData.ion);
end

function test_certified_default_controls(tc)
ion = tc.TestData.ion;
eopts = invzt_ewald_defaults(ion);
alpha0 = sqrt(pi)/ion.Vc^(1/3);
verifyEqual(tc, eopts.alpha, alpha0, 'AbsTol', 0);
verifyEqual(tc, eopts.r_cut, 5.5/alpha0, 'AbsTol', 0);
verifyEqual(tc, eopts.g_cut, 11*alpha0, 'AbsTol', 0);
verifyEqual(tc, eopts.boundary, 'conducting_k0_omitted');
verifyError(tc, @() invzt_ewald_defaults(struct()), 'invzt:ewaldDefaultsIon');
end

function test_direct_primitive_assembly_and_gamma_convention(tc)
ion = tc.TestData.ion;
eopts = tc.TestData.eopts;
qvec = [0 0 0; 0.17 -0.23 0.11];
lat = invzt_jq_tensor(ion, qvec, struct('dipole', 'ewald', ...
    'ewald', eopts, 'cache', false));

C = invz_const();
expected = complex(zeros(12,12,size(qvec,1)));
for iq = 1:size(qvec,1)
    dip = invz_dipole_ewald(qvec(iq,:), ion.a, ion.tau, eopts);
    ex = exchange(qvec(iq,:), abs(ion.J12), ion.a, ion.tau);
    for mu = 1:3
        for nu = 1:3
            blk = -C.gfac*squeeze(dip(mu,nu,:,:)) ...
                + sign(ion.J12)*squeeze(ex(mu,nu,:,:));
            if mu == nu, blk = (blk + blk')/2; end
            expected(mu:3:12,nu:3:12,iq) = blk;
        end
    end
end

verifyEqual(tc, lat.Jt, expected, 'AbsTol', 2e-13);
verifyEqual(tc, lat.JtGamma, expected(:,:,1), 'AbsTol', 2e-13);
for iq = 1:size(qvec,1)
    verifyLessThan(tc, norm(lat.Jt(:,:,iq) - lat.Jt(:,:,iq)', 'fro'), 2e-12);
end
verifyEqual(tc, lat.info.dipole.backend, 'ewald');
verifyEqual(tc, lat.info.dipole.ewald, eopts);
verifyTrue(tc, isnan(lat.info.dpRng));

% The Ewald primitive is already the regularized Gamma tensor.  Its uniform
% cc projection therefore receives no wrapper +lorz in either JtGamma or info.
v = ones(4,1)/2;
cc = 3:3:12;
dip0 = invz_dipole_ewald([0 0 0], ion.a, ion.tau, eopts);
jcc0d = -C.gfac*squeeze(dip0(3,3,:,:));
jcc0d = (jcc0d + jcc0d')/2;
verifyEqual(tc, lat.info.Jcc0_dipole, real(v'*jcc0d*v), 'AbsTol', 2e-13);
verifyEqual(tc, lat.info.Jcc0, real(v'*lat.JtGamma(cc,cc)*v), 'AbsTol', 2e-13);
end

function test_parts_difference_is_exchange_only_for_ewald(tc)
ion = tc.TestData.ion;
eopts = tc.TestData.eopts;
qvec = [0 0 0; 0.09 0.14 -0.21];
full = invzt_jq_tensor(ion, qvec, struct('dipole', 'ewald', ...
    'ewald', eopts, 'cache', false));
diponly = invzt_jq_tensor(ion, qvec, struct('dipole', 'ewald', ...
    'ewald', eopts, 'parts', 'dipole', 'cache', false));

expected = complex(zeros(size(full.Jt)));
for iq = 1:size(qvec,1)
    ex = exchange(qvec(iq,:), abs(ion.J12), ion.a, ion.tau);
    for mu = 1:3
        for nu = 1:3
            blk = sign(ion.J12)*squeeze(ex(mu,nu,:,:));
            if mu == nu, blk = (blk + blk')/2; end
            expected(mu:3:12,nu:3:12,iq) = blk;
        end
    end
end
verifyEqual(tc, full.Jt - diponly.Jt, expected, 'AbsTol', 2e-13);
end

function test_reciprocal_shift_preserves_full_page_spectrum(tc)
ion = tc.TestData.ion;
eopts = tc.TestData.eopts;
q = [0.13 -0.22 0.19];
lat = invzt_jq_tensor(ion, [q; q + [1 -1 2]], struct('dipole', 'ewald', ...
    'ewald', eopts, 'cache', false));
e1 = sort(real(eig(lat.Jt(:,:,1))));
e2 = sort(real(eig(lat.Jt(:,:,2))));
verifyEqual(tc, e2, e1, 'AbsTol', 2e-11);
end

function test_backend_option_validation(tc)
ion = tc.TestData.ion;
eopts = tc.TestData.eopts;
verifyError(tc, @() invzt_jq_tensor(ion, [0 0 0], ...
    struct('dipole', 'unknown', 'cache', false)), 'invzt:jqTensorBackend');
verifyError(tc, @() invzt_jq_tensor(ion, [0 0 0], ...
    struct('dipole', 'bruteforce', 'ewald', eopts, 'cache', false)), ...
    'invzt:jqTensorEwaldOptsUnexpected');
verifyError(tc, @() invzt_jq_tensor(ion, [0 0 0], ...
    struct('dpRng', 6, 'cache', false)), 'invzt:jqTensorDpRngUnexpected');
bad = eopts; bad.boundary = 'vacuum';
verifyError(tc, @() invzt_jq_tensor(ion, [0 0 0], ...
    struct('dipole', 'ewald', 'ewald', bad, 'cache', false)), ...
    'invzt:jqTensorEwaldBoundary');
end

function test_absent_backend_is_certified_ewald_and_brute_is_explicit(tc)
ion = tc.TestData.ion;
qvec = [0 0 0; 0.12 -0.08 0.19];
implicit = invzt_jq_tensor(ion, qvec, struct('cache', false));
explicit = invzt_jq_tensor(ion, qvec, struct('dipole', 'ewald', ...
    'ewald', invzt_ewald_defaults(ion), 'cache', false));
verifyEqual(tc, implicit.Jt, explicit.Jt, 'AbsTol', 0);
verifyEqual(tc, implicit.JtGamma, explicit.JtGamma, 'AbsTol', 0);
verifyEqual(tc, implicit.info.dipole.backend, 'ewald');
verifyEqual(tc, implicit.info.dipole.ewald, invzt_ewald_defaults(ion));

brute = invzt_jq_tensor(ion, qvec, struct('dipole', 'bruteforce', ...
    'dpRng', 6, 'cache', false));
verifyEqual(tc, brute.info.dipole.backend, 'bruteforce');
verifyEqual(tc, brute.info.dpRng, 6);
verifyEqual(tc, brute.info.dipole.ewald, ...
    struct('alpha', [], 'r_cut', [], 'g_cut', [], 'boundary', ''));
end

function test_cache_miss_hit_and_corruption_recompute(tc)
ion = tc.TestData.ion;
qvec = [0.123456789 -0.234567891 0.345678912];
cacheDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'cache');
before = dir(fullfile(cacheDir, 'jqt2_ewald_*.mat'));
beforeNames = {before.name};
cleanup = onCleanup(@() remove_new_cache_files(cacheDir, beforeNames));

opts = struct('dipole', 'ewald', 'ewald', tc.TestData.eopts, 'cache', true);
fresh = invzt_jq_tensor(ion, qvec, opts);
after = dir(fullfile(cacheDir, 'jqt2_ewald_*.mat'));
newNames = setdiff({after.name}, beforeNames);
verifyEqual(tc, numel(newNames), 1);
cacheFile = fullfile(cacheDir, newNames{1});

hit = invzt_jq_tensor(ion, qvec, opts);
verifyEqual(tc, hit, fresh);

S = load(cacheFile);
S.lat.Jt = S.lat.Jt(:,:,[]);
save(cacheFile, '-struct', 'S');
repaired = invzt_jq_tensor(ion, qvec, opts);
verifyEqual(tc, repaired, fresh);
S2 = load(cacheFile);
verifyEqual(tc, size(S2.lat.Jt), size(fresh.Jt));
end

function test_realaxis_explicit_q_preserves_point_backend(tc)
ion = tc.TestData.ion;
eopts = tc.TestData.eopts;
g = invzt_qgrid(4, 'halfopen');
lat = invzt_jq_tensor(ion, g, struct('dipole', 'ewald', ...
    'ewald', eopts, 'cache', false));
pt = invzt_solve_point(ion, 0.1, [6 0 0], lat, struct());
verifyTrue(tc, pt.converged);
out = invzt_chi_realaxis(ion, 0.1, [6 0 0], pt, [0.005; 0.01], ...
    struct('qsel', [0.1 0 0], 'cache', false, 'force_sigma0', true));
verifyEqual(tc, out.dipole.backend, 'ewald');
verifyEqual(tc, out.explicit_q_dipole.backend, 'ewald');
verifyEqual(tc, out.explicit_q_dipole.ewald, eopts);
end

function test_realaxis_legacy_point_falls_back_to_bruteforce(tc)
ion = tc.TestData.ion;
g = invzt_qgrid(4, 'halfopen');
lat = invzt_jq_tensor(ion, g, struct('dipole', 'bruteforce', ...
    'dpRng', 6, 'cache', false));
pt = invzt_solve_point(ion, 0.1, [6 0 0], lat, struct());
verifyTrue(tc, pt.converged);
pt.lat.info = rmfield(pt.lat.info, 'dipole');
out = invzt_chi_realaxis(ion, 0.1, [6 0 0], pt, [0.005; 0.01], ...
    struct('qsel', [0.1 0 0], 'dpRng', 6, 'cache', false, 'force_sigma0', true));
verifyEqual(tc, out.explicit_q_dipole.backend, 'bruteforce');
verifyEqual(tc, out.explicit_q_dipole.ewald, ...
    struct('alpha', [], 'r_cut', [], 'g_cut', [], 'boundary', ''));
end

function remove_new_cache_files(cacheDir, beforeNames)
now = dir(fullfile(cacheDir, 'jqt2_ewald_*.mat'));
newNames = setdiff({now.name}, beforeNames);
for k = 1:numel(newNames)
    delete(fullfile(cacheDir, newNames{k}));
end
end
