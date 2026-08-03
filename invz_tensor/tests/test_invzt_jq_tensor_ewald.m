function tests = test_invzt_jq_tensor_ewald
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
tc.TestData.ion = invz_ion();
tc.TestData.eopts = struct('alpha', 0.3, 'r_cut', 16, 'g_cut', 3, ...
    'boundary', 'conducting_k0_omitted');
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
    struct('ewald', eopts, 'cache', false)), 'invzt:jqTensorEwaldOptsUnexpected');
verifyError(tc, @() invzt_jq_tensor(ion, [0 0 0], ...
    struct('dipole', 'ewald', 'cache', false)), 'invzt:jqTensorEwaldOptsFields');
bad = eopts; bad.boundary = 'vacuum';
verifyError(tc, @() invzt_jq_tensor(ion, [0 0 0], ...
    struct('dipole', 'ewald', 'ewald', bad, 'cache', false)), ...
    'invzt:jqTensorEwaldBoundary');
end

function test_absent_backend_is_bit_exact_explicit_bruteforce(tc)
ion = tc.TestData.ion;
qvec = [0 0 0; 0.12 -0.08 0.19];
a = invzt_jq_tensor(ion, qvec, struct('dpRng', 6, 'cache', false));
b = invzt_jq_tensor(ion, qvec, struct('dipole', 'bruteforce', ...
    'dpRng', 6, 'cache', false));
verifyEqual(tc, b.Jt, a.Jt, 'AbsTol', 0);
verifyEqual(tc, b.JtGamma, a.JtGamma, 'AbsTol', 0);
verifyEqual(tc, a.info.dipole.backend, 'bruteforce');
verifyEqual(tc, a.info.dipole.ewald, ...
    struct('alpha', [], 'r_cut', [], 'g_cut', [], 'boundary', ''));
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
