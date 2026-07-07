function test_dipole_batch()
% Test and benchmark MF_dipole (looped) vs dipole_batch (vectorized)
% Usage: run in MATLAB: test_dipole_batch

% Lattice (cubic)
a0 = 5.0;
a = [a0 0 0; 0 a0 0; 0 0 a0];

% Basis
tau = [
    0.0 0.0 0.0;
    0.5 0.5 0.0;
    0.5 0.0 0.5;
    0.0 0.5 0.5
];

% Dipole range (keep small for test)
N = 2;

% q points
Nq = 60;  % modest for quick test
qvec = rand(Nq,3) - 0.5;

% Reference loop
dip_ref = zeros(3,3,size(tau,1),size(tau,1),Nq);
t_ref = tic;
for iq = 1:Nq
    dip_ref(:,:,:,:,iq) = MF_dipole(qvec(iq,:), N, a, tau);
end
t_ref = toc(t_ref);

% Batched
t_bat = tic;
dip_all = dipole_batch(qvec, N, a, tau);
t_bat = toc(t_bat);

% Compare
abs_diff = abs(dip_all - dip_ref);
max_abs = max(abs_diff(:));
rel_diff = abs_diff ./ max(abs(dip_ref(:)), eps);
max_rel = max(rel_diff(:));

fprintf('dipole vs dipole\\_batch: Nq=%d, N=%d\n', Nq, N);
fprintf('  max |abs diff| = %.3e\n', max_abs);
fprintf('  max |rel diff| = %.3e\n', max_rel);
fprintf('  time looped    = %.3fs\n', t_ref);
fprintf('  time batched   = %.3fs\n', t_bat);
fprintf('  speedup        = %.2fx\n', t_ref / max(t_bat, eps));

assert(max_abs < 1e-10 || max_rel < 1e-10, 'dipole_batch mismatch exceeds tolerance');

end

