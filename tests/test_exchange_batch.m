function test_exchange_batch()
% Test and benchmark exchange.m vs exchange_batch.m
% Usage: run in MATLAB: test_exchange_batch

% Lattice (cubic, a0 = 5 Angstrom)
a0 = 5.0;
a = [a0 0 0; 0 a0 0; 0 0 a0];

% 4-site basis (fractional coordinates)
tau = [
    0.0 0.0 0.0;
    0.5 0.5 0.0;
    0.5 0.0 0.5;
    0.0 0.5 0.5
];

% Exchange strength (meV)
Jex = 0.1;

% Random q-points in reduced units (Miller indices)
Nq = 200;  % keep moderate for quick test
qvec = rand(Nq,3) - 0.5;  % uniform in [-0.5, 0.5]

% Reference: looped exchange
ex_ref = zeros(3,3,size(tau,1),size(tau,1),Nq);
t_ref = tic;
for iq = 1:Nq
    ex_ref(:,:,:,:,iq) = exchange(qvec(iq,:), Jex, a, tau);
end
t_ref = toc(t_ref);

% Batched: exchange_batch
t_bat = tic;
ex_all = exchange_batch(qvec, Jex, a, tau);
t_bat = toc(t_bat);

% Agreement metrics
abs_diff = abs(ex_all - ex_ref);
max_abs = max(abs_diff(:));
rel_diff = abs_diff ./ max(abs(ex_ref(:)), eps);
max_rel = max(rel_diff(:));

fprintf('exchange vs exchange\\_batch: Nq=%d\n', Nq);
fprintf('  max |abs diff| = %.3e\n', max_abs);
fprintf('  max |rel diff| = %.3e\n', max_rel);
fprintf('  time looped    = %.3fs\n', t_ref);
fprintf('  time batched   = %.3fs\n', t_bat);
fprintf('  speedup        = %.2fx\n', t_ref / max(t_bat, eps));

% Sanity threshold (tolerate roundoff)
assert(max_abs < 1e-10 || max_rel < 1e-10, 'exchange_batch mismatch exceeds tolerance');

end

