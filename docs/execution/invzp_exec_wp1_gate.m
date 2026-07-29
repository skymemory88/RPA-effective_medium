function R = invzp_exec_wp1_gate(save_path)
%INVZP_EXEC_WP1_GATE Identity gates for the 136-state exact local source oracle (WP1).
%
% Execution packet S7 of docs/execution/invzp_plan_execution_diary.md.
% invzp_convg_fix.md WP1 requires the local generator's cumulants to be verified against
% "KMS, permutation, Hermiticity, static derivative, and high-frequency identities", with
% the gate "analytic/source derivatives match high-precision finite differences and direct
% Lehmann sums". These are those checks, run on all 136 electronuclear states through
% invzf_electronuclear_cumulant.
%
%   P1 HERMITICITY. Jz in the eigenbasis must be Hermitian to rounding.
%   P2 STATIC DERIVATIVE. C2(0) from the cumulant engine must equal dm/dhz obtained by a
%      high-order central finite difference of <Jz> with Richardson refinement. This is
%      the identity that connects the source and the correlator, and it is the one that
%      breaks first if a convention or a sign is wrong.
%   P3 SECOND STATIC DERIVATIVE. The fully static C3(0,0,0) must equal
%      d2m/dhz2 by the same finite-difference route. This grades rank 3, which had no
%      direct 136-state route at all before this wiring.
%   P4 PERMUTATION SYMMETRY. C3 and C4 must be invariant under permutation of their
%      Matsubara labels.
%   P5 KMS / REALITY. C_k(-n) must equal conj(C_k(n)); with a real operator and a
%      frequency multiset invariant under n -> -n the cumulant must be real.
%   P6 HIGH-FREQUENCY TAIL. wn^2 * C2(n) must approach a constant; the gate checks that
%      the ratio between two large indices tends to 1.
%   P7 RANK LADDER. C2(0) is recomputed with the local rank truncated; the engine's
%      discarded_boltzmann_weight is reported alongside the ACTUAL change, because the
%      former bounds only omitted thermal weight, not virtual-state truncation.
if nargin < 1, save_path = ''; end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);
addpath(fullfile(root, 'invz_functional'));
addpath(fullfile(root, 'invz_common'));

ion = invz_ion();
T = 0.30;  B = [4.0 0 0];  hz = 0.02;      % warm enough that several levels are populated
% dense_beta_span_limit is set BELOW beta*energy_span deliberately, to select the
% exponential-action backend. This is not tolerance tuning: the engine's default dense
% backend returns status = 'nonfinite' (NaN) for rank-2 cumulants at nonzero Matsubara
% labels at this temperature -- see the measured-defect note in
% invzf_electronuclear_cumulant.m. P0 below cross-checks the two backends where both work.
base = struct('hyp', true, 'transverse_mf', 'none', 'dense_beta_span_limit', 500);
fprintf('=== WP1 exact local source oracle gates (136 electronuclear states) ===\n');
fprintf('fixture: T = %.2f K, B = [%.2f 0 0] T, hz = %.4g meV, transverse_mf = none\n\n', ...
    T, B(1), hz);

% ---- P0 backend cross-check, at a temperature where BOTH backends run ----------------------
Tx = 0.60;
bx_dense  = struct('hyp', true, 'transverse_mf', 'none', 'dense_beta_span_limit', 5000);
bx_action = struct('hyp', true, 'transverse_mf', 'none', 'dense_beta_span_limit', 100);
xd = invzf_electronuclear_cumulant(ion, Tx, B, hz, [3 -3], bx_dense);
xa = invzf_electronuclear_cumulant(ion, Tx, B, hz, [3 -3], bx_action);
P0_rel = abs(xd.value - xa.value)/max(abs(xd.value), realmin);
P0 = P0_rel < 1e-10;
fprintf('P0 backend cross-check at T = %.2f K (beta*span %.4g), labels [3 -3]:\n', ...
    Tx, xd.cumulant.beta_energy_span);
fprintf('   dense_block_exponential      %.14g\n   dense_generator_exp_action   %.14g\n', ...
    xd.value, xa.value);
fprintf('   relative difference %.3e -> %s\n\n', P0_rel, pf(P0));

c2 = invzf_electronuclear_cumulant(ion, T, B, hz, [0 0], base);
beta = c2.beta;
n_levels = c2.nlevels;
fprintf('levels retained: %d   beta = %.6g meV^-1   m = <Jz> = %.12g\n', ...
    n_levels, beta, c2.m);

% ---- P1 Hermiticity -----------------------------------------------------------------------
Mz = c2.si.Mz;
herm = norm(Mz - Mz', 'fro')/max(norm(Mz, 'fro'), realmin);
P1 = herm < 1e-13;
fprintf('P1 Hermiticity of Jz in the eigenbasis: ||M-M''||_F/||M||_F = %.3e -> %s\n', herm, pf(P1));

% ---- P2 static derivative ------------------------------------------------------------------
mfun = @(x) local_m(ion, T, B, x, base);
[dm, dm_err] = richardson_derivative(mfun, hz, 1, 2e-3);
C2_0 = real(c2.value);
P2_rel = abs(C2_0 - dm)/max(abs(C2_0), realmin);
P2 = P2_rel < 1e-7 && dm_err < 1e-6*max(abs(dm),1);
fprintf('P2 static derivative: C2(0) = %.12g   dm/dhz = %.12g (fd err %.2e)   rel diff %.3e -> %s\n', ...
    C2_0, dm, dm_err, P2_rel, pf(P2));

% ---- P3 second static derivative ------------------------------------------------------------
c3 = invzf_electronuclear_cumulant(ion, T, B, hz, [0 0 0], base);
% NORMALISATION (corrected 2026-07-29 after the gate caught it): with the engine's
% convention X(tau) = beta^-1 sum_n X_n exp(-i w_n tau), the fully static rank-3 cumulant
% is C3(0,0,0) = d2m/dhz2 DIRECTLY. The first version of this gate divided by beta and
% failed by a factor of exactly beta (38.68), which is what identified the error.
C3_static = real(c3.value);
[d2m, d2m_err] = richardson_derivative(mfun, hz, 2, 4e-3);
P3_rel = abs(C3_static - d2m)/max(abs(d2m), realmin);
P3 = P3_rel < 1e-4;
fprintf('P3 second derivative: C3(0,0,0) = %.12g   d2m/dhz2 = %.12g (fd err %.2e)   rel diff %.3e -> %s\n', ...
    C3_static, d2m, d2m_err, P3_rel, pf(P3));

% ---- P4 permutation symmetry ---------------------------------------------------------------
L3 = [2 -5 3];
p3set = perms(L3);
v3 = zeros(size(p3set,1),1);
for k = 1:size(p3set,1)
    v3(k) = invzf_electronuclear_cumulant(ion, T, B, hz, p3set(k,:), base).value;
end
P4a = max(abs(v3 - v3(1)))/max(abs(v3(1)), realmin);
L4 = [1 4 -2 -3];
p4set = perms(L4);
v4 = zeros(size(p4set,1),1);
for k = 1:size(p4set,1)
    v4(k) = invzf_electronuclear_cumulant(ion, T, B, hz, p4set(k,:), base).value;
end
P4b = max(abs(v4 - v4(1)))/max(abs(v4(1)), realmin);
P4 = P4a < 1e-9 && P4b < 1e-9;
fprintf('P4 permutation symmetry: C3 over %d orderings spread %.3e; C4 over %d orderings spread %.3e -> %s\n', ...
    size(p3set,1), P4a, size(p4set,1), P4b, pf(P4));

% ---- P5 KMS / reality -----------------------------------------------------------------------
vp = invzf_electronuclear_cumulant(ion, T, B, hz, [3 -3], base).value;
vm = invzf_electronuclear_cumulant(ion, T, B, hz, [-3 3], base).value;
P5a = abs(vp - conj(vm))/max(abs(vp), realmin);
P5b = abs(imag(vp))/max(abs(vp), realmin);
v3p = invzf_electronuclear_cumulant(ion, T, B, hz, [2 -5 3], base).value;
v3m = invzf_electronuclear_cumulant(ion, T, B, hz, [-2 5 -3], base).value;
P5c = abs(v3p - conj(v3m))/max(abs(v3p), realmin);
P5 = P5a < 1e-10 && P5b < 1e-10 && P5c < 1e-10;
fprintf('P5 KMS/reality: |C2(n)-conj(C2(-n))| %.3e; |Im C2| %.3e; C3 conjugate pair %.3e -> %s\n', ...
    P5a, P5b, P5c, pf(P5));

% ---- P6 high-frequency tail ------------------------------------------------------------------
ns = [40 80 160 320];
tail = zeros(size(ns));
for k = 1:numel(ns)
    v = invzf_electronuclear_cumulant(ion, T, B, hz, [ns(k) -ns(k)], base).value;
    w = 2*pi*ns(k)/beta;
    tail(k) = real(v)*w^2;
end
ratio = tail(2:end)./tail(1:end-1);
P6 = all(isfinite(ratio)) && abs(ratio(end) - 1) < 5e-2;
fprintf('P6 high-frequency tail wn^2*C2(n): %s\n', mat2str(tail, 8));
fprintf('   successive ratios %s (last should approach 1) -> %s\n', mat2str(ratio, 6), pf(P6));

% ---- P7 rank ladder ---------------------------------------------------------------------------
ranks = [n_levels, round(n_levels/2), round(n_levels/4)];
c2r = zeros(size(ranks));  dw = zeros(size(ranks));
for k = 1:numel(ranks)
    o = base;  o.local_rank = ranks(k);
    r = invzf_electronuclear_cumulant(ion, T, B, hz, [0 0], o);
    c2r(k) = real(r.value);  dw(k) = r.discarded_boltzmann_weight;
end
rel_change = abs(c2r - c2r(1))/max(abs(c2r(1)), realmin);
P7 = isfinite(c2r(1));
fprintf('P7 rank ladder on C2(0):\n');
for k = 1:numel(ranks)
    fprintf('   local_rank %4d -> C2(0) = %.12g   rel change %.3e   discarded weight %.3e\n', ...
        ranks(k), c2r(k), rel_change(k), dw(k));
end
fprintf('   NOTE: discarded weight bounds omitted THERMAL weight only, never virtual-state\n');
fprintf('   truncation -- compare it with the actual rel change above before trusting a truncation.\n');

all_pass = P0 && P1 && P2 && P3 && P4 && P5 && P6;
fprintf('\nOVERALL: %s  (P0 %s, P1 %s, P2 %s, P3 %s, P4 %s, P5 %s, P6 %s; P7 informational)\n', ...
    pf(all_pass), pf(P0), pf(P1), pf(P2), pf(P3), pf(P4), pf(P5), pf(P6));

R = struct('P0', struct('pass', P0, 'rel', P0_rel, 'dense', xd.value, 'action', xa.value, 'T', Tx), ...
    'T', T, 'B', B, 'hz', hz, 'n_levels', n_levels, 'beta', beta, 'm', c2.m, ...
    'P1_herm', herm, 'P2', struct('C2_0', C2_0, 'dm', dm, 'fd_err', dm_err, 'rel', P2_rel), ...
    'P3', struct('C3_static', C3_static, 'd2m', d2m, 'fd_err', d2m_err, 'rel', P3_rel), ...
    'P4', struct('C3_spread', P4a, 'C4_spread', P4b), ...
    'P5', struct('c2_conj', P5a, 'c2_imag', P5b, 'c3_conj', P5c), ...
    'P6', struct('n', ns, 'tail', tail, 'ratio', ratio), ...
    'P7', struct('ranks', ranks, 'C2_0', c2r, 'rel_change', rel_change, 'discarded', dw), ...
    'pass', [P0 P1 P2 P3 P4 P5 P6], 'all_pass', all_pass);
if ~isempty(save_path), save(save_path, 'R'); fprintf('saved: %s\n', save_path); end
end

function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end

function m = local_m(ion, T, B, hz, base)
o = base;  o.hz_fixed = hz;
si = invz_single_ion(ion, T, B, o);
m = si.Jexp(3);
end

function [d, err] = richardson_derivative(f, x, order, h0)
%RICHARDSON_DERIVATIVE Central difference of the given order with one Richardson step.
% Returns the extrapolated value and |extrapolated - coarse| as the measured error.
switch order
    case 1
        D = @(h) (f(x+h) - f(x-h))/(2*h);
    case 2
        D = @(h) (f(x+h) - 2*f(x) + f(x-h))/h^2;
    otherwise
        error('order must be 1 or 2');
end
d1 = D(h0);  d2 = D(h0/2);  d3 = D(h0/4);
r1 = (4*d2 - d1)/3;         % O(h^4) after one Richardson step
r2 = (4*d3 - d2)/3;
d = (16*r2 - r1)/15;        % O(h^6)
err = abs(d - r2);
end
