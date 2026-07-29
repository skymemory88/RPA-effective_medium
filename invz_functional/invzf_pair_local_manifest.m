function man = invzf_pair_local_manifest(Delta, M, h, beta, opts)
%INVZF_PAIR_LOCAL_MANIFEST Coupling-order manifest for the two-site bond, from LOCAL data only.
%
%   MAN = INVZF_PAIR_LOCAL_MANIFEST(DELTA, M, H, BETA, OPTS) predicts the Taylor
%   coefficients of the exact two-site cluster free energy
%
%       F_pair(J) = sum_k a_k J^k,        H = sum_i [D|1><1|_i - H X_i] - J X_1 X_2,
%
%   using ONLY single-site quantities: the local moment m = <X> and the local
%   connected Matsubara correlator C_n at J = 0. Nothing from the two-site
%   diagonalisation enters. This is the "compute it independently" half of
%   invzp_convg_fix.md Sec. 6 -- the side that a candidate functional's diagram
%   generator must eventually reproduce, and that the exact cluster grades.
%
%   DERIVATION (linked-cluster / cumulant expansion; recorded here because the
%   coefficients are the load-bearing content, not the code):
%     ln Z - ln Z_0 = sum_k J^k kappa_k / k!,   A = int_0^beta X_1(tau) X_2(tau) dtau,
%   with kappa_k the k-th time-ordered cumulant of A in the DECOUPLED ensemble.
%   Writing the per-site two-point function as <T x(t1) x(t2)> = m^2 + G(t1-t2) and
%   using C_n = int_0^beta G(tau) exp(i w_n tau) dtau, and that the two sites are
%   statistically independent:
%     kappa_1 = beta m^2
%     kappa_2 = 2 beta m^2 C_0 + sum_n C_n^2
%   and F = -(ln Z)/beta gives
%     a_1 = -m^2
%     a_2 = -m^2 C_0 - (1/(2 beta)) sum_n C_n^2.
%
%   TOPOLOGY PROVENANCE (what each retained monomial IS -- Sec. 5.3's classification
%   obligation, at the orders derived here):
%     a_1                        one lattice line joining a C1 vertex on each site.
%     a_2 term "C1C1|C2"         two lattice lines; one site carries a connected C2,
%                                the other carries two disconnected C1 vertices;
%                                2 site-orientations / 2! = symmetry factor 1.
%     a_2 term "ring C2|C2"      two lattice lines; a connected C2 on each site;
%                                symmetry factor 1/2. This is the single-ring
%                                (RPA-like) term.
%
%   ORDER 3 IS DERIVED ONLY AT H = 0. There the local Hamiltonian is diagonal in the
%   basis in which X is purely off-diagonal, so every odd local correlator vanishes:
%   m = 0 and C3 = 0, hence a_1 = a_3 = 0 exactly. More strongly, the local unitary
%   that flips X on one site maps J -> -J at H = 0, so F_pair is EVEN in J and every
%   odd coefficient vanishes. At H ~= 0, a_3 requires the frequency-resolved C3
%   contractions (m^3 gamma_3, m^2 C_0^2, m*<G C3> and <C3 C3>) and is NOT derived
%   here; those slots are returned with derived = false so the harness reports them
%   as ungraded rather than silently passing.
%
%   OPTS: .ncut (4000) Matsubara index cutoff for sum_n C_n^2; .tail_check (true)
%   re-evaluates the sum at ncut/2 and reports the difference as the truncation error.
%
%   MAN fields: .a (predicted coefficients, NaN where underived), .derived (logical),
%   .terms (per-order cell of topology records), .m, .C0, .sumC2, .sumC2_trunc_err,
%   .parity_zero (orders forced to zero by the H = 0 symmetry), .provenance.
if nargin < 5, opts = struct(); end
ncut = getf_local(opts, 'ncut', 4000);
tail_check = getf_local(opts, 'tail_check', true);
validateattributes(beta, {'numeric'}, {'real','scalar','finite','positive'});

% --- single-site exact data (one site, zero coupling) --------------------------------------
wn = (-ncut:ncut).';
loc = invzf_cluster_exact(Delta, M, h, 0, beta, wn, 1);
m  = loc.m_site(1);
Cn = loc.C_local(:, 1);
i0 = find(wn == 0, 1);
C0 = Cn(i0);
sumC2 = sum(Cn.^2);

sumC2_err = NaN;
if tail_check
    half = abs(wn) <= floor(ncut/2);
    sumC2_err = abs(sumC2 - sum(Cn(half).^2));
end

is_zero_field = (h == 0);

a = nan(5,1);
derived = false(5,1);
a(1) = loc.F * 2;                       % a_0: two independent sites  (index 1 <-> order 0)
derived(1) = true;
a(2) = -m^2;                            derived(2) = true;   % order 1
a(3) = -m^2*C0 - sumC2/(2*beta);        derived(3) = true;   % order 2
if is_zero_field
    a(4) = 0;  derived(4) = true;                             % order 3, by parity
end
% order 4 is not derived at any field here.

terms = cell(5,1);
terms{1} = struct('order', 0, 'label', 'decoupled', ...
    'topology', 'two independent local partition functions, no lattice line', ...
    'symmetry_factor', 1, 'value', a(1));
terms{2} = {struct('order', 1, 'label', 'C1|C1', ...
    'topology', 'one lattice line joining a C1 (moment) vertex on each site', ...
    'symmetry_factor', 1, 'value', -m^2)};
terms{3} = {struct('order', 2, 'label', 'C1C1|C2', ...
    'topology', ['two lattice lines; connected C2 on one site, two disconnected ' ...
                 'C1 vertices on the other; 2 orientations / 2!'], ...
    'symmetry_factor', 1, 'value', -m^2*C0), ...
            struct('order', 2, 'label', 'ring C2|C2', ...
    'topology', 'two lattice lines; connected C2 on each site (single ring)', ...
    'symmetry_factor', 0.5, 'value', -sumC2/(2*beta))};
if is_zero_field
    terms{4} = {struct('order', 3, 'label', 'parity-forbidden', ...
        'topology', ['all order-3 topologies carry an odd number of legs on each ' ...
                     'site; every odd local correlator vanishes at h = 0'], ...
        'symmetry_factor', 0, 'value', 0)};
else
    terms{4} = {};
end
terms{5} = {};

parity_zero = [];
if is_zero_field, parity_zero = 1:2:4; end     % orders 1 and 3 (and all odd orders)

man = struct('a', a, 'derived', derived, 'terms', {terms}, ...
    'm', m, 'C0', C0, 'sumC2', sumC2, 'sumC2_trunc_err', sumC2_err, ...
    'ncut', ncut, 'is_zero_field', is_zero_field, 'parity_zero', parity_zero, ...
    'fixture', struct('Delta', Delta, 'M', M, 'h', h, 'beta', beta), ...
    'provenance', struct('schema', 'invzf_pair_local_manifest/v1', ...
        'derivation', 'linked-cluster cumulant expansion of a single X-X bond', ...
        'orders_derived_any_field', [0 1 2], ...
        'orders_derived_zero_field_only', 3, ...
        'orders_not_derived', 4, ...
        'source_of_local_data', 'invzf_cluster_exact with N = 1, J = 0'));
end

function v = getf_local(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
