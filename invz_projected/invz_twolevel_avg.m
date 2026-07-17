function [tla, avg] = invz_twolevel_avg(ion, T, Bx, C, opts)
%INVZ_TWOLEVEL_AVG Gauss-Hermite-dressed doublet over quenched transverse-field disorder (T3.2).
% [tla, avg] = invz_twolevel_avg(ion, T, Bx, C, opts)
%
% Averages the electronic two-level response over a quenched Gaussian distribution of
% transverse fields (h_a, h_b) ~ N(0, C) — C the 2x2 internal-field covariance in meV^2
% from invz_odd_fieldvar (T3.1) — then fits EFFECTIVE two-level parameters
% (Delta_eff, n01_eff, M2_eff) so the Sigma machinery keeps its two-level algebra.
%
% *** FLAG (plan T3.2, verbatim in spirit): this quenched-Gaussian, STATIC dressing of the
% *** doublet is the LEAST RIGOROUS step of the whole plan — the thermal/quantum field
% *** fluctuations are treated as a frozen classical distribution, and the averaged
% *** response is re-compressed into a single pole. Treat downstream numbers accordingly.
%
% Node fields enter the single-ion problem as a FIELD SHIFT (no single-ion code changes):
%   B_node = invz_field_vec(Bx) + [ha, hb, 0]/(ion.gL*C0.muB)
% so the Zeeman term -gL*muB*B.J carries exactly -ha*Ja - hb*Jb ((ha,hb) in meV,
% conjugate to J_a/J_b).
%
% C = zeros(2) (all-zero) SHORT-CIRCUITS: tla = invz_twolevel(ion, T, Bx, opts_pass)
% bitwise, no quadrature (avg.nodes = 1; avg.G0/chi0cc0 still honored when opts.G0 = true,
% so the Tier-2 consumer sees a uniform contract across the C = 0 limit).
%
% opts:
%   .ngh  (default 7)          Gauss-Hermite nodes per axis (2-D tensor product).
%   .avg  (default 'response') 'response': average g(i*wn) node-by-node (the physically
%         preferred disorder-averaged propagator, plan design decision) and FIT
%         (Delta_eff, n01_eff) by matching EXACTLY (i) gbar(0) and (ii) the leading
%         wn^-2 tail coefficient lim wn^2*gbar = sum_w w*2*n01*Delta (analytic per node);
%         'params': average (Delta, M2, n01) directly (inequivalent variant, kept for the
%         plan's one-shot comparison).
%   .wn   Matsubara grid from invz_matsubara (REQUIRED in response mode: invz:tlavgArgs).
%   .G0   (default FALSE) opt-in: also compute avg.G0/avg.chi0cc0, the node-averaged FULL
%         electronuclear cc propagator (hyp = true, elastic on — the disorder-averaged G0
%         the Tier-2 outer loop swaps in; Task 10 passes G0 = true). Requires opts.wn
%         (invz:tlavgArgs otherwise). AMENDED 2026-07-17 (Task-9 adjudication): opt-in,
%         not a side effect of opts.wn — the 49 electronuclear diagonalizations per call
%         were burning fast-suite budget in tests that only need the two-level fit.
%   .Jxx0, .transverse_mf      forwarded to invz_twolevel / invz_single_ion when given.
%
% Fit algebra (response mode): with gbar0 = sum_w w*g0_node and tail = sum_w w*2*n01*Delta,
%   Delta_eff = sqrt(tail/gbar0),  n01_eff = gbar0*Delta_eff/2
% match conditions (i)+(ii) exactly. M2_eff = sum_w(w*M2*g0)/gbar0 — RESPONSE-weighted:
% M^2 sits OUTSIDE g in this codebase (chi = M2*g), so the g(0)+tail fit pins (Delta, n01)
% only and M2 needs its own weighting; weighting each node's M2 by its static response g0
% is the disorder-averaged-propagator-consistent choice (chi_bar(0) = M2_eff*gbar0).
% The third two-level identity, the sum rule (1/beta)*sum_n wts.*g = n01*coth(beta*Delta/2)
% = tanh*coth = 1 EXACTLY per node (analytic, infinite sum; the truncated grid gives
% ~ 1 - 2*n01*Delta/(pi*Ecut)), is NOT imposed by the fit: avg.fit_resid reports the max
% relative mismatch of the three conditions, of which g(0) and the tail vanish by
% construction and the truncated-grid sum-rule mismatch |sumrule_fit - sumrule_avg| is the
% honest fit residual (second order in the node-parameter spread: the fitted (Delta_eff,
% n01_eff) pair no longer satisfies n01 = tanh(beta*Delta/2) exactly).
%
% Degenerate nodes (T3.4 groundwork): nodes where invz_twolevel throws
% invz:degenerateDoublet (Delta_node < 1e-4 meV, e.g. the exact-origin node at Bx = 0) are
% evaluated in their h -> 0+ limit (other errors, incl. invz:orderedPhase, are rethrown):
%   g0_node = beta          (2*tanh(beta*D/2)/D -> beta as D -> 0)
%   tail_node = 2*n01*D -> 0,  n01_node -> 0,  m_node -> 0
%   g_node(i*wn) = beta*(wn == 0)   (elastic Curie weight; per-node sum rule = 1 exactly)
%   M2_node = |Mz(1,2)|^2 + |Mz(1,1)|^2  from a direct invz_single_ion(hyp=false) call.
% NOTE the M2 form (controller-approved, Task-9 adjudication): on an exactly degenerate
% doublet eig returns an ARBITRARY basis of the 2-dim subspace, so abs(Mz(1,2))^2 alone is
% basis-dependent (anything in [0, m^2]); the sum |Mz11|^2 + |Mz12|^2 = (Mz^2)(1,1) is
% basis-INVARIANT — in the Jz-DIAGONAL (Ising) basis of the doublet it equals
% Mz(1,1)^2 = m^2, the doublet's Curie weight — and equals the h -> 0+ limit of the
% invz_twolevel construction (split states carry Mz11 -> 0, |Mz12| -> m). Such nodes are
% counted in avg.n_degenerate.
%
% Rank-deficient C: eigenvalues of C are clamped at 0 (tiny negatives from roundoff), so a
% zero eigenvalue collapses that quadrature axis to duplicated nodes at h = 0 — weights
% still sum to 1 and the average is exact for the surviving axis; no further special-casing.
%
% tla carries the SAME field names as invz_twolevel (Delta, M2, m, n01, g0, transverse_mf)
% plus tla.gh = true (absent on the C = 0 short-circuit, which returns the literal
% invz_twolevel struct). tla.m is the w-averaged node moment (~0 in the PM phase);
% tla.g0 = 2*n01_eff/Delta_eff (consistency). Delta_eff < 1e-4 meV errors with
% invz:degenerateDoublet exactly like the bare constructor.
%
% avg: mode, nodes, ngh, Delta_eff, M2_eff, n01_eff, fit_resid, sumrule_avg, node_Delta
% (per-node splittings, spread diagnostic), n_degenerate, gbar0, tail; response average
% gbar (when opts.wn given) and G0 [nwn x 1] + chi0cc0 (when opts.G0 = true).

if nargin < 5, opts = struct(); end
C0 = invz_const();
beta = 1/(C0.kB*T);

% ---- validate C ------------------------------------------------------------------------
if ~isnumeric(C) || ~isreal(C) || ~isequal(size(C), [2 2]) || ~all(isfinite(C(:)))
    error('invz:tlavgArgs', 'C must be a real, finite 2x2 covariance matrix (meV^2).');
end
if max(abs(C - C.'), [], 'all') > 1e-8*max(abs(C(:)))
    error('invz:tlavgArgs', 'C must be symmetric (max asymmetry %.3g meV^2).', ...
        max(abs(C - C.'), [], 'all'));
end

% opts_pass: forward ONLY the single-ion knobs invz_twolevel understands, when given.
op = struct();
if isfield(opts, 'Jxx0'),          op.Jxx0          = opts.Jxx0;          end
if isfield(opts, 'transverse_mf'), op.transverse_mf = opts.transverse_mf; end
tmf = getf(opts, 'transverse_mf', 'legacy_x');

mode = getf(opts, 'avg', 'response');
if ~(ischar(mode) || isstring(mode)) || ~any(strcmp(mode, {'response', 'params'}))
    error('invz:tlavgArgs', 'opts.avg must be ''response'' or ''params''.');
end
mode = char(mode);

has_wn = isfield(opts, 'wn') && ~isempty(opts.wn);
if has_wn
    wn = opts.wn(:);
    if ~isreal(wn) || any(wn < 0) || wn(1) ~= 0 || any(diff(wn) <= 0)
        error('invz:tlavgArgs', ...
            'opts.wn must be the invz_matsubara grid: real, ascending, wn(1) = 0.');
    end
    wts = 1 + (wn > 0);            % invz_matsubara doubling weights (n >= 0 grid)
end
do_G0 = getf(opts, 'G0', false);
if do_G0 && ~has_wn
    error('invz:tlavgArgs', 'opts.G0 = true requires opts.wn (from invz_matsubara).');
end

% ---- C = 0 short-circuit (bitwise reproduction, plan T3.2 acceptance) -------------------
if norm(C) == 0
    tla = invz_twolevel(ion, T, Bx, op);
    avg = struct('mode', mode, 'nodes', 1, 'ngh', getf(opts, 'ngh', 7), ...
        'Delta_eff', tla.Delta, 'M2_eff', tla.M2, 'n01_eff', tla.n01, ...
        'fit_resid', 0, 'node_Delta', tla.Delta, 'n_degenerate', 0, ...
        'gbar0', tla.g0, 'tail', 2*tla.n01*tla.Delta);
    if has_wn
        avg.gbar = real(invz_g(tla, 1i*wn));
        avg.sumrule_avg = sum(wts.*avg.gbar)/beta;
        if do_G0
            [avg.G0, avg.chi0cc0] = local_g0(ion, T, invz_field_vec(Bx), wn, op);
        end
    else
        avg.sumrule_avg = NaN;
    end
    return
end

if strcmp(mode, 'response') && ~has_wn
    error('invz:tlavgArgs', 'response averaging requires opts.wn (from invz_matsubara).');
end

% ---- 2-D Gauss-Hermite nodes -------------------------------------------------------------
ngh = getf(opts, 'ngh', 7);
[Vc, s2] = eig((C + C.')/2, 'vector');       % C = Vc*diag(s2)*Vc', s2 ascending real
if any(s2 < -1e-10*max(s2))
    error('invz:tlavgArgs', 'C must be PSD (eigenvalues %.3g, %.3g meV^2).', s2(1), s2(2));
end
s2 = max(s2, 0);                             % clamp roundoff negatives (rank-deficient C)
[x, wgh] = gh_nodes(ngh);
[xi, xj] = ndgrid(x, x);
[wi, wj] = ndgrid(wgh, wgh);
w = wi(:).*wj(:)/pi;                          % normalized: sum = (sqrt(pi))^2/pi = 1
if abs(sum(w) - 1) > 1e-12
    error('invz:tlavgWeights', 'GH weights sum to 1%+.3e, need |.| < 1e-12.', sum(w) - 1);
end
hab = Vc*[sqrt(2*s2(1))*xi(:).'; sqrt(2*s2(2))*xj(:).'];   % [ha; hb] per node, meV
nn = numel(w);

% ---- per-node loop -----------------------------------------------------------------------
Bbase = invz_field_vec(Bx);
scal  = ion.gL*C0.muB;                        % meV/T
gbar0 = 0;  tail = 0;  M2g0 = 0;  m_acc = 0;
Dp = 0;  M2p = 0;  n01p = 0;                  % 'params'-mode plain averages
node_Delta = zeros(nn, 1);
n_deg = 0;
if has_wn
    gbar = zeros(numel(wn), 1);
end
if do_G0
    G0   = zeros(numel(wn), 1);
    chi0cc0 = 0;
end
sop = op;  sop.hyp = false;                   % degenerate-branch single-ion opts
for k = 1:nn
    B_node = Bbase + [hab(1,k), hab(2,k), 0]/scal;
    try
        tln = invz_twolevel(ion, T, B_node, op);
        Dn = tln.Delta;  n01n = tln.n01;  M2n = tln.M2;  mn = tln.m;  g0n = tln.g0;
        tailn = 2*n01n*Dn;
        if has_wn
            gn = real(invz_g(tln, 1i*wn));
            sn = sum(wts.*gn)/beta;           % ~ 1 - 2*n01*Delta/(pi*Ecut) (truncation)
        end
    catch ME
        if ~strcmp(ME.identifier, 'invz:degenerateDoublet'), rethrow(ME); end
        % h -> 0+ limit of the two-level response (header note; T3.4 groundwork)
        si = invz_single_ion(ion, T, B_node, sop);
        Dn = si.E(2) - si.E(1);               % raw tiny splitting, recorded for the spread
        M2n = abs(si.Mz(1,2))^2 + abs(si.Mz(1,1))^2;   % basis-invariant (Mz^2)(1,1) = m^2
        mn = 0;  n01n = 0;  g0n = beta;  tailn = 0;
        if has_wn
            gn = beta*(wn == 0);              % elastic Curie weight
            sn = 1;                           % exact: wts(wn==0) = 1
        end
        n_deg = n_deg + 1;
    end
    node_Delta(k) = Dn;
    gbar0 = gbar0 + w(k)*g0n;
    tail  = tail  + w(k)*tailn;
    M2g0  = M2g0  + w(k)*M2n*g0n;
    m_acc = m_acc + w(k)*mn;
    Dp = Dp + w(k)*Dn;  M2p = M2p + w(k)*M2n;  n01p = n01p + w(k)*n01n;
    if has_wn
        if abs(sn - 1) >= 2e-2               % per-node truncated sum rule (analytic = 1)
            error('invz:tlavgSumrule', ...
                'node %d truncated sum rule %.6g departs from 1 by >= 2e-2 (Ecut too low?).', k, sn);
        end
        gbar = gbar + w(k)*gn;
    end
    if do_G0
        [G0n, c0n] = local_g0(ion, T, B_node, wn, op);
        G0 = G0 + w(k)*G0n;
        chi0cc0 = chi0cc0 + w(k)*c0n;
    end
end

% ---- effective two-level fit -------------------------------------------------------------
if strcmp(mode, 'response')
    Delta_eff = sqrt(tail/gbar0);             % matches gbar(0) and the wn^-2 tail exactly
    n01_eff   = gbar0*Delta_eff/2;
    M2_eff    = M2g0/gbar0;                   % response(g0)-weighted (header note)
else
    Delta_eff = Dp;  n01_eff = n01p;  M2_eff = M2p;   % plain parameter averages
end
if ~(isfinite(Delta_eff) && Delta_eff >= 1e-4)
    error('invz:degenerateDoublet', ...
        'Dressed splitting %.2e meV too small: Bx=0 limit needs the closed-form Sigma_c (invz_sigma_crit).', Delta_eff);
end

avg = struct('mode', mode, 'nodes', nn, 'ngh', ngh, ...
    'Delta_eff', Delta_eff, 'M2_eff', M2_eff, 'n01_eff', n01_eff, ...
    'node_Delta', node_Delta, 'n_degenerate', n_deg, 'gbar0', gbar0, 'tail', tail);
if has_wn
    % fit residual: max relative mismatch of the three matched conditions when reproduced
    % by the fitted form (g(0) and tail are exact by construction in response mode; the
    % truncated sum rule is the honest residual — in params mode all three are the
    % inequivalence diagnostic, not a fit quality).
    tlf = struct('Delta', Delta_eff, 'n01', n01_eff);
    sumrule_fit = sum(wts.*real(invz_g(tlf, 1i*wn)))/beta;
    avg.sumrule_avg = sum(wts.*gbar)/beta;    % == sum_w w*sn by linearity
    r1 = abs(2*n01_eff/Delta_eff - gbar0)/abs(gbar0);
    r2 = abs(2*n01_eff*Delta_eff - tail)/abs(tail);
    r3 = abs(sumrule_fit - avg.sumrule_avg)/abs(avg.sumrule_avg);
    avg.fit_resid = max([r1, r2, r3]);
    avg.gbar = gbar;
    if do_G0
        avg.G0 = G0;
        avg.chi0cc0 = chi0cc0;
    end
else
    avg.sumrule_avg = NaN;
    avg.fit_resid = NaN;
end

tla.Delta = Delta_eff;
tla.M2    = M2_eff;
tla.m     = m_acc;
tla.n01   = n01_eff;
tla.g0    = 2*n01_eff/Delta_eff;
tla.transverse_mf = tmf;
tla.gh    = true;
end

% -------------------------------------------------------------------------------------------
function [G0, c0cc0] = local_g0(ion, T, B_node, wn, op)
%LOCAL_G0 Full electronuclear cc propagator at one node (hyp = true, elastic on).
% G0 = -real(chi0_cc(i*wn)) (codebase convention chi = -G); c0cc0 = static cc value at
% z = 0 (wn(1) = 0 entry carries the elastic beta-term), so c0cc0 = -G0(1).
soph = op;  soph.hyp = true;
si_h = invz_single_ion(ion, T, B_node, soph);
c0 = invz_chi0z(si_h, T, 1i*wn, struct('elastic', true));
G0 = -real(squeeze(c0(3,3,:)));
c0cc0 = real(c0(3,3,1));
end

function [x, w] = gh_nodes(n)
%GH_NODES Gauss-Hermite nodes/weights (physicists' weight exp(-x^2)), Golub-Welsch:
% eig of the symmetric tridiagonal Jacobi matrix (diag 0, off-diag sqrt(k/2)); nodes are
% the ascending eigenvalues, weights sqrt(pi)*(first eigenvector component)^2.
% Validation values: n = 2: x = -/+ 1/sqrt(2) = -/+ 0.707106781, w = sqrt(pi)/2 each;
%                    n = 3: x = 0, -/+ sqrt(3/2) = -/+ 1.224744871,
%                           w = 2*sqrt(pi)/3 (center), sqrt(pi)/6 (outer).
if n == 1
    x = 0;  w = sqrt(pi);  return
end
b = sqrt((1:n-1)/2);
J = diag(b, 1) + diag(b, -1);                % symmetric tridiagonal Jacobi matrix
[V, d] = eig(J, 'vector');                   % symmetric input: ascending real eigenvalues
[x, ix] = sort(real(d));                     % defensive sort (eig already ascending)
V = V(:, ix);
w = sqrt(pi)*(V(1, :).').^2;
end
