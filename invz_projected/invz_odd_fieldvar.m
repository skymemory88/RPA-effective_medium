function [C, info] = invz_odd_fieldvar(ion, pt, S, T, opts)
%INVZ_ODD_FIELDVAR Tier-2 internal transverse-field covariance C (E3), meV^2.
%   [C, info] = INVZ_ODD_FIELDVAR(ion, pt, S, T, opts) evaluates the ODD plan
%   eq (E3) at a CONVERGED 1/z ODD point solve (plan T3.1, Dollberg's
%   variable-moments mechanism): the covariance of the internal transverse
%   field h_i^alpha = sum_j J_ij^{alpha c} Jz_j seen by the doublet,
%
%     C_ab = (1/(4*Nq)) * sum_q tr[ V_ac(q) * Scc(q) * V_cb(q) ]      (E3)
%
%   with V_ac(q) = Vca(:,:,iq)' (J^{ac}(q) = J^{ca}(q)', blocks are NOT
%   Hermitian individually), V_cb(q) = Vcb(:,:,iq), and the equal-time cc
%   structure factor assembled from the converged EMT lattice propagator
%   (framework eq 7, same form as invz_emt_scalar):
%
%     Gq_nu(i w_n) = G ./ (1 + (Jnu_odd(q,nu) - K).*G),  G = pt.G, K = pt.K
%     S_nu(q)      = -(1/beta) * sum_n wts_n * Gq_nu(i w_n)
%     Scc(q)       = sum_nu u_nu(q) * S_nu(q) * u_nu(q)'
%
%   SIGN CONVENTION (mirrors the solver's sum rule): invz_solve_point checks
%   sum(wts.*pt.G)/beta = -si.JzJz_fluct, i.e. the equal-time fluctuation
%   correlator is MINUS the weighted Matsubara sum of the propagator (chi =
%   -G). S_nu(q) above is therefore >= 0 for physical modes at a converged
%   paramagnetic point, which makes every per-(q,nu) term of C PSD and C
%   itself real symmetric PSD after the full-mesh q-sum.
%
%   BASIS CHOICE (plan T3.1 "document the choice"): the contraction is done in
%   the MODE basis. With y_alpha(nu) = Vc_alpha(:,:,iq)' * u_nu (alpha = a,b),
%   cyclic invariance gives per q
%     tr[V_ac Scc V_cb] = sum_nu S_nu(q) * (y_b(nu)' * y_a(nu))
%   -- note the CONJUGATE ORDER: the (a,b) element carries y_b' * y_a, so
%   C_ba(q) = conj(C_ab(q)) and each per-q contribution is Hermitian (exactly)
%   and PSD (when S_nu >= 0). This needs only two 4x4 products (Vca'*U,
%   Vcb'*U) per q on top of the eig that is recomputed anyway, and never forms
%   Scc explicitly; the sublattice-basis route (explicit Scc + traces) is the
%   brute-force cross-check in test_invz_odd_fieldvar. The eigendecomposition
%   [U, ev] = eig(Vcc + dJ hermitized) is sorted ascending with U permuted in
%   lockstep (the T2.1 w_odd pattern in invz_solve_point). dJ is rebuilt from
%   pt.odd.Xp via invz_odd_deltaJ, so C is evaluated at the SAME converged
%   state as the point solve (deterministic: same Xp -> same dJ exactly).
%
%   Imaginary parts of the off-diagonal cancel in the q-sum by the +/-q
%   pairing of a full uniform mesh (J(-q) = conj(J(q))); the assembled C is
%   asserted real symmetric to 1e-10 relative ('invz:fieldvarHerm') and then
%   defensively symmetrized to real((C + C')/2).
%
%   INPUTS
%     ion  : invz_ion() struct (gL used for the field conversion below).
%     pt   : CONVERGED output of invz_solve_point with opts.odd = true --
%            needs pt.G, pt.K (local EMT propagator/coupling on the Matsubara
%            grid), pt.odd.Xp, pt.converged. Error 'invz:fieldvarArgs' if the
%            ODD fields are missing or pt.converged is false (a non-converged
%            pt gives garbage C -- error, never warn).
%     S    : struct('Vca','Vcb','Vcc',...) -- the SAME blocks struct passed to
%            the solve (full uniform BZ mesh, invz_odd_blocks).
%     T    : temperature (K) of the solve (sets beta and the Matsubara grid).
%     opts : .Ecut (default 40) -- MUST equal the Ecut of the point solve: the
%            grid is rebuilt here and numel(pt.G) == numel(wn) is asserted
%            ('invz:fieldvarGrid').
%            .static_approx (default false) -- gates the classical shortcut
%            S_nu(q) ~ kB*T*chi_nu(q,0) = -Gq_nu(0)/beta (n = 0 term only).
%            NEVER silent: info.static_approx reports it and the caller logs
%            the full-vs-static difference once (plan T3.1).
%            .debug (default false) -- info.Cq [2,2,nq] per-q UNNORMALIZED E3
%            integrands (test hook for the two-basis cross-check).
%
%   OUTPUTS
%     C    [2,2] real symmetric PSD (meV^2): internal-field covariance in the
%          (a,b) transverse plane. C_aa = C_bb at Bx = 0 (C4).
%     info struct:
%       .heq_T         sqrt(max(eig(C))) / (ion.gL * C.muB): equivalent
%                      transverse-field scale in Tesla -- the qualitative
%                      comparator against Dollberg's distribution width
%                      h ~ 0.4 T (their Fig. 4), log-only, never asserted.
%       .tail_share    |last-frequency contribution|/|total| (Frobenius): the
%                      w_n^-2 truncation-tail diagnostic (NaN when
%                      static_approx -- single-term sum, no tail).
%       .static_approx logical, see opts.static_approx.
%       .d             E5 uniform reduction from the dJ rebuild (meV).
%       .Snu_min       min over (q,nu) of S_nu(q) (>= -rounding at a physical
%                      PM point; a real negative value means a sign/state bug).
%       .herm_resid    max |Craw - real(sym(Craw))| pre-symmetrization (meV^2).
%       .Cq            [2,2,nq] only when opts.debug (see above).
%
%   See also INVZ_SOLVE_POINT, INVZ_ODD_DELTAJ, INVZ_ODD_BLOCKS,
%   INVZ_EMT_SCALAR, INVZ_MATSUBARA.

if nargin < 5, opts = struct(); end
staticApprox = getf(opts, 'static_approx', false);
staticApprox = ~isempty(staticApprox) && ~isequal(staticApprox, false);
debugOn = getf(opts, 'debug', false);
debugOn = ~isempty(debugOn) && ~isequal(debugOn, false);
Ecut = getf(opts, 'Ecut', 40);

% --- guards: converged ODD pt + matching blocks + matching Matsubara grid ---
if ~(isstruct(pt) && isscalar(pt) && all(isfield(pt, {'G', 'K', 'converged'})) ...
        && isfield(pt, 'odd') && isfield(pt.odd, 'Xp'))
    error('invz:fieldvarArgs', ['pt must be an invz_solve_point output with ' ...
        'opts.odd = true (needs pt.G, pt.K, pt.odd.Xp): the covariance is ' ...
        'evaluated at the converged ODD state, not at a flag-off solve.']);
end
if ~pt.converged
    error('invz:fieldvarArgs', ['pt.converged is false: a non-converged point ' ...
        'has no paramagnetic EMT fixed point and its G/K give garbage C.']);
end
if ~(isstruct(S) && isscalar(S) && all(isfield(S, {'Vca', 'Vcb', 'Vcc'})))
    error('invz:fieldvarArgs', ['S must be the invz_odd_blocks struct(''Vca'',' ...
        '''Vcb'',''Vcc'',...) the point was solved with.']);
end
[~, wts, beta] = invz_matsubara(T, Ecut);
if numel(pt.G) ~= numel(wts)
    error('invz:fieldvarGrid', ['Matsubara grid mismatch: numel(pt.G) = %d vs ' ...
        '%d for (T = %g K, Ecut = %g meV) -- pt must have been solved at this ' ...
        'T with the SAME Ecut (opts.Ecut, default 40).'], ...
        numel(pt.G), numel(wts), T, Ecut);
end

% --- dJ at the SAME converged state (E1/E4/E5 from the stored chi_perp) ---
[dJ, d] = invz_odd_deltaJ(S.Vca, S.Vcb, pt.odd.Xp);

nq = size(S.Vcc, 3);
G = pt.G(:);  K = pt.K(:);
Craw  = zeros(2, 2);
Ctail = zeros(2, 2);
SnuMin = Inf;
if debugOn, Cq_dbg = zeros(2, 2, nq); end
for iq = 1:nq
    % Converged ODD cc modes at this q (Task-7 pattern: sort ev ascending,
    % permute the eigenvectors in lockstep).
    M = S.Vcc(:,:,iq) + dJ(:,:,iq);
    M = (M + M')/2;                        % both terms Hermitian; cleans rounding
    [U, ev] = eig(M, 'vector');
    [ev, isrt] = sort(real(ev));  U = U(:, isrt);
    % EMT lattice propagator per mode on the full grid (framework eq 7):
    % implicit expansion, G/K [nwn,1] x ev.' [1,4] -> Gq [nwn,4].
    Gq = G ./ (1 + (ev.' - K).*G);
    if staticApprox
        Snu = -Gq(1,:)/beta;               % kB*T*chi_nu(q,0): n = 0 term only
        SnuTail = zeros(1, 4);
    else
        Snu = -(wts.' * Gq)/beta;          % equal-time S_nu(q), sum-rule sign
        SnuTail = -(wts(end)*Gq(end,:))/beta;
    end
    SnuMin = min(SnuMin, min(Snu));
    % Mode-basis contraction (see header): columns of Ya/Yb are y_a/y_b per nu.
    Ya = S.Vca(:,:,iq)' * U;
    Yb = S.Vcb(:,:,iq)' * U;
    naa = sum(abs(Ya).^2, 1);              % y_a'*y_a  -> C(1,1)
    nbb = sum(abs(Yb).^2, 1);              % y_b'*y_b  -> C(2,2)
    nba = sum(conj(Yb).*Ya, 1);            % y_b'*y_a  -> C(1,2); conj -> C(2,1)
    Cq = [Snu*naa.', Snu*nba.'; Snu*nba', Snu*nbb.'];
    Craw = Craw + Cq;
    Ctail = Ctail + [SnuTail*naa.', SnuTail*nba.'; SnuTail*nba', SnuTail*nbb.'];
    if debugOn, Cq_dbg(:,:,iq) = Cq; end
end
Craw  = Craw / (4*nq);
Ctail = Ctail / (4*nq);

% --- +/-q imaginary cancellation: assert real symmetric, then symmetrize ---
Csym = (Craw + Craw')/2;
resid = max(abs(Craw - real(Csym)), [], 'all');
scaleC = max(abs(Craw(:)));
if resid > 1e-10*max(scaleC, 1e-30)
    error('invz:fieldvarHerm', ['assembled C not real symmetric to 1e-10 ' ...
        'relative (resid %.3g vs scale %.3g meV^2): the +/-q pairing of a ' ...
        'full uniform mesh should cancel the imaginary parts -- check the ' ...
        'qvec behind S and the conjugate order of the contraction.'], ...
        resid, scaleC);
end
C = real(Csym);

% --- diagnostics ---
C0 = invz_const();
info.heq_T = sqrt(max(eig(C))) / (ion.gL * C0.muB);
if staticApprox
    info.tail_share = NaN;                 % single-term sum: no truncation tail
else
    info.tail_share = norm(Ctail, 'fro') / max(norm(Craw, 'fro'), realmin);
end
info.static_approx = staticApprox;
info.d = d;
info.Snu_min = SnuMin;
info.herm_resid = resid;
if debugOn, info.Cq = Cq_dbg; end
end
