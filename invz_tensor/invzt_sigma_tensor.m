function out = invzt_sigma_tensor(si, T, lat_eff, wn, beta, opts)
%INVZT_SIGMA_TENSOR  A3 genuine tensor 1/z self-energy from the exact four-point vertex.
%
%   out = INVZT_SIGMA_TENSOR(si, T, lat_eff, wn, beta, opts) assembles the
%   component-labelled tensor self-energy correction Vmat_{mu nu}(i omega_n)
%   [3,3,nwn] for the three-state single-ion toy si (INVZT_THREESTATE), resums it
%   into a renormalized local susceptibility, closes it against the lattice through
%   the A2 matrix effective medium, and iterates the whole map to self-consistency.
%   This is the A3 rung of docs/superpowers/plans/2026-07-17-invz-tensor-full.md
%   (Task 12); its central validation is the framework SS11.8 EMERGENCE gate (A3's
%   Gaussian truncation reproduces A1/E1), exercised by INVZT_SOLVE_POINT mode 'a3'.
%
%   INPUTS
%     si      : three-state toy struct (E, P, Mx, My, Mz, Jexp; INVZT_THREESTATE).
%     T       : temperature (K).
%     lat_eff : the (odd-handled) LOCKED lattice struct (INVZT_JQ_TENSOR); its ODD
%               (c<->a,b) blocks are already zeroed when odd = false.
%     wn      : Matsubara frequencies (numel = nwn); external bosonic indices are 0..nwn-1.
%     beta    : inverse temperature (meV^-1).
%     opts    : mix_outer (0.5), tol_outer (1e-9), max_outer (150), rank_tol (1e-12),
%               Lmax (nwn-1 -- the internal l-sum then spans the SAME grid the scalar
%               invz_lambdas sums, so the two-level limit is exact), Vmat_seed [].
%
%   THE A3 MAP (one outer iteration; the fixed-point variable is Vmat).
%     (0) es = (si.E, si.P); ops.{a,b,c} = CENTERED si.{Mx,My,Mz} (delta-J: the vertex
%         path sum + Wick-pairing subtraction assume zero-mean Hermitian operators, so
%         the elastic Curie term is handled automatically and matches -invz_chi0z).
%         G0 = -invz_chi0z(si, T, i*wn) [3,3,nwn] is the BARE local tensor propagator.
%     (1) chi_til = -G0til with the LOCKED SYMMETRIC BRACKET resummation (constraint 2 --
%         carry the additive correction, never a ratio; never inv(G0)):
%             G0til = G0 * ((G0 + Vmat) \ G0)          on the ACTIVE SUBSPACE P
%         (scalar reduction: G0til = G0^2/(G0+V) = G0/(1+Sigma) with V = G0*Sigma, so
%         chi_til_cc = chi0_cc/(1+Sigma) -- exactly A1's cdom/(1+Sigma) factor). P is
%         the frequency-consistent union-of-ranges active subspace of G0 (shared policy
%         with INVZT_EMT_MATRIX): with rho -> 0 the transverse channels are null and P
%         collapses to cc, so the cc sector decouples EXACTLY into the scalar chain.
%     (2) [Kmat, chi_bar] = INVZT_EMT_MATRIX(chi_til, lat_eff): the lattice RPA of
%         chi_til + the DIRECT matrix effective-medium closure -> the internal-channel
%         kernel K_{rho sigma}(i omega_l) [3,3,nwn] (the FULL non-Hermitian K -- its
%         gyrotropic imaginary-antisymmetric part is physical; consumed AS-IS, never
%         symmetrized, hard-math constraint 9).
%     (3) Vmat = (1/2beta) sum_l sum_{rho sigma} K_{rho sigma}(l) Gamma4_{mu nu;rho sigma}(n,l),
%         the exact connected four-point contraction (INVZT_VERTEX4). Gamma4 is
%         Kmat-INDEPENDENT, so it is PRECOMPUTED ONCE (via stage 'Gamma' over the full
%         (n,l) grid) and only the cheap K-contraction repeats each outer iteration --
%         bit-identical to a stage-'V' call every iteration, ~30x faster. Negative l is
%         reconstructed by the LOCKED transpose K_{rho sigma}(-l)=K_{sigma rho}(+l).
%     Damped outer mix on Vmat; converged when max|Vmat_new - Vmat| < tol_outer.
%
%   OUTPUTS (struct out)
%     Vmat        [3,3,nwn] : the converged tensor self-energy correction (V = G0.Sigma).
%     chi_til     [3,3,nwn] : Dyson (symmetric-bracket) renormalized local chi.
%     chi_til_add [3,3,nwn] : ADDITIVE resummation -(G0 - Vmat) (constraint 8: the
%                             Cartesian Dyson-vs-additive spread is the O(1/z^2) method
%                             error bar; the caller reports crit(dyson)-crit(additive)).
%     Kmat        [3,3,nwn] : internal-channel medium kernel (for pt.Kmat / diagnostics).
%     chi_bar     [3,3,nwn] : effective-medium response (chi_bar_cc -> Gloc = -real chi_bar_cc).
%     Gloc        [nwn,1]   : -real(chi_bar_cc), the site-local cc effective medium.
%     G0          [3,3,nwn] : the bare local propagator (so the caller forms Sigma_cc_equiv).
%     emtinfo, converged, outer_iters, active_rank, pruned_bound.
%
%   See also INVZT_THREESTATE, INVZT_VERTEX4, INVZT_EMT_MATRIX, INVZT_SOLVE_POINT.
if nargin < 6 || isempty(opts), opts = struct(); end
mixo  = getf(opts, 'mix_outer', 0.5);
tolo  = getf(opts, 'tol_outer', 1e-9);
maxo  = getf(opts, 'max_outer', 400);
rtol  = getf(opts, 'rank_tol', 1e-12);
nwn   = numel(wn);
Lmax  = getf(opts, 'Lmax', nwn - 1);
% dress: 'full' (default) dresses ALL Cartesian channels via the four-point vertex;
% 'dominant' is E1 (Jensen's "dominant renormalized, weak at RPA" rule) -- the vertex
% dresses ONLY the DOMINANT TRANSITION: the ground doublet cc transition (states with
% E < Esplit), with SPECTATOR-normalized populations (constraint 3), leaving the transverse
% (a,b) spectator AND the non-dominant (|3>) cc paths BARE; the caller (INVZT_SOLVE_POINT)
% additionally builds the crit from E1's dominant-renormalized/rest-bare ctil0. The lattice
% ca/cb mediation (Jt -> emt K, carrying the E1 dJpre = chi_perp*(Vca Vca' + Vcb Vcb') with
% the BARE transverse chi) is untouched. This is the MATCHED truncation for the framework
% 11.8 ODD-sector emergence gate (A3-dominant reduces to A1/E1 up to the constraint-8
% O(1/z^2) resummation-scheme difference); the FULL-vs-dominant crit-shift difference is the
% genuine beyond-E1 transverse-spectator dressing (~11% at the T=2 K/Bx=0.5 T PM anchor).
dress = getf(opts, 'dress', 'full');
if ~ismember(char(dress), {'full','dominant'})
    error('invzt:dress', 'opts.dress must be ''full'' or ''dominant''; got %s.', char(dress));
end
dom = strcmp(char(dress), 'dominant');

% --- (0) centred operators, bare local propagator -------------------------------
% The single-ion STATE space has dimension N = numel(si.E) (3 for the toy/e3, 6 for
% e6, 24 for e3xI8, ... up to 136); the centering identity must be N x N -- NOT the
% external Cartesian 3 x 3 (which is the (mu,nu) tensor index, a separate axis). A
% hardcoded eye(3) here is correct only at the toy dimension and crashes for any
% larger rung (mirrors the dominant branch's eye(numel(di)) below).
N    = numel(si.E);
p    = si.P(:);
Jexp = si.Jexp(:);
es   = struct('E', si.E(:), 'p', p);
ops  = struct('a', si.Mx - Jexp(1)*eye(N), ...
              'b', si.My - Jexp(2)*eye(N), ...
              'c', si.Mz - Jexp(3)*eye(N));
G0 = -invz_chi0z(si, T, 1i*wn, struct('elastic', true));        % [3,3,nwn] (Cartesian)

% --- active-subspace projector of G0 (union over frequencies; shared A2 policy) --
P = active_projector(G0, rtol);
r = size(P, 2);

% --- (3-precompute) connected Gamma4[npair,nc,nc,nwn,nl], Kmat-INDEPENDENT --------
labs  = {'a', 'b', 'c'};  nc = 3;
idx   = struct('a', 1, 'b', 2, 'c', 3);
ext   = {{'a','a'},{'a','b'},{'a','c'},{'b','a'},{'b','b'},{'b','c'},{'c','a'},{'c','b'},{'c','c'}};
npair = numel(ext);
nl    = 2*Lmax + 1;
next  = (0:nwn-1).';
lvals = (-Lmax:Lmax).';
mu_i = zeros(npair, 1);  nu_i = zeros(npair, 1);
for ip = 1:npair, mu_i(ip) = idx.(ext{ip}{1});  nu_i(ip) = idx.(ext{ip}{2}); end
% Gamma4[npair,nc,nc,nwn,nl] -- vectorized over the (n,l) grid (bit-identical to
% invzt_vertex4 stage 'Gamma', gated by test_invzt_a3_threestate). Cached across
% solves at the same toy: the same (T,B,toy) drives every odd-on/off and lambda-scaled
% A3 solve (Gamma4 depends on the single ion only, NOT the lattice).
if dom
    % DOMINANT-TRANSITION dressing (E1): restrict the cc self-energy to the dominant
    % transition -- the ground doublet (states with E < Esplit) -- with SPECTATOR-
    % normalized populations p~ = p/(p_dom sum) (constraint 3). This is the two-level
    % Gamma4 (the |3> multi-level cc paths are dropped, exactly as E1's dominant-only
    % rule does), placed in the cc;cc slot of a full-shaped G4; the transverse a,b
    % spectator stays bare. Reproduces invz_sigma's two-level self-energy at the medium
    % K_cc, i.e. A1/E1, so the framework 11.8 emergence gate is a MATCHED truncation.
    Esplit = getf(opts, 'Esplit', 0.4653);
    di = find(si.E(:) < Esplit);
    if numel(di) < 2
        error('invzt:domGroup', 'dominant group (E < %.4g) has %d < 2 states.', Esplit, numel(di));
    end
    Ed = si.E(di) - si.E(di(1));
    pd = si.P(di);  pd = pd / sum(pd);                 % spectator normalization
    Mzd = si.Mz(di, di);
    ops_d = struct('c', Mzd - real(sum(pd .* diag(Mzd)))*eye(numel(di)));
    G4dom = gamma4_cached(struct('E', Ed, 'p', pd), ops_d, {'c'}, next, lvals, beta);
    G4 = complex(zeros(npair, nc, nc, nwn, nl));
    ipcc = find(mu_i == 3 & nu_i == 3);                % cc;cc slot
    G4(ipcc, 3, 3, :, :) = reshape(G4dom, [1, 1, 1, nwn, nl]);
else
    G4 = gamma4_cached(es, ops, labs, next, lvals, beta);
end

% --- (1)-(3) outer self-consistent loop on Vmat ----------------------------------
Vmat = complex(zeros(3, 3, nwn));
if isfield(opts, 'Vmat_seed') && isequal(size(opts.Vmat_seed), [3 3 nwn])
    Vmat = opts.Vmat_seed;
end
converged = false;
chi_til = complex(zeros(3,3,nwn));  Kmat = complex(zeros(3,3,nwn));
chi_bar = complex(zeros(3,3,nwn));  emtinfo = struct();  outer = 0;
% Outer mixing on Vmat. DEFAULT depth 1 (plain damped mix): the dominant cost is the
% one-time Gamma4 precompute (now vectorized + cached), NOT the iteration count, so each
% outer iteration is cheap (one INVZT_EMT_MATRIX). Anderson (depth > 1) is OPT-IN and
% OFF by default because it was observed to converge to a SPURIOUS A3 self-energy fixed
% point at the exact rho->0 proxy (0.19 vs the correct 0.36 that plain mixing + the
% scalar chain both give) -- the unsafeguarded least-squares step jumps basins. Plain
% damped mixing converges monotonically to the correct fixed point.
mAA = getf(opts, 'anderson_depth', 1);
wstate = warning('off', 'MATLAB:rankDeficientMatrix');
cleaner = onCleanup(@() warning(wstate));  %#ok<NASGU>
Fhist = cell(1, 0);  Xhist = cell(1, 0);
Vnew = Vmat;
for outer = 1:maxo
    chi_til = resum_dyson(G0, Vmat, P);                         % symmetric bracket
    [Kmat, chi_bar, emtinfo] = invzt_emt_matrix(chi_til, lat_eff, struct('rank_tol', rtol));
    Vnew = contract_vertex(G4, Kmat, mu_i, nu_i, Lmax, nwn, nc, beta, dom);
    f  = Vnew - Vmat;
    dV = max(abs(f(:)));
    if dV < tolo, converged = true; break; end
    fv = f(:);  xv = Vnew(:);
    Fhist{end+1} = fv;  Xhist{end+1} = xv;   %#ok<AGROW>
    if numel(Fhist) > mAA, Fhist(1) = [];  Xhist(1) = []; end
    mk = numel(Fhist);
    if mk == 1
        Vmat = Vmat + mixo*f;                                   % damped mix
    else
        DF = zeros(numel(fv), mk-1);  DX = zeros(numel(fv), mk-1);
        for j = 1:mk-1, DF(:,j) = Fhist{j+1} - Fhist{j};  DX(:,j) = Xhist{j+1} - Xhist{j}; end
        gcoef = DF \ fv;                                        % least-squares min residual
        Vv = xv - DX*gcoef;
        if any(~isfinite(Vv)), Vv = Vmat(:) + mixo*fv; end      % safeguard -> damped mix
        Vmat = reshape(Vv, 3, 3, nwn);
    end
end
% At exit chi_til/Kmat/chi_bar are consistent with Vmat to within tol; adopt Vnew as
% the reported self-energy (the vertex output at the converged medium).
Vmat = Vnew;

chi_til_add = -(G0 - Vmat);                                    % additive resummation
Gloc = -real(squeeze(chi_bar(3, 3, :)));

out = struct();
out.Vmat = Vmat;  out.chi_til = chi_til;  out.chi_til_add = chi_til_add;
out.Kmat = Kmat;  out.chi_bar = chi_bar;  out.Gloc = Gloc(:);  out.G0 = G0;
out.emtinfo = emtinfo;  out.converged = converged;  out.outer_iters = outer;
out.active_rank = r;  out.pruned_bound = 0;                     % no matrix-element pruning (exact)
out.dress = char(dress);
end

% ======================================================================= %
%  Symmetric-bracket resummation on the active subspace P:
%     chi_til = -G0 * ((G0 + Vmat) \ G0)   restricted to range(P), 0 on the complement.
% ======================================================================= %
function chi_til = resum_dyson(G0, Vmat, P)
nwn = size(G0, 3);  r = size(P, 2);
chi_til = complex(zeros(3, 3, nwn));
if r == 0, return; end
for k = 1:nwn
    G0p = P' * G0(:,:,k) * P;                   % [r,r]
    Vp  = P' * Vmat(:,:,k) * P;
    G0til_p = G0p * ((G0p + Vp) \ G0p);         % LOCKED symmetric bracket
    chi_til(:,:,k) = -P * G0til_p * P';
end
end

% ----- cheap K-contraction of the precomputed Gamma4 (bit-identical to vertex4 'V') -----
function Vmat = contract_vertex(G4, Kmat, mu_i, nu_i, Lmax, nwn, nc, beta, dom)
nl = 2*Lmax + 1;
% Build Karr[nc,nc,nl] indexed l = -Lmax..Lmax (l>=0 direct; l<0 transpose relation
% K_{rho sigma}(-l) = K_{sigma rho}(+l)) -- Kmat(:,:,m) is Matsubara index m-1 = |l|.
Karr = complex(zeros(nc, nc, nl));
for li = 1:nl
    l = li - Lmax - 1;                          % lvals(li)
    if l >= 0
        Karr(:, :, li) = Kmat(:, :, l + 1);
    else
        Karr(:, :, li) = Kmat(:, :, -l + 1).';  % plain transpose (constraint 9)
    end
end
if dom
    % DOMINANT-ONLY dressing (E1): the vertex self-energy lives in the cc channel ALONE.
    % Restrict BOTH the internal contraction (Karr -> cc component only) AND the external
    % legs (compute V_cc only); every transverse/cross Vmat component stays 0 (bare
    % spectator). V_cc = (1/2beta) sum_l K_cc(l) Gamma4_{cc;cc}(n,l).
    Kd = complex(zeros(nc, nc, nl));  Kd(3,3,:) = Karr(3,3,:);  Karr = Kd;
    iplist = find(mu_i == 3 & nu_i == 3);       % external cc pair only
else
    iplist = 1:numel(mu_i);
end
Kb = reshape(Karr, [nc, nc, 1, nl]);            % broadcast over the frequency axis
Vmat = complex(zeros(3, 3, nwn));
for ip = iplist(:).'
    G4p   = reshape(G4(ip, :, :, :, :), [nc, nc, nwn, nl]);
    contr = sum(sum(sum(G4p .* Kb, 1), 2), 4);  % [1,1,nwn,1]
    Vmat(mu_i(ip), nu_i(ip), :) = reshape(contr, [1, 1, nwn]) / (2*beta);
end
end

% ----- frequency-consistent active-subspace projector (union of ranges) -----
function P = active_projector(X, rtol)
nz = size(X, 3);
G3 = zeros(3);
for k = 1:nz
    Ck = (X(:,:,k) + X(:,:,k)')/2;              % Hermitized channel content
    G3 = G3 + Ck'*Ck;
end
G3 = (G3 + G3')/2;
[V, dv] = eig(G3, 'vector');
dv = real(dv);
mag = sqrt(max(dv, 0));
smax = max(mag);  if smax <= 0, smax = 1; end
active = mag > rtol*smax;
P = V(:, active);
end

% ----- session cache of Gamma4 (keyed+verified by the toy signature) -----
function G4 = gamma4_cached(es, ops, comps, next, lvals, beta)
% Gamma4 depends only on the single ion (es, ops) and the (n,l) grid -- NOT on the
% lattice -- so every odd-on/off and lambda-scaled A3 solve at the SAME (T,B,toy)
% reuses one precompute. Session-persistent, signature-verified (isequal) on hit.
persistent CACHE
if isempty(CACHE), CACHE = struct('key', {}, 'sig', {}, 'G4', {}); end
nc = numel(comps);
ext = cell(1, nc*nc);  c = 0;
for i = 1:nc, for j = 1:nc, c = c + 1; ext{c} = {comps{i}, comps{j}}; end, end
sig = [es.E(:); es.p(:); beta; next(:); lvals(:)];
for i = 1:nc, O = ops.(comps{i}); sig = [sig; real(O(:)); imag(O(:))]; end %#ok<AGROW>
key = sprintf('%dv_%08x', numel(sig), typecast(single(sum(sig.*(1:numel(sig))')), 'uint32'));
for i = 1:numel(CACHE)
    if strcmp(CACHE(i).key, key) && isequal(CACHE(i).sig, sig)
        G4 = CACHE(i).G4;  return;
    end
end
G4 = invzt_gamma4(es, ops, ext, comps, next, lvals, beta);
CACHE(end+1) = struct('key', key, 'sig', sig, 'G4', G4);  %#ok<AGROW>
if numel(CACHE) > 6, CACHE(1) = []; end                    % bounded (each ~18 MB)
end
