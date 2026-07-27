function out = invzt_vertex4(es, ops, Kmat, wn, beta, opts)
%INVZT_VERTEX4  Component-labelled four-point vertex engine (A3 tensor 1/z).
%
%   out = INVZT_VERTEX4(es, ops, Kmat, wn, beta, opts)
%
%   The DENSE ordered path-sum reference-of-record for the tensor 1/z self-energy
%   correction V_{mu nu}(i omega_n).  Reproduces the mpmath oracle of
%   ../verify_tensor_vertex.py (Task 10) to ~1e-14 on every
%   fixture system.  The engine has one PRODUCTION entry (stage 'V', the contracted
%   correction) plus three intermediate stages ('F','Gamma','C') that expose the
%   un-contracted four-/two-point objects the oracle stores row-by-row.
%
%   INPUTS
%     es    : struct, es.E (N x 1 real, ground-shifted energies), es.p (N x 1
%             Boltzmann populations, sum 1).  logp = log(es.p) is formed once and
%             folded into every kernel (hard constraint 4, KMS boundary layers).
%     ops   : struct of named centred Hermitian operators (fields are channel
%             labels, e.g. ops.a, ops.b, ops.c), each N x N.
%     Kmat  : internal-channel kernel K_{rho sigma}(i omega_l)  (stage 'V' only):
%               * function handle  @(ri,si,l) -> scalar  (1-based ri,si; any l), OR
%               * numeric array [nc,nc,nl] indexed l = 0..nl-1 for l >= 0; negative
%                 l is reconstructed by the LOCKED transpose relation (constraint 9)
%                 K_{rho sigma}(-i omega_l) = K_{sigma rho}(+i omega_l).
%             Pass [] for stages F/Gamma/C.
%     wn    : vector of external bosonic indices n (stage 'V').  Ignored for
%             F/Gamma/C (which read (n,l) from opts.nl).
%     beta  : inverse temperature (meV^-1).
%     opts  : struct (all optional)
%             .impl  'dense' (default) | 'factored'.  'factored' is DISABLED --
%                    Task-10 check 10 recorded factored_ok = false (naive O(N^3)
%                    resolvent chain drops the phi/KMS + initial-state structure);
%                    selecting it errors invzt:factoredUnproven.
%             .stage 'V' (default) | 'F' | 'Gamma' | 'C'
%             .quad  F/Gamma: cellstr {mu,nu,rho,sigma}; C: {rho,sigma} (labels in ops)
%             .nl    F/Gamma: M x 2 array of (n,l); C: M x 1 vector of n
%             .ext   V: cell of {mu,nu} label pairs (default: all ordered pairs)
%             .comps V: cellstr of internal channel labels (default: all ops fields)
%             .Lmax  V: internal Matsubara cutoff for the l-sum (default 200)
%             .mtol  matrix-element pruning tolerance (default 0 = no pruning);
%                    dropped |element| < mtol, with the discarded weight bounded in
%                    out.pruned_bound.
%
%   OUTPUT (struct)
%     stage 'F' / 'Gamma' : out.val = M x 1 complex (one per opts.nl row); out.tags
%     stage 'C'           : out.val = M x 1 complex two-point C_{rho sigma}(i w_n)
%     stage 'V'           : out.val = [npairs x numel(wn)] complex; out.pairs (cell),
%                           out.wn, out.pruned_bound
%
%   Frequency assignment (LOCKED oracle convention): the four external legs of
%   F_{mu nu; rho sigma}(n,l) carry  O^mu:+z_n  O^nu:-z_n  O^rho:+z_l  O^sigma:-z_l,
%   z_k = 2*pi*i*k/beta.  O^nu is pinned at tau=0 and the other three
%   (operator,frequency) pairs are summed over their six descending-time orderings,
%   each an ordered simplex == the divided-difference kernel I3 (invzt_kernels).
%
%     F_{mu nu;rho sigma}(n,l) = sum_{perm in S3} sum_{r,s,t,u}
%          p_r (O^{pi1})_{rs}(O^{pi2})_{st}(O^{pi3})_{tu}(O^nu)_{ur}
%          * I3( E_r-E_s+xi1, E_s-E_t+xi2, E_t-E_u+xi3 )
%
%   Connected cumulant (Gamma4 = F - Wick pairings):
%     Gamma4 = F - beta C_{mu nu}(z_n) C_{rho sigma}(z_l)
%                - beta delta_{n,-l} C_{mu rho}(z_n) C_{sigma nu}(z_n)
%                - beta delta_{n, l} C_{mu sigma}(z_n) C_{rho nu}(z_n).
%   Contraction (constraint 2, additive correction 𝒱 = G0.Sigma; never a ratio):
%     V_{mu nu}(n) = (1/(2 beta)) sum_l sum_{rho sigma} K_{rho sigma}(l) Gamma4(n,l).
%
%   See also INVZT_KERNELS, INVZT_VERTEX3.

if nargin < 6 || isempty(opts), opts = struct(); end
impl  = getfield_default(opts, 'impl',  'dense');
stage = getfield_default(opts, 'stage', 'V');
mtol  = getfield_default(opts, 'mtol',  0);

if strcmp(impl, 'factored')
    % Task-10 verdict: FACTORIZATION NOT ESTABLISHED (fixture factored.factored_ok
    % = false).  The O(N^3) chained-resolvent form is an unproven conjecture that
    % drops the initial-state-referenced KMS structure of I3; the dense path-sum is
    % the sole production engine.
    error('invzt:factoredUnproven', ...
        ['opts.impl=''factored'' is disabled: Task-10 check 10 recorded ', ...
         'factored_ok = false (naive O(N^3) resolvent chain drops the phi/KMS ', ...
         'boundary and initial-state structure). Use the dense engine.']);
end
if ~strcmp(impl, 'dense')
    error('invzt:badImpl', 'Unknown opts.impl ''%s'' (use ''dense'').', impl);
end

k    = invzt_kernels();
E    = es.E(:);
p    = es.p(:);
logp = log(p);

switch stage
    case {'F', 'Gamma'}
        q = opts.quad;
        Om = ops.(q{1}); On = ops.(q{2}); Or = ops.(q{3}); Os = ops.(q{4});
        nl = opts.nl;
        m  = size(nl, 1);
        val = complex(zeros(m, 1));
        for i = 1:m
            n = nl(i, 1); l = nl(i, 2);
            val(i) = fourpoint(k, E, logp, Om, On, Or, Os, beta, n, l, mtol);
            if strcmp(stage, 'Gamma')
                val(i) = val(i) - pairings(E, p, Om, On, Or, Os, beta, n, l);
            end
        end
        out.val  = val;
        out.tags = opts.quad;

    case 'C'
        q  = opts.quad;
        Oa = ops.(q{1}); Ob = ops.(q{2});
        nn = opts.nl(:);
        val = complex(zeros(numel(nn), 1));
        for i = 1:numel(nn)
            val(i) = twopoint(E, p, Oa, Ob, beta, nn(i));
        end
        out.val = val;

    case 'V'
        labels = fieldnames(ops);
        comps  = getfield_default(opts, 'comps', labels);
        ext    = getfield_default(opts, 'ext',   allpairs(labels));
        Lmax   = getfield_default(opts, 'Lmax',  200);
        nc     = numel(comps);
        Ccell  = cell(nc, 1);
        for i = 1:nc, Ccell{i} = ops.(comps{i}); end
        Kget = make_Kget(Kmat);
        npair = numel(ext);
        nn    = wn(:).';
        val   = complex(zeros(npair, numel(nn)));
        pruned = 0;
        for ip = 1:npair
            Om = ops.(ext{ip}{1});
            On = ops.(ext{ip}{2});
            for iw = 1:numel(nn)
                n = nn(iw);
                acc = complex(0);
                for l = -Lmax:Lmax
                    for ri = 1:nc
                        for si = 1:nc
                            Kv = Kget(ri, si, l);
                            if Kv == 0, continue; end
                            [F, pb] = fourpoint(k, E, logp, Om, On, Ccell{ri}, Ccell{si}, beta, n, l, mtol);
                            g = F - pairings(E, p, Om, On, Ccell{ri}, Ccell{si}, beta, n, l);
                            acc = acc + Kv * g;
                            pruned = pruned + abs(Kv) * pb;
                        end
                    end
                end
                val(ip, iw) = acc / (2 * beta);
            end
        end
        out.val   = val;
        out.pairs = ext;
        out.wn    = nn;
        out.pruned_bound = pruned / (2 * beta);

    otherwise
        error('invzt:badStage', 'Unknown opts.stage ''%s''.', stage);
end
end

% ======================================================================= %
%  Four-point path sum F_{mu nu; rho sigma}(n,l)  (dense, KMS-folded, O(N^4)).
%  Vectorised over the innermost loop index u.
% ======================================================================= %
function [F, pruned] = fourpoint(k, E, logp, Om, On, Or, Os, beta, n, l, mtol)
N  = numel(E);
zn = 1i * 2 * pi * n / beta;
zl = 1i * 2 * pi * l / beta;
triples = {Om, zn; Or, zl; Os, -zl};       % the three permuted (operator,freq) legs
perms = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];
Eu = E(:).';                                % 1 x N row over u
F  = complex(0);
pruned = 0;
for ipm = 1:6
    A  = triples{perms(ipm,1), 1}; x0 = triples{perms(ipm,1), 2};
    B  = triples{perms(ipm,2), 1}; x1 = triples{perms(ipm,2), 2};
    Cm = triples{perms(ipm,3), 1}; x2 = triples{perms(ipm,3), 2};
    for r = 1:N
        lpr = logp(r);
        Onr = On(:, r).';                   % 1 x N row: O^nu(u,r)
        for s = 1:N
            Ars = A(r, s);
            if Ars == 0, continue; end
            if mtol > 0 && abs(Ars) < mtol, pruned = pruned + abs(Ars); continue; end
            a1 = E(r) - E(s) + x0;
            for t = 1:N
                Bst = B(s, t);
                if Bst == 0, continue; end
                if mtol > 0 && abs(Bst) < mtol, pruned = pruned + abs(Ars*Bst); continue; end
                a2 = E(s) - E(t) + x1;
                wu = (Ars * Bst) * (Cm(t, :) .* Onr);      % 1 x N over u
                if mtol > 0
                    drop = abs(wu) < mtol * abs(Ars * Bst);
                    if any(drop), pruned = pruned + sum(abs(wu(drop))); wu(drop) = 0; end
                end
                nz = wu ~= 0;
                if ~any(nz), continue; end
                a3 = E(t) - Eu(nz) + x2;                   % 1 x nnz over u
                ker = k.I3s(lpr, a1, a2, a3, beta);
                F = F + sum(wu(nz) .* ker);
            end
        end
    end
end
end

% ----- two-point Lehmann C_{Oa Ob}(i omega_n) -----
function C = twopoint(E, p, Oa, Ob, beta, n)
N  = numel(E);
wn = 2 * pi * n / beta;
C  = complex(0);
for r = 1:N
    for s = 1:N
        a = Oa(r, s); b = Ob(s, r);
        if a == 0 || b == 0, continue; end
        d = E(r) - E(s);
        if abs(d) < 1e-12 && n == 0
            C = C + beta * p(r) * a * b;
        elseif abs(d) > 1e-12
            C = C + (p(s) - p(r)) * a * b / (1i * wn + d);
        end
    end
end
end

% ----- the three Wick pairings (subtracted from F to give the connected Gamma4) -----
function t = pairings(E, p, Om, On, Or, Os, beta, n, l)
t = beta * twopoint(E, p, Om, On, beta, n) * twopoint(E, p, Or, Os, beta, l);
if n == -l
    t = t + beta * twopoint(E, p, Om, Or, beta, n) * twopoint(E, p, Os, On, beta, n);
end
if n == l
    t = t + beta * twopoint(E, p, Om, Os, beta, n) * twopoint(E, p, Or, On, beta, n);
end
end

% ----- internal-channel kernel accessor (function handle OR array + transpose) -----
function Kget = make_Kget(Kmat)
if isa(Kmat, 'function_handle')
    Kget = Kmat;
elseif isnumeric(Kmat) && ~isempty(Kmat)
    % [nc,nc,nl], l = 0..nl-1; negative l via transpose relation (constraint 9).
    Kget = @(ri, si, l) arr_K(Kmat, ri, si, l);
else
    error('invzt:badKmat', 'Kmat must be a function handle or a [nc,nc,nl] array.');
end
end

function v = arr_K(Kmat, ri, si, l)
if l >= 0
    v = Kmat(ri, si, l + 1);
else
    v = Kmat(si, ri, -l + 1);      % K_{rho sigma}(-l) = K_{sigma rho}(+l)
end
end

% ----- misc helpers -----
function pr = allpairs(labels)
n = numel(labels);
pr = cell(1, n*n);
c = 0;
for i = 1:n
    for j = 1:n
        c = c + 1;
        pr{c} = {labels{i}, labels{j}};
    end
end
end

function v = getfield_default(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
