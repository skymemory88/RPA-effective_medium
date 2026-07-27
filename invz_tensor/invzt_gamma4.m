function G4 = invzt_gamma4(es, ops, ext, comps, next, lvals, beta)
%INVZT_GAMMA4  Vectorized connected four-point cumulant Gamma4 over the (n,l) grid.
%
%   G4 = INVZT_GAMMA4(es, ops, ext, comps, next, lvals, beta) returns the connected
%   cumulant Gamma4_{mu nu; rho sigma}(n, l) as a [npair, nc, nc, nwn, nl] array, for
%   every external pair ext{ip} = {mu,nu}, internal channels (comps{ri}, comps{si}),
%   external bosonic index n in next, and internal index l in lvals.  This is the
%   A3 self-energy's Kmat-INDEPENDENT ingredient (INVZT_SIGMA_TENSOR precomputes it
%   ONCE and only the cheap K-contraction repeats each outer iteration).
%
%   It is bit-identical (to ~1e-12) to calling INVZT_VERTEX4 stage 'Gamma' row by row,
%   but VECTORIZED over the whole (n,l) grid: the matrix-element path sums (the O(N^4)
%   part) are (n,l)-INDEPENDENT -- only the divided-difference kernel arguments a1/a2/a3
%   carry the external/internal Matsubara shifts z_n, z_l -- so the path is walked once
%   and its Hermite kernel evaluated across the grid in one shot.  This is the ~100x
%   speedup A3 needs when the toy operators are dense (rho != 0).  The kernels reuse the
%   VERIFIED exprel-stable phis/dphis/d2phis handles of INVZT_KERNELS (already array-safe);
%   only the I2/I3 divided-difference combinators are re-expressed array-safe here (the
%   sole change vs invzt_kernels is out(small)=ser(SMALL), needed because x,y are now
%   arrays too).  Equivalence with INVZT_KERNELS is the invariant to preserve here.
%
%   Frequency assignment (LOCKED oracle convention, mirrors INVZT_VERTEX4): the four
%   legs carry O^mu:+z_n, O^nu:-z_n(pinned tau=0), O^rho:+z_l, O^sigma:-z_l; O^nu is
%   pinned and the other three (operator,frequency) legs are summed over their six
%   descending-time orderings (each an ordered simplex == the divided difference I3).
%   Gamma4 = F - beta C_{mu nu}(z_n) C_{rho sigma}(z_l) - delta_{n,-l} beta C_{mu rho} C_{sigma nu}
%              - delta_{n,l} beta C_{mu sigma} C_{rho nu}.
%
%   See also INVZT_VERTEX4, INVZT_KERNELS, INVZT_SIGMA_TENSOR.
k    = invzt_kernels();
E    = es.E(:);
p    = es.p(:);
logp = log(p);
N    = numel(E);
nc   = numel(comps);
npair = numel(ext);
nwn  = numel(next);
nl   = numel(lvals);

% operator lookup by channel label
opmap = struct();
for i = 1:nc, opmap.(comps{i}) = ops.(comps{i}); end

% grid frequency shifts (column vectors over the flattened [nwn x nl] grid, n-major)
[NN, LL] = ndgrid(next(:), lvals(:));
ng = nwn*nl;
zn = 1i*2*pi*NN(:)/beta;              % +z_n  [ng]
zl = 1i*2*pi*LL(:)/beta;              % +z_l  [ng]
frq = {zn, zl, -zl};                  % leg frequencies: {Om:+zn, Or:+zl, Os:-zl}
perms = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];

G4 = complex(zeros(npair, nc, nc, nwn, nl));
for ip = 1:npair
    Om = opmap.(ext{ip}{1});
    On = opmap.(ext{ip}{2});
    for ri = 1:nc
        Or = opmap.(comps{ri});
        for si = 1:nc
            Os = opmap.(comps{si});
            F = fourpoint_grid(k, E, logp, Om, On, Or, Os, beta, frq, perms, N, ng);
            F = F - pairings_grid(E, p, Om, On, Or, Os, beta, NN(:), LL(:));
            G4(ip, ri, si, :, :) = reshape(F, [1, 1, 1, nwn, nl]);
        end
    end
end
end

% ======================================================================= %
%  Four-point path sum F over the (n,l) grid (vectorised over the grid; the
%  O(N^4) matrix-element walk is grid-independent).  triples = {Om,Or,Os} with
%  leg frequencies {zn, zl, -zl}; On is pinned.
% ======================================================================= %
function F = fourpoint_grid(k, E, logp, Om, On, Or, Os, beta, frq, perms, N, ng)
op3 = {Om, Or, Os};
F = complex(zeros(ng, 1));
for ipm = 1:6
    A  = op3{perms(ipm,1)};  x0 = frq{perms(ipm,1)};
    B  = op3{perms(ipm,2)};  x1 = frq{perms(ipm,2)};
    Cm = op3{perms(ipm,3)};  x2 = frq{perms(ipm,3)};
    for r = 1:N
        lpr = logp(r);
        for s = 1:N
            Ars = A(r, s);
            if Ars == 0, continue; end
            d1 = E(r) - E(s);
            a1 = d1 + x0;                       % [ng]
            for t = 1:N
                Bst = B(s, t);
                if Bst == 0, continue; end
                d2 = E(s) - E(t);
                a2 = d2 + x1;                   % [ng]
                AB = Ars * Bst;
                for u = 1:N
                    w = AB * Cm(t, u) * On(u, r);
                    if w == 0, continue; end
                    a3 = (E(t) - E(u)) + x2;    % [ng]
                    F = F + w * I3s_a(k, lpr, a1, a2, a3, beta);
                end
            end
        end
    end
end
end

% ----- Wick pairings over the (n,l) grid (subtracted from F -> connected Gamma4) -----
function P = pairings_grid(E, p, Om, On, Or, Os, beta, nvec, lvec)
% main pairing beta*C_{mu nu}(n) C_{rho sigma}(l): factorises over n and l.
un = unique(nvec);  ul = unique(lvec);
Cmn_u = twopoint_vec(E, p, Om, On, beta, un);
Crs_u = twopoint_vec(E, p, Or, Os, beta, ul);
[~, in] = ismember(nvec, un);
[~, il] = ismember(lvec, ul);
P = beta * Cmn_u(in) .* Crs_u(il);
% delta_{n,-l} beta C_{mu rho}(n) C_{sigma nu}(n)
m1 = (nvec == -lvec);
if any(m1)
    nn = nvec(m1);
    P(m1) = P(m1) + beta * twopoint_vec(E, p, Om, Or, beta, nn) .* twopoint_vec(E, p, Os, On, beta, nn);
end
% delta_{n,l} beta C_{mu sigma}(n) C_{rho nu}(n)
m2 = (nvec == lvec);
if any(m2)
    nn = nvec(m2);
    P(m2) = P(m2) + beta * twopoint_vec(E, p, Om, Os, beta, nn) .* twopoint_vec(E, p, Or, On, beta, nn);
end
end

% ----- two-point Lehmann C_{Oa Ob}(i omega_m) for an ARRAY of indices m -----
function C = twopoint_vec(E, p, Oa, Ob, beta, marr)
marr = marr(:);
wm = 2*pi*marr/beta;
C = complex(zeros(numel(marr), 1));
N = numel(E);
is0 = (marr == 0);
for r = 1:N
    for s = 1:N
        a = Oa(r, s);  b = Ob(s, r);
        if a == 0 || b == 0, continue; end
        d = E(r) - E(s);
        if abs(d) < 1e-12
            C(is0) = C(is0) + beta * p(r) * a * b;      % elastic, only at m = 0
        else
            C = C + (p(s) - p(r)) * a * b ./ (1i*wm + d);
        end
    end
end
end

% ======================================================================= %
%  Array-safe divided-difference kernels.  phis/dphis/d2phis (already array-safe)
%  come from the VERIFIED invzt_kernels handles; only the I2/I3 combinators are
%  re-expressed here so x AND y may be arrays (out(small)=ser(SMALL), the sole change
%  vs invzt_kernels, which assumes scalar x,y).  Return p_r * I3 (logpr folded).
% ======================================================================= %
function out = I3s_a(k, logpr, x, y, z, beta)
out = (I2h_a(k, logpr, x, y + z, beta) - I2h_a(k, logpr, x, y, beta)) ./ z;
serz = dI2y_a(k, logpr, x, y, beta);
small = abs(z) <= 1e-10;
out(small) = serz(small);
end

function out = I2h_a(k, logpr, x, Y, beta)
out = (k.phis(logpr, x + Y, beta) - k.phis(logpr, x, beta)) ./ Y;
ser = k.dphis(logpr, x, beta) + Y .* k.d2phis(logpr, x, beta) ./ 2;
small = abs(Y) <= 1e-10;
out(small) = ser(small);
end

function out = dI2y_a(k, logpr, x, Y, beta)
out = (k.dphis(logpr, x + Y, beta) ...
       - (k.phis(logpr, x + Y, beta) - k.phis(logpr, x, beta)) ./ Y) ./ Y;
ser = k.d2phis(logpr, x, beta) ./ 2;
small = abs(Y) <= 1e-10;
out(small) = ser(small);
end
