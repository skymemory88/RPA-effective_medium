function [dm, T3] = invzt_vertex3(es, ops, wn, beta, opts)
%INVZT_VERTEX3  Component-labelled three-point vertex (sum-rule route, A3 tensor).
%
%   [dm, T3] = INVZT_VERTEX3(es, ops, wn, beta, opts)
%
%   The three-point cumulant that carries the generalized sum-rule / condition-space
%   leg of the tensor 1/z expansion (three_level_1z_extension: three-point L_n / T^3_l).
%   Same divided-difference / KMS kernel machinery as INVZT_VERTEX4, but the two
%   internal legs meet at ONE intermediate index so the ordered simplex reduces to
%   the FIRST divided difference I2 (two-node) rather than I3.
%
%   For external operator O^mu and internal pair (O^rho, O^sigma) the frequency-
%   resolved three-point is (KMS-scaled, folding the initial log-weight logpr):
%
%     T3(l) = sum_{(e1,e2) in {(+z_l,-z_l),(-z_l,+z_l)}} sum_{r,s,t}
%                 p_r (O^rho)_{rs} (O^sigma)_{st} (O^mu)_{tr}
%                 * I2( E_r-E_s+e1, E_s-E_t+e2 ),          z_l = 2*pi*i*l/beta
%
%   (the two etas are the two orderings of the internal legs; O^mu closes the loop
%   at zero external frequency -- the response point).  The equal-time moment
%
%     dm = sum_r p_r (O^rho O^sigma O^mu)_{rr} = <O^rho O^sigma O^mu>_0
%
%   is the contact term the l-sum obeys as a sum rule (a consistency diagnostic,
%   never a closure -- hard constraint 6).
%
%   INPUTS
%     es   : struct es.E (N x 1), es.p (N x 1)  (as INVZT_VERTEX4).
%     ops  : struct of named centred operators.
%     wn   : vector of internal bosonic indices l at which to evaluate T3.
%     beta : inverse temperature.
%     opts : .quad = {mu, rho, sigma} channel labels (external, internal, internal).
%
%   OUTPUTS
%     dm : scalar equal-time three-point moment <O^rho O^sigma O^mu>_0.
%     T3 : numel(wn) x 1 complex, the three-point vertex at each l in wn.
%
%   Regular at exact degeneracy (Bx=0 doublet) by the exprel Hermite limits of I2;
%   validated in test_invzt_vertex against independent Gauss-Legendre simplex
%   quadrature (the same methodology the oracle uses for the four-point F).
%
%   See also INVZT_VERTEX4, INVZT_KERNELS.

if nargin < 5 || isempty(opts), opts = struct(); end
k    = invzt_kernels();
E    = es.E(:);
p    = es.p(:);
logp = log(p);
N    = numel(E);
q    = opts.quad;
Omu  = ops.(q{1});           % external (closes loop)
Orho = ops.(q{2});           % internal leg 1
Osig = ops.(q{3});           % internal leg 2

% equal-time moment  dm = <O^rho O^sigma O^mu>_0 = sum_r p_r (Orho Osig Omu)_{rr}
dm = sum(p .* diag(Orho * Osig * Omu));

lvals = wn(:);
T3 = complex(zeros(numel(lvals), 1));
Et = E(:).';                  % 1 x N row over t
for il = 1:numel(lvals)
    zl = 1i * 2 * pi * lvals(il) / beta;
    etas = [zl, -zl; -zl, zl];
    tot = complex(0);
    for ie = 1:2
        e1 = etas(ie, 1); e2 = etas(ie, 2);
        for r = 1:N
            lpr = logp(r);
            Omur = Omu(:, r).';                 % 1 x N: (O^mu)_{t r}
            for s = 1:N
                Ors = Orho(r, s);
                if Ors == 0, continue; end
                x  = E(r) - E(s) + e1;
                wt = Ors * (Osig(s, :) .* Omur); % 1 x N over t
                nz = wt ~= 0;
                if ~any(nz), continue; end
                yt  = E(s) - Et(nz) + e2;        % 1 x nnz over t
                ker = k.I2s(lpr, x, yt, beta);
                tot = tot + sum(wt(nz) .* ker);
            end
        end
    end
    T3(il) = tot;
end
end
