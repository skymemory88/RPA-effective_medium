function out = invzf_taylor_extract(f, K, opts)
%INVZF_TAYLOR_EXTRACT Certified Taylor coefficients of a real-analytic function of one coupling.
%
%   OUT = INVZF_TAYLOR_EXTRACT(F, K, OPTS) returns the coefficients a_0..a_K of
%   F(J) = sum_k a_k J^k, together with a MEASURED extraction error for each
%   coefficient. F must be a handle accepting a real scalar J and returning a real
%   scalar; it is assumed analytic in a disc containing the sampling radius.
%
%   This is the differentiation half of the exact-cluster coefficient gate
%   (invzp_convg_fix.md Sec. 6): "compute its coefficient by exact finite-cluster
%   differentiation". The cluster oracle it is applied to is exact at every J, so
%   the ONLY error here is the extraction, and it must be measured rather than
%   assumed -- a coefficient comparison is worthless if the comparator's own error
%   is unknown.
%
%   Method. F is sampled at Chebyshev-Lobatto nodes on [-r, r], a degree-(K+PAD)
%   polynomial is fitted in the scaled variable t = J/r, and a_k = b_k / r^k. The
%   whole extraction is repeated at several radii r0, r0/2, r0/4, ... The three
%   estimates of a_k differ by (i) truncation of the fit, which shrinks with r, and
%   (ii) roundoff amplified by r^-k, which grows as r shrinks. Their SPREAD is
%   therefore an honest upper bound on the extraction error, and it is reported
%   per coefficient rather than as a single global figure.
%
%   OPTS fields (defaults):
%     r0        (0.05)   largest sampling radius in J
%     n_radii   (3)      number of successive halvings
%     pad       (4)      extra fitted degrees above K, to absorb the tail
%     oversample(2)      node count = oversample*(K+pad)+1
%
%   OUT:
%     .a          (K+1)x1 coefficients (median over radii)
%     .err        (K+1)x1 measured spread over radii -- the extraction error bar
%     .a_by_radius(K+1)x n_radii, every estimate, unaveraged
%     .radii      the radii used
%     .fit_resid  max relative residual of the polynomial fits
%     .cond       max condition number of the fitted Vandermonde systems
%     .n_evals    number of F evaluations
if nargin < 3, opts = struct(); end
r0    = getf_local(opts, 'r0', 0.05);
nrad  = getf_local(opts, 'n_radii', 3);
pad   = getf_local(opts, 'pad', 4);
osamp = getf_local(opts, 'oversample', 2);
validateattributes(K, {'numeric'}, {'scalar','integer','>=',0});
validateattributes(r0, {'numeric'}, {'scalar','real','finite','positive'});

deg = K + pad;
nn  = osamp*deg;                       % nn+1 Chebyshev-Lobatto nodes
tt  = cos(pi*(0:nn).'/nn);             % in [-1, 1]
V   = tt.^(0:deg);                     % Vandermonde in the SCALED variable
radii = r0 * 2.^-(0:nrad-1);

A = nan(K+1, nrad);
fit_resid = 0;  cnd = 0;  nev = 0;
for ir = 1:nrad
    r = radii(ir);
    y = zeros(nn+1, 1);
    for j = 1:nn+1
        y(j) = f(r*tt(j));
        nev = nev + 1;
    end
    if ~all(isfinite(y))
        error('invzf:taylorExtract', ...
            'F returned a non-finite value at radius %g; the sampling disc is outside its domain.', r);
    end
    b = V\y;
    cnd = max(cnd, cond(V));
    yscale = max(abs(y));
    fit_resid = max(fit_resid, norm(V*b - y, inf)/max(yscale, realmin));
    A(:, ir) = b(1:K+1) ./ (r.^(0:K).');
end

a   = median(A, 2);
err = max(abs(A - a), [], 2);

out = struct('a', a, 'err', err, 'a_by_radius', A, 'radii', radii, ...
    'fit_resid', fit_resid, 'cond', cnd, 'n_evals', nev, ...
    'degree_fitted', deg, 'n_nodes', nn+1);
end

function v = getf_local(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
