function med = invz_emt_scalar(G0, Sigma, Jnu_flat, opts)
%INVZ_EMT_SCALAR Scalar effective-medium fixed point at fixed self-energy.
% G  = G0./(1 + Sigma + K.*G0)            (Dyson, R eq 9)
% Gq = G ./ (1 + (J_nu - K).*G)           (R eq 7)
% K  = (1/N) sum_{q,nu} J_nu.*Gq ./ G     (R eq 8)
%
% Closed form: K cancels out of Gq, so the fixed point of R eqs 7-9 at fixed
% Sigma is solved directly instead of iterated. With D = 1 + Sigma,
%   Gq = G0./(D + J.*G0)                       (K-free: substitute G into Gq)
%   A  = (1/N) sum_q  J ./ (D + J.*G0)         (K-independent)
%   K  = A.*D ./ (1 - A.*G0)                   (linear in K -> exact)
%   G  = G0./(D + K.*G0)
% At this fixed point R eq 8 also gives mean_q Gq = G exactly, so the closure
% diagnostic is machine-zero when a solution exists.
if nargin < 4, opts = struct(); end
% tol/max_iter/mix/K0 are accepted for backward compatibility but unused: the
% solve is direct, so there is nothing to iterate, seed, or damp.
blk = 4096; if isfield(opts,'freq_block') && ~isempty(opts.freq_block), blk = opts.freq_block; end
G0 = G0(:);  Sigma = Sigma(:);  Jf = Jnu_flat(:);
D = 1 + Sigma;  nw = numel(D);

% A(w) = mean_q J./(D + J.*G0), evaluated in frequency blocks so the [nJ x nw]
% implicit-expansion temporary never materializes at full width.
A = zeros(nw,1);
for i0 = 1:blk:nw
    idx = i0:min(i0+blk-1, nw);
    A(idx) = mean(Jf ./ (D(idx).' + Jf*G0(idx).'), 1).';
end
K = A .* D ./ (1 - A .* G0);
G = G0 ./ (D + K .* G0);

med.G = G;  med.K = K;
% Closure diagnostic (R eq 8: q-averaged Gq must equal the local G). The closed
% form makes this machine-zero whenever a solution exists, so it never affects
% the returned G/K -- but it costs a second full-width [nJ x nw] pass, doubling
% the solve. Compute it only when opts.debug is set (the tests do); production
% detects a missing solution through med.converged below, not through closure.
if isfield(opts,'debug') && ~isempty(opts.debug) && opts.debug
    cl = zeros(nw,1);
    for i0 = 1:blk:nw
        idx = i0:min(i0+blk-1, nw);
        Gqb = G(idx).' ./ (1 + (Jf - K(idx).').*G(idx).');    % [nJ, numel(idx)]
        cl(idx) = abs(mean(Gqb,1).' - G(idx)) ./ max(abs(G(idx)), 1e-14);
    end
    med.closure = max(cl);
else
    med.closure = NaN;                                        % not computed unless opts.debug
end
% Singular medium (1 - A.*G0 -> 0, i.e. a vanishing RPA denominator) leaves K
% or G non-finite: report that as not-converged so invz_solve_point treats the
% point as having no paramagnetic solution. The old iteration masked this
% because max() ignores the NaN residual it produced.
med.converged = all(isfinite(G)) && all(isfinite(K));
med.iters = 1;
end
