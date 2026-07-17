function med = invz_emt_scalar(G0, Sigma, Jnu_flat, opts)
%INVZ_EMT_SCALAR Scalar effective-medium fixed point at fixed self-energy.
% Solves the scalar effective-medium equations (R eqs 7-9):
%   G  = G0./(1 + Sigma + K.*G0)            (Dyson, eq 9)
%   Gq = G ./ (1 + (J_nu - K).*G)           (eq 7)
%   K  = (1/N) sum_{q,nu} J_nu.*Gq ./ G     (eq 8)
% K cancels out of Gq algebraically, so the fixed point is solved in closed
% form (see code) rather than iterated; consequently R eq 8's closure
% (mean_q Gq = G) is machine-zero whenever a solution exists.
%
% Jnu_flat: [nJ,1] static mode spectrum (one BZ average shared by every
% frequency, the original interface), OR [nJ,nw] per-frequency spectra (T2.1
% retardation): column n is the mode spectrum entering frequency n's BZ
% average, paired with the caller's Matsubara grid in invz_matsubara order
% (wn(1) = 0 first). Any vector input takes the original code path unchanged;
% a matrix with identical columns is bitwise-equal to the vector path (each
% element sees the same multiply/divide, no accumulation differences --
% regression-gated by test_invz_odd_retarded/test_emt_matrix_column_constant_bitwise).
if nargin < 4, opts = struct(); end
% tol/max_iter/mix/K0 are accepted for backward compatibility but unused: the
% solve is direct, so there is nothing to iterate, seed, or damp.
blk = 4096; if isfield(opts,'freq_block') && ~isempty(opts.freq_block), blk = opts.freq_block; end
G0 = G0(:);  Sigma = Sigma(:);
D = 1 + Sigma;  nw = numel(D);
retMat = ~isvector(Jnu_flat);          % [nJ,nw] per-frequency branch (T2.1)
if retMat
    if size(Jnu_flat, 2) ~= nw
        error('invz:emtJnu', ...
            'Jnu_flat matrix must be [nJ, nw]: got %d columns for nw = %d frequencies.', ...
            size(Jnu_flat, 2), nw);
    end
    Jf = Jnu_flat;
else
    Jf = Jnu_flat(:);
end

% A(w) = mean_q J./(D + J.*G0), evaluated in frequency blocks so the [nJ x nw]
% implicit-expansion temporary never materializes at full width.
A = zeros(nw,1);
if retMat
    for i0 = 1:blk:nw
        idx = i0:min(i0+blk-1, nw);
        A(idx) = mean(Jf(:,idx) ./ (D(idx).' + Jf(:,idx) .* G0(idx).'), 1).';
    end
else
for i0 = 1:blk:nw
    idx = i0:min(i0+blk-1, nw);
    A(idx) = mean(Jf ./ (D(idx).' + Jf*G0(idx).'), 1).';
end
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
    if retMat
        for i0 = 1:blk:nw
            idx = i0:min(i0+blk-1, nw);
            Gqb = G(idx).' ./ (1 + (Jf(:,idx) - K(idx).').*G(idx).');    % [nJ, numel(idx)]
            cl(idx) = abs(mean(Gqb,1).' - G(idx)) ./ max(abs(G(idx)), 1e-14);
        end
    else
    for i0 = 1:blk:nw
        idx = i0:min(i0+blk-1, nw);
        Gqb = G(idx).' ./ (1 + (Jf - K(idx).').*G(idx).');    % [nJ, numel(idx)]
        cl(idx) = abs(mean(Gqb,1).' - G(idx)) ./ max(abs(G(idx)), 1e-14);
    end
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
