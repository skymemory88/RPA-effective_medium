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
% element sees the same multiply/divide, no accumulation differences -- that
% column-constant bitwise equality is the invariant to preserve).
%
% opts.static_medium ('resummed' default): 'strict_1z_dyson_ref' / 'strict_1z_bare_ref' replace
% the omega_n = 0 slot with the one-shot moment closure K0 = Jbar - mu2*Gref (spec SS4.2). Set it
% ONCE via the public opts.static_medium and let invz_check_static_medium stamp it here -- never
% by hand, or the two legs can diverge in truncation order. opts.Jmom (optional) is the
% invz_coupling_moments struct for the static column; derived here when absent. opts.ref_margin
% (1e-6) is the reference-denominator floor. med.medium/med.medium_status report the reference
% and closure outcome; a domain event leaves K(1)/G(1) = NaN and med.converged false, WITH a
% status that says why. med.dynamic_converged checks slots 2:end only for ordered callers that
% replace slot 1; PM callers continue to use med.converged.
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

% --- strict-O(1/z) static slot: omega_n = 0 ONLY (spec SS4.2, SSB) -------------------------
% The closed-form solve above resums mean_q 1/(D + J*G0) to all orders. That resummation's
% feedback into K exceeds retained order and carries the q-denominator that dies below Bc, so
% under a strict scheme slot 1 is replaced by the one-shot moment closure. K(2:end)/G(2:end)
% and both Jnu branches above are untouched: this is a single-slot substitution, and the
% resulting O(1/z^2) mismatch between K(1) and K(2) is a DOCUMENTED artifact measured by the
% G7 gate, not an assumption.
smid = getf(opts, 'static_medium', 'resummed');
medium = struct('scheme', smid, 'status', 'not_applicable', 'ref', [], 'closure', []);
if ~strcmp(smid, 'resummed')
    if isfield(opts, 'Jmom') && ~isempty(opts.Jmom)
        mom_all = opts.Jmom;                   % hot-path: threaded once per resolved point
    else
        mom_all = invz_coupling_moments(Jf);   % per-column for [nJ,nw]
    end
    mom = local_static_mom(mom_all);            % omega_n=0 is column/index 1, never flatten
    [Gref, refi] = invz_static_medium_reference(G0(1), Sigma(1), smid, ...
        struct('ref_margin', getf(opts, 'ref_margin', 1e-6)));
    [K0s, clo] = invz_medium_moment_closure(Gref, mom, smid);
    medium.ref = refi;  medium.closure = clo;
    if strcmp(refi.status, 'ok') && strcmp(clo.status, 'ok')
        medium.status = 'ok';
        K(1) = K0s;
        G(1) = G0(1)/(D(1) + K0s*G0(1));       % Dyson at the strict medium, same form as :52
    else
        if strcmp(refi.status, 'ok'), medium.status = clo.status;
        else,                         medium.status = refi.status;  end
        K(1) = NaN;  G(1) = NaN;               % domain event: reported, never regularised
    end
end

med.G = G;  med.K = K;
med.medium = medium;  med.medium_status = medium.status;
% Ordered callers replace slot 1 with the elastic-hybrid static sector. They must be able to
% validate the dynamic slots without letting the discarded PM slot vote on node acceptance.
med.dynamic_converged = all(isfinite(G(2:end))) && all(isfinite(K(2:end)));
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

% ---------------------------------------------------------------------------------------------
function mom0 = local_static_mom(mom_all)
%LOCAL_STATIC_MOM Static (omega_n = 0) scalar-field moment struct: index/column 1 of each field
% of an invz_coupling_moments struct, copied explicitly (never structfun, so the field set and
% shape are visible at the call site, not left implicit in iteration order). A missing or empty
% field means the caller did not pass an invz_coupling_moments-shaped struct -- a wiring error,
% not a domain event -- so this throws rather than returning a status.
req = {'Jbar', 'mu2', 'mu3', 'mu4', 'n'};
for k = 1:numel(req)
    f = req{k};
    if ~isfield(mom_all, f) || isempty(mom_all.(f))
        error('invz:staticMedium', ['local_static_mom: mom.%s is missing or empty; ' ...
            'expected an invz_coupling_moments struct.'], f);
    end
end
mom0 = struct('Jbar', mom_all.Jbar(1), 'mu2', mom_all.mu2(1), 'mu3', mom_all.mu3(1), ...
              'mu4', mom_all.mu4(1), 'n', mom_all.n(1));
end
