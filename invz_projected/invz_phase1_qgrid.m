function g = invz_phase1_qgrid(ion, N, offsetFlags, convention, gammaPolicy)
%INVZ_PHASE1_QGRID Cubic BZ grid with explicit endpoint, offset and Gamma policy.
% offsetFlags selects an unshifted or half-step-shifted axis independently.
% Returned q rows are wrapped into [-1/2,1/2); weights always sum to one.

if ~(isnumeric(N) && isscalar(N) && isfinite(N) && N == round(N) && N >= 2)
    error('invz:phase1Config', 'N must be a scalar integer >= 2; got %s.', mat2str(N));
end
offsetFlags = logical(offsetFlags(:).');
if ~isequal(size(offsetFlags), [1 3])
    error('invz:phase1Config', ...
        'offsetFlags must be a 1x3 logical/0-1 vector; got size %s.', ...
        mat2str(size(offsetFlags)));
end
if ~(ischar(convention) || isstring(convention)) || ...
        ~ismember(convention, {'halfopen','legacy_inclusive'})
    error('invz:phase1Config', ...
        'convention must be ''halfopen'' or ''legacy_inclusive''.');
end
if ~(ischar(gammaPolicy) || isstring(gammaPolicy)) || ...
        ~ismember(gammaPolicy, {'P_complete','P_drop'})
    error('invz:phase1Config', ...
        'gammaPolicy must be ''P_complete'' or ''P_drop''.');
end

endpointFlag = strcmp(convention, 'legacy_inclusive');
qax = qVec_generator(ion.a, 'mode', 'grid', 'grid', [N 1 1], ...
    'range', [-0.5 0.5], 'endpoint', endpointFlag, 'verbose', false);
ax0 = qax(:,1).';
ax1 = ax0 + (ax0(2) - ax0(1))/2;
axesByOffset = {ax0, ax1};
axh = axesByOffset{double(offsetFlags(1)) + 1};
axk = axesByOffset{double(offsetFlags(2)) + 1};
axl = axesByOffset{double(offsetFlags(3)) + 1};

[QH, QK, QL] = meshgrid(axh, axk, axl);
qvec = mod([QH(:), QK(:), QL(:)] + 0.5, 1) - 0.5;
nominal = size(qvec, 1);
isGamma = false(nominal, 1);
for k = 1:nominal
    isGamma(k) = invz_is_gamma_equiv(qvec(k,:), ion.tau);
end
nGamma = nnz(isGamma);

if strcmp(gammaPolicy, 'P_drop')
    qvec = qvec(~isGamma, :);
    if isempty(qvec)
        error('invz:phase1Config', ...
            'P_drop removed every point from the N=%d grid.', N);
    end
end
w = ones(size(qvec,1), 1)/size(qvec,1);

g = struct('qvec', qvec, 'w', w, 'n_gamma', nGamma, ...
    'nominal', nominal, 'N', N, 'offsetFlags', offsetFlags, ...
    'convention', char(convention), 'gammaPolicy', char(gammaPolicy));
g.note = sprintf('%s N=%d offset=[%d %d %d], %s: %d/%d rows retained.', ...
    g.convention, N, offsetFlags, g.gammaPolicy, size(qvec,1), nominal);
end
