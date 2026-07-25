function mom = invz_coupling_moments(Jnu_flat)
%INVZ_COUPLING_MOMENTS Population central moments of the coupling multiset (strict-order
% static medium, spec SS4.1). G = -chi (meV^-1), ferromagnetic positive J.
%   mom.Jbar = mean_q J,  mom.mu_n = mean((J - Jbar).^n) for n = 2,3,4,  mom.n = count
% NORMALIZATION IS POPULATION (divide by N), matching the BZ average mean_q that
% invz_emt_scalar / invz_emt_static_ordered actually take. MATLAB's default var() is
% sample-normalized (N-1) and is NOT interchangeable: the difference is 6e-5 relative at
% N = 16384 but 4% at the N = 24 synthetic test fixtures, i.e. largest exactly where it would
% go unnoticed.
%
% Jnu_flat: [nJ,1] column vector -> scalar fields.
%           [nJ,nw] with nw > 1 (including nJ = 1, i.e. a 1 x nw row) -> 1 x nw row-vector
%           fields, ONE moment set PER COLUMN (T2.1 retardation interface). Static callers
%           use index 1 only; flattening columns into one static multiset is never a valid
%           interpretation (spec SS4.3).
% The exact q/branch weighting of the input is preserved: every entry counts once.
if ~isnumeric(Jnu_flat) || isempty(Jnu_flat) || ~isreal(Jnu_flat) || ~all(isfinite(Jnu_flat(:))) ...
        || ndims(Jnu_flat) > 2
    error('invz:couplingMoments', ...
        'Jnu_flat must be a nonempty real finite 2-D numeric array; got %s (%d elements).', ...
        class(Jnu_flat), numel(Jnu_flat));
end
J = Jnu_flat;
Jbar = mean(J, 1);
d    = J - Jbar;                       % implicit expansion over columns
mom  = struct('Jbar', Jbar, ...
              'mu2',  mean(d.^2, 1), ...
              'mu3',  mean(d.^3, 1), ...
              'mu4',  mean(d.^4, 1), ...
              'n',    repmat(size(J, 1), 1, size(J, 2)));
end
