function out = invzf_multilevel_cumulant(E,O,beta,labels,opts)
%INVZF_MULTILEVEL_CUMULANT Exact finite-level connected Matsubara cumulant.
%
%   OUT=INVZF_MULTILEVEL_CUMULANT(E,O,BETA,LABELS,OPTS) evaluates a
%   connected rank-2, rank-3, or rank-4 local cumulant.  E are energy
%   eigenvalues and O is the Hermitian operator in the same eigenbasis.
%   Integer bosonic LABELS must sum to zero.  The convention is
%
%     X(tau)=beta^(-1)*sum_n X_n*exp(-i*omega_n*tau).
%
%   Unlike invzf_twolevel_cumulant, this routine does not enumerate state
%   sequences.  Each imaginary-time ordering is one block-triangular
%   matrix exponential of dimension rank*capped_local_rank.  This is exact
%   for the retained finite Hilbert space and includes Hermite limits
%   automatically when energies or cumulative frequencies coincide.  A
%   block-diagonal similarity scaling keeps the insertion blocks O(1) in
%   the exponential; its known corner factor is removed analytically.  If
%   the retained energy span is too wide for a stable dense scaling-and-
%   squaring exponential, a dense-generator exponential action propagates
%   all required right-hand sides without forming the full exponential.
%
%   opts.local_rank (default numel(E)) keeps the lowest-energy states.
%   opts.dense_beta_span_limit (default 2000) selects the dense-exponential
%   versus exponential-action backend.  opts.degeneracy_tolerance sets the
%   rank-cut multiplet guard.  opts.max_corner_amplification may declare a
%   caller-side similarity-scaling budget; if omitted, that numerical gate
%   remains explicitly ungraded.
%
%   OUT.discarded_boltzmann_weight reports only omitted thermal weight; it
%   is not a bound on virtual intermediate-state truncation.  A rank ladder
%   is therefore mandatory before a truncated result is accepted.
%
%   This is an isolated local oracle.  It does not build a lattice
%   functional or dispatch into production.

if nargin < 5 || isempty(opts), opts = struct(); end
if ~(isstruct(opts) && isscalar(opts))
    error('invzf:multilevelCumulantOpts', ...
        'opts must be a scalar struct.');
end
validateattributes(E,{'numeric'}, ...
    {'real','vector','finite','nonempty'},mfilename,'E');
validateattributes(O,{'numeric'}, ...
    {'2d','finite','nonempty'},mfilename,'O');
validateattributes(beta,{'numeric'}, ...
    {'real','scalar','finite','positive'},mfilename,'beta');
validateattributes(labels,{'numeric'}, ...
    {'real','vector','finite','integer'},mfilename,'labels');

E = E(:);
labels = labels(:).';
nfull = numel(E);
if ~isequal(size(O),[nfull,nfull])
    error('invzf:multilevelCumulantShape', ...
        'O must be square with size numel(E).');
end
rankCumulant = numel(labels);
if rankCumulant < 2 || rankCumulant > 4
    error('invzf:multilevelCumulantRank', ...
        'labels must contain two, three, or four Matsubara indices.');
end
if sum(labels) ~= 0
    error('invzf:multilevelCumulantFrequency', ...
        'The Matsubara labels must sum exactly to zero.');
end

hermiticityResidual = norm(O-O','fro');
hermiticityScale = max(1,norm(O,'fro'));
hermiticityTolerance = 4096*eps(hermiticityScale);
if hermiticityResidual > hermiticityTolerance
    error('invzf:multilevelCumulantHermiticity', ...
        'O is not Hermitian within the roundoff-scaled tolerance.');
end
O = (O+O')/2;

[E,order] = sort(E);
O = O(order,order);
localRank = get_opt(opts,'local_rank',nfull);
validateattributes(localRank,{'numeric'}, ...
    {'real','scalar','finite','integer','positive','<=',nfull}, ...
    mfilename,'opts.local_rank');
denseSpanLimit = get_opt(opts,'dense_beta_span_limit',2000);
validateattributes(denseSpanLimit,{'numeric'}, ...
    {'real','scalar','finite','positive'}, ...
    mfilename,'opts.dense_beta_span_limit');
degeneracyScale = max(1,max(abs(E)));
degeneracyTolerance = get_opt(opts,'degeneracy_tolerance', ...
    4096*eps(degeneracyScale));
validateattributes(degeneracyTolerance,{'numeric'}, ...
    {'real','scalar','finite','nonnegative'}, ...
    mfilename,'opts.degeneracy_tolerance');
maxCornerAmplification = get_opt(opts,'max_corner_amplification',[]);
if ~isempty(maxCornerAmplification)
    validateattributes(maxCornerAmplification,{'numeric'}, ...
        {'real','scalar','finite','positive'}, ...
        mfilename,'opts.max_corner_amplification');
end

if localRank < nfull
    rankCutGap = E(localRank+1)-E(localRank);
else
    rankCutGap = Inf;
end
splitMultiplet = localRank < nfull && rankCutGap <= degeneracyTolerance;
if splitMultiplet
    error('invzf:multilevelCumulantSplitMultiplet', ...
        ['local_rank=%d cuts a degenerate/near-degenerate multiplet: ' ...
         'gap %.3e <= tolerance %.3e.'], ...
        localRank,rankCutGap,degeneracyTolerance);
end

EshiftFull = E-min(E);
logWeightsFull = -beta*EshiftFull;
logZFull = logsumexp(logWeightsFull);
probabilityFull = exp(logWeightsFull-logZFull);
discardedBoltzmannWeight = sum(probabilityFull(localRank+1:end));

E = E(1:localRank);
O = O(1:localRank,1:localRank);
E = E-min(E);
logWeights = -beta*E;
logZ = logsumexp(logWeights);
probability = exp(logWeights-logZ);
meanO = real(sum(probability.*diag(O)));
A = O-meanO*eye(localRank);
energySpan = E(end)-E(1);

[full,numericalRecords] = ...
    ordered_full_block(E,A,logZ,beta,labels,denseSpanLimit);
connected = full;
disconnected = complex(0);
if rankCumulant == 4
    pairing = [1 2 3 4;1 3 2 4;1 4 2 3];
    for k = 1:3
        left = pairing(k,1:2);
        right = pairing(k,3:4);
        if sum(labels(left)) == 0 && sum(labels(right)) == 0
            [Cleft,leftRecord] = ordered_full_block( ...
                E,A,logZ,beta,labels(left),denseSpanLimit);
            [Cright,rightRecord] = ordered_full_block( ...
                E,A,logZ,beta,labels(right),denseSpanLimit);
            numericalRecords(end+1) = leftRecord; %#ok<AGROW>
            numericalRecords(end+1) = rightRecord; %#ok<AGROW>
            disconnected = disconnected+beta*Cleft*Cright;
        end
    end
    connected = full-disconnected;
end

valueScale = max(1,abs(connected));
if abs(imag(connected)) <= 2048*eps(valueScale)
    connected = real(connected);
end
fullScale = max(1,abs(full));
if abs(imag(full)) <= 2048*eps(fullScale)
    full = real(full);
end
if abs(imag(disconnected)) <= 2048*eps(max(1,abs(disconnected)))
    disconnected = real(disconnected);
end

backends = unique({numericalRecords.backend},'stable');
if isscalar(backends)
    evaluationBackend = backends{1};
else
    evaluationBackend = 'mixed_dense_and_exponential_action';
end
cornerAmplifications = [numericalRecords.corner_amplification];
maxObservedCornerAmplification = max(cornerAmplifications);
scalingBudgetDeclared = ~isempty(maxCornerAmplification);
scalingBudgetPass = scalingBudgetDeclared && ...
    maxObservedCornerAmplification <= maxCornerAmplification;

out = struct( ...
    'schema','invzf_multilevel_cumulant/v1', ...
    'status','ok', ...
    'rank',rankCumulant, ...
    'labels',labels, ...
    'full',full, ...
    'disconnected',disconnected, ...
    'connected',connected, ...
    'mean',meanO, ...
    'local_rank',localRank, ...
    'full_local_rank',nfull, ...
    'rank_truncated',localRank < nfull, ...
    'rank_cut_gap',rankCutGap, ...
    'degeneracy_tolerance',degeneracyTolerance, ...
    'split_multiplet',splitMultiplet, ...
    'discarded_boltzmann_weight',discardedBoltzmannWeight, ...
    'virtual_truncation_certified',localRank == nfull, ...
    'rank_ladder_accepted',localRank == nfull, ...
    'retained_energy_span',energySpan, ...
    'beta_energy_span',beta*energySpan, ...
    'block_dimension',rankCumulant*localRank, ...
    'dense_beta_span_limit',denseSpanLimit, ...
    'evaluation_backend',evaluationBackend, ...
    'evaluation_backends',{backends}, ...
    'numerical_contributions',{numericalRecords}, ...
    'max_corner_amplification',maxObservedCornerAmplification, ...
    'scaling_budget_declared',scalingBudgetDeclared, ...
    'scaling_budget',{maxCornerAmplification}, ...
    'scaling_budget_pass',scalingBudgetPass, ...
    'functional_use_authorized',false, ...
    'retained_energies',E, ...
    'retained_probabilities',probability, ...
    'operator_centered',A, ...
    'operator_hermiticity_residual',hermiticityResidual, ...
    'operator_hermiticity_tolerance',hermiticityTolerance, ...
    'convention','X(tau)=beta^-1 sum_n X_n exp(-i*omega_n*tau)');
if any(~isfinite([real(full),imag(full),real(connected), ...
        imag(connected),real(disconnected),imag(disconnected)]))
    out.status = 'nonfinite';
end
end

function [value,record] = ordered_full_block( ...
        E,A,logZ,beta,labels,denseSpanLimit)
rankCumulant = numel(labels);
timed = rankCumulant-1;
orderings = perms(1:timed);
omega = 2*pi*labels/beta;
n = numel(E);
base = diag(-E);
insertionScale = 1/max(1,beta*norm(A,1));
cornerScale = insertionScale^(rankCumulant-1);
if ~(isfinite(cornerScale) && cornerScale > 0)
    error('invzf:multilevelCumulantScaling', ...
        'The insertion similarity corner scale underflowed or is nonfinite.');
end
useExponentialAction = ...
    beta*max(E) > denseSpanLimit && any(labels ~= 0);
if useExponentialAction
    backend = 'dense_generator_exponential_action';
    actionTolerance = 2^-53;
else
    backend = 'dense_block_exponential';
    actionTolerance = NaN;
end
record = struct( ...
    'labels',labels, ...
    'backend',backend, ...
    'block_dimension',rankCumulant*n, ...
    'beta_energy_span',beta*max(E), ...
    'dense_beta_span_limit',denseSpanLimit, ...
    'insertion_scale',insertionScale, ...
    'corner_scale',cornerScale, ...
    'corner_amplification',1/cornerScale, ...
    'expmv_default_tolerance',actionTolerance);
value = complex(0);
for ip = 1:size(orderings,1)
    cumulativeOmega = [0,cumsum(omega(orderings(ip,:)))];
    generator = complex(zeros(rankCumulant*n));
    for k = 1:rankCumulant
        index = (k-1)*n+(1:n);
        generator(index,index) = ...
            base+1i*cumulativeOmega(k)*eye(n);
        if k < rankCumulant
            next = k*n+(1:n);
            generator(index,next) = insertionScale*A;
        end
    end
    if useExponentialAction
        blockDimension = rankCumulant*n;
        rightHandSides = complex(zeros(blockDimension,n));
        rightHandSides((rankCumulant-1)*n+(1:n),:) = A;
        action = @(flag,x) block_rhs_action( ...
            flag,x,generator,blockDimension,n);
        propagated = expmv(action,rightHandSides(:),beta);
        propagated = reshape(propagated,blockDimension,n);
        tracedCorner = trace(propagated(1:n,:));
    else
        propagated = expm(beta*generator);
        corner = propagated(1:n,(rankCumulant-1)*n+(1:n));
        tracedCorner = trace(corner*A);
    end
    value = value+exp(-logZ)*tracedCorner/cornerScale;
end
end

function y = block_rhs_action(flag,x,generator,blockDimension,nrhs)
if strcmp(flag,'real')
    y = false;
    return
end
nprobe = size(x,2);
y = complex(zeros(size(x)));
for k = 1:nprobe
    slab = reshape(x(:,k),blockDimension,nrhs);
    if strcmp(flag,'notransp')
        product = generator*slab;
    elseif strcmp(flag,'transp')
        product = generator'*slab;
    else
        error('invzf:multilevelCumulantAction', ...
            'Unexpected exponential-action flag %s.',flag);
    end
    y(:,k) = product(:);
end
end

function z = logsumexp(x)
xmax = max(x);
z = xmax+log(sum(exp(x-xmax)));
end

function value = get_opt(opts,name,default)
if isfield(opts,name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = default;
end
end
