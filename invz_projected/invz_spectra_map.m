function S = invz_spectra_map(ion, T, fields, w, opts)
%INVZ_SPECTRA_MAP Uniform-mode chi''_cc over a transverse-field sweep.
% The 1/z leg uses the Jensen-consistent ordered/PM dispatch. The RPA leg is
% phased independently from its own bare-PM mass and bare ordered mean field;
% it never reuses the 1/z-selected state.
% S.phase / S.phase_rpa use 1=ordered, 2=paramagnet, 0=masked.
% S.Bc (alias S.Bc_1z) and S.Bc_rpa are separate sweep-midpoint boundaries;
% S.rpa_mass_pm is 1-Jcc0*chi0cc0 on the independently computed PM state.
%
% opts:
%   hyp       true
%   grid      [16 16 16]
%   dpRng     30
%   cache     false
%   eta       5e-3 meV
%   parallel  false
%   peak_wmin 0 meV
%   verbose   true
%   solve_opts struct passed to the point solvers
%     hmf_integral_mode='missing_area_ensemble' is an opt-in production
%     approximation.  hmf_missing_area_factors must contain shape factors
%     >= 0.5, including 1; the map evaluates every member and exports
%     uncertainty ranges.  The default solve remains strict full_profile.
%     hmf_adjacent_retry=false optionally retries a cold-pass approximation
%     mask from independently accepted ordered states on both sides;
%     hmf_adjacent_retry_max_span=0.30 limits their total field separation.
%     Both results must pass the unchanged gates and agree before admission.
%     hmf_ordered_boundary_retry=false optionally retries the last masked
%     ordered-side field from two independently accepted lower-field states.
%     It additionally requires an independently converged negative PM mass,
%     a PM field above the target, and explicit D_uni/F' margins.  Cold
%     phase==2 points and recovered states are never eligible.
%   Jnu/info optional precomputed coupling pair
% Coupling backend/grid fields are forwarded to invz_bz_couplings.
if nargin < 5, opts = struct(); end

hyp      = getf(opts, 'hyp', true);
eta      = getf(opts, 'eta', 5e-3);
parallel = getf(opts, 'parallel', false);
wmin     = getf(opts, 'peak_wmin', 0);
verbose  = getf(opts, 'verbose', true);
sxtra    = getf(opts, 'solve_opts', struct());
if strcmp(getf(sxtra,'hmf_integral_mode','full_profile'), ...
        'missing_area_ensemble')
    S = missing_area_ensemble_map(ion,T,fields,w,opts,sxtra);
    return;
end
invz_check_solve_opts(sxtra);

if ~(isnumeric(fields) && isreal(fields) && all(isfinite(fields(:))) && all(fields(:) >= 0))
    error('invz:fields', 'fields must contain finite nonnegative transverse fields.');
end
fields = fields(:).';
w = w(:);
nB = numel(fields);
nw = numel(w);

hasJnu = isfield(opts, 'Jnu');
hasInfo = isfield(opts, 'info');
if xor(hasJnu, hasInfo)
    error('invz:couplings', 'opts.Jnu and opts.info must be supplied together.');
end
if hasJnu
    Jnu = opts.Jnu(:);
    info = opts.info;
else
    [Jnu, info] = invz_bz_couplings(ion, opts);
end
Jcc0 = info.Jcc0;
Jaa0 = ion.Jxx0;
if isfield(info, 'Jaa0'), Jaa0 = info.Jaa0; end
Jshape = 0;
if isfield(info, 'Jshape_cc'), Jshape = info.Jshape_cc; end

sopts = sxtra;
sopts.hyp = hyp;
sopts.J0eff = Jcc0;
sopts.Jxx0 = Jaa0;
retryEnabled = getf(sopts,'hmf_adjacent_retry',false);
if ~(isscalar(retryEnabled) && (islogical(retryEnabled) || ...
        (isnumeric(retryEnabled) && isfinite(retryEnabled) && ...
        any(retryEnabled == [0 1]))))
    error('invz:hmfAdjacentRetry', ...
        'hmf_adjacent_retry must be a finite logical scalar.');
end
retryEnabled = logical(retryEnabled);
retryMaxSpan = getf(sopts,'hmf_adjacent_retry_max_span',0.30);
if ~(isnumeric(retryMaxSpan) && isreal(retryMaxSpan) && ...
        isscalar(retryMaxSpan) && isfinite(retryMaxSpan) && retryMaxSpan > 0)
    error('invz:hmfAdjacentRetry', ...
        'hmf_adjacent_retry_max_span must be a finite positive scalar in T.');
end
boundaryRetryEnabled = getf(sopts,'hmf_ordered_boundary_retry',false);
if ~(isscalar(boundaryRetryEnabled) && (islogical(boundaryRetryEnabled) || ...
        (isnumeric(boundaryRetryEnabled) && isfinite(boundaryRetryEnabled) && ...
        any(boundaryRetryEnabled == [0 1]))))
    error('invz:hmfOrderedBoundaryRetry', ...
        'hmf_ordered_boundary_retry must be a finite logical scalar.');
end
boundaryRetryEnabled = logical(boundaryRetryEnabled);
boundaryRetryMaxSpan = getf(sopts, ...
    'hmf_ordered_boundary_retry_max_span',0.20);
boundaryRetryMinD = getf(sopts,'hmf_ordered_boundary_retry_min_D_uni',1e-3);
boundaryRetryMinFprime = getf(sopts, ...
    'hmf_ordered_boundary_retry_min_Fprime',1e-3);
boundaryPmMix = getf(sopts,'hmf_ordered_boundary_retry_pm_mix',0.25);
boundaryPmMax = getf(sopts,'hmf_ordered_boundary_retry_pm_max_outer',1000);
if ~(isnumeric(boundaryRetryMaxSpan) && isreal(boundaryRetryMaxSpan) && ...
        isscalar(boundaryRetryMaxSpan) && isfinite(boundaryRetryMaxSpan) && ...
        boundaryRetryMaxSpan > 0)
    error('invz:hmfOrderedBoundaryRetry', ...
        'hmf_ordered_boundary_retry_max_span must be finite and positive.');
end
if ~(isnumeric(boundaryRetryMinD) && isreal(boundaryRetryMinD) && ...
        isscalar(boundaryRetryMinD) && isfinite(boundaryRetryMinD) && ...
        boundaryRetryMinD > 0 && isnumeric(boundaryRetryMinFprime) && ...
        isreal(boundaryRetryMinFprime) && isscalar(boundaryRetryMinFprime) && ...
        isfinite(boundaryRetryMinFprime) && boundaryRetryMinFprime > 0)
    error('invz:hmfOrderedBoundaryRetry', ...
        'Boundary-retry D_uni and Fprime margins must be finite and positive.');
end
if ~(isnumeric(boundaryPmMix) && isreal(boundaryPmMix) && ...
        isscalar(boundaryPmMix) && isfinite(boundaryPmMix) && ...
        boundaryPmMix > 0 && boundaryPmMix <= 1 && ...
        isnumeric(boundaryPmMax) && isreal(boundaryPmMax) && ...
        isscalar(boundaryPmMax) && isfinite(boundaryPmMax) && ...
        boundaryPmMax >= 1 && boundaryPmMax == floor(boundaryPmMax))
    error('invz:hmfOrderedBoundaryRetry', ...
        'Boundary-retry PM mix/max_outer settings are invalid.');
end

chiz = nan(nw, nB);
chirpa = nan(nw, nB);
Sigma0 = nan(1, nB);
phase = zeros(1, nB);
phaseRpa = zeros(1, nB);
rpaMass = nan(1, nB);
moment = nan(1, nB);
Dord = nan(1, nB);
hmf = nan(1,nB);
missingArea = nan(1,nB);
componentEdge = nan(1,nB);
pointState = cell(1,nB);

nWorkers = 0;
if parallel && ~isempty(ver('parallel')), nWorkers = Inf; end
parfor (k = 1:nB, nWorkers)
    [chiz(:,k), chirpa(:,k), Sigma0(k), phase(k), phaseRpa(k), ...
        rpaMass(k), moment(k), Dord(k), hmf(k), missingArea(k), ...
        componentEdge(k), pointState{k}] = ...
        one_field(ion, T, fields(k), Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts);
    if verbose
        labels = {'masked', 'ordered', 'paramagnet'};
        fprintf('  B = %5.2f T : 1/z %-11s RPA %-11s Sigma0 = %s\n', ...
            fields(k), labels{phase(k)+1}, labels{phaseRpa(k)+1}, num2str(Sigma0(k)));
    end
end

% Both recovery policies are defined against this immutable cold pass.
% A recovered state is output-only and can never become another seed.
coldPhase = phase;
coldPointState = pointState;

retryAttempted = false(1,nB);
retryUsed = false(1,nB);
retrySources = nan(nB,2);
retryHmfDelta = nan(1,nB);
retrySigmaDelta = nan(1,nB);
retryEdgeDelta = nan(1,nB);
if retryEnabled && strcmp(getf(sopts,'hmf_integral_mode','full_profile'), ...
        'missing_area_approx')
    % Freeze the cold-pass labels: recovered points never seed another retry.
    % This makes the result independent of retry traversal order.
    retryTargets = zeros(1,0);
    retrySeedIndex = zeros(0,2);
    for k = find(coldPhase == 0)
        orderedIndex = find(coldPhase == 1);
        leftIndex = orderedIndex(fields(orderedIndex) < fields(k));
        rightIndex = orderedIndex(fields(orderedIndex) > fields(k));
        il = [];
        ir = [];
        if ~isempty(leftIndex)
            [~,nearest] = max(fields(leftIndex));
            il = leftIndex(nearest);
        end
        if ~isempty(rightIndex)
            [~,nearest] = min(fields(rightIndex));
            ir = rightIndex(nearest);
        end
        if isempty(il) || isempty(ir) || ...
                fields(ir)-fields(il) > retryMaxSpan || ...
                isempty(coldPointState{il}) || isempty(coldPointState{ir})
            continue;
        end
        retryAttempted(k) = true;
        retrySources(k,:) = [fields(il) fields(ir)];
        retryTargets(end+1) = k; %#ok<AGROW>
        retrySeedIndex(end+1,:) = [il ir]; %#ok<AGROW>
    end
    nRetry = numel(retryTargets);
    nJobs = 2*nRetry;
    jobTarget = repelem(retryTargets,2);
    jobField = fields(jobTarget);
    jobSeed = cell(1,nJobs);
    for job = 1:nJobs
        source = coldPointState{retrySeedIndex(ceil(job/2),mod(job-1,2)+1)};
        jobSeed{job} = struct('Sigma',source.Sigma, ...
            'lambda',source.lambda,'K0',source.K(1));
    end
    trial = cell(1,nJobs);
    roptsBase = remove_fields(sopts,{'hmf_adjacent_retry', ...
        'hmf_adjacent_retry_max_span'});
    parfor (job = 1:nJobs,nWorkers)
        ropts = roptsBase;
        ropts.hmf_profile_state_seed = jobSeed{job};
        ropts.hmf_profile_sweep_direction = 'descending';
        try
            trial{job} = invz_solve_point_ordered( ...
                ion,T,jobField(job),Jnu,ropts);
        catch err
            if ~strncmp(err.identifier,'invz:',5), rethrow(err); end
            trial{job} = [];
        end
    end
    for retryIndex = 1:nRetry
        k = retryTargets(retryIndex);
        left = trial{2*retryIndex-1};
        right = trial{2*retryIndex};
        if isempty(left) || isempty(right)
            continue;
        end
        if ~(left.converged && right.converged && left.is_ordered && ...
                right.is_ordered && isfinite(left.hmf) && isfinite(right.hmf))
            continue;
        end
        leftEdge = left.hmf_prof.missing_area_integral.component_edge;
        rightEdge = right.hmf_prof.missing_area_integral.component_edge;
        retryHmfDelta(k) = abs(left.hmf-right.hmf);
        retrySigmaDelta(k) = max(abs(left.Sigma-right.Sigma));
        retryEdgeDelta(k) = abs(leftEdge-rightEdge);
        if retryHmfDelta(k) >= 1e-6 || retrySigmaDelta(k) >= 1e-6 || ...
                retryEdgeDelta(k) >= 1e-9
            continue;
        end
        selected = left; % independently checked against the right-seeded state
        copts = struct('Jsel',Jcc0,'eta',eta,'Jxx0',Jaa0, ...
            'Jshape',Jshape,'hyp',hyp,'si',selected.si);
        try
            out = invz_chi_realaxis(ion,T,fields(k),selected,w,copts);
        catch err
            if ~strncmp(err.identifier,'invz:',5), rethrow(err); end
            continue;
        end
        chiz(:,k) = imag(out.chi_cc_q(1,:)).';
        Sigma0(k) = selected.Sigma0;
        phase(k) = 1;
        moment(k) = selected.m0;
        Dord(k) = selected.D_uni;
        hmf(k) = selected.hmf;
        missingArea(k) = selected.hmf_prof.missing_area;
        componentEdge(k) = leftEdge;
        pointState{k} = selected;
        retryUsed(k) = true;
    end
end

boundaryAttempted = false(1,nB);
boundaryUsed = false(1,nB);
boundarySources = nan(nB,2);
boundaryPmConverged = false(1,nB);
boundaryPmMass = nan(1,nB);
boundaryHmfDelta = nan(1,nB);
boundarySigmaDelta = nan(1,nB);
boundaryEdgeDelta = nan(1,nB);
boundaryFprimeMin = nan(1,nB);
boundaryDuniMin = nan(1,nB);
boundaryStatus = repmat("not_eligible",1,nB);
if boundaryRetryEnabled && ...
        strcmp(getf(sopts,'hmf_integral_mode','full_profile'), ...
        'missing_area_approx')
    boundaryTargets = zeros(1,0);
    boundarySeedIndex = zeros(0,2);
    orderedIndex = find(coldPhase == 1);
    for k = find(coldPhase == 0)
        lower = orderedIndex(fields(orderedIndex) < fields(k));
        upper = orderedIndex(fields(orderedIndex) > fields(k));
        rightKnown = find(fields > fields(k) & coldPhase ~= 0,1,'first');
        if numel(lower) < 2 || ~isempty(upper) || isempty(rightKnown) || ...
                coldPhase(rightKnown) ~= 2
            continue;
        end
        sourceIndex = lower(end-1:end);
        if fields(k)-fields(sourceIndex(1)) > boundaryRetryMaxSpan || ...
                any(cellfun(@isempty,coldPointState(sourceIndex)))
            continue;
        end
        boundaryAttempted(k) = true;
        boundaryStatus(k) = "attempted";
        boundarySources(k,:) = fields(sourceIndex);
        boundaryTargets(end+1) = k; %#ok<AGROW>
        boundarySeedIndex(end+1,:) = sourceIndex; %#ok<AGROW>
    end

    nBoundary = numel(boundaryTargets);
    boundaryTrials = cell(1,nBoundary);
    boundaryPm = cell(1,nBoundary);
    boundaryFields = fields(boundaryTargets);
    boundarySeedStates = cell(nBoundary,2);
    for retryIndex = 1:nBoundary
        for sourceNumber = 1:2
            boundarySeedStates{retryIndex,sourceNumber} = coldPointState{ ...
                boundarySeedIndex(retryIndex,sourceNumber)};
        end
    end
    boundaryRemove = {'hmf_adjacent_retry', ...
        'hmf_adjacent_retry_max_span','hmf_ordered_boundary_retry', ...
        'hmf_ordered_boundary_retry_max_span', ...
        'hmf_ordered_boundary_retry_min_D_uni', ...
        'hmf_ordered_boundary_retry_min_Fprime', ...
        'hmf_ordered_boundary_retry_pm_mix', ...
        'hmf_ordered_boundary_retry_pm_max_outer'};
    boundaryOptsBase = remove_fields(sopts,boundaryRemove);
    parfor (retryIndex = 1:nBoundary,nWorkers)
        targetField = boundaryFields(retryIndex);
        pmopts = boundaryOptsBase;
        pmopts.mix_outer = boundaryPmMix;
        pmopts.max_outer = boundaryPmMax;
        try
            boundaryPm{retryIndex} = invz_solve_point( ...
                ion,T,targetField,Jnu,pmopts);
        catch err
            if ~strncmp(err.identifier,'invz:',5), rethrow(err); end
            boundaryPm{retryIndex} = [];
        end
        pair = cell(1,2);
        for sourceNumber = 1:2
            source = boundarySeedStates{retryIndex,sourceNumber};
            ropts = boundaryOptsBase;
            ropts.hmf_profile_state_seed = struct('Sigma',source.Sigma, ...
                'lambda',source.lambda,'K0',source.K(1));
            ropts.hmf_profile_sweep_direction = 'descending';
            try
                pair{sourceNumber} = invz_solve_point_ordered( ...
                    ion,T,targetField,Jnu,ropts);
            catch err
                if ~strncmp(err.identifier,'invz:',5), rethrow(err); end
                pair{sourceNumber} = [];
            end
        end
        boundaryTrials{retryIndex} = pair;
    end

    for retryIndex = 1:nBoundary
        k = boundaryTargets(retryIndex);
        pm = boundaryPm{retryIndex};
        if ~isempty(pm)
            boundaryPmConverged(k) = pm.converged;
            boundaryPmMass(k) = pm.crit;
        end
        if isempty(pm) || ~pm.converged || ~isfinite(pm.crit) || pm.crit >= 0
            boundaryStatus(k) = "pm_gate_failed";
            continue;
        end
        pair = boundaryTrials{retryIndex};
        left = pair{1};
        right = pair{2};
        [leftGate,leftFprime] = boundary_point_gate( ...
            left,Jnu,Jcc0,boundaryRetryMinD,boundaryRetryMinFprime);
        [rightGate,rightFprime] = boundary_point_gate( ...
            right,Jnu,Jcc0,boundaryRetryMinD,boundaryRetryMinFprime);
        if ~(leftGate && rightGate)
            boundaryStatus(k) = "ordered_gate_failed";
            continue;
        end
        leftEdge = left.hmf_prof.missing_area_integral.component_edge;
        rightEdge = right.hmf_prof.missing_area_integral.component_edge;
        boundaryHmfDelta(k) = abs(left.hmf-right.hmf);
        boundarySigmaDelta(k) = max(abs(left.Sigma-right.Sigma));
        boundaryEdgeDelta(k) = abs(leftEdge-rightEdge);
        boundaryFprimeMin(k) = min(leftFprime,rightFprime);
        boundaryDuniMin(k) = min(left.D_uni,right.D_uni);
        sameComponent = ...
            left.hmf_prof.missing_area_integral.node_count == ...
            right.hmf_prof.missing_area_integral.node_count;
        if boundaryHmfDelta(k) >= 1e-6 || ...
                boundarySigmaDelta(k) >= 1e-6 || ...
                boundaryEdgeDelta(k) >= 1e-9 || ~sameComponent
            boundaryStatus(k) = "agreement_failed";
            continue;
        end
        selected = left;
        copts = struct('Jsel',Jcc0,'eta',eta,'Jxx0',Jaa0, ...
            'Jshape',Jshape,'hyp',hyp,'si',selected.si);
        try
            out = invz_chi_realaxis(ion,T,fields(k),selected,w,copts);
        catch err
            if ~strncmp(err.identifier,'invz:',5), rethrow(err); end
            boundaryStatus(k) = "realaxis_failed";
            continue;
        end
        candidateSpectrum = imag(out.chi_cc_q(1,:)).';
        if any(~isfinite(candidateSpectrum))
            boundaryStatus(k) = "realaxis_nonfinite";
            continue;
        end
        chiz(:,k) = candidateSpectrum;
        Sigma0(k) = selected.Sigma0;
        phase(k) = 1;
        moment(k) = selected.m0;
        Dord(k) = selected.D_uni;
        hmf(k) = selected.hmf;
        missingArea(k) = selected.hmf_prof.missing_area;
        componentEdge(k) = leftEdge;
        pointState{k} = selected;
        boundaryUsed(k) = true;
        boundaryStatus(k) = "used";
    end
end

S = struct('fields', fields, 'w', w, 'T', T, 'info', info, ...
    'chiz', chiz, 'chirpa', chirpa, 'Sigma0', Sigma0, ...
    'phase', phase, 'phase_rpa', phaseRpa, 'rpa_mass_pm', rpaMass, ...
    'moment', moment, 'D_ord', Dord,'hmf',hmf, ...
    'missing_area',missingArea,'component_edge',componentEdge);
S.hmf_integral_mode = getf(sxtra,'hmf_integral_mode','full_profile');
S.visual_only = ismember(S.hmf_integral_mode, ...
    {'endpoint_trapezoid_visual','filtered_profile_visual'});
S.production_approximation = strcmp(S.hmf_integral_mode, ...
    'missing_area_approx');
S.approximation_only = S.visual_only || S.production_approximation;
S.missing_area_factor = getf(sxtra,'hmf_missing_area_factor',NaN);
S.approximation_branch = getf(sxtra,'hmf_approx_branch', ...
    'not_applicable');
S.adjacent_retry = struct('enabled',retryEnabled, ...
    'max_seed_span_T',retryMaxSpan,'attempted',retryAttempted, ...
    'used',retryUsed,'source_fields',retrySources, ...
    'hmf_delta',retryHmfDelta,'Sigma_delta',retrySigmaDelta, ...
    'component_edge_delta',retryEdgeDelta, ...
    'interpretation',['A retry is accepted only when independently seeded ' ...
    'ordered states from both sides pass every original gate and agree. ' ...
    'Cold-pass labels are frozen, so retries cannot seed other retries.']);
S.ordered_boundary_retry = struct('enabled',boundaryRetryEnabled, ...
    'max_seed_span_T',boundaryRetryMaxSpan, ...
    'min_D_uni',boundaryRetryMinD,'min_Fprime',boundaryRetryMinFprime, ...
    'pm_mix_outer',boundaryPmMix,'pm_max_outer',boundaryPmMax, ...
    'attempted',boundaryAttempted,'used',boundaryUsed, ...
    'source_fields',boundarySources,'pm_converged',boundaryPmConverged, ...
    'pm_mass',boundaryPmMass,'hmf_delta',boundaryHmfDelta, ...
    'Sigma_delta',boundarySigmaDelta, ...
    'component_edge_delta',boundaryEdgeDelta, ...
    'Fprime_min',boundaryFprimeMin,'D_uni_min',boundaryDuniMin, ...
    'status',boundaryStatus, ...
    'interpretation',['A boundary retry is attempted only for a frozen ' ...
    'cold phase==0 point above every accepted ordered cold point and below ' ...
    'an accepted PM cold point. Two untouched lower ordered states must ' ...
    'independently recover the same stable root, and an independent PM ' ...
    'solve must converge with negative mass.']);
S.Bc = invz_boundary_field(fields, phase == 1, phase == 2);
S.Bc_1z = S.Bc;
S.Bc_rpa = invz_boundary_field(fields, phaseRpa == 1, phaseRpa == 2);
S.Epeak = invz_peak_energy(chiz, w, wmin);
S.Epeak_rpa = invz_peak_energy(chirpa, w, wmin);
end

function [chiz, chirpa, Sigma0, phase, phaseRpa, rpaMass, moment, Dord, ...
    hmf, missingArea, componentEdge, pointState] = ...
    one_field(ion, T, B, Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts)
nw = numel(w);
chiz = nan(nw,1);
chirpa = nan(nw,1);
Sigma0 = NaN;
rpaMass = NaN;
moment = NaN;
Dord = NaN;
hmf = NaN;
missingArea = NaN;
componentEdge = NaN;
pointState = [];

[pt, phase] = invz_solve_auto(ion, T, B, Jnu, sopts);
if phase ~= 0 && ~isempty(pt)
    copts = struct('Jsel', Jcc0, 'eta', eta, 'Jxx0', Jaa0, ...
        'Jshape', Jshape, 'hyp', hyp, 'si', pt.si);
    try
        out = invz_chi_realaxis(ion, T, B, pt, w, copts);
        chiz = imag(out.chi_cc_q(1,:)).';
        pointState = pt;
        Sigma0 = pt.Sigma0;
        if phase == 1
            moment = pt.m0;
            Dord = pt.D_uni;
            hmf = pt.hmf;
            if isfield(pt,'hmf_prof') && ...
                    isfield(pt.hmf_prof,'missing_area')
                missingArea = pt.hmf_prof.missing_area;
            end
            if isfield(pt,'hmf_prof') && ...
                    isfield(pt.hmf_prof,'missing_area_integral') && ...
                    isfield(pt.hmf_prof.missing_area_integral,'component_edge')
                componentEdge = ...
                    pt.hmf_prof.missing_area_integral.component_edge;
            end
        end
    catch err
        if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
        phase = 0;
        pointState = [];
        Sigma0 = NaN;
        moment = NaN;
        Dord = NaN;
    end
end

% Independent bare-RPA leg. J0eff selects/constructs the MF state; the exact same
% scalar is Jsel in the denominator. This remains available when the 1/z column is masked.
ropts = struct('Jsel', Jcc0, 'J0eff', Jcc0, 'bare_rpa', true, ...
    'eta', eta, 'Jxx0', Jaa0, 'Jshape', Jshape, 'hyp', hyp, 'npass', 1);
try
    out0 = invz_chi_realaxis(ion, T, B, [], w, ropts);
    chirpa = imag(out0.chi_cc_q(1,:)).';
    phaseRpa = out0.phase_rpa;
    rpaMass = out0.rpa_mass_pm;
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    phaseRpa = 0;
end
end

function S = missing_area_ensemble_map(ion,T,fields,w,opts,sxtra)
% Evaluate explicit area members through the scalar approximation mode.
factors = getf(sxtra,'hmf_missing_area_factors',[0.75 1 1.5]);
if ~(isnumeric(factors) && isreal(factors) && isvector(factors) && ...
        all(isfinite(factors(:))) && all(factors(:) >= 0.5) && ...
        numel(unique(factors(:))) == numel(factors))
    error('invz:hmfMissingAreaFactors', ...
        ['hmf_missing_area_factors must be distinct finite values >= 0.5 ' ...
         'for nonnegative linear completions.']);
end
factors = factors(:).';
central = find(abs(factors-1) <= 16*eps,1);
if isempty(central)
    error('invz:hmfMissingAreaFactors', ...
        'hmf_missing_area_factors must include factor 1 as the central member.');
end
branch = getf(sxtra,'hmf_approx_branch', ...
    'picard_attracting_contiguous_high_h_component');
if ~strcmp(branch,'picard_attracting_contiguous_high_h_component')
    error('invz:hmfApproxBranch', ...
        'The ensemble supports only the declared Picard-attracting branch.');
end

base_opts = opts;
if xor(isfield(base_opts,'Jnu'),isfield(base_opts,'info'))
    error('invz:couplings', ...
        'opts.Jnu and opts.info must be supplied together.');
end
if ~isfield(base_opts,'Jnu')
    [base_opts.Jnu,base_opts.info] = invz_bz_couplings(ion,opts);
end
base_verbose = getf(opts,'verbose',true);
base_opts.verbose = false;
members = cell(1,numel(factors));
for k = 1:numel(factors)
    member_opts = sxtra;
    member_opts = remove_fields(member_opts, ...
        {'hmf_missing_area_factors'});
    member_opts.hmf_integral_mode = 'missing_area_approx';
    member_opts.hmf_missing_area_factor = factors(k);
    member_opts.hmf_approx_branch = branch;
    base_opts.solve_opts = member_opts;
    members{k} = invz_spectra_map(ion,T,fields,w,base_opts);
    if base_verbose
        fprintf(['missing-area ensemble factor %.6g: ordered %d, PM %d, ' ...
            'masked %d\n'],factors(k),nnz(members{k}.phase == 1), ...
            nnz(members{k}.phase == 2),nnz(members{k}.phase == 0));
    end
end

S = members{central};
nmember = numel(members);
nfield = numel(fields);
phase = nan(nmember,nfield);
hmf = nan(nmember,nfield);
area = nan(nmember,nfield);
edge = nan(nmember,nfield);
sigma = nan(nmember,nfield);
moment = nan(nmember,nfield);
Dord = nan(nmember,nfield);
peak = nan(nmember,nfield);
chiz = nan(numel(w),nfield,nmember);
retryAttempted = false(nmember,nfield);
retryUsed = false(nmember,nfield);
retrySources = nan(nmember,nfield,2);
retryHmfDelta = nan(nmember,nfield);
retrySigmaDelta = nan(nmember,nfield);
retryEdgeDelta = nan(nmember,nfield);
boundaryAttempted = false(nmember,nfield);
boundaryUsed = false(nmember,nfield);
boundarySources = nan(nmember,nfield,2);
boundaryPmConverged = false(nmember,nfield);
boundaryPmMass = nan(nmember,nfield);
boundaryHmfDelta = nan(nmember,nfield);
boundarySigmaDelta = nan(nmember,nfield);
boundaryEdgeDelta = nan(nmember,nfield);
boundaryFprimeMin = nan(nmember,nfield);
boundaryDuniMin = nan(nmember,nfield);
boundaryStatus = strings(nmember,nfield);
for k = 1:nmember
    phase(k,:) = members{k}.phase;
    hmf(k,:) = members{k}.hmf;
    area(k,:) = members{k}.missing_area;
    edge(k,:) = members{k}.component_edge;
    sigma(k,:) = members{k}.Sigma0;
    moment(k,:) = members{k}.moment;
    Dord(k,:) = members{k}.D_ord;
    peak(k,:) = members{k}.Epeak;
    chiz(:,:,k) = members{k}.chiz;
    retryAttempted(k,:) = members{k}.adjacent_retry.attempted;
    retryUsed(k,:) = members{k}.adjacent_retry.used;
    retrySources(k,:,:) = reshape( ...
        members{k}.adjacent_retry.source_fields,[1 nfield 2]);
    retryHmfDelta(k,:) = members{k}.adjacent_retry.hmf_delta;
    retrySigmaDelta(k,:) = members{k}.adjacent_retry.Sigma_delta;
    retryEdgeDelta(k,:) = members{k}.adjacent_retry.component_edge_delta;
    boundaryAttempted(k,:) = members{k}.ordered_boundary_retry.attempted;
    boundaryUsed(k,:) = members{k}.ordered_boundary_retry.used;
    boundarySources(k,:,:) = reshape( ...
        members{k}.ordered_boundary_retry.source_fields,[1 nfield 2]);
    boundaryPmConverged(k,:) = ...
        members{k}.ordered_boundary_retry.pm_converged;
    boundaryPmMass(k,:) = members{k}.ordered_boundary_retry.pm_mass;
    boundaryHmfDelta(k,:) = members{k}.ordered_boundary_retry.hmf_delta;
    boundarySigmaDelta(k,:) = members{k}.ordered_boundary_retry.Sigma_delta;
    boundaryEdgeDelta(k,:) = ...
        members{k}.ordered_boundary_retry.component_edge_delta;
    boundaryFprimeMin(k,:) = ...
        members{k}.ordered_boundary_retry.Fprime_min;
    boundaryDuniMin(k,:) = members{k}.ordered_boundary_retry.D_uni_min;
    boundaryStatus(k,:) = members{k}.ordered_boundary_retry.status;
end

[hmf_min,hmf_max,hmf_complete] = finite_range(hmf,1);
[sigma_min,sigma_max,sigma_complete] = finite_range(sigma,1);
[moment_min,moment_max,moment_complete] = finite_range(moment,1);
[Dord_min,Dord_max,Dord_complete] = finite_range(Dord,1);
[peak_min,peak_max,peak_complete] = finite_range(peak,1);
[chiz_min,chiz_max,chiz_complete] = finite_range(chiz,3);
phase_consistent = all(phase == phase(central,:),1) & phase(central,:) ~= 0;

ensemble = struct('factors',factors, ...
    'origin_edge_ratios',2*factors-1,'central_index',central, ...
    'central_factor',factors(central),'branch_assumption',branch, ...
    'profile_node_count',getf(sxtra,'nH',33), ...
    'resolution_qualification',['The area band does not include profile-' ...
    'resolution uncertainty; validate the certified edge on a node ladder.'], ...
    'phase',phase,'phase_consistent',phase_consistent, ...
    'hmf',hmf,'hmf_interval',[hmf_min(:) hmf_max(:)], ...
    'hmf_complete',hmf_complete,'missing_area',area, ...
    'component_edge',edge,'Sigma0',sigma, ...
    'Sigma0_interval',[sigma_min(:) sigma_max(:)], ...
    'Sigma0_complete',sigma_complete,'moment',moment, ...
    'moment_interval',[moment_min(:) moment_max(:)], ...
    'moment_complete',moment_complete,'D_ord',Dord, ...
    'D_ord_interval',[Dord_min(:) Dord_max(:)], ...
    'D_ord_complete',Dord_complete,'Epeak',peak, ...
    'Epeak_interval',[peak_min(:) peak_max(:)], ...
    'Epeak_complete',peak_complete,'chiz',chiz, ...
    'chiz_min',chiz_min,'chiz_max',chiz_max, ...
    'chiz_complete',chiz_complete, ...
    'adjacent_retry',struct('attempted',retryAttempted, ...
        'used',retryUsed,'source_fields',retrySources, ...
        'hmf_delta',retryHmfDelta,'Sigma_delta',retrySigmaDelta, ...
        'component_edge_delta',retryEdgeDelta), ...
    'ordered_boundary_retry',struct('attempted',boundaryAttempted, ...
        'used',boundaryUsed,'source_fields',boundarySources, ...
        'pm_converged',boundaryPmConverged,'pm_mass',boundaryPmMass, ...
        'hmf_delta',boundaryHmfDelta,'Sigma_delta',boundarySigmaDelta, ...
        'component_edge_delta',boundaryEdgeDelta, ...
        'Fprime_min',boundaryFprimeMin,'D_uni_min',boundaryDuniMin, ...
        'status',boundaryStatus), ...
    'member_maps',{members}, ...
    'interpretation',['Ranges span explicit positive missing-area shape ' ...
    'completions on one declared numerical branch. They are approximation ' ...
    'sensitivity bands, not confidence intervals or equilibrium selection.']);
S.hmf_integral_mode = 'missing_area_ensemble';
S.visual_only = false;
S.production_approximation = true;
S.approximation_only = true;
S.missing_area_factor = factors(central);
S.approximation_branch = branch;
S.approximation_ensemble = ensemble;
end

function [pass,Fprime] = boundary_point_gate(point,Jnu,J0eff,minD,minFprime)
pass = false;
Fprime = NaN;
if isempty(point) || ~point.converged || ~point.is_ordered || ...
        ~isfinite(point.hmf) || ~isfinite(point.D_uni) || ...
        point.D_uni < minD || ~isfinite(point.final_resid) || ...
        point.final_resid >= 1e-8
    return;
end
p = point.hmf_prof;
if ~(isfield(p,'r_star') && isfinite(p.r_star) && p.r_star > 0 && ...
        isfield(p,'Gstat_star') && isfinite(p.Gstat_star) && ...
        isfield(p,'root_bracket_count') && p.root_bracket_count == 1 && ...
        isfield(p,'root_bracket_bridged') && ~p.root_bracket_bridged && ...
        isfield(p,'missing_area_integral') && ...
        p.missing_area_integral.node_count >= 2 && ...
        isfinite(p.missing_area_integral.component_edge) && ...
        isfinite(point.K(1)))
    return;
end
Gstat = p.Gstat_star;
K0 = point.K(1);
Gtilde = Gstat/(1-K0*Gstat);
Fprime = p.r_star*(1+J0eff*Gtilde);
supremumMass = 1+J0eff*Gtilde;
if isvector(Jnu), Jstatic = Jnu(:); else, Jstatic = Jnu(:,1); end
meshXMass = min(1+Jstatic*Gtilde);
meshMediumMass = min(1+(Jstatic-K0)*Gstat);
pass = isfinite(Fprime) && Fprime >= minFprime && ...
    isfinite(supremumMass) && supremumMass > 0 && ...
    isfinite(meshXMass) && meshXMass > 0 && ...
    isfinite(meshMediumMass) && meshMediumMass > 0;
end

function s = remove_fields(s,names)
for k = 1:numel(names)
    if isfield(s,names{k}), s = rmfield(s,names{k}); end
end
end

function [lo,hi,complete] = finite_range(values,dim)
complete = all(isfinite(values),dim);
lo = min(values,[],dim);
hi = max(values,[],dim);
lo(~complete) = NaN;
hi(~complete) = NaN;
end
