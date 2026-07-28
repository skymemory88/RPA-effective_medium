function result = invzp_select_smooth_r_oracle(paths, spec)
%INVZP_SELECT_SMOOTH_R_ORACLE Rank complete paths by frozen smooth-r rules.
%
%   RESULT = INVZP_SELECT_SMOOTH_R_ORACLE(PATHS,SPEC) is a pure diagnostic
%   selector. PATHS must already be enumerated, branch-tracked, complete,
%   and independently audited. This function does not solve, continue,
%   differentiate, interpolate, smooth, or repair a path.

invalidId = 'invzp:SmoothSelector:InvalidInput';
validateSpec(spec,invalidId);

emptyIds = strings(0,1);
emptyMetrics = table(emptyIds,zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    false(0,1),false(0,1),false(0,1),false(0,1), ...
    'VariableNames',{'id','hstar','derivative_energy', ...
    'normalized_max_slope_change','normalized_r_curvature', ...
    'raw_max_slope_change','raw_integrated_r_curvature', ...
    'max_state_curvature','qcp_distance','trace_agreement_available', ...
    'trace_agreement_ok','shape_metric_defined','gate_ok'});
emptyRejected = table(emptyIds,emptyIds, ...
    'VariableNames',{'id','reasons'});
result = struct( ...
    'schema','invzp_smooth_r_selector_result/v2', ...
    'status','branch_unresolved', ...
    'status_detail','no candidate paths were supplied', ...
    'selected_id',emptyIds, ...
    'selected_index',[], ...
    'selected_by','', ...
    'candidate_ids',emptyIds, ...
    'survivors',emptySurvivors(), ...
    'metrics',emptyMetrics, ...
    'rejected',emptyRejected, ...
    'spec_version',string(spec.version));

if ~isstruct(paths)
    fail(invalidId,'paths must be a struct array.');
end
if isempty(paths)
    return
end

required = { ...
    'id','hstar','x','r','drdx','d2rdx2', ...
    'r_jump_free','complete','admissible','endpoint_ok','qcp_ok', ...
    'max_state_curvature','qcp_distance', ...
    'trace_agreement_available','trace_agreement_ok'};
for k = 1:numel(required)
    if ~isfield(paths,required{k})
        fail(invalidId,'Every path must contain field "%s".',required{k});
    end
end

npath = numel(paths);
ids = strings(npath,1);
hstar = zeros(npath,1);
derivativeEnergy = zeros(npath,1);
normalizedSlopeChange = nan(npath,1);
normalizedCurvature = nan(npath,1);
rawSlopeChange = zeros(npath,1);
rawCurvature = zeros(npath,1);
maxStateCurvature = nan(npath,1);
qcpDistance = nan(npath,1);
traceAgreement = false(npath,1);
traceAvailable = false(npath,1);
shapeMetricDefined = false(npath,1);
gateOk = false(npath,1);
reasons = strings(npath,1);
nx = numel(spec.x);

for k = 1:npath
    path = paths(k);
    ids(k) = normalizeId(path.id,invalidId);
    hstar(k) = finitePositiveScalar(path.hstar,'path.hstar',invalidId);
    validatePathVector(path.x,'path.x',nx,invalidId);
    if ~isequal(path.x(:),spec.x(:))
        fail(invalidId, ...
            'Every path.x must exactly equal spec.x; resampling is prohibited.');
    end
    validatePathVector(path.r,'path.r',nx,invalidId);
    validatePathVector(path.drdx,'path.drdx',nx,invalidId);
    validatePathVector(path.d2rdx2,'path.d2rdx2',nx,invalidId);

    jumpFree = logicalScalar(path.r_jump_free,'path.r_jump_free',invalidId);
    complete = logicalScalar(path.complete,'path.complete',invalidId);
    admissible = logicalScalar(path.admissible,'path.admissible',invalidId);
    endpointOk = logicalScalar(path.endpoint_ok,'path.endpoint_ok',invalidId);
    qcpOk = logicalScalar(path.qcp_ok,'path.qcp_ok',invalidId);
    traceAvailable(k) = logicalScalar( ...
        path.trace_agreement_available,'path.trace_agreement_available',invalidId);
    traceAgreement(k) = logicalScalar( ...
        path.trace_agreement_ok,'path.trace_agreement_ok',invalidId);
    if ~traceAvailable(k) && traceAgreement(k)
        fail(invalidId, ...
            'path.trace_agreement_ok cannot be true when its audit is unavailable.');
    end
    maxStateCurvature(k) = optionalNonnegativeScalar( ...
        path.max_state_curvature,'path.max_state_curvature',invalidId);
    qcpDistance(k) = optionalNonnegativeScalar( ...
        path.qcp_distance,'path.qcp_distance',invalidId);

    rawSlopeChange(k) = max(abs(diff(path.drdx(:))));
    derivativeEnergy(k) = trapz(spec.x(:),abs(path.drdx(:)).^2);
    rawCurvature(k) = trapz(spec.x(:),abs(path.d2rdx2(:)).^2);
    if ~isfinite(rawSlopeChange(k)) || ~isfinite(derivativeEnergy(k)) || ...
            ~isfinite(rawCurvature(k))
        fail(invalidId,'Path "%s" overflows a smoothness metric.',ids(k));
    end
    if derivativeEnergy(k) <= spec.degenerate.derivative_energy_max
        if rawCurvature(k) <= spec.degenerate.curvature_energy_max
            normalizedSlopeChange(k) = 0;
            normalizedCurvature(k) = 0;
            shapeMetricDefined(k) = true;
        end
    else
        normalizedSlopeChange(k) = rawSlopeChange(k)/sqrt(derivativeEnergy(k));
        normalizedCurvature(k) = rawCurvature(k)/derivativeEnergy(k);
        shapeMetricDefined(k) = isfinite(normalizedSlopeChange(k)) && ...
            isfinite(normalizedCurvature(k));
    end
    gateOk(k) = jumpFree && complete && admissible && endpointOk && qcpOk && ...
        traceAvailable(k) && traceAgreement(k) && shapeMetricDefined(k);

    failed = strings(0,1);
    if ~jumpFree, failed(end+1) = "r_jump_free=false"; end %#ok<AGROW>
    if ~complete, failed(end+1) = "complete=false"; end %#ok<AGROW>
    if ~admissible, failed(end+1) = "admissible=false"; end %#ok<AGROW>
    if ~endpointOk, failed(end+1) = "endpoint_ok=false"; end %#ok<AGROW>
    if ~qcpOk, failed(end+1) = "qcp_ok=false"; end %#ok<AGROW>
    if ~traceAvailable(k), failed(end+1) = "trace_agreement_available=false"; end %#ok<AGROW>
    if traceAvailable(k) && ~traceAgreement(k)
        failed(end+1) = "trace_agreement_ok=false"; %#ok<AGROW>
    end
    if ~shapeMetricDefined(k)
        failed(end+1) = "shape_metric_unresolved"; %#ok<AGROW>
    end
    reasons(k) = strjoin(failed,'; ');
end

if numel(unique(ids)) ~= npath
    fail(invalidId,'path.id values must be unique.');
end

result.metrics = table(ids,hstar,derivativeEnergy,normalizedSlopeChange, ...
    normalizedCurvature,rawSlopeChange,rawCurvature,maxStateCurvature, ...
    qcpDistance,traceAvailable,traceAgreement,shapeMetricDefined,gateOk, ...
    'VariableNames',{'id','hstar','derivative_energy', ...
    'normalized_max_slope_change','normalized_r_curvature', ...
    'raw_max_slope_change','raw_integrated_r_curvature', ...
    'max_state_curvature','qcp_distance','trace_agreement_available', ...
    'trace_agreement_ok','shape_metric_defined','gate_ok'});
rejected = ~gateOk;
result.rejected = table(ids(rejected),reasons(rejected), ...
    'VariableNames',{'id','reasons'});

survivors = find(gateOk);
result.candidate_ids = ids(survivors);
result.survivors.after_admissibility = ids(survivors);
if isempty(survivors)
    result.status_detail = ...
        ['no path passed all continuity, completeness, endpoint, QCP, ', ...
        'reconstruction, and shape-definition gates'];
    return
end
if isscalar(survivors)
    result = select(result,survivors,'sole_admissible',ids);
    result.survivors = fillSurvivors(result.survivors,ids(survivors));
    return
end

survivors = retainMinimum(survivors,normalizedSlopeChange, ...
    spec.tie.slope_abs,spec.tie.slope_rel);
result.survivors.after_normalized_slope = ids(survivors);
if isscalar(survivors)
    result = select(result,survivors,'normalized_max_slope_change',ids);
    result.survivors = fillSurvivors(result.survivors,ids(survivors));
    return
end

survivors = retainMinimum(survivors,normalizedCurvature, ...
    spec.tie.curvature_abs,spec.tie.curvature_rel);
result.survivors.after_normalized_curvature = ids(survivors);
if isscalar(survivors)
    result = select(result,survivors,'normalized_r_curvature',ids);
    result.survivors = fillSurvivors(result.survivors,ids(survivors));
    return
end

if any(~isfinite(maxStateCurvature(survivors)))
    result.status = 'branch_ambiguous';
    result.status_detail = ...
        'state-curvature tie breaker is unavailable for a tied path';
    result.survivors = fillSurvivors(result.survivors,ids(survivors));
    return
end
survivors = retainMinimum(survivors,maxStateCurvature, ...
    spec.tie.state_abs,spec.tie.state_rel);
result.survivors.after_state_curvature = ids(survivors);
if isscalar(survivors)
    result = select(result,survivors,'max_state_curvature',ids);
    result.survivors = fillSurvivors(result.survivors,ids(survivors));
    return
end

if any(~isfinite(qcpDistance(survivors)))
    result.status = 'branch_ambiguous';
    result.status_detail = ...
        'QCP-distance tie breaker is unavailable for a tied path';
    result.survivors = fillSurvivors(result.survivors,ids(survivors));
    return
end
survivors = retainMinimum(survivors,qcpDistance, ...
    spec.tie.qcp_abs,spec.tie.qcp_rel);
result.survivors.after_qcp_distance = ids(survivors);
if isscalar(survivors)
    result = select(result,survivors,'qcp_distance',ids);
    result.survivors = fillSurvivors(result.survivors,ids(survivors));
    return
end

result.status = 'branch_ambiguous';
result.status_detail = ...
    'multiple paths remain indistinguishable under every registered criterion';
result.survivors = fillSurvivors(result.survivors,ids(survivors));
end

function validateSpec(spec,invalidId)
if ~isstruct(spec) || ~isscalar(spec)
    fail(invalidId,'spec must be a scalar struct.');
end
if ~isfield(spec,'schema') || ~isfield(spec,'version') || ...
        ~isfield(spec,'x') || ~isfield(spec,'tie') || ...
        ~isfield(spec,'degenerate')
    fail(invalidId,'spec requires schema, version, x, tie, and degenerate fields.');
end
schema = normalizeId(spec.schema,invalidId);
if schema ~= "invzp_smooth_r_selector/v2"
    fail(invalidId,'spec.schema must be "invzp_smooth_r_selector/v2".');
end
normalizeId(spec.version,invalidId);
x = spec.x;
if ~isnumeric(x) || ~isreal(x) || ~isvector(x) || numel(x) < 2 || ...
        any(~isfinite(x),'all')
    fail(invalidId,'spec.x must be a finite real vector with at least two nodes.');
end
x = x(:);
if x(1) ~= 0 || x(end) ~= 1 || any(diff(x) <= 0)
    fail(invalidId,'spec.x must be strictly increasing from exactly 0 to 1.');
end
if ~isstruct(spec.tie) || ~isscalar(spec.tie)
    fail(invalidId,'spec.tie must be a scalar struct.');
end
names = {'slope_abs','slope_rel','curvature_abs','curvature_rel', ...
    'state_abs','state_rel','qcp_abs','qcp_rel'};
for k = 1:numel(names)
    name = names{k};
    if ~isfield(spec.tie,name)
        fail(invalidId,'spec.tie.%s is required.',name);
    end
    value = spec.tie.(name);
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
            ~isfinite(value) || value < 0
        fail(invalidId,'spec.tie.%s must be a finite nonnegative scalar.',name);
    end
end
if ~isstruct(spec.degenerate) || ~isscalar(spec.degenerate)
    fail(invalidId,'spec.degenerate must be a scalar struct.');
end
names = {'derivative_energy_max','curvature_energy_max'};
for k = 1:numel(names)
    name = names{k};
    if ~isfield(spec.degenerate,name)
        fail(invalidId,'spec.degenerate.%s is required.',name);
    end
    value = spec.degenerate.(name);
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
            ~isfinite(value) || value < 0
        fail(invalidId, ...
            'spec.degenerate.%s must be a finite nonnegative scalar.',name);
    end
end
end

function id = normalizeId(value,invalidId)
if ischar(value)
    ok = isrow(value);
elseif isstring(value)
    ok = isscalar(value) && ~ismissing(value);
else
    ok = false;
end
if ~ok
    fail(invalidId, ...
        'Text identifiers and spec schema/version must be character rows or string scalars.');
end
id = string(value);
if strlength(id) == 0
    fail(invalidId,'Text identifiers and spec schema/version must be nonempty.');
end
end

function value = finitePositiveScalar(value,name,invalidId)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value <= 0
    fail(invalidId,'%s must be a finite positive scalar.',name);
end
end

function validatePathVector(value,name,nx,invalidId)
if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || ...
        numel(value) ~= nx || any(~isfinite(value),'all')
    fail(invalidId,'%s must be a finite real vector matching spec.x.',name);
end
end

function value = logicalScalar(value,name,invalidId)
if ~islogical(value) || ~isscalar(value)
    fail(invalidId,'%s must be a scalar logical certificate.',name);
end
end

function value = optionalNonnegativeScalar(value,name,invalidId)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        (~isnan(value) && (~isfinite(value) || value < 0))
    fail(invalidId,'%s must be NaN or a finite nonnegative scalar.',name);
end
end

function survivors = retainMinimum(survivors,values,absTol,relTol)
best = min(values(survivors));
delta = abs(values(survivors)-best);
scale = max(abs(values(survivors)),abs(best));
tied = delta <= absTol + relTol.*scale;
survivors = survivors(tied);
end

function result = select(result,index,criterion,ids)
result.status = 'selected';
result.status_detail = ['unique path selected by ',criterion];
result.selected_id = ids(index);
result.selected_index = index;
result.selected_by = criterion;
end

function survivors = emptySurvivors()
empty = strings(0,1);
survivors = struct( ...
    'after_admissibility',empty, ...
    'after_normalized_slope',empty, ...
    'after_normalized_curvature',empty, ...
    'after_state_curvature',empty, ...
    'after_qcp_distance',empty);
end

function survivors = fillSurvivors(survivors,ids)
names = fieldnames(survivors);
for k = 1:numel(names)
    if isempty(survivors.(names{k}))
        survivors.(names{k}) = ids;
    end
end
end

function fail(identifier,message,varargin)
error(identifier,message,varargin{:});
end
