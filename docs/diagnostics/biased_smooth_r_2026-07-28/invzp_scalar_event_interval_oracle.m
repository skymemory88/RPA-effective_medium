function result = invzp_scalar_event_interval_oracle( ...
        Jq,enclosure,spec)
%INVZP_SCALAR_EVENT_INTERVAL_ORACLE Certify invZ events from a w enclosure.
%
%   RESULT = INVZP_SCALAR_EVENT_INTERVAL_ORACLE(JQ,ENCLOSURE,SPEC)
%   bounds the pole and reciprocal-mean event margins over a caller-certified
%   interval for
%
%       w(h) = z(h)-K0(h).
%
%   For fixed real JQ, every static denominator is w+Jq and
%
%       Gbar(w) = mean(1./(w+Jq)).
%
%   Gbar is strictly decreasing between its poles because
%   dGbar/dw = -mean(1./(w+Jq).^2). Thus the event calculation is scalar
%   once ENCLOSURE proves that the complete branch edge stays inside its
%   declared w interval.
%
%   The static closure gives a second scalar identity:
%
%       Jloc = 1/Gbar-w,  z = w+K0 = 1/Gbar+(K0-Jloc).
%
%   A separately validated normalization enclosure may turn the lower bound
%   on abs(Gbar) into the required upper bound on abs(z). A caller-supplied
%   number is not sufficient: ENCLOSURE.normalization must carry its own
%   validated flag and nonempty proof payload.
%
%   Arithmetic is outward-rounded in binary64 by stepping every elementary
%   interval result to the adjacent floating-point numbers. A passing
%   result is returned in the exact event-evidence schema consumed by
%   INVZP_FIXED_H_TRACE_GRAPH_INPUTS.
%
%   This function does NOT construct or validate either enclosure.
%   ENCLOSURE.validated and ENCLOSURE.normalization.validated must be
%   caller-owned certificates backed by nonempty proof payloads. Without
%   both, event_bracket_ok is always false. Absolute pole distance and the
%   reciprocal-average interval remain useful conditional diagnostics.

invalidId = 'invzp:ScalarEventInterval:InvalidInput';
[Jq,enclosure,cfg] = validateInputs(Jq,enclosure,spec,invalidId);

Jscale = max(abs(Jq));
[denLower,denUpper,arithmeticFinite] = denominatorIntervals( ...
    enclosure.w_interval,Jq);

polePossible = ~arithmeticFinite || any( ...
    denLower <= 0 & denUpper >= 0);
poleDistanceLower = 0;
reciprocalInterval = [-Inf Inf];
if ~polePossible
    distance = zeros(size(Jq));
    positive = denLower > 0;
    distance(positive) = denLower(positive);
    distance(~positive) = -denUpper(~positive);
    poleDistanceLower = max(0,min(distance));
    reciprocalInterval = reciprocalAverageInterval(denLower,denUpper);
end

poleMarginLower = 0;
meanMarginLower = 0;
scaleUpper = Inf;
if ~polePossible && all(isfinite(reciprocalInterval))
    meanAbsoluteLower = intervalAbsLower(reciprocalInterval);
    meanMarginLower = max(0,nextDownScalar(meanAbsoluteLower*Jscale));
    if meanAbsoluteLower > 0 && enclosure.normalization.validated
        zAbsUpper = nextUpScalar( ...
            nextUpScalar(1/meanAbsoluteLower)+ ...
            enclosure.normalization.static_closure_abs_bound);
        scaleUpper = max(zAbsUpper,Jscale);
        poleMarginLower = max(0,nextDownScalar( ...
            poleDistanceLower/scaleUpper));
    end
end

poleEventLower = nextDownScalar( ...
    poleMarginLower-cfg.pole_margin_min);
meanEventLower = nextDownScalar( ...
    meanMarginLower-cfg.mean_margin_min);
meanZeroPossible = reciprocalInterval(1) <= 0 && ...
    reciprocalInterval(2) >= 0;
boundsPass = arithmeticFinite && ~polePossible && ~meanZeroPossible && ...
    poleEventLower > 0 && meanEventLower > 0;
conditionalBoundsClear = arithmeticFinite && ~polePossible && ...
    ~meanZeroPossible && poleDistanceLower > 0 && meanEventLower > 0;
eventPass = enclosure.validated && ...
    enclosure.normalization.validated && boundsPass;

if ~enclosure.validated
    status = "enclosure_unvalidated";
    detail = "the scalar bounds are conditional on an unvalidated branch enclosure";
elseif ~enclosure.normalization.validated
    status = "normalization_unvalidated";
    detail = "absolute scalar clearances exist, but the normalized pole margin lacks a validated z bound";
elseif ~arithmeticFinite
    status = "arithmetic_inconclusive";
    detail = "nonfinite interval arithmetic prevents certification";
elseif polePossible
    status = "pole_possible";
    detail = "the w enclosure intersects at least one lattice pole";
elseif meanZeroPossible
    status = "mean_zero_possible";
    detail = "the reciprocal-average interval contains zero";
elseif ~boundsPass
    status = "margin_too_small";
    detail = "a full-edge lower bound does not exceed its registered event floor";
else
    status = "certified_clear";
    detail = "the validated scalar enclosure excludes both registered events";
end

marginNames = ["pole","mean"];
lowerBounds = [poleEventLower meanEventLower];
proof = struct( ...
    'schema','invzp_scalar_stieltjes_interval_proof/v1', ...
    'identity','Gbar=mean(1/(w+Jq)); dGbar/dw=-mean(1/(w+Jq)^2)', ...
    'arithmetic','outward_binary64_adjacent/v1', ...
    'enclosure',enclosure, ...
    'normalization',enclosure.normalization, ...
    'Jscale',Jscale,'scale_upper',scaleUpper, ...
    'denominator_lower',denLower, ...
    'denominator_upper',denUpper, ...
    'pole_distance_lower',poleDistanceLower, ...
    'pole_margin_lower',poleMarginLower, ...
    'reciprocal_average_interval',reciprocalInterval, ...
    'mean_margin_lower',meanMarginLower, ...
    'pole_event_margin_lower',poleEventLower, ...
    'mean_event_margin_lower',meanEventLower);
payload = struct( ...
    'schema',"invzp_signed_event_edge/v1", ...
    'method',"scalar_stieltjes_interval", ...
    'source_id',enclosure.source_id, ...
    'from_record_index',enclosure.from_record_index, ...
    'to_record_index',enclosure.to_record_index, ...
    'margin_names',marginNames, ...
    'lower_bounds',lowerBounds, ...
    'proof',proof);
sourceRecordId = makeEventId( ...
    enclosure.source_id,enclosure.from_record_index, ...
    enclosure.to_record_index);
eventEvidence = struct( ...
    'schema',cfg.certificate_schema, ...
    'from_record_index',enclosure.from_record_index, ...
    'to_record_index',enclosure.to_record_index, ...
    'event_bracket_ok',logical(eventPass), ...
    'source_record_id',sourceRecordId, ...
    'payload',payload);

result = struct( ...
    'schema','invzp_scalar_event_interval_oracle/v1', ...
    'status',char(status),'status_detail',char(detail), ...
    'spec_version',cfg.version, ...
    'source_id',enclosure.source_id, ...
    'from_record_index',enclosure.from_record_index, ...
    'to_record_index',enclosure.to_record_index, ...
    'Jscale',Jscale,'scale_upper',scaleUpper, ...
    'w_interval',enclosure.w_interval, ...
    'pole_possible',logical(polePossible), ...
    'mean_zero_possible',logical(meanZeroPossible), ...
    'pole_distance_lower',poleDistanceLower, ...
    'pole_margin_lower',poleMarginLower, ...
    'reciprocal_average_interval',reciprocalInterval, ...
    'mean_margin_lower',meanMarginLower, ...
    'event_margin_lower_bounds',lowerBounds, ...
    'enclosure_validated',enclosure.validated, ...
    'normalization_validated',enclosure.normalization.validated, ...
    'conditional_bounds_clear',logical(conditionalBoundsClear), ...
    'event_bracket_ok',logical(eventPass), ...
    'event_evidence',eventEvidence, ...
    'complete_graph_edge_claim',false);
end

function [Jq,enclosure,cfg] = validateInputs( ...
        Jq,enclosure,spec,invalidId)
if ~isa(Jq,'double') || ~isreal(Jq) || ~isvector(Jq) || ...
        isempty(Jq) || any(~isfinite(Jq),'all')
    fail(invalidId,'Jq must be a nonempty finite real double vector.');
end
Jq = Jq(:);
Jscale = max(abs(Jq));
if ~(isfinite(Jscale) && Jscale > 0)
    fail(invalidId,'Jq must have a positive finite absolute scale.');
end

if ~isstruct(spec) || ~isscalar(spec)
    fail(invalidId,'spec must be a scalar struct.');
end
requiredSpec = {'version','certificate_schema','pole_margin_min', ...
    'mean_margin_min','allowed_enclosure_methods'};
requireFields(spec,requiredSpec,'spec',invalidId);
cfg = struct( ...
    'version',textScalar(spec.version,'spec.version',invalidId), ...
    'certificate_schema',textScalar( ...
        spec.certificate_schema,'spec.certificate_schema',invalidId), ...
    'pole_margin_min',positiveScalar( ...
        spec.pole_margin_min,'spec.pole_margin_min',invalidId), ...
    'mean_margin_min',positiveScalar( ...
        spec.mean_margin_min,'spec.mean_margin_min',invalidId), ...
    'allowed_enclosure_methods',textVector( ...
        spec.allowed_enclosure_methods, ...
        'spec.allowed_enclosure_methods',invalidId));
if numel(unique(cfg.allowed_enclosure_methods)) ~= ...
        numel(cfg.allowed_enclosure_methods)
    fail(invalidId,'spec.allowed_enclosure_methods must be unique.');
end

if ~isstruct(enclosure) || ~isscalar(enclosure)
    fail(invalidId,'enclosure must be a scalar struct.');
end
requiredEnclosure = {'schema','source_id','from_record_index', ...
    'to_record_index','method','validated','w_interval', ...
    'proof','normalization'};
requireFields(enclosure,requiredEnclosure,'enclosure',invalidId);
if textScalar(enclosure.schema,'enclosure.schema',invalidId) ~= ...
        "invzp_scalar_w_enclosure/v1"
    fail(invalidId, ...
        'enclosure must use schema invzp_scalar_w_enclosure/v1.');
end
enclosure.source_id = textScalar( ...
    enclosure.source_id,'enclosure.source_id',invalidId);
enclosure.from_record_index = positiveInteger( ...
    enclosure.from_record_index, ...
    'enclosure.from_record_index',invalidId);
enclosure.to_record_index = positiveInteger( ...
    enclosure.to_record_index,'enclosure.to_record_index',invalidId);
if enclosure.to_record_index ~= enclosure.from_record_index+1
    fail(invalidId,'enclosure must bind consecutive record indices.');
end
enclosure.method = textScalar( ...
    enclosure.method,'enclosure.method',invalidId);
if ~ismember(enclosure.method,cfg.allowed_enclosure_methods)
    fail(invalidId, ...
        'Enclosure method "%s" is not preregistered.',enclosure.method);
end
enclosure.validated = logicalScalar( ...
    enclosure.validated,'enclosure.validated',invalidId);
if ~isa(enclosure.w_interval,'double') || ...
        ~isreal(enclosure.w_interval) || ...
        ~isequal(size(enclosure.w_interval),[1 2]) || ...
        any(~isfinite(enclosure.w_interval),'all') || ...
        enclosure.w_interval(1) > enclosure.w_interval(2)
    fail(invalidId, ...
        'enclosure.w_interval must be a finite increasing double row.');
end
if ~isstruct(enclosure.proof) || ~isscalar(enclosure.proof)
    fail(invalidId,'enclosure.proof must be a scalar struct.');
elseif enclosure.validated && isempty(fieldnames(enclosure.proof))
    fail(invalidId, ...
        'A validated enclosure requires a nonempty proof payload.');
end
if ~isstruct(enclosure.normalization) || ...
        ~isscalar(enclosure.normalization)
    fail(invalidId,'enclosure.normalization must be a scalar struct.');
end
requireFields(enclosure.normalization, ...
    {'schema','validated','static_closure_abs_bound','proof'}, ...
    'enclosure.normalization',invalidId);
if textScalar(enclosure.normalization.schema, ...
        'enclosure.normalization.schema',invalidId) ~= ...
        "invzp_scalar_normalization_enclosure/v1"
    fail(invalidId, ...
        ['enclosure.normalization must use schema ', ...
        'invzp_scalar_normalization_enclosure/v1.']);
end
enclosure.normalization.validated = logicalScalar( ...
    enclosure.normalization.validated, ...
    'enclosure.normalization.validated',invalidId);
if enclosure.normalization.validated
    enclosure.normalization.static_closure_abs_bound = ...
        nonnegativeScalar( ...
        enclosure.normalization.static_closure_abs_bound, ...
        'enclosure.normalization.static_closure_abs_bound',invalidId);
    if ~isstruct(enclosure.normalization.proof) || ...
            ~isscalar(enclosure.normalization.proof) || ...
            isempty(fieldnames(enclosure.normalization.proof))
        fail(invalidId, ...
            'A validated normalization requires a nonempty proof payload.');
    end
else
    if ~isa(enclosure.normalization.static_closure_abs_bound,'double') || ...
            ~isreal(enclosure.normalization.static_closure_abs_bound) || ...
            ~isscalar(enclosure.normalization.static_closure_abs_bound) || ...
            ~isnan(enclosure.normalization.static_closure_abs_bound)
        fail(invalidId, ...
            ['An unvalidated normalization must use NaN rather than an ', ...
            'unproved static-closure bound.']);
    end
    if ~isstruct(enclosure.normalization.proof) || ...
            ~isscalar(enclosure.normalization.proof)
        fail(invalidId, ...
            'enclosure.normalization.proof must be a scalar struct.');
    end
end
end

function [lower,upper,finite] = denominatorIntervals(wInterval,Jq)
lower = zeros(size(Jq));
upper = zeros(size(Jq));
for k = 1:numel(Jq)
    lower(k) = nextDownScalar(wInterval(1)+Jq(k));
    upper(k) = nextUpScalar(wInterval(2)+Jq(k));
end
finite = all(isfinite(lower),'all') && all(isfinite(upper),'all');
end

function interval = reciprocalAverageInterval(denLower,denUpper)
sumLower = 0;
sumUpper = 0;
for k = 1:numel(denLower)
    reciprocalLower = nextDownScalar(1/denUpper(k));
    reciprocalUpper = nextUpScalar(1/denLower(k));
    sumLower = nextDownScalar(sumLower+reciprocalLower);
    sumUpper = nextUpScalar(sumUpper+reciprocalUpper);
end
n = numel(denLower);
interval = [nextDownScalar(sumLower/n),nextUpScalar(sumUpper/n)];
end

function value = intervalAbsLower(interval)
if interval(1) > 0
    value = interval(1);
elseif interval(2) < 0
    value = -interval(2);
else
    value = 0;
end
end

function y = nextUpScalar(x)
if isnan(x) || x == Inf
    y = x;
elseif x == -Inf
    y = -realmax('double');
elseif x == 0
    y = typecast(uint64(1),'double');
else
    bits = typecast(x,'uint64');
    if x > 0
        bits = bits+uint64(1);
    else
        bits = bits-uint64(1);
    end
    y = typecast(bits,'double');
end
end

function y = nextDownScalar(x)
if isnan(x) || x == -Inf
    y = x;
elseif x == Inf
    y = realmax('double');
elseif x == 0
    y = -typecast(uint64(1),'double');
else
    bits = typecast(x,'uint64');
    if x > 0
        bits = bits-uint64(1);
    else
        bits = bits+uint64(1);
    end
    y = typecast(bits,'double');
end
end

function value = positiveScalar(value,name,invalidId)
if ~isa(value,'double') || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value <= 0
    fail(invalidId,'%s must be a finite positive double scalar.',name);
end
end

function value = nonnegativeScalar(value,name,invalidId)
if ~isa(value,'double') || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value < 0
    fail(invalidId,'%s must be a finite nonnegative double scalar.',name);
end
end

function value = positiveInteger(value,name,invalidId)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value < 1 || value ~= floor(value)
    fail(invalidId,'%s must be a positive integer.',name);
end
end

function value = logicalScalar(value,name,invalidId)
if ~islogical(value) || ~isscalar(value)
    fail(invalidId,'%s must be a scalar logical.',name);
end
end

function value = textScalar(value,name,invalidId)
if ischar(value) && isrow(value)
    value = string(value);
elseif isstring(value) && isscalar(value) && ~ismissing(value)
    value = string(value);
else
    fail(invalidId,'%s must be a nonmissing text scalar.',name);
end
if strlength(value) == 0
    fail(invalidId,'%s must not be empty.',name);
end
end

function values = textVector(value,name,invalidId)
if ischar(value) && isrow(value)
    values = string(value);
elseif isstring(value) && isvector(value) && ...
        ~any(ismissing(value),'all')
    values = value(:);
else
    fail(invalidId,'%s must be a nonmissing text vector.',name);
end
if isempty(values) || any(strlength(values) == 0)
    fail(invalidId,'%s must contain nonempty entries.',name);
end
end

function requireFields(value,names,label,invalidId)
for k = 1:numel(names)
    if ~isfield(value,names{k})
        fail(invalidId,'%s.%s is required.',label,names{k});
    end
end
end

function id = makeEventId(sourceId,first,second)
id = sourceId+"#event-"+compose("%06d",first)+"-"+compose("%06d",second);
end

function fail(identifier,message,varargin)
error(identifier,message,varargin{:});
end
