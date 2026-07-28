function result = invzp_assemble_audited_graph(vertices,segments,spec)
%INVZP_ASSEMBLE_AUDITED_GRAPH Assemble only explicitly certified branch evidence.
%
%   RESULT = INVZP_ASSEMBLE_AUDITED_GRAPH(VERTICES,SEGMENTS,SPEC) is an
%   evidence-normalization layer. It does not solve, audit, cluster, merge,
%   interpolate, extrapolate, split folds, construct Jensen sections, or
%   select a branch. Every returned edge corresponds to one supplied
%   adjacency whose registered continuity gates all pass. Missing evidence
%   is retained as a blocker and never inferred from geometric proximity.

invalidId = 'invzp:AuditedGraph:InvalidInput';
cfg = validateSpec(spec,invalidId);
vertices = normalizeVertices(vertices,cfg,invalidId);
segments = normalizeSegments(segments,vertices,cfg,invalidId);

result = struct( ...
    'schema','invzp_audited_graph/v1', ...
    'status','no_admissible_vertices', ...
    'status_detail','no accepted vertices were supplied', ...
    'spec_version',cfg.version, ...
    'vertices',vertices, ...
    'segments',segments, ...
    'directed_edges',repmat(blankEdge(),1,0), ...
    'undirected_edges',repmat(blankUndirectedEdge(),1,0), ...
    'components',repmat(blankComponent(),1,0), ...
    'component_index',zeros(0,1), ...
    'isolated_vertex_ids',strings(0,1), ...
    'missing_required_vertex_ids',strings(0,1), ...
    'adjacency_blockers',repmat(blankAdjacencyBlocker(),1,0), ...
    'termination_blockers',repmat(blankTerminationBlocker(),1,0), ...
    'complete_candidate_claim',false);

if isempty(vertices)
    if ~isempty(segments)
        fail(invalidId, ...
            'segments cannot be supplied without accepted vertices.');
    end
    result.missing_required_vertex_ids = cfg.required_vertex_ids;
    return
end

vertexIds = reshape(string({vertices.id}),[],1);
edges = repmat(blankEdge(),1,0);
adjacencyBlockers = repmat(blankAdjacencyBlocker(),1,0);
terminationBlockers = repmat(blankTerminationBlocker(),1,0);

for k = 1:numel(segments)
    segment = segments(k);
    ids = segment.vertex_ids;
    for j = 1:numel(segment.adjacency)
        certificate = segment.adjacency(j);
        fromIndex = find(vertexIds == ids(j),1);
        toIndex = find(vertexIds == ids(j+1),1);
        distance = stateDistance( ...
            vertices(fromIndex),vertices(toIndex), ...
            cfg.sigma_scale,cfg.k0_scale,cfg.m_scale);
        reasons = adjacencyFailureReasons(certificate,distance,cfg);
        if isempty(reasons)
            edge = blankEdge();
            edge.id = segment.id+"#"+compose("%06d",j);
            edge.segment_id = segment.id;
            edge.adjacency_index = j;
            edge.from_id = ids(j);
            edge.to_id = ids(j+1);
            edge.from_index = fromIndex;
            edge.to_index = toIndex;
            edge.state_distance = distance;
            edge.predictor_distance = certificate.predictor_distance;
            edge.tangent_required = certificate.tangent_required;
            edge.same_h = certificate.same_h;
            edge.tangent_line_error = certificate.tangent_line_error;
            edge.source_record_id = certificate.source_record_id;
            edge.evidence_payload = certificate.evidence_payload;
            edges(end+1) = edge; %#ok<AGROW>
        else
            blocker = blankAdjacencyBlocker();
            blocker.segment_id = segment.id;
            blocker.adjacency_index = j;
            blocker.from_id = ids(j);
            blocker.to_id = ids(j+1);
            blocker.reasons = reasons;
            blocker.source_record_id = certificate.source_record_id;
            adjacencyBlockers(end+1) = blocker; %#ok<AGROW>
        end
    end
    if ~segment.termination.resolved
        blocker = blankTerminationBlocker();
        blocker.segment_id = segment.id;
        blocker.status = segment.termination.status;
        blocker.reason = segment.termination.reason;
        blocker.source_record_id = segment.termination.source_record_id;
        terminationBlockers(end+1) = blocker; %#ok<AGROW>
    end
end

[components,componentIndex,isolatedIds] = ...
    connectedComponents(vertices,edges);
undirectedEdges = collapseUndirected(edges);
missing = cfg.required_vertex_ids( ...
    ~ismember(cfg.required_vertex_ids,vertexIds));

result.directed_edges = edges;
result.undirected_edges = undirectedEdges;
result.components = components;
result.component_index = componentIndex;
result.isolated_vertex_ids = isolatedIds;
result.missing_required_vertex_ids = missing;
result.adjacency_blockers = adjacencyBlockers;
result.termination_blockers = terminationBlockers;

incomplete = ~isempty(adjacencyBlockers) || ...
    ~isempty(terminationBlockers) || ~isempty(isolatedIds) || ...
    ~isempty(missing);
if incomplete
    result.status = 'evidence_incomplete';
    result.status_detail = sprintf( ...
        ['%d accepted vertices, %d certified directed edges, ', ...
        '%d adjacency blockers, %d termination blockers, ', ...
        '%d isolated vertices, %d missing required vertices'], ...
        numel(vertices),numel(edges),numel(adjacencyBlockers), ...
        numel(terminationBlockers),numel(isolatedIds),numel(missing));
else
    result.status = 'ok';
    result.status_detail = sprintf( ...
        '%d accepted vertices and %d certified directed edges assembled', ...
        numel(vertices),numel(edges));
end
end

function cfg = validateSpec(spec,invalidId)
if ~isstruct(spec) || ~isscalar(spec)
    fail(invalidId,'spec must be a scalar struct.');
end
required = {'version','certificate_schema','allowed_source_schemas', ...
    'sigma_scale','k0_scale','m_scale','state_distance_max', ...
    'predictor_tube_max','tangent_line_error_max', ...
    'required_event_margins','required_vertex_ids', ...
    'resolved_termination_statuses'};
for k = 1:numel(required)
    if ~isfield(spec,required{k})
        fail(invalidId,'spec.%s is required.',required{k});
    end
end
version = normalizeId(spec.version,'spec.version',invalidId);
certificateSchema = normalizeId( ...
    spec.certificate_schema,'spec.certificate_schema',invalidId);
allowedSchemas = normalizeIdVector( ...
    spec.allowed_source_schemas,'spec.allowed_source_schemas', ...
    false,invalidId);
if numel(unique(allowedSchemas)) ~= numel(allowedSchemas)
    fail(invalidId,'spec.allowed_source_schemas must be unique.');
end
sigmaScale = positiveVector(spec.sigma_scale,'spec.sigma_scale',invalidId);
k0Scale = positiveScalar(spec.k0_scale,'spec.k0_scale',invalidId);
mScale = positiveScalar(spec.m_scale,'spec.m_scale',invalidId);
stateDistanceMax = positiveScalar( ...
    spec.state_distance_max,'spec.state_distance_max',invalidId);
predictorTubeMax = positiveScalar( ...
    spec.predictor_tube_max,'spec.predictor_tube_max',invalidId);
tangentLineErrorMax = positiveScalar( ...
    spec.tangent_line_error_max,'spec.tangent_line_error_max',invalidId);
marginNames = normalizeIdVector( ...
    spec.required_event_margins,'spec.required_event_margins', ...
    true,invalidId);
for k = 1:numel(marginNames)
    if ~isvarname(char(marginNames(k)))
        fail(invalidId, ...
            'spec.required_event_margins entries must be MATLAB field names.');
    end
end
if numel(unique(marginNames)) ~= numel(marginNames)
    fail(invalidId,'spec.required_event_margins must be unique.');
end
requiredVertexIds = normalizeIdVector( ...
    spec.required_vertex_ids,'spec.required_vertex_ids',true,invalidId);
if numel(unique(requiredVertexIds)) ~= numel(requiredVertexIds)
    fail(invalidId,'spec.required_vertex_ids must be unique.');
end
resolvedTerminationStatuses = normalizeIdVector( ...
    spec.resolved_termination_statuses, ...
    'spec.resolved_termination_statuses',false,invalidId);
if numel(unique(resolvedTerminationStatuses)) ~= ...
        numel(resolvedTerminationStatuses)
    fail(invalidId,'spec.resolved_termination_statuses must be unique.');
end
cfg = struct( ...
    'version',version, ...
    'certificate_schema',certificateSchema, ...
    'allowed_source_schemas',allowedSchemas, ...
    'sigma_scale',sigmaScale, ...
    'k0_scale',k0Scale, ...
    'm_scale',mScale, ...
    'state_distance_max',stateDistanceMax, ...
    'predictor_tube_max',predictorTubeMax, ...
    'tangent_line_error_max',tangentLineErrorMax, ...
    'required_event_margins',marginNames, ...
    'required_vertex_ids',requiredVertexIds, ...
    'resolved_termination_statuses',resolvedTerminationStatuses);
end

function vertices = normalizeVertices(vertices,cfg,invalidId)
if ~isstruct(vertices)
    fail(invalidId,'vertices must be a struct array.');
end
if isempty(vertices)
    vertices = repmat(blankVertex(),1,0);
    return
end
required = {'id','h','Sigma','K0s','m','source_schema','source_id', ...
    'certificate'};
for k = 1:numel(required)
    if ~isfield(vertices,required{k})
        fail(invalidId,'Every vertex must contain field "%s".',required{k});
    end
end
normalized = repmat(blankVertex(),1,numel(vertices));
ids = strings(numel(vertices),1);
for k = 1:numel(vertices)
    vertex = vertices(k);
    ids(k) = normalizeId(vertex.id,'vertex.id',invalidId);
    h = nonnegativeScalar(vertex.h,'vertex.h',invalidId);
    Sigma = finiteVector(vertex.Sigma,'vertex.Sigma',invalidId);
    if numel(Sigma) ~= numel(cfg.sigma_scale)
        fail(invalidId, ...
            'vertex.Sigma has %d entries; spec.sigma_scale has %d.', ...
            numel(Sigma),numel(cfg.sigma_scale));
    end
    K0s = finiteScalar(vertex.K0s,'vertex.K0s',invalidId);
    m = finiteScalar(vertex.m,'vertex.m',invalidId);
    sourceSchema = normalizeId( ...
        vertex.source_schema,'vertex.source_schema',invalidId);
    if ~ismember(sourceSchema,cfg.allowed_source_schemas)
        fail(invalidId, ...
            'vertex.source_schema "%s" is not registered.',sourceSchema);
    end
    sourceId = normalizeId(vertex.source_id,'vertex.source_id',invalidId);
    certificate = normalizeVertexCertificate( ...
        vertex.certificate,cfg,invalidId);
    normalized(k) = struct( ...
        'id',ids(k),'h',h,'Sigma',Sigma,'K0s',K0s,'m',m, ...
        'source_schema',sourceSchema,'source_id',sourceId, ...
        'certificate',certificate);
end
if numel(unique(ids)) ~= numel(ids)
    fail(invalidId,'vertex.id values must be unique.');
end
[~,order] = sort(ids);
vertices = normalized(order);
end

function certificate = normalizeVertexCertificate(certificate,cfg,invalidId)
if ~isstruct(certificate) || ~isscalar(certificate)
    fail(invalidId,'vertex.certificate must be a scalar struct.');
end
required = {'schema','residual_inf','residual_tolerance', ...
    'audit_accepted','audit_payload','domain_ok','event_margins'};
for k = 1:numel(required)
    if ~isfield(certificate,required{k})
        fail(invalidId, ...
            'vertex.certificate.%s is required.',required{k});
    end
end
schema = normalizeId( ...
    certificate.schema,'vertex.certificate.schema',invalidId);
if schema ~= cfg.certificate_schema
    fail(invalidId, ...
        'vertex certificate schema "%s" does not match spec.',schema);
end
residual = nonnegativeScalar( ...
    certificate.residual_inf,'vertex.certificate.residual_inf',invalidId);
tolerance = positiveScalar( ...
    certificate.residual_tolerance, ...
    'vertex.certificate.residual_tolerance',invalidId);
auditAccepted = logicalScalar( ...
    certificate.audit_accepted, ...
    'vertex.certificate.audit_accepted',invalidId);
if ~isstruct(certificate.audit_payload) || ...
        ~isscalar(certificate.audit_payload)
    fail(invalidId, ...
        'vertex.certificate.audit_payload must be a scalar struct.');
end
domainOk = logicalScalar( ...
    certificate.domain_ok,'vertex.certificate.domain_ok',invalidId);
if ~(residual < tolerance && auditAccepted && domainOk)
    fail(invalidId, ...
        ['A graph vertex must carry a strict residual pass, an accepted ', ...
        'A-D audit, and a passing domain certificate.']);
end
if ~isstruct(certificate.event_margins) || ...
        ~isscalar(certificate.event_margins)
    fail(invalidId, ...
        'vertex.certificate.event_margins must be a scalar struct.');
end
margins = certificate.event_margins;
for k = 1:numel(cfg.required_event_margins)
    name = char(cfg.required_event_margins(k));
    if ~isfield(margins,name)
        fail(invalidId, ...
            'vertex.certificate.event_margins.%s is required.',name);
    end
    positiveScalar(margins.(name), ...
        ['vertex.certificate.event_margins.',name],invalidId);
end
certificate = struct( ...
    'schema',schema,'residual_inf',residual, ...
    'residual_tolerance',tolerance, ...
    'audit_accepted',auditAccepted, ...
    'audit_payload',certificate.audit_payload, ...
    'domain_ok',domainOk,'event_margins',margins);
end

function segments = normalizeSegments(segments,vertices,cfg,invalidId)
if ~isstruct(segments)
    fail(invalidId,'segments must be a struct array.');
end
if isempty(segments)
    segments = repmat(blankSegment(),1,0);
    return
end
required = {'id','vertex_ids','source_schema','source_id', ...
    'adjacency','termination'};
for k = 1:numel(required)
    if ~isfield(segments,required{k})
        fail(invalidId,'Every segment must contain field "%s".',required{k});
    end
end
vertexIds = reshape(string({vertices.id}),[],1);
normalized = repmat(blankSegment(),1,numel(segments));
ids = strings(numel(segments),1);
for k = 1:numel(segments)
    segment = segments(k);
    ids(k) = normalizeId(segment.id,'segment.id',invalidId);
    orderedIds = normalizeIdVector( ...
        segment.vertex_ids,'segment.vertex_ids',false,invalidId);
    if numel(orderedIds) < 2
        fail(invalidId,'segment.vertex_ids must contain at least two IDs.');
    elseif numel(unique(orderedIds)) ~= numel(orderedIds)
        fail(invalidId, ...
            'segment.vertex_ids must not repeat a vertex.');
    elseif any(~ismember(orderedIds,vertexIds))
        missing = orderedIds(~ismember(orderedIds,vertexIds));
        fail(invalidId,'segment references unknown vertices: %s.', ...
            strjoin(missing,', '));
    end
    sourceSchema = normalizeId( ...
        segment.source_schema,'segment.source_schema',invalidId);
    if ~ismember(sourceSchema,cfg.allowed_source_schemas)
        fail(invalidId, ...
            'segment.source_schema "%s" is not registered.',sourceSchema);
    end
    sourceId = normalizeId(segment.source_id,'segment.source_id',invalidId);
    adjacency = normalizeAdjacency( ...
        segment.adjacency,orderedIds,cfg,invalidId);
    termination = normalizeTermination(segment.termination,cfg,invalidId);
    h = zeros(numel(orderedIds),1);
    for j = 1:numel(orderedIds)
        h(j) = vertices(vertexIds == orderedIds(j)).h;
    end
    for j = 1:numel(adjacency)
        adjacency(j).same_h = h(j) == h(j+1);
    end
    dh = diff(h);
    monotoneH = all(dh > 0) || all(dh < 0);
    normalized(k) = struct( ...
        'id',ids(k),'vertex_ids',orderedIds, ...
        'source_schema',sourceSchema,'source_id',sourceId, ...
        'adjacency',adjacency,'termination',termination, ...
        'monotone_h',monotoneH);
end
if numel(unique(ids)) ~= numel(ids)
    fail(invalidId,'segment.id values must be unique.');
end
[~,order] = sort(ids);
segments = normalized(order);
end

function adjacency = normalizeAdjacency(adjacency,orderedIds,cfg,invalidId)
if ~isstruct(adjacency) || ...
        numel(adjacency) ~= numel(orderedIds)-1
    fail(invalidId, ...
        'segment.adjacency must contain one struct per adjacent vertex pair.');
end
required = {'schema','from_id','to_id','local_continuity_ok', ...
    'predictor_distance','predictor_tube_ok','tangent_required', ...
    'tangent_line_error','tangent_line_ok','event_bracket_ok', ...
    'source_record_id','evidence_payload'};
for k = 1:numel(required)
    if ~isfield(adjacency,required{k})
        fail(invalidId,'Every adjacency must contain field "%s".',required{k});
    end
end
normalized = repmat(blankAdjacency(),1,numel(adjacency));
for k = 1:numel(adjacency)
    item = adjacency(k);
    schema = normalizeId(item.schema,'adjacency.schema',invalidId);
    if schema ~= cfg.certificate_schema
        fail(invalidId, ...
            'adjacency schema "%s" does not match spec.',schema);
    end
    fromId = normalizeId(item.from_id,'adjacency.from_id',invalidId);
    toId = normalizeId(item.to_id,'adjacency.to_id',invalidId);
    if fromId ~= orderedIds(k) || toId ~= orderedIds(k+1)
        fail(invalidId, ...
            'adjacency IDs must exactly match the ordered segment vertices.');
    end
    localOk = logicalScalar( ...
        item.local_continuity_ok,'adjacency.local_continuity_ok',invalidId);
    predictorDistance = nonnegativeScalar( ...
        item.predictor_distance,'adjacency.predictor_distance',invalidId);
    predictorOk = logicalScalar( ...
        item.predictor_tube_ok,'adjacency.predictor_tube_ok',invalidId);
    tangentRequired = logicalScalar( ...
        item.tangent_required,'adjacency.tangent_required',invalidId);
    tangentError = optionalNonnegativeScalar( ...
        item.tangent_line_error,'adjacency.tangent_line_error',invalidId);
    tangentOk = logicalScalar( ...
        item.tangent_line_ok,'adjacency.tangent_line_ok',invalidId);
    if tangentRequired && ~isfinite(tangentError)
        fail(invalidId, ...
            'A tangent-required adjacency needs finite tangent_line_error.');
    end
    eventOk = logicalScalar( ...
        item.event_bracket_ok,'adjacency.event_bracket_ok',invalidId);
    sourceRecordId = normalizeId( ...
        item.source_record_id,'adjacency.source_record_id',invalidId);
    if ~isstruct(item.evidence_payload) || ~isscalar(item.evidence_payload)
        fail(invalidId, ...
            'adjacency.evidence_payload must be a scalar struct.');
    end
    normalized(k) = struct( ...
        'schema',schema,'from_id',fromId,'to_id',toId, ...
        'local_continuity_ok',localOk, ...
        'predictor_distance',predictorDistance, ...
        'predictor_tube_ok',predictorOk, ...
        'tangent_required',tangentRequired, ...
        'tangent_line_error',tangentError, ...
        'tangent_line_ok',tangentOk, ...
        'event_bracket_ok',eventOk, ...
        'source_record_id',sourceRecordId, ...
        'evidence_payload',item.evidence_payload, ...
        'same_h',false);
end
adjacency = normalized;
end

function termination = normalizeTermination(termination,cfg,invalidId)
if ~isstruct(termination) || ~isscalar(termination)
    fail(invalidId,'segment.termination must be a scalar struct.');
end
required = {'schema','resolved','status','reason','source_record_id','payload'};
for k = 1:numel(required)
    if ~isfield(termination,required{k})
        fail(invalidId,'segment.termination.%s is required.',required{k});
    end
end
schema = normalizeId( ...
    termination.schema,'segment.termination.schema',invalidId);
if schema ~= cfg.certificate_schema
    fail(invalidId, ...
        'termination schema "%s" does not match spec.',schema);
end
resolved = logicalScalar( ...
    termination.resolved,'segment.termination.resolved',invalidId);
status = normalizeId( ...
    termination.status,'segment.termination.status',invalidId);
registeredResolved = ismember(status,cfg.resolved_termination_statuses);
if resolved ~= registeredResolved
    fail(invalidId, ...
        ['segment.termination.resolved must agree with the registered ', ...
        'resolved termination statuses.']);
end
reason = normalizeText( ...
    termination.reason,'segment.termination.reason',invalidId,resolved);
sourceRecordId = normalizeId( ...
    termination.source_record_id, ...
    'segment.termination.source_record_id',invalidId);
if ~isstruct(termination.payload) || ~isscalar(termination.payload)
    fail(invalidId,'segment.termination.payload must be a scalar struct.');
end
termination = struct('schema',schema,'resolved',resolved,'status',status, ...
    'reason',reason,'source_record_id',sourceRecordId, ...
    'payload',termination.payload);
end

function reasons = adjacencyFailureReasons(certificate,distance,cfg)
reasons = strings(0,1);
if ~certificate.local_continuity_ok
    reasons(end+1) = "local_continuity_not_demonstrated";
end
if distance > cfg.state_distance_max
    reasons(end+1) = "state_distance_exceeded";
end
if ~certificate.predictor_tube_ok || ...
        certificate.predictor_distance > cfg.predictor_tube_max
    reasons(end+1) = "predictor_tube_failed";
end
if certificate.tangent_required && ...
        (~certificate.tangent_line_ok || ...
        certificate.tangent_line_error > cfg.tangent_line_error_max)
    reasons(end+1) = "tangent_line_failed";
end
if certificate.same_h && ~certificate.tangent_required
    reasons(end+1) = "same_h_tangent_evidence_missing";
end
if ~certificate.event_bracket_ok
    reasons(end+1) = "signed_event_bracket_missing";
end
end

function distance = stateDistance(first,second,sigmaScale,k0Scale,mScale)
dSigma = (second.Sigma-first.Sigma)./sigmaScale;
dK0 = (second.K0s-first.K0s)/k0Scale;
dM = (second.m-first.m)/mScale;
distance = sqrt(mean(abs(dSigma).^2)+abs(dK0).^2+abs(dM).^2);
end

function [components,componentIndex,isolatedIds] = ...
        connectedComponents(vertices,edges)
nvertex = numel(vertices);
parent = (1:nvertex).';
degree = zeros(nvertex,1);
for k = 1:numel(edges)
    first = edges(k).from_index;
    second = edges(k).to_index;
    degree(first) = degree(first)+1;
    degree(second) = degree(second)+1;
    parent = unite(parent,first,second);
end
for k = 1:nvertex
    parent(k) = rootOf(parent,k);
end
roots = unique(parent,'stable');
componentIndex = zeros(nvertex,1);
components = repmat(blankComponent(),1,numel(roots));
vertexIds = reshape(string({vertices.id}),[],1);
edgeIds = reshape(string({edges.id}),[],1);
for k = 1:numel(roots)
    members = find(parent == roots(k));
    componentIndex(members) = k;
    componentEdge = false(numel(edges),1);
    for j = 1:numel(edges)
        componentEdge(j) = any(edges(j).from_index == members);
    end
    h = [vertices(members).h].';
    components(k) = struct( ...
        'id',"component_"+compose("%06d",k), ...
        'vertex_ids',vertexIds(members), ...
        'directed_edge_ids',edgeIds(componentEdge), ...
        'h_min',min(h),'h_max',max(h));
end
isolatedIds = vertexIds(degree == 0);
end

function parent = unite(parent,first,second)
rootFirst = rootOf(parent,first);
rootSecond = rootOf(parent,second);
if rootFirst == rootSecond
    return
end
if rootFirst < rootSecond
    parent(rootSecond) = rootFirst;
else
    parent(rootFirst) = rootSecond;
end
end

function root = rootOf(parent,index)
root = index;
while parent(root) ~= root
    root = parent(root);
end
end

function undirected = collapseUndirected(edges)
undirected = repmat(blankUndirectedEdge(),1,0);
for k = 1:numel(edges)
    pair = sort([edges(k).from_id;edges(k).to_id]);
    found = 0;
    for j = 1:numel(undirected)
        if isequal(undirected(j).vertex_ids,pair)
            found = j;
            break
        end
    end
    if found == 0
        item = blankUndirectedEdge();
        item.vertex_ids = pair;
        item.directed_edge_ids = edges(k).id;
        undirected(end+1) = item; %#ok<AGROW>
    else
        undirected(found).directed_edge_ids(end+1,1) = edges(k).id;
    end
end
end

function value = normalizeId(value,name,invalidId)
if ischar(value)
    valid = isrow(value);
elseif isstring(value)
    valid = isscalar(value) && ~ismissing(value);
else
    valid = false;
end
if ~valid
    fail(invalidId,'%s must be a nonmissing character row or string scalar.',name);
end
value = string(value);
if strlength(value) == 0
    fail(invalidId,'%s must be nonempty.',name);
end
end

function values = normalizeIdVector(value,name,allowEmpty,invalidId)
if ischar(value)
    if isempty(value) && allowEmpty
        values = strings(0,1);
        return
    end
    if ~isrow(value)
        fail(invalidId,'%s must contain text IDs.',name);
    end
    values = string({value});
elseif iscellstr(value)
    values = string(value(:));
elseif isstring(value)
    values = value(:);
else
    fail(invalidId,'%s must be a string, character row, or cellstr.',name);
end
if isempty(values)
    if allowEmpty
        return
    end
    fail(invalidId,'%s must be nonempty.',name);
end
if any(ismissing(values)) || any(strlength(values) == 0)
    fail(invalidId,'%s contains a missing or empty ID.',name);
end
end

function value = normalizeText(value,name,invalidId,allowEmpty)
if ischar(value)
    valid = isempty(value) || isrow(value);
elseif isstring(value)
    valid = isscalar(value) && ~ismissing(value);
else
    valid = false;
end
if ~valid
    fail(invalidId,'%s must be a nonmissing text scalar.',name);
end
value = string(value);
if ~allowEmpty && strlength(value) == 0
    fail(invalidId,'%s must be nonempty for unresolved termination.',name);
end
end

function value = finiteScalar(value,name,invalidId)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value)
    fail(invalidId,'%s must be a finite real scalar.',name);
end
end

function value = nonnegativeScalar(value,name,invalidId)
value = finiteScalar(value,name,invalidId);
if value < 0
    fail(invalidId,'%s must be nonnegative.',name);
end
end

function value = positiveScalar(value,name,invalidId)
value = finiteScalar(value,name,invalidId);
if value <= 0
    fail(invalidId,'%s must be positive.',name);
end
end

function value = optionalNonnegativeScalar(value,name,invalidId)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        (~isnan(value) && (~isfinite(value) || value < 0))
    fail(invalidId,'%s must be NaN or a finite nonnegative scalar.',name);
end
end

function value = finiteVector(value,name,invalidId)
if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || ...
        isempty(value) || any(~isfinite(value),'all')
    fail(invalidId,'%s must be a nonempty finite real vector.',name);
end
value = value(:);
end

function value = positiveVector(value,name,invalidId)
value = finiteVector(value,name,invalidId);
if any(value <= 0)
    fail(invalidId,'%s must contain only positive values.',name);
end
end

function value = logicalScalar(value,name,invalidId)
if ~islogical(value) || ~isscalar(value)
    fail(invalidId,'%s must be a scalar logical.',name);
end
end

function item = blankVertex()
item = struct('id',"",'h',NaN,'Sigma',zeros(0,1),'K0s',NaN,'m',NaN, ...
    'source_schema',"",'source_id',"",'certificate',struct());
end

function item = blankSegment()
item = struct('id',"",'vertex_ids',strings(0,1), ...
    'source_schema',"",'source_id',"", ...
    'adjacency',repmat(blankAdjacency(),1,0), ...
    'termination',struct(),'monotone_h',false);
end

function item = blankAdjacency()
item = struct('schema',"",'from_id',"",'to_id',"", ...
    'local_continuity_ok',false,'predictor_distance',NaN, ...
    'predictor_tube_ok',false,'tangent_required',false, ...
    'tangent_line_error',NaN,'tangent_line_ok',false, ...
    'event_bracket_ok',false,'source_record_id',"", ...
    'evidence_payload',struct(),'same_h',false);
end

function item = blankEdge()
item = struct('id',"",'segment_id',"",'adjacency_index',0, ...
    'from_id',"",'to_id',"",'from_index',0,'to_index',0, ...
    'state_distance',NaN,'predictor_distance',NaN, ...
    'tangent_required',false,'same_h',false,'tangent_line_error',NaN, ...
    'source_record_id',"",'evidence_payload',struct());
end

function item = blankUndirectedEdge()
item = struct('vertex_ids',strings(0,1), ...
    'directed_edge_ids',strings(0,1));
end

function item = blankComponent()
item = struct('id',"",'vertex_ids',strings(0,1), ...
    'directed_edge_ids',strings(0,1),'h_min',NaN,'h_max',NaN);
end

function item = blankAdjacencyBlocker()
item = struct('segment_id',"",'adjacency_index',0, ...
    'from_id',"",'to_id',"",'reasons',strings(0,1), ...
    'source_record_id',"");
end

function item = blankTerminationBlocker()
item = struct('segment_id',"",'status',"",'reason',"", ...
    'source_record_id',"");
end

function fail(identifier,message,varargin)
error(identifier,message,varargin{:});
end
