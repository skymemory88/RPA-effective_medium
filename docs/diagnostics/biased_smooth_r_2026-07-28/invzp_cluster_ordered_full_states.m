function result = invzp_cluster_ordered_full_states(candidates, spec)
%INVZP_CLUSTER_ORDERED_FULL_STATES Pure independent-state root equivalence oracle.
%
%   CANDIDATES requires unique id and the complete independent coordinates
%   state={Sigma,K0s}; K and lambda are derived. The relation has an
%   explicit uncertainty band:
%     d <= same_tol       same-root evidence
%     d > distinct_tol   distinct-root evidence
%     otherwise          unresolved

invalidId = 'invzp:RootCluster:InvalidInput';
validateSpec(spec,invalidId);
emptyIds = strings(0,1);
result = struct( ...
    'schema','invzp_ordered_full_state_cluster/v1', ...
    'status','no_admissible_roots', ...
    'status_detail','no accepted candidates were supplied', ...
    'candidate_ids',emptyIds, ...
    'distance_matrix',zeros(0), ...
    'pairwise',emptyPairwise(), ...
    'clusters',repmat(blankCluster(),1,0), ...
    'ambiguous_components',repmat(blankAmbiguous(),1,0), ...
    'spec',spec);
if ~isstruct(candidates)
    error(invalidId,'candidates must be a struct array.');
end
if isempty(candidates), return; end
if ~isfield(candidates,'id') || ~isfield(candidates,'state')
    error(invalidId,'Every candidate requires id and state fields.');
end

ncandidate = numel(candidates);
ids = strings(ncandidate,1);
nw = [];
for k = 1:ncandidate
    ids(k) = normalizeId(candidates(k).id,invalidId);
    state = candidates(k).state;
    if ~isstruct(state) || ~isscalar(state) || ...
            ~isfield(state,'Sigma') || ~isfield(state,'K0s')
        error(invalidId,'candidate.state requires Sigma and K0s.');
    end
    if isempty(nw)
        nw = numel(state.Sigma);
    end
    if ~isnumeric(state.Sigma) || ~isreal(state.Sigma) || ...
            ~isvector(state.Sigma) || isempty(state.Sigma) || ...
            numel(state.Sigma) ~= nw || ...
            any(~isfinite(state.Sigma),'all')
        error(invalidId, ...
            'candidate.state.Sigma must be a nonempty finite real common-length vector.');
    end
    if ~isnumeric(state.K0s) || ~isreal(state.K0s) || ...
            ~isscalar(state.K0s) || ~isfinite(state.K0s)
        error(invalidId,'candidate.state.K0s must be a finite real scalar.');
    end
end
if numel(unique(ids)) ~= ncandidate
    error(invalidId,'candidate.id values must be unique.');
end
if ~isscalar(spec.sigma_scale) && numel(spec.sigma_scale) ~= nw
    error(invalidId, ...
        'spec.sigma_scale must be scalar or match candidate Sigma length.');
end
[ids,order] = sort(ids);
candidates = candidates(order);
result.candidate_ids = ids;

distance = zeros(ncandidate);
pairI = strings(0,1);
pairJ = strings(0,1);
pairD = zeros(0,1);
pairRelation = strings(0,1);
relation = strings(ncandidate);
for i = 1:ncandidate
    relation(i,i) = "same_evidence";
    for j = i+1:ncandidate
        dsigma = (candidates(i).state.Sigma(:)- ...
            candidates(j).state.Sigma(:))./spec.sigma_scale(:);
        dK0 = (candidates(i).state.K0s-candidates(j).state.K0s)/spec.k0_scale;
        dij = sqrt(mean(abs(dsigma).^2)+abs(dK0).^2);
        distance(i,j) = dij;
        distance(j,i) = dij;
        if dij <= spec.same_tol
            label = "same_evidence";
        elseif dij > spec.distinct_tol
            label = "distinct_evidence";
        else
            label = "unresolved_relation";
        end
        relation(i,j) = label;
        relation(j,i) = label;
        pairI(end+1,1) = ids(i); %#ok<AGROW>
        pairJ(end+1,1) = ids(j); %#ok<AGROW>
        pairD(end+1,1) = dij; %#ok<AGROW>
        pairRelation(end+1,1) = label; %#ok<AGROW>
    end
end
result.distance_matrix = distance;
result.pairwise = table(pairI,pairJ,pairD,pairRelation, ...
    'VariableNames',{'id_i','id_j','distance','relation'});

adjacency = relation ~= "distinct_evidence";
components = connectedComponents(adjacency);
validClusters = repmat(blankCluster(),1,0);
ambiguous = repmat(blankAmbiguous(),1,0);
for ic = 1:numel(components)
    members = components{ic};
    memberRelation = relation(members,members);
    allSame = all(memberRelation == "same_evidence",'all');
    if ~allSame
        [memberIds,memberOrder] = sort(ids(members));
        members = members(memberOrder);
        item = blankAmbiguous();
        item.member_ids = memberIds;
        item.relations = relation(members,members);
        item.distance_matrix = distance(members,members);
        ambiguous(end+1) = item; %#ok<AGROW>
        continue
    end

    sums = sum(distance(members,members).^2,2);
    best = min(sums);
    tiedMembers = members(sums == best);
    if numel(tiedMembers) > 1
        [~,order] = sort(ids(tiedMembers));
        tiedMembers = tiedMembers(order);
    end
    memberIds = sort(ids(members));
    cluster = blankCluster();
    cluster.id = "cluster:"+strjoin(memberIds,"+");
    cluster.member_ids = memberIds;
    cluster.medoid_id = ids(tiedMembers(1));
    cluster.diameter = max(distance(members,members),[],'all');
    validClusters(end+1) = cluster; %#ok<AGROW>
end

[~,order] = sort([validClusters.id]);
result.clusters = validClusters(order);
if ~isempty(ambiguous)
    result.status = 'cluster_ambiguous';
    result.status_detail = ...
        'at least one tolerance component contains an unresolved or distinct pair';
    result.ambiguous_components = ambiguous;
    return
end
result.status = 'ok';
result.status_detail = sprintf('%d distinct accepted root cluster(s)', ...
    numel(result.clusters));
end

function validateSpec(spec,invalidId)
required = {'sigma_scale','k0_scale','same_tol','distinct_tol'};
if ~isstruct(spec) || ~isscalar(spec)
    error(invalidId,'spec must be a scalar struct.');
end
for k = 1:numel(required)
    if ~isfield(spec,required{k})
        error(invalidId,'spec.%s is required.',required{k});
    end
end
if ~isnumeric(spec.sigma_scale) || ~isreal(spec.sigma_scale) || ...
        ~isvector(spec.sigma_scale) || isempty(spec.sigma_scale) || ...
        any(~isfinite(spec.sigma_scale),'all') || any(spec.sigma_scale <= 0)
    error(invalidId,'spec.sigma_scale must be a positive finite vector.');
end
positive(spec.k0_scale,'k0_scale',invalidId);
nonnegative(spec.same_tol,'same_tol',invalidId);
nonnegative(spec.distinct_tol,'distinct_tol',invalidId);
if spec.distinct_tol < spec.same_tol
    error(invalidId,'spec.distinct_tol must be at least spec.same_tol.');
end
end

function components = connectedComponents(adjacency)
n = size(adjacency,1);
seen = false(n,1);
components = {};
for start = 1:n
    if seen(start), continue; end
    queue = start;
    seen(start) = true;
    members = zeros(0,1);
    while ~isempty(queue)
        current = queue(1);
        queue(1) = [];
        members(end+1,1) = current; %#ok<AGROW>
        neighbours = find(adjacency(current,:)).';
        unseen = neighbours(~seen(neighbours));
        seen(unseen) = true;
        queue = [queue;unseen]; %#ok<AGROW>
    end
    components{end+1,1} = sort(members); %#ok<AGROW>
end
end

function tableOut = emptyPairwise()
tableOut = table(strings(0,1),strings(0,1),zeros(0,1),strings(0,1), ...
    'VariableNames',{'id_i','id_j','distance','relation'});
end

function cluster = blankCluster()
cluster = struct('id',"",'member_ids',strings(0,1),'medoid_id',"", ...
    'diameter',NaN);
end

function item = blankAmbiguous()
item = struct('member_ids',strings(0,1), ...
    'relations',strings(0),'distance_matrix',zeros(0));
end

function id = normalizeId(value,invalidId)
if ischar(value)
    valid = isrow(value) && ~isempty(value);
elseif isstring(value)
    valid = isscalar(value) && ~ismissing(value) && strlength(value) > 0;
else
    valid = false;
end
if ~valid
    error(invalidId,'candidate.id must be a nonempty character row or string scalar.');
end
id = string(value);
end

function value = positive(value,name,invalidId)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value <= 0
    error(invalidId,'spec.%s must be a finite positive scalar.',name);
end
end

function value = nonnegative(value,name,invalidId)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value < 0
    error(invalidId,'spec.%s must be a finite nonnegative scalar.',name);
end
end
