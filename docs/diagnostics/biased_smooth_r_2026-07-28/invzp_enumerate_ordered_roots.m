function report = invzp_enumerate_ordered_roots(ctx, h, seeds, spec)
%INVZP_ENUMERATE_ORDERED_ROOTS Run explicit fixed-h seeds and cluster roots.
%
%   This diagnostic driver never invents a seed. It runs seeds in lexical
%   ID order, preserves every submitted attempt, and exports no roots when
%   equivalence in the complete independent [Sigma;K0] coordinates is
%   ambiguous. K and lambda are derived and are not separate coordinates.

invalidId = 'invzp:RootEnumerator:InvalidInput';
validateInputs(ctx,h,seeds,spec,invalidId);
nw = numel(ctx.wn);
ids = strings(numel(seeds),1);
submitted = (1:numel(seeds)).';
for k = 1:numel(seeds)
    ids(k) = normalizeId(seeds(k).id,invalidId);
    state = seeds(k).state;
    if ~isstruct(state) || ~isscalar(state) || ...
            ~isfield(state,'Sigma') || ~isfield(state,'K0s') || ...
            ~isnumeric(state.Sigma) || ~isreal(state.Sigma) || ...
            ~isvector(state.Sigma) || numel(state.Sigma) ~= nw || ...
            any(~isfinite(state.Sigma),'all') || ...
            ~isnumeric(state.K0s) || ~isreal(state.K0s) || ...
            ~isscalar(state.K0s) || ~isfinite(state.K0s)
        error(invalidId, ...
            'Every seed.state requires finite real nw-vector Sigma and scalar K0s.');
    end
end
if numel(unique(ids)) ~= numel(ids)
    error(invalidId,'seed.id values must be unique.');
end
[sortedIds,solveOrder] = sort(ids);

report = struct( ...
    'schema','invzp_fixed_h_root_enumerator/v2', ...
    'status','not_run', ...
    'status_detail','', ...
    'h',h,'node_meta',struct(),'spec',spec, ...
    'submitted_seed_ids',ids, ...
    'solve_seed_ids',sortedIds, ...
    'attempts',repmat(blankAttempt(),1,0), ...
    'accepted_attempt_ids',strings(0,1), ...
    'pairwise',emptyPairwise(), ...
    'clusters',repmat(blankCluster(),1,0), ...
    'ambiguous_components',repmat(blankAmbiguous(),1,0), ...
    'roots',repmat(blankRoot(),1,0));

try
    [node,nodeMeta] = invz_ordered_make_node(ctx,h);
catch ME
    report.node_meta = exceptionInfo(ME);
    for k = 1:numel(solveOrder)
        index = solveOrder(k);
        attempt = makeNotRunAttempt( ...
            seeds(index),ids(index),submitted(index),k, ...
            "node_exception:"+string(ME.identifier));
        report.attempts(end+1) = attempt;
    end
    report.status = 'node_construction_failed';
    report.status_detail = char("node_exception:"+string(ME.identifier));
    return
end
report.node_meta = nodeMeta;
if isempty(node)
    for k = 1:numel(solveOrder)
        index = solveOrder(k);
        attempt = makeNotRunAttempt( ...
            seeds(index),ids(index),submitted(index),k,nodeMeta.status);
        report.attempts(end+1) = attempt;
    end
    report.status = 'node_unavailable';
    report.status_detail = nodeMeta.status;
    return
end

newtonOpts = getf(spec,'newton',struct());
for k = 1:numel(solveOrder)
    index = solveOrder(k);
    seedState = seeds(index).state;
    try
        [state,info] = invz_ordered_node_newton(node,seedState,newtonOpts);
    catch ME
        attempt = blankAttempt();
        attempt.seed_id = ids(index);
        attempt.submitted_index = submitted(index);
        attempt.solve_order = k;
        attempt.seed_state = seedState;
        attempt.reason = "solver_exception:"+string(ME.identifier);
        attempt.returned_state = seedState;
        attempt.info = struct('accepted',false,'reason','solver_exception', ...
            'exception_identifier',ME.identifier, ...
            'exception_message',ME.message);
        report.attempts(end+1) = attempt;
        continue
    end
    attempt = blankAttempt();
    attempt.seed_id = ids(index);
    attempt.submitted_index = submitted(index);
    attempt.solve_order = k;
    attempt.seed_state = seedState;
    attempt.returned_state = state;
    attempt.info = info;
    if ~validSolverOutput(state,info,nw)
        attempt.reason = "invalid_solver_output";
        report.attempts(end+1) = attempt;
        continue
    end
    if info.accepted
        attempt.accepted = info.audit.accepted;
    end
    attempt.reason = string(info.reason);
    report.attempts(end+1) = attempt;
end

accepted = find([report.attempts.accepted]);
report.accepted_attempt_ids = string({report.attempts(accepted).seed_id}).';
if isempty(accepted)
    report.status = 'no_admissible_roots';
    report.status_detail = 'no explicit seed produced an A-D-accepted root';
    return
end

candidates = repmat(struct('id',"",'state',struct()),1,numel(accepted));
for k = 1:numel(accepted)
    candidates(k).id = report.attempts(accepted(k)).seed_id;
    candidates(k).state = report.attempts(accepted(k)).returned_state;
end
clustered = invzp_cluster_ordered_full_states(candidates,spec.cluster);
report.pairwise = clustered.pairwise;
report.clusters = clustered.clusters;
report.ambiguous_components = clustered.ambiguous_components;
if ~strcmp(clustered.status,'ok')
    report.status = clustered.status;
    report.status_detail = clustered.status_detail;
    return
end

for k = 1:numel(clustered.clusters)
    cluster = clustered.clusters(k);
    medoidAttempt = find([report.attempts.accepted] & ...
        string({report.attempts.seed_id}) == cluster.medoid_id,1);
    root = blankRoot();
    root.id = cluster.id;
    root.member_seed_ids = cluster.member_ids;
    root.medoid_seed_id = cluster.medoid_id;
    root.cluster_diameter = cluster.diameter;
    root.state = report.attempts(medoidAttempt).returned_state;
    root.info = report.attempts(medoidAttempt).info;
    report.roots(end+1) = root;
end
report.status = 'ok';
report.status_detail = sprintf('%d distinct accepted root(s)',numel(report.roots));
end

function validateInputs(ctx,h,seeds,spec,invalidId)
if ~isstruct(ctx) || ~isscalar(ctx) || ~isfield(ctx,'schema') || ...
        ~strcmp(ctx.schema,'invzp_ordered_node_context/v1')
    error(invalidId,'ctx must come from invz_ordered_node_context.');
end
if ~isnumeric(h) || ~isreal(h) || ~isscalar(h) || ~isfinite(h)
    error(invalidId,'h must be a finite real scalar.');
end
if ~isstruct(seeds) || isempty(seeds) || ...
        ~isfield(seeds,'id') || ~isfield(seeds,'state')
    error(invalidId,'seeds must be a nonempty struct array with id and state.');
end
if ~isstruct(spec) || ~isscalar(spec) || ~isfield(spec,'cluster')
    error(invalidId,'spec.cluster is required.');
end
invzp_cluster_ordered_full_states(struct([]),spec.cluster);
nw = numel(ctx.wn);
if ~isscalar(spec.cluster.sigma_scale) && ...
        numel(spec.cluster.sigma_scale) ~= nw
    error(invalidId, ...
        'spec.cluster.sigma_scale must be scalar or match ctx.wn length.');
end
if isfield(spec,'newton') && ...
        (~isstruct(spec.newton) || ~isscalar(spec.newton))
    error(invalidId,'spec.newton must be a scalar struct.');
end
end

function attempt = makeNotRunAttempt(seed,id,submittedIndex,solveOrder,reason)
attempt = blankAttempt();
attempt.seed_id = id;
attempt.submitted_index = submittedIndex;
attempt.solve_order = solveOrder;
attempt.seed_state = seed.state;
attempt.reason = "not_run:"+string(reason);
end

function attempt = blankAttempt()
state = struct('Sigma',[],'K',[],'lam',[],'K0s',NaN);
attempt = struct('seed_id',"",'submitted_index',NaN,'solve_order',NaN, ...
    'seed_state',struct('Sigma',[],'K0s',NaN), ...
    'accepted',false,'reason',"",'returned_state',state,'info',struct());
end

function root = blankRoot()
state = struct('Sigma',[],'K',[],'lam',[],'K0s',NaN);
root = struct('id',"",'member_seed_ids',strings(0,1), ...
    'medoid_seed_id',"",'cluster_diameter',NaN, ...
    'state',state,'info',struct());
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

function valid = validSolverOutput(state,info,nw)
stateValid = isstruct(state) && isscalar(state) && ...
    isfield(state,'Sigma') && isnumeric(state.Sigma) && ...
    isreal(state.Sigma) && isvector(state.Sigma) && ...
    numel(state.Sigma) == nw && all(isfinite(state.Sigma),'all') && ...
    isfield(state,'K0s') && isnumeric(state.K0s) && ...
    isreal(state.K0s) && isscalar(state.K0s) && isfinite(state.K0s);
baseValid = stateValid && isstruct(info) && isscalar(info) && ...
    isfield(info,'accepted') && islogical(info.accepted) && ...
    isscalar(info.accepted) && isfield(info,'audit') && ...
    isfield(info,'reason') && validText(info.reason);
if ~baseValid
    valid = false;
    return
end
auditValid = isstruct(info.audit) && isscalar(info.audit) && ...
    isfield(info.audit,'accepted') && islogical(info.audit.accepted) && ...
    isscalar(info.audit.accepted);
valid = auditValid || (~info.accepted && isempty(info.audit));
end

function valid = validText(value)
valid = (ischar(value) && (isempty(value) || isrow(value))) || ...
    (isstring(value) && isscalar(value) && ~ismissing(value));
end

function info = exceptionInfo(exception)
info = struct('status','exception', ...
    'exception_identifier',exception.identifier, ...
    'exception_message',exception.message);
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
    error(invalidId,'seed.id must be a nonempty character row or string scalar.');
end
id = string(value);
end

function value = getf(s,name,default)
if isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
