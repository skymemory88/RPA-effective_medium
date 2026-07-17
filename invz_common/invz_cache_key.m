function key = invz_cache_key(prefix, dpRng, qmeta, pkey)
%INVZ_CACHE_KEY Weak-hash cache filename, reproducing the invz_projected scheme.
%   key = INVZ_CACHE_KEY(prefix, dpRng, qmeta, pkey) returns
%   '<prefix>_<dpRng>_<hash(qmeta)>_<hash(pkey)>.mat', using the SAME
%   sprintf('%dv_%08x', numel(v), typecast(single(sum(v.*(1:numel(v))')), 'uint32'))
%   weak hash as the local hash_vec in invz_jq_modes.m / invz_odd_blocks.m.
%   qmeta is the q-mesh descriptor with any convention-distinguishing scalar
%   folded in, e.g. qmeta = [qvec(:); convcode], so distinct tensor
%   conventions at the same q-mesh never alias to the same key; pkey is the
%   physical/numerical parameter vector (ion geometry, couplings, schema
%   version, ...).
%
%   WEAK HASH: collisions are possible (32-bit truncation of a single-
%   precision, index-weighted sum). A key match is only a CACHE-HIT
%   CANDIDATE -- callers MUST isequal-verify the loaded qmeta/pkey against
%   the ones just hashed before trusting any other field in the cache file
%   (mirror the isequal(S.pkey,pkey) && isequal(S.qvec,qvec) guard in
%   invz_jq_modes.m).
%
%   Tensor-only: invz_projected keeps its own local hash_vec unchanged (pure
%   move); this helper is for invz_tensor and must never be retrofitted there.
key = sprintf('%s_%d_%s_%s.mat', prefix, dpRng, hash_vec(qmeta(:)), hash_vec(pkey(:)));
end

function h = hash_vec(v)
h = sprintf('%dv_%08x', numel(v), ...
    typecast(single(sum(v.*(1:numel(v))')), 'uint32'));
end
