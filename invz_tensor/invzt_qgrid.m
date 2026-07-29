function g = invzt_qgrid(n, conv)
%INVZT_QGRID Gamma-excluded q-grid contract (r.l.u.) for invz_tensor (LOCKED).
%   g = INVZT_QGRID(n, conv) returns a q-mesh descriptor struct:
%     g.qvec [nq,3] : q points in reduced coordinates (r.l.u.); the Gamma
%                     point ([0 0 0], to 1e-12) is always dropped.
%     g.w    [nq,1] : quadrature weights, UNIFORM over the kept points,
%                     sum(g.w) == 1 (the dropped-Gamma point simply is not
%                     represented -- see g.note).
%     g.conv        : echoes the input, 'halfopen' | 'legacy_inclusive'.
%     g.hash        : weak-hash grid-identity tag (n, convention, qvec
%                     content) -- distinct across n/conv/content, used for
%                     cache keys and lat/report provenance.
%     g.note        : human-readable string documenting the dropped-Gamma
%                     convention (both) and the duplicated reciprocal-
%                     periodic face convention (legacy_inclusive only).
%
%   CONVENTIONS (historical tensor Global Constraints, LOCKED -- never mix
%   conventions inside one
%   comparison):
%     'halfopen' (default for new tensor work)
%       ndgrid((0:n-1)/n) replicated over the 3 axes -- n DISTINCT points per
%       axis, no duplicate reciprocal-periodic boundary face. The Gamma point
%       ([0 0 0], occurring exactly once at the (0,0,0) grid node) is dropped.
%     'legacy_inclusive' (parity with projected/legacy anchors ONLY)
%       EXACTLY qVec_generator(ion.a, 'mode', 'grid', 'grid', [n n n],
%       'range', [-0.5 0.5]) -- i.e. the historical endpoint-inclusive
%       linspace(-0.5, 0.5, n) per axis, spacing 1/(n-1) -- then the SAME
%       Gamma-row filter as 'halfopen'. This range/endpoint combination
%       duplicates the +-0.5 reciprocal-periodic faces (an N^3 grid then has
%       only (N-1)^3 distinct points); NOT deduplicated here by design (see
%       g.note) -- 'legacy_inclusive' exists solely to reproduce the
%       projected production convention bit-for-bit, duplicate faces and all.
%
%   invz_ion() supplies the lattice (ion.a) used to build the legacy grid;
%   qVec_generator's own console diagnostics are captured and discarded
%   (evalc) so this function stays silent like every other invzt_ function.
%
%   See also INVZT_JQ_TENSOR, QVEC_GENERATOR.
if nargin < 2, conv = 'halfopen'; end

switch conv
    case 'halfopen'
        ax = (0:n-1)/n;
        [QX, QY, QZ] = ndgrid(ax);
        qvec = [QX(:), QY(:), QZ(:)];
        note = ['halfopen: ndgrid((0:n-1)/n) over 3 axes, n distinct points ' ...
                'per axis, no duplicate reciprocal-periodic face; Gamma point ' ...
                '[0 0 0] dropped (occurs exactly once on this grid). Weights ' ...
                'are uniform (1/nq) over the KEPT points only -- sum(w) == 1 ' ...
                'renormalizes away the dropped Gamma point''s measure rather ' ...
                'than carrying it as an explicit zero-weighted row.'];
    case 'legacy_inclusive'
        ion = invz_ion();
        qvec = evalc_call1(@() qVec_generator(ion.a, 'mode', 'grid', ...
            'grid', [n n n], 'range', [-0.5 0.5], 'verbose', false));
        note = ['legacy_inclusive: qVec_generator endpoint-inclusive linspace ' ...
                '(-0.5:0.5, spacing 1/(n-1)) -- the +-0.5 faces are DUPLICATE ' ...
                'reciprocal-periodic points, NOT deduplicated (an n^3 grid has ' ...
                'only (n-1)^3 distinct points; the redundant faces each still ' ...
                'get their own uniform weight below). Gamma point dropped where ' ...
                'it occurs exactly (all components ~0); weights are uniform ' ...
                '(1/nq) over the KEPT points only, renormalizing away the ' ...
                'dropped Gamma point''s measure (sum(w) == 1), same convention ' ...
                'as halfopen.'];
    otherwise
        error('invzt:qgridConv', ...
            'conv must be ''halfopen'' or ''legacy_inclusive''; got %s.', ...
            invzt_str(conv));
end

keep = any(abs(qvec) > 1e-12, 2);     % drop the Gamma row (every component ~0)
qvec = qvec(keep, :);
nq   = size(qvec, 1);

g.qvec = qvec;
g.w    = ones(nq, 1) / nq;
g.conv = conv;
g.note = note;
g.hash = qgrid_hash(n, conv_code(conv), qvec);
end

% ------------------------------- local helpers ------------------------------

function out = evalc_call1(fn)
%EVALC_CALL1 Run fn() (single output) with stdout captured/discarded.
[~, out] = evalc('fn()');
end

function c = conv_code(convstr)
%CONV_CODE Numeric tag folded into cache/hash keys so distinct grid
% conventions at the same q-mesh content never alias (mirrored in
% invzt_jq_tensor's own conv_code for its cache key; kept as independent
% per-file copies, matching the codebase's existing local-hash-helper
% convention -- see invz_cache_key.m header).
switch convstr
    case 'explicit',         c = 0;
    case 'halfopen',         c = 1;
    case 'legacy_inclusive', c = 2;
    otherwise
        error('invzt:conv', 'Unknown grid convention %s.', invzt_str(convstr));
end
end

function h = qgrid_hash(n, convcode, qvec)
% Same weak-hash formula as invz_cache_key/invz_jq_modes/invz_odd_blocks'
% local hash_vec, applied to [n; convcode; qvec(:)] so grid size, convention,
% and content all participate.
v = [n; convcode; qvec(:)];
h = sprintf('%dv_%08x', numel(v), typecast(single(sum(v.*(1:numel(v))')), 'uint32'));
end
