function [sm, eopts, eso] = invz_check_static_medium(opts, eopts, eso)
%INVZ_CHECK_STATIC_MEDIUM Sole public authority for the omega_n = 0 static-medium scheme
% (spec SS4.2), following the shared-validator idiom of invz_check_coupling_opts.m /
% invz_check_solve_opts.m / invz_common/invz_check_transverse_mf.m.
%
%   opts.static_medium : 'resummed'             (DEFAULT -- legacy, bit-identical)
%                        'strict_1z_dyson_ref'  (selected strict candidate, spec SS0.3)
%                        'strict_1z_bare_ref'   (systematic comparator)
%   opts.ref_margin    : reference-denominator floor (default 1e-6, FROZEN in
%                        docs/invzp_strict_medium_prereg.md SS2)
%
% Returns sm = struct('scheme','is_strict','ref_margin') and STAMPS the resolved scheme plus
% ref_margin into BOTH internal leg option structs -- eopts (invz_emt_scalar, the PM static
% slot) and eso (invz_emt_static_ordered, the ordered static sector). Stamping from one
% authority is what makes it impossible for the two legs to run different truncation orders,
% which spec SS0.2 requires: a formally O(1/z^2) PM/FM mismatch is dangerous at a continuous
% boundary, where the target mass is exactly zero.
%
% Setting the scheme directly on one leg (opts.emt.static_medium or
% opts.emt_static.static_medium) is a CONFLICT, not an override: it would split the sectors.
% 'strict_1z_selfconsistent' is deliberately NOT selectable -- the rejected self-consistent
% quadratic needs different inputs (G0, D, a branch choice) and lives in a separately named
% diagnostic comparator (spec SSB).
% This validator does not infer coupling shape from option names. The PM leaf supports an
% [nJ,nw] matrix by using column 1 for the strict static slot. The ordered public solver,
% after resolving ODD/retardation, rejects a non-vector Jnu_flat under strict mode rather
% than allowing invz_emt_static_ordered.m:43 to flatten all columns.
if nargin < 1 || isempty(opts), opts = struct(); end
if nargin < 2 || isempty(eopts), eopts = struct(); end
if nargin < 3 || isempty(eso),   eso   = struct(); end

valid = {'resummed', 'strict_1z_dyson_ref', 'strict_1z_bare_ref'};
scheme = getf(opts, 'static_medium', 'resummed');
if ~(ischar(scheme) || isstring(scheme)) || ~any(strcmp(char(scheme), valid))
    error('invz:staticMedium', ['opts.static_medium must be one of %s; got %s. ' ...
        '(''strict_1z_selfconsistent'' is a diagnostic comparator, not a selectable ' ...
        'production scheme -- see the spec SSB.)'], strjoin(valid, ' | '), ...
        local_describe(scheme));
end
scheme = char(scheme);

% A per-leg value that AGREES with the resolved scheme is this function's own stamp being
% re-read, so validation is IDEMPOTENT: invz_solve_point_ordered forwards its full numerical
% context (opts.emt / opts.emt_static included) into invz_hmf_ordered, which validates again.
% A value that DISAGREES is a genuine conflict -- it would split the two sectors across
% different truncation orders, which spec SS0.2 forbids.
for f = {'emt', 'emt_static'}
    if isfield(opts, f{1}) && isstruct(opts.(f{1})) && isfield(opts.(f{1}), 'static_medium')
        inner = opts.(f{1}).static_medium;
        if ~(ischar(inner) || isstring(inner)) || ~strcmp(char(inner), scheme)
            error('invz:staticMedium', ['opts.%s.static_medium = %s conflicts with the ' ...
                'resolved scheme ''%s''. Set the scheme ONCE via opts.static_medium and let ' ...
                'this function stamp both legs (spec SS4.2); a per-leg override would split ' ...
                'the sectors across different truncation orders.'], f{1}, ...
                local_describe(inner), scheme);
        end
    end
end

% Also validate the explicit internal structs passed as arguments; otherwise a disagreeing
% eopts/eso value would simply be overwritten and the conflict hidden.
for pair = {{'eopts', eopts}, {'eso', eso}}
    label = pair{1}{1}; innerStruct = pair{1}{2};
    if isfield(innerStruct, 'static_medium')
        inner = innerStruct.static_medium;
        if ~(ischar(inner) || isstring(inner)) || ~strcmp(char(inner), scheme)
            error('invz:staticMedium', '%s.static_medium conflicts with resolved scheme ''%s''.', ...
                  label, scheme);
        end
    end
end

is_strict = ~strcmp(scheme, 'resummed');
ref_margin = getf(opts, 'ref_margin', 1e-6);
if ~(isscalar(ref_margin) && isfinite(ref_margin) && ref_margin > 0)
    error('invz:staticMedium', 'opts.ref_margin must be a positive finite scalar.');
end
sm = struct('scheme', scheme, 'is_strict', is_strict, 'ref_margin', ref_margin);
eopts.static_medium = sm.scheme;   eopts.ref_margin = sm.ref_margin;
eso.static_medium   = sm.scheme;   eso.ref_margin   = sm.ref_margin;
end

% ---------------------------------------------------------------------------------------------
function s = local_describe(v)
%LOCAL_DESCRIBE Readable rendering of a rejected value (mirrors invz_check_coupling_opts.m).
if ischar(v) || isstring(v)
    s = sprintf('''%s''', char(v));
elseif isnumeric(v) && isscalar(v)
    s = sprintf('%g (%s)', v, class(v));
else
    s = sprintf('a %s', class(v));
end
end
