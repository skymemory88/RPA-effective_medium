function O = invzp_area_rule_oracle(prof, J0eff)
%INVZP_AREA_RULE_ORACLE Compare the two Jensen-path free-energy constructions.
%
% The production profile stores
%   F_direct(h) = h0(h) - J0eff*m(h),  h0 = integral r(h') dh'.
% Since dm/dh = -G0bare and crit = r + J0eff*G0bare, the same continuous
% quantity obeys
%   F_area(h) = integral crit(h') dh'.
%
% This diagnostic uses the profile's own nodes and trapezoidal rule for both
% routes.  Their difference therefore measures finite-grid quadrature and
% state consistency; it is not a replacement for a missing/nonconverged node.
required = {'hgrid','h0','m','F','crit','slope0'};
if ~isstruct(prof) || ~isscalar(prof) || ~all(isfield(prof, required))
    error('invzp:areaRule', ...
        'prof must be a scalar invz_hmf_ordered profile with fields %s.', ...
        strjoin(required, ', '));
end
if ~(isscalar(J0eff) && isfinite(J0eff))
    error('invzp:areaRule', 'J0eff must be a finite scalar.');
end

h = reshape(prof.hgrid, 1, []);
h0 = reshape(prof.h0, 1, []);
m = reshape(prof.m, 1, []);
Fd = reshape(prof.F, 1, []);
crit = reshape(prof.crit, 1, []);
if any([numel(h0),numel(m),numel(Fd),numel(crit)] ~= numel(h)) || ...
        any(~isfinite(h)) || any(diff(h) <= 0)
    error('invzp:areaRule', ...
        'Profile arrays must match a finite, strictly increasing hgrid.');
end

Fa = cumtrapz([0 h], [prof.slope0 crit]);
Fa = Fa(2:end);
identity_resid = Fd-(h0-J0eff*m);
route_diff = Fd-Fa;
subtraction_scale = abs(h0)+abs(J0eff*m);
condition_direct = subtraction_scale./abs(Fd);
condition_direct(Fd == 0 & subtraction_scale == 0) = 1;

O = struct('schema', 'invzp_area_rule_oracle/v1', ...
    'hgrid', h, 'F_direct', Fd, 'F_area', Fa, ...
    'route_diff', route_diff, ...
    'direct_identity_resid', identity_resid, ...
    'subtraction_scale', subtraction_scale, ...
    'condition_direct', condition_direct, ...
    'max_abs_route_diff', finite_max(abs(route_diff)), ...
    'max_scaled_route_diff', ...
        finite_max(abs(route_diff)./max(subtraction_scale, eps)), ...
    'max_abs_direct_identity_resid', finite_max(abs(identity_resid)), ...
    'direct_zero_linear', last_increasing_zero(h, Fd), ...
    'area_zero_linear', last_increasing_zero(h, Fa), ...
    'profile_status', getf(prof, 'status', 'unknown'), ...
    'all_nodes_accepted', isfield(prof, 'node_conv') && ...
                          all(logical(prof.node_conv)));
end

function x0 = last_increasing_zero(x, y)
idx = find(isfinite(y(1:end-1)) & isfinite(y(2:end)) & ...
           y(1:end-1) < 0 & y(2:end) >= 0, 1, 'last');
if isempty(idx)
    x0 = NaN;
else
    x0 = x(idx)-y(idx)*(x(idx+1)-x(idx))/(y(idx+1)-y(idx));
end
end

function v = finite_max(x)
x = x(isfinite(x));
if isempty(x), v = NaN; else, v = max(x); end
end

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
