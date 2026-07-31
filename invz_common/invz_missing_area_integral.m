function [h0, out] = invz_missing_area_integral(h, r, certified, area)
%INVZ_MISSING_AREA_INTEGRAL Integrate one contiguous certified high-h component.
%   [h0,out] = invz_missing_area_integral(h,r,certified,area) selects the
%   terminal contiguous block for which CERTIFIED is true and R is finite.
%   If h_e is that block's lower edge, the returned values are
%
%       h0(h) = area + integral_{h_e}^{h} r(s) ds.
%
%   Values below h_e are NaN.  They are represented only by the explicit
%   nonnegative AREA; no rejected/nonfinite node is bridged or assigned a
%   zero contribution.  This helper performs quadrature only.  It does not
%   select a thermodynamic branch or claim that AREA is physically derived.

if nargin < 4
    error('invz:missingAreaIntegral', ...
        'h, r, certified, and area are all required.');
end
if ~(isnumeric(h) && isreal(h) && isvector(h) && ~isempty(h) && ...
        all(isfinite(h(:))) && all(h(:) > 0))
    error('invz:missingAreaIntegral', ...
        'h must be a nonempty finite real positive vector.');
end
if ~(isnumeric(r) && isreal(r) && isvector(r) && numel(r) == numel(h))
    error('invz:missingAreaIntegral', ...
        'r must be a real vector with the same number of elements as h.');
end
if ~(islogical(certified) || isnumeric(certified)) || ...
        ~isreal(certified) || ~isvector(certified) || ...
        numel(certified) ~= numel(h) || any(~isfinite(certified(:)))
    error('invz:missingAreaIntegral', ...
        'certified must be a finite real/logical vector matching h.');
end
if ~(isnumeric(area) && isreal(area) && isscalar(area) && ...
        isfinite(area) && area >= 0)
    error('invz:missingAreaIntegral', ...
        'area must be a finite real nonnegative scalar.');
end

was_row = isrow(h);
hcol = h(:);
rcol = r(:);
eligible = logical(certified(:)) & isfinite(rcol);
if any(diff(hcol) <= 0)
    error('invz:missingAreaIntegral', 'h must be strictly increasing.');
end

last_gap = find(~eligible, 1, 'last');
if isempty(last_gap)
    edge_index = 1;
else
    edge_index = last_gap + 1;
end
selected = false(size(hcol));
if edge_index <= numel(hcol)
    selected(edge_index:end) = true;
end

h0col = nan(size(hcol));
if any(selected)
    idx = find(selected);
    hkeep = hcol(idx);
    rkeep = rcol(idx);
    h0col(idx) = area + cumtrapz(hkeep, rkeep);
    edge_h = hkeep(1);
    edge_r = rkeep(1);
else
    idx = zeros(0,1);
    edge_index = NaN;
    edge_h = NaN;
    edge_r = NaN;
end

out = struct('selected_mask',reshape(selected,size(h)), ...
    'selected_indices',idx,'node_count',numel(idx), ...
    'edge_index',edge_index,'component_edge',edge_h,'edge_r',edge_r, ...
    'missing_area',area,'quadrature_rule',"trapezoid", ...
    'branch_selector',"external", ...
    'interior_bridge_count',0,'uses_uncertified_nodes',false, ...
    'certified_below_edge_count',nnz(eligible & ~selected), ...
    'unresolved_below_edge_count',nnz(~eligible & ~selected), ...
    'approximation_only',true);

if was_row
    h0 = h0col.';
else
    h0 = h0col;
end
end
