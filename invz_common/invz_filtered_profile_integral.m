function [h0, out] = invz_filtered_profile_integral(h, r, eligible, anchor_r)
%INVZ_FILTERED_PROFILE_INTEGRAL Visual-only trapezoid over eligible finite nodes.
%   [h0,out] = invz_filtered_profile_integral(h,r,eligible,anchor_r)
%   removes every ineligible node or node with a nonfinite integrand, prepends
%   the finite visual anchor (0,anchor_r), and applies cumtrapz to the retained
%   sequence. h0 is returned on the original grid, with NaN at omitted nodes.
%   The lower missing interval is represented by one straight trapezoid from
%   the external endpoint to the first retained positive-h node.
%
%   If retained nodes were separated by omitted grid nodes, cumtrapz connects
%   their endpoints with one straight trapezoid. out records every such
%   bridge. This helper deliberately does not claim that the resulting
%   integral is a certified realization of equation (45).

if nargin < 4
    error('invz:filteredProfileIntegral', ...
        'h, r, eligible, and anchor_r are all required.');
end
if ~(isnumeric(h) && isreal(h) && isvector(h) && ~isempty(h) && ...
        all(isfinite(h(:))) && all(h(:) >= 0))
    error('invz:filteredProfileIntegral', ...
        'h must be a nonempty finite real nonnegative vector.');
end
if ~(isnumeric(r) && isreal(r) && isvector(r) && numel(r) == numel(h))
    error('invz:filteredProfileIntegral', ...
        'r must be a real vector with the same number of elements as h.');
end
if ~(islogical(eligible) || isnumeric(eligible)) || ...
        ~isreal(eligible) || ~isvector(eligible) || numel(eligible) ~= numel(h) || ...
        any(~isfinite(eligible(:)))
    error('invz:filteredProfileIntegral', ...
        'eligible must be a finite real/logical vector matching h.');
end
if ~(isnumeric(anchor_r) && isreal(anchor_r) && isscalar(anchor_r) && ...
        isfinite(anchor_r))
    error('invz:filteredProfileIntegral', ...
        'anchor_r must be a finite real scalar.');
end

was_row = isrow(h);
hcol = h(:);
rcol = r(:);
ecol = logical(eligible(:));
if any(diff(hcol) <= 0)
    error('invz:filteredProfileIntegral', 'h must be strictly increasing.');
end

used = ecol & isfinite(rcol);
idx = find(used);
h0col = nan(size(hcol));
if isempty(idx)
    first_h = NaN;
    panel_widths = zeros(0,1);
    missing_between = zeros(0,1);
    bridge_mask = false(0,1);
else
    hkeep = hcol(idx);
    rkeep = rcol(idx);
    anchored = cumtrapz([0;hkeep],[anchor_r;rkeep]);
    h0col(idx) = anchored(2:end);
    first_h = hkeep(1);
    panel_widths = diff(hkeep);
    missing_between = diff(idx)-1;
    bridge_mask = missing_between > 0;
end

bridge_widths = panel_widths(bridge_mask);
if isempty(bridge_widths)
    max_bridge_width = 0;
else
    max_bridge_width = max(bridge_widths);
end

out = struct();
out.used_mask = reshape(used, size(h));
out.used_indices = idx;
out.node_count = numel(idx);
out.omitted_node_count = numel(hcol)-numel(idx);
out.ineligible_node_count = nnz(~ecol);
out.nonfinite_eligible_node_count = nnz(ecol & ~isfinite(rcol));
out.first_retained_h = first_h;
out.lower_omitted_width = 0;
out.lower_bridge_width = first_h;
out.anchor = "external_h0_r";
out.anchor_r = anchor_r;
out.bridge_panel_mask = bridge_mask;
out.bridge_missing_node_count = missing_between(bridge_mask);
out.bridge_widths = bridge_widths;
out.bridge_panel_count = nnz(bridge_mask);
out.bridged_missing_node_count = sum(missing_between(bridge_mask));
out.bridged_width_total = sum(bridge_widths);
out.max_bridge_width = max_bridge_width;
out.visual_only = true;

if was_row
    h0 = h0col.';
else
    h0 = h0col;
end
end
