function out = invzt_relaxed_longitudinal_response(chi, Jrelax, relaxed)
%INVZT_RELAXED_LONGITUDINAL_RESPONSE Longitudinal response on a relaxed-MF path.
%   OUT = INVZT_RELAXED_LONGITUDINAL_RESPONSE(CHI,JRELAX,RELAXED) eliminates
%   the transverse internal-field coordinates RELAXED (indices drawn from
%   [1 2]) while coordinate 3 is varied. With dm=CHI*dh and
%       dh_r = JRELAX*dm_r,
%   the path derivative is
%       dh_r/dh_3 = (I-JRELAX*CHI_rr)^-1 JRELAX*CHI_r3,
%       dm/dh_3   = CHI_:3 + CHI_:r dh_r/dh_3.
%   OUT.chi_eff is the third component of dm/dh_3. This is the Schur-
%   complement response required when the applied transverse field is fixed
%   but its molecular field relaxes.
if nargin < 3 || isempty(relaxed), relaxed = 1; end
if ~(isnumeric(chi) && isequal(size(chi),[3 3]) ...
        && all(isfinite(chi(:))))
    error('invzt:relaxedResponseChi','chi must be a finite 3x3 matrix.');
end
if ~(isnumeric(relaxed) && isvector(relaxed) ...
        && all(ismember(relaxed,[1 2])) && numel(unique(relaxed))==numel(relaxed))
    error('invzt:relaxedResponseIndices', ...
        'relaxed must contain unique indices drawn from [1 2].');
end
relaxed = relaxed(:).'; nr = numel(relaxed);
if isscalar(Jrelax), Jrelax=Jrelax*eye(nr); end
if ~(isnumeric(Jrelax) && isequal(size(Jrelax),[nr nr]) ...
        && all(isfinite(Jrelax(:))))
    error('invzt:relaxedResponseCoupling', ...
        'Jrelax must be a finite scalar or numel(relaxed)-square matrix.');
end

chi = real((chi+chi')/2);
A = eye(nr)-Jrelax*chi(relaxed,relaxed);
if rcond(A)<1e-12
    error('invzt:relaxedResponseSingular', ...
        'The transverse molecular-field relaxation matrix is singular.');
end
field_slope = A \ (Jrelax*chi(relaxed,3));
moment_slope = chi(:,3)+chi(:,relaxed)*field_slope;
residual = A*field_slope-Jrelax*chi(relaxed,3);

out = struct('schema','invzt_relaxed_longitudinal_response/v1', ...
    'relaxed_indices',relaxed,'Jrelax',Jrelax,'chi',chi, ...
    'relaxation_matrix',A,'relaxation_rcond',rcond(A), ...
    'field_slope',field_slope,'moment_slope',moment_slope, ...
    'chi_eff',moment_slope(3),'closure_residual_inf',max(abs(residual)), ...
    'unrelaxed_chi_cc',chi(3,3), ...
    'feedback_correction',moment_slope(3)-chi(3,3));
end
