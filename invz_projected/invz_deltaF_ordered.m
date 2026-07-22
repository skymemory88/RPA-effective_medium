function [dF, out] = invz_deltaF_ordered(ion, T, Bx, Jnu_flat, opts)
%INVZ_DELTAF_ORDERED Free-energy correction dF(m=0), route A (framework SS9.4, J 2.34):
%   dF = -int_0^{M0} (h0(m') - hmf(m')) d<Jz>'
% evaluated on the invz_hmf_ordered profile as -trapz(m, h0 - hgrid), with the profile
% extended toward saturation (hmax_fac default 4: fluctuations quench as the splitting
% grows, dh -> 0). out.tail_est = |dh(end)|*(M0 - m(end)) estimates the cut tail and must
% be reported alongside dF whenever dF is quoted. Validation-only -- no production path
% depends on this function.
if nargin < 5, opts = struct(); end
o = opts;  o.hmax_fac = getf(opts, 'hmax_fac', 4);  o.nH = getf(opts, 'nH', 65);
[~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu_flat, o);
% round-4 P1-C: integrate ONLY trusted profiles -- an 'unresolved' or 'node_failed'
% profile can carry a nonempty grid and finite arrays, but its contract forbids use.
trusted = strcmp(prof.status, 'ok') && ~isempty(prof.hgrid) && all(prof.node_conv);
if ~trusted
    dF = NaN;  out = struct('dF_routeA', NaN, 'tail_est', NaN, 'prof', prof);  return;
end
dh = prof.h0 - prof.hgrid;
dF = -trapz(prof.m, dh);
M0 = ion.J;                                          % max <Jz> = J exactly
out = struct('dF_routeA', dF, 'tail_est', abs(dh(end))*(M0 - prof.m(end)), 'prof', prof);
end
