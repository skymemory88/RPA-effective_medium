function [dF, out] = invz_deltaF_ordered(ion, T, Bx, Jnu_flat, opts)
%INVZ_DELTAF_ORDERED PARTIAL HYBRID DIAGNOSTIC (plan SS7b): finite-domain Legendre line
% integral over the H_MF profile. NOT dF(m=0) -- the hybrid's saturation-normalized
% absolute free energy is outside its validated domain (route-A divergence recorded in
% the task-6 execution record); cross-state comparison is meaningful only at a COMMON
% cutoff (use opts.hmax_abs). The J 2.34 two-route identity is validated in the closed
% 2x2 model (test_invz_deltaF).
%   dF = -int_0^{hmax} (h0(h') - hmf(h')) d<Jz>'
% evaluated on the invz_hmf_ordered profile as -trapz(m, h0 - hgrid), truncated at the
% profile's cutoff (opts.hmax_fac, or the exact opts.hmax_abs override -- see
% invz_hmf_ordered). out.endpoint_dh = |dh(end)| is a NON-CLAIMING endpoint diagnostic
% (the magnitude of the un-integrated tail's driving term only) -- it is NEVER a bound
% or a saturation remainder, and NEVER describes a continuation to ion.J. out.hmax_used
% reports the actual cutoff used (NaN when the profile is empty or untrusted).
% Validation-only -- no production path depends on this function.
if nargin < 5, opts = struct(); end
o = opts;  o.hmax_fac = getf(opts, 'hmax_fac', 4);  o.nH = getf(opts, 'nH', 65);
[~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu_flat, o);
% round-4 P1-C: integrate ONLY trusted profiles -- an 'unresolved' or 'node_failed'
% profile can carry a nonempty grid and finite arrays, but its contract forbids use.
trusted = strcmp(prof.status, 'ok') && ~isempty(prof.hgrid) && all(prof.node_conv);
if ~trusted
    dF = NaN;
    out = struct('dF_partial', NaN, 'endpoint_dh', NaN, 'hmax_used', NaN, 'prof', prof);
    return;
end
dh = prof.h0 - prof.hgrid;
dF = -trapz(prof.m, dh);
out = struct('dF_partial', dF, 'endpoint_dh', abs(dh(end)), ...
             'hmax_used', max(prof.hgrid), 'prof', prof);
end
