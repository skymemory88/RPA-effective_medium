function [Jnu, info, Jaa0] = invz_anchor_couplings(ion, opts)
%INVZ_ANCHOR_COUPLINGS Thin, presence-preserving Gate-E coupling-injection helper (Ewald Step-5
% Task 5; docs/invzp_ewald_prereg.md Sec.7 "Frozen Gate-E -- physics anchors", docs/invzp_ewald_
% design.md item 6 "Physics-anchor wrapper API", docs/invzp_ewald_integration_map.md Sec.5.4(B),
% docs/superpowers/plans/2026-07-24-ewald-step5-integration.md Task 5).
%
% A pure pass-through wrapper over invz_bz_couplings(ion, opts2) -- the SAME [Jnu, info, Jaa0]
% outputs, computed nowhere else. This file computes NO physics of its own: no Sigma_c, Tc, Bc,
% Jensen field, or anchor inequality -- those belong to the Step-7 Gate-E anchor-wrapper test
% files that will call this helper; this task only builds and parity-tests the injection point.
%
% ONLY two defaults exist, applied ONLY when the caller's opts does not already have the field:
%   opts.dpRng (30)    real-space dipole cutoff (bruteforce backend only; ignored by Ewald).
%   opts.cache (true)  invz_jq_modes file cache.
% Every other field -- opts.grid, opts.dipole, opts.ewald, opts.gridConvention, opts.gridOffset,
% opts.gammaPolicy -- is forwarded to invz_bz_couplings completely UNCHANGED, BY PRESENCE: a
% field this caller's opts does not have stays absent from the forwarded call too (never
% synthesized, never given a value here, regardless of what invz_bz_couplings' own default for
% that field happens to be). In particular, a caller who supplies none of the three grid-policy
% fields and no backend field reproduces invz_bz_couplings' own unmodified legacy
% qVec_generator+Gamma-drop route exactly (info.grid absent, 'bruteforce' backend) -- i.e. is
% legacy-identical by construction (prereg Sec.7 / design item 6: "the helper omits optional
% grid/backend fields when absent, so its default path is genuinely legacy-identical").
%
% See invz_bz_couplings.m for the full opts.grid/dpRng/cache/dipole/ewald/gridConvention/
% gridOffset/gammaPolicy contract, accepted values, and info.grid/info.dipole provenance shape --
% this helper does not repeat, re-validate, or re-default any of it beyond dpRng/cache above.
if nargin < 2, opts = struct(); end
opts2 = opts;                                      % exact copy: every field already present in
                                                    % opts (grid, dipole, ewald, gridConvention,
                                                    % gridOffset, gammaPolicy, ...) rides through
                                                    % untouched, unrenamed, un-re-defaulted.
if ~isfield(opts2, 'dpRng'), opts2.dpRng = 30;   end
if ~isfield(opts2, 'cache'), opts2.cache = true; end
[Jnu, info, Jaa0] = invz_bz_couplings(ion, opts2);
end
