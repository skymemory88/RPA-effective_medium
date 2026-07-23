function c = invz_phase1_couplings(ion, g, dpRng)
%INVZ_PHASE1_COUPLINGS Phase-1 coupling evaluator: invz_jq_modes over a Phase-1 grid (stage-2c
% coupling-only audit, ADDITIVE). Calls invz_jq_modes -- the SAME production coupling function
% invz_bz_couplings.m calls -- and does NOT reimplement any physics (no MF_dipole/exchange logic
% here). NO field argument anywhere: Phase 1 is field-independent (prereg "Scope").
%
% g        struct from invz_phase1_qgrid (needs g.qvec, g.w).
% dpRng    real-space dipole cutoff (30/40/50 in the frozen ladder).
%
% OUTPUT c (struct):
%   c.Jnu_unflat  [Npts x 4] sorted eigenvalue branches, invz_jq_modes' own first output,
%                 unmodified.
%   c.Jnu_flat    Jnu_unflat(:) -- column-major flatten, matching invz_bz_couplings.m:17's own
%                 convention exactly -- the "flat coupling multiset" item 5 summarizes.
%   c.w_flat      g.w repeated once per branch (repmat(g.w,4,1)), row-aligned with c.Jnu_flat
%                 under the SAME column-major convention (branch b's Npts entries all carry that
%                 q-point's own weight). Exported for provenance; NOT used by invz_phase1_checks'
%                 item-5 summaries, which are UNWEIGHTED statistics of c.Jnu_flat -- matching
%                 invz_emt_scalar.m/invz_emt_static_ordered.m's own unweighted `mean` over
%                 Jnu_flat with no weight argument (prereg "Weights").
%   c.info        invz_jq_modes' own info struct (Jcc0, Jaa0, dpRng, ...), unmodified.
%   c.J0eff       info.Jcc0 -- production's driver-facing name for this scalar (invz_hmf_ordered.m
%                 / invz_task2_couplings_shifted_grid.m opts.J0eff); LITERALLY the same number as
%                 c.Jcc0 below, exported under both names because prereg item 5 lists both
%                 "J0eff, Jcc0" as separate named energy scalars.
%   c.Jcc0        info.Jcc0 under its invz_jq_modes-native field name.
%   c.maxJnu      max(c.Jnu_flat) -- the largest branch value anywhere on THIS finite grid. Equals
%                 Jcc0 when Gamma is present (the analytic maximum sits exactly at Gamma); may be
%                 strictly less than Jcc0 under P_drop, where Gamma has been removed.
if nargin < 3 || isempty(dpRng), dpRng = 30; end
[Jnu_unflat, info] = invz_jq_modes(ion, g.qvec, struct('dpRng', dpRng, 'cache', true));
c.Jnu_unflat = Jnu_unflat;
c.Jnu_flat   = Jnu_unflat(:);
c.w_flat     = repmat(g.w, 4, 1);
c.info       = info;
c.J0eff      = info.Jcc0;
c.Jcc0       = info.Jcc0;
c.maxJnu     = max(c.Jnu_flat);
end
