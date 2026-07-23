function [Jnu_flat, J0eff, Jxx0, qc, Jnu_unflat, latinfo, n_gamma_dropped] = ...
    invz_task2_couplings_shifted_grid(ion, grid, dpRng, shift_mode)
%INVZ_TASK2_COUPLINGS_SHIFTED_GRID Real-lattice couplings on an (optionally) shifted q-grid
% (stage-2c task 2b-driver; prereg docs/invzp_task2_prereg.md Sec. D / brief G6 grid-offset
% ladder). ADDITIVE, READ-ONLY w.r.t. production: calls qVec_generator/invz_jq_modes (the SAME
% two helpers invz_bz_couplings.m itself calls) but does NOT modify invz_bz_couplings.m.
%
% ============================================================================================
% THE DELICATE PART -- read before changing this file (brief: "if you cannot implement a
% correct, tested dipolar-Gamma offset cleanly, STOP" -- this is that implementation):
%
% invz_bz_couplings.m:15 drops the Gamma point from its q-list via a HARDCODED zero check,
% `qc = qc(any(abs(qc) > 1e-12, 2), :)`. For the STANDARD unshifted grid
% (qVec_generator(..., 'range', [-0.5 0.5]), grid=[16 16 16]) this is a NO-OP: linspace(-0.5,
% 0.5, 16) has 15 intervals (an ODD count), so NO axis value is ever exactly 0 (nearest values
% are +/- half a step, ~0.0333) -- confirmed numerically here and independently in this task's
% report; the real lattice's own nJ = 16384 = 4*16^3 (Task-2a report) already proves no point
% was ever dropped. So "drop q==0" and "drop the Gamma-equivalent point(s)" happen to COINCIDE
% for THIS specific grid, and the hardcoded check has never actually fired in production.
%
% A "naive q-translation" (the brief's own words) breaks this coincidence in exactly the way
% prereg Sec. D warns about: shifting the WHOLE grid by half of its own axis step (the
% "half-step body-centred" offset below) moves grid index i=7 (0-based) on EVERY axis from
% -1/30 to EXACTLY 0 -- so the shifted grid DOES contain the literal Gamma point (verified
% numerically: exactly 1 of the 4096 shifted points is (0,0,0), see
% test_invz_task2_matrix.m's test_shifted_grid_half_step_drops_exactly_gamma). If the
% shift were applied AFTER a hardcoded "drop q==0" filter that ran on the UNSHIFTED grid (i.e.
% naively reusing invz_bz_couplings.m's own filter-then-nothing order), the filter would be a
% no-op (nothing to drop pre-shift) and the reintroduced Gamma point would NEVER be excluded
% post-shift -- silently leaving an extra, redundant Dq entry (Gamma's own top branch equals
% J0eff exactly, i.e. Dq(Gamma) = D_uni identically) in the flattened fluctuation-mode array,
% which would also silently perturb the finite-lattice EMT average (invz_emt_static_ordered's
% mean(Gq)/mean(Jf.*Gq) closure) by an O(1/nq) amount at exactly the largest-coupling point --
% precisely the kind of artifact that could produce a misleading Sec. D mesh-convergence
% verdict (a "shifted grid disagrees with unshifted" result that is actually just an artifact
% of inconsistent Gamma bookkeeping, not a genuine physical/lattice-convergence signal).
%
% This function avoids that bug by construction: it (1) shifts qc FIRST, then (2) filters
% using invz_is_gamma_equiv -- the SAME general Gamma-equivalence test invz_jq_modes.m:85
% itself uses to decide Lorentz-cavity placement -- evaluated on the SHIFTED q-values, never
% the pre-shift ones. This is provably general (not merely "happens to work"): for THIS ion's
% 4-site tau basis, the Gamma-equivalence conditions (tau*q' all-integer) reduce, within any
% q-box no larger than roughly [-0.5,0.5]^3 (worked out by hand in this task's report and
% confirmed numerically over the full 3D grid), to the single solution q=(0,0,0) -- so any
% shift that reintroduces a Gamma-equivalent point reintroduces EXACTLY (0,0,0), and the
% invz_is_gamma_equiv filter catches it unconditionally, independent of which grid index it
% lands on. invz_jq_modes' OWN per-point Lorentz-term addition (its internal
% invz_is_gamma_equiv check, unconditional on any pre-filtering) is unaffected by any of this
% either way -- the fix here is specifically about keeping the OUTER "exclude Gamma from the
% flattened fluctuation array" convention consistent (general, post-shift) across every rung
% of the offset ladder, not about the Lorentz term itself (which invz_jq_modes always gets
% right, regardless of what the caller pre-filters).
% ============================================================================================
%
% ion       invz_ion()-shaped struct (needs .a, .tau -- SAME fields invz_jq_modes reads).
% grid      [1x3] BZ grid size, e.g. [16 16 16] (SAME convention invz_bz_couplings uses).
% dpRng     real-space dipole cutoff (SAME convention invz_bz_couplings uses).
% shift_mode 'unshifted' (delta = 0; reproduces invz_bz_couplings bit-for-bit, see
%             test_invz_task2_matrix.m's test_shifted_grid_unshifted_matches_bz_couplings) |
%            'half_step' (delta = HALF of the existing per-axis linspace step, added
%             UNIFORMLY to h, k, AND l -- a body-centred offset of the whole grid: every new
%             point sits at the geometric centre of a cube of 8 neighbouring unshifted-grid
%             points, in direct analogy to a simple-cubic -> body-centred-cubic secondary
%             lattice. Since grid/range are identical and isotropic across h,k,l here, the
%             per-axis step is identical on all three axes, so ONE scalar delta correctly
%             implements a per-axis half-step shift on each of them.)
%
% Returns the SAME shape invz_bz_couplings/invz_jq_modes jointly produce: Jnu_flat (column-
% major flatten of Jnu_unflat, matching invz_bz_couplings.m:17's OWN convention exactly, so
% invz_ordered_trace_resolve's flat-index inversion formula applies unchanged), J0eff (=
% latinfo.Jcc0), Jxx0 (= latinfo.Jaa0), qc ([nq x 3], POST-shift, POST-Gamma-filter -- row-
% aligned with Jnu_unflat, the q/branch provenance this task's driver needs since
% invz_task2_run_config's synthetic-Jnu_flat injection path cannot carry qc/Jnu_unflat through
% (see invz_task2_matrix.m's header)), Jnu_unflat ([nq x 4]), latinfo (invz_jq_modes' own info
% struct), n_gamma_dropped (0 or 1 for a single uniform shift of this grid; report-only, so a
% human/2b-report can see whether the general filter actually fired on this particular rung,
% or was a no-op like the unshifted baseline).
if nargin < 3 || isempty(dpRng), dpRng = 30; end
if nargin < 4 || isempty(shift_mode), shift_mode = 'unshifted'; end

range = [-0.5 0.5];
qc0 = qVec_generator(ion.a, 'mode', 'grid', 'grid', grid, 'range', range, 'verbose', false);

switch shift_mode
    case 'unshifted'
        delta = 0;
    case 'half_step'
        step  = (range(2) - range(1)) / (grid(1) - 1);   % SAME step invz_bz_couplings' own
                                                           % grid uses (grid assumed isotropic
                                                           % across axes, as invz_bz_couplings
                                                           % itself assumes via one scalar range).
        delta = step / 2;
    otherwise
        error('invz:task2Config', 'invz_task2_couplings_shifted_grid: unknown shift_mode ''%s'' (expected ''unshifted'' or ''half_step'').', shift_mode);
end
qc_shifted = qc0 + delta;   % uniform body-centred shift of ALL three reduced components

% --- explicit, tested Gamma handling (the delicate step -- see header) --------------------
is_g = false(size(qc_shifted, 1), 1);
for i = 1:size(qc_shifted, 1)
    is_g(i) = invz_is_gamma_equiv(qc_shifted(i, :), ion.tau);
end
n_gamma_dropped = nnz(is_g);
qc = qc_shifted(~is_g, :);

[Jnu_unflat, latinfo] = invz_jq_modes(ion, qc, struct('dpRng', dpRng, 'cache', true));
Jnu_flat = Jnu_unflat(:);
J0eff = latinfo.Jcc0;
Jxx0  = latinfo.Jaa0;
end
