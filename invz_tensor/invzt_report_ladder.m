function txt = invzt_report_ladder(out)
%INVZT_REPORT_LADDER  Serialize a completed A4 ladder run for ODD-LOG SS-A4 pasting.
%
%   txt = INVZT_REPORT_LADDER(out) formats the struct returned by INVZT_RUN_LADDER into a
%   plain-text report block (per-rung table + budget verdict + monitors + cross-validation
%   comparators) and RETURNS it as a char array. With no output argument it also prints to
%   stdout. It WRITES NO FILES: the controller pastes the returned text into
%   docs/ODD-LOG.md (Section A4) -- this serializer never touches ODD-LOG or any tracked
%   file (Task-13 data-only discipline).
%
%   The table columns:
%     rung / N / dim   basis label, vertex Hilbert dimension, actual (multiplet-complete) dim
%     mc               multiplet_complete flag
%     crit_shift_odd   crit(a3-full, odd on) - crit(odd off); + = ODD lowers Tc (projected
%                      DeltaTc>0 direction)
%     rf / rd          full-A3 / A1 and dominant-dress / A1 odd-shift ratios (rf = beyond-E1
%                      transverse-spectator dressing; rd = matched E1 truncation)
%     collapse         |rd-1|/|rf-1| (transverse-dressing collapse at the ladder point)
%     eps_el/eps_cr    elastic-sector control; cross-Cartesian leakage (constraints 7/9)
%     sumrl/resum      sum-rule residual; Dyson-vs-additive crit spread (constraints 6/8)
%     vdef             chi0 virtual-completeness deficit vs the full-136 chi0 (DIAGNOSTIC)
%     conv / t[s]      converged flag; wall time
%
%   See also INVZT_RUN_LADDER.
if nargin < 1 || ~isstruct(out), error('invzt:report', 'invzt_report_ladder needs the invzt_run_ladder struct.'); end
m = out.meta;
nr = numel(out.rungs);
L = {};
p = @(varargin) local_push(varargin{:});

L = p(L, '=== A4 basis-defined state-space ladder (invzt_run_ladder) ===');
L = p(L, 'point: T = %.3f K, B = [%.3g %.3g %.3g] T | grid %d^3 %s, dpRng %d, Ecut %d meV | production=%d', ...
    m.T, m.B(1), m.B(2), m.B(3), m.ngrid, m.conv, m.dpRng, m.Ecut, m.production);
L = p(L, 'lattice: Jcc0 = %.6g meV, Jaa0 = %.6g meV, qhash %s | git %s | %s', ...
    m.Jcc0, m.Jaa0, m.qhash, m.git, m.date);
L = p(L, 'full-136 chi0(T,B): cc = %.4f, perp = %.4f meV^-1 (virtual-deficit reference)', ...
    m.chi0_full136_ccdiag, m.chi0_full136_perpdiag);
L = p(L, 'budget: %g h/solve; %d rungs run, %d refused.', m.budget_hours, m.nrun, m.nskip);
L = p(L, '');

% --- per-rung table --------------------------------------------------------------------
hdr = sprintf('%-8s %4s %4s %3s %14s %7s %7s %8s %7s %8s %7s %9s %7s %6s %6s', ...
    'rung', 'N', 'dim', 'mc', 'crit_shift_odd', 'rf', 'rd', 'collapse', ...
    'eps_el', 'eps_cr', 'sumrl', 'resum', 'vdef', 'conv', 't[s]');
L = p(L, '%s', hdr);
L = p(L, '%s', repmat('-', 1, numel(hdr)));
for i = 1:nr
    L = p(L, '%-8s %4d %4d %3d %+14.5e %7.3f %7.3f %8.3f %7.3f %8.4f %7.3f %+9.2e %7.3f %6d %6.0f', ...
        out.rungs{i}, out.N(i), out.dim_actual(i), out.multiplet_complete(i), ...
        out.crit_shift_odd(i), out.rf(i), out.rd(i), out.matched_collapse(i), ...
        out.eps_el(i), out.eps_cross(i), out.sumrule_rel(i), out.resum_spread(i), ...
        out.chi0_virtual_deficit(i), out.converged(i), out.t_sec(i));
end
L = p(L, '');

% --- raw crits + A1 odd shift (context for the ratios) ---------------------------------
L = p(L, 'raw crits (a3-full): odd-off / odd-on, and the A1 odd shift d1 (ratio denominator):');
for i = 1:nr
    L = p(L, '  %-8s crit_off=%+.5e crit_on=%+.5e | d1(A1)=%+.5e | npaths=%d proj=%.3g h', ...
        out.rungs{i}, out.crit_oddoff(i), out.crit_oddon(i), out.a1_odd_shift(i), ...
        out.npaths(i), out.projected_hours(i));
end
L = p(L, '');

% --- small-Bx-proxy Tc (if computed) ---------------------------------------------------
if m.tc
    L = p(L, 'small-Bx-proxy Tc (Bx = %.3g T; tensor A3 true-B=0 Tc DEFERRED -- proxy for EVERY rung):', m.tc_bx);
    for i = 1:nr
        if isfinite(out.tc_proxy(i))
            L = p(L, '  %-8s Tc_proxy = %.4f K  (window T = [%.4f %.4f]) vs proj 1.509 / no-ODD 1.743 / exp 1.53 K', ...
                out.rungs{i}, out.tc_proxy(i), out.tc_used(i,1), out.tc_used(i,2));
        else
            L = p(L, '  %-8s Tc_proxy = NaN (no PM extrapolation window on the grid)', out.rungs{i});
        end
    end
    L = p(L, '');
end

% --- budget-refusal verdict ------------------------------------------------------------
L = p(L, 'budget verdict (LOCKED T11 dense-vertex scaling; refuse projected > %g h):', m.budget_hours);
if isempty(out.skipped_rungs)
    L = p(L, '  no rung refused.');
else
    for i = 1:numel(out.skipped_rungs)
        s = out.skipped_rungs(i);
        L = p(L, '  REFUSE %-8s N=%-3d dim=%-3d : %s', s.name, s.N, s.dim, s.reason);
    end
end
L = p(L, '');

% --- beyond-E1 (beyond-Gaussian) share per rung, vs the projected Tier-2 share ---------
L = p(L, 'beyond-E1 (beyond-Gaussian) dressing share = rf - 1 per rung (vs projected Tier-2 ~%.1f%%):', ...
    100*m.proj_Tier2_share);
for i = 1:nr
    L = p(L, '  %-8s rf-1 = %+.3f  (%+.1f%%)   [rd-1 matched = %+.3f]', ...
        out.rungs{i}, out.rf(i) - 1, 100*(out.rf(i) - 1), out.rd(i) - 1);
end
L = p(L, '');

% --- cross-validation comparators + regime note (REPORT, never tune) -------------------
L = p(L, 'cross-validation (REPORT): projected Tier-1+2 Tc = %.3f K, DeltaTc = %.4f K, Tier-2 share ~%.1f%%.', ...
    m.proj_Tier1p2_K, m.proj_DeltaTc_K, 100*m.proj_Tier2_share);
pm = nr > 0 && all(out.crit_oddoff > 0);
if pm
    L = p(L, ['regime: this point is a STABLE PM anchor (crit(odd-off) > 0, single-root with the ', ...
        'seed continuity), so rf/rd/collapse are the CLEAN emergence numbers -- rf-1 is the ', ...
        'genuine beyond-E1 transverse-spectator dressing, the tensor analogue of the projected ', ...
        '~2.8%% Tier-2 share; rd->1 is the matched E1 truncation (transverse dressing collapsed).']);
else
    L = p(L, ['regime: this point sits BELOW the tensor Tc (crit(odd-off) < 0, metastable-PM branch ', ...
        'reached by the Anderson map with seed continuity), so crit_shift_odd is a well-defined ', ...
        'criticality-margin difference but the odd-shift ratios rf/rd reflect that near/sub-critical ', ...
        'regime, NOT the clean emergence. Re-run at the stable anchor (T=2.0 K, B=(0.5 0 0) T) for ', ...
        'the emergence table; T12 recorded rf(1)=1.113 (beyond-E1 ~11%%), rd(1)=1.016 there.']);
end

txt = strjoin(L, newline);
if nargout == 0
    fprintf('%s\n', txt);
    clear txt;
end
end

function L = local_push(L, fmt, varargin)
L{end+1} = sprintf(fmt, varargin{:});
end
