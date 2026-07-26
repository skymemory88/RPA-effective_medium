function v = invz_pm_verdict(crit_pm, pm_ok, crit_band)
%INVZ_PM_VERDICT Three-way PM-probe verdict for the 1/z phase dispatcher (spec SS4.4).
%   'pm'                converged finite PM with crit_pm >  crit_band
%   'ordered_eligible'  converged finite PM with crit_pm < -crit_band
%   'unknown'           PM non-convergence / non-finite / recoverable error, OR
%                       |crit_pm| <= crit_band (boundary-indeterminate)
%
% A FAILED PM PROBE IS NOT EVIDENCE FOR ORDER. The historical dispatch ran the ordered solver
% whenever the PM probe was not valid and then labelled a converging jensen result 'ordered' --
% which makes solver availability the phase criterion. 'unknown' may run the ordered solver for
% DIAGNOSTICS but must never emit phase_1z = 1 without a separately validated free-energy /
% branch-selection rule.
% crit_band is the frozen crit_tol. Field-resolution uncertainty is represented separately
% by S.Bc_1z_interval; a parfor point never owns dcrit/dB.
if ~pm_ok || ~isfinite(crit_pm)
    v = 'unknown';
elseif crit_pm > crit_band
    v = 'pm';
elseif crit_pm < -crit_band
    v = 'ordered_eligible';
else
    v = 'unknown';
end
end
