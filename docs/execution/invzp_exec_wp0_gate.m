function G = invzp_exec_wp0_gate(save_path, dpRng_ladder)
%INVZP_EXEC_WP0_GATE WP0 gate (invzp_convg_fix.md Sec.9, L389-401): "two independent
% constructors reproduce all frozen microscopic anchors and the same serialized coupling
% set." Builds the frozen contract (invzf_wp0_contract) and checks:
%   G0a reproducibility -- a SECOND independent invz_bz_couplings call, per backend,
%                          reproduces the contract's exact-byte digest.
%   G0b independent constructors agree at finite q -- bruteforce vs ewald, mode by mode, at
%                          every q the grid retains (Gamma excluded -- the exclusion is
%                          VERIFIED via invz_is_gamma_equiv, not assumed). The bruteforce leg
%                          is swept over a dpRng_ladder (default [30 45 60 90], the production
%                          default dpRng=30 always included as the first rung) against the
%                          SAME single ewald reference, because bruteforce truncates a
%                          conditionally-convergent real-space sum at a finite cutoff while
%                          ewald is convergent -- a shrinking bf/ewald gap with dpRng is a
%                          truncation artifact, not a constructor disagreement; a plateau
%                          would be. 2026-07-29 measured ladder (this grid/ewald opts): RMS
%                          diff over all compared pairs falls monotonically at every rung
%                          (3.32e-6 -> 9.82e-7 -> 8.81e-7 -> 4.00e-7 meV, dpRng=30/45/60/90),
%                          a power-law fit of ln(rms) vs ln(dpRng) gives exponent ~-1.83 with
%                          no sign of a nonzero floor; max|diff| falls ~9x overall (1.14e-4 ->
%                          1.27e-5 meV) but is not strictly monotone (dpRng=60 ticks up from
%                          dpRng=45), consistent with shell-commensurability noise on a
%                          decaying envelope, not a fixed convention offset. dpRng=120 was not
%                          included in the default ladder: dpRng=90 alone took ~23 min on this
%                          grid (~64x dpRng=30's ~1 min; cost scales close to dpRng^2.4-3), so
%                          a fifth rung is a separate, explicit opt-in (pass it in
%                          dpRng_ladder), not part of the routine gate.
%   G0c Gamma-point difference is the documented one -- quantified against the Lorentz term.
%   G0d limit order -- RECORDED (not tested), transcribed from invzp_convg_fix.md Sec.9/2.5/4.2.
%
% Chosen (not repo-frozen) tolerances, stated plainly so they are auditable, never tuned
% after seeing the numbers below:
%   G0b: max|branch diff| < 1% of the bruteforce Jcc0 scale (physical coupling energy scale;
%        an eigenvalue-wise relative error is not used because Jcc(q) branches cross zero).
%        Evaluated at the BEST-TESTED (largest) dpRng rung -- the question G0b answers is
%        whether the two constructors agree in the converged limit; production_pass (at
%        dpRng=30 specifically) is reported alongside as a separate, not-swept-under-the-rug
%        finding about the CURRENT production default's own convergence adequacy.
%   G0c: |Jcc0 diff| and max|Jgamma_cc diff| both < 5% of the Lorentz term itself (a mismatch
%        of ORDER the Lorentz term would mean the correction was not applied to one backend).
if nargin < 1, save_path = ''; end
if nargin < 2 || isempty(dpRng_ladder), dpRng_ladder = [30 45 60 90]; end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); addpath(fullfile(root,'invz_common')); addpath(fullfile(root,'invz_projected'));
addpath(fullfile(root,'invz_functional'));

contract = invzf_wp0_contract();
ion = contract.ion;

fprintf('=== WP0 gate: freeze the microscopic contract ===\n');

% ---- G0a: reproducibility (second independent call per backend) -------------------------
backends = {'bruteforce','ewald'};
G0a = struct();
for i = 1:2
    b = backends{i};
    rec = contract.backends.(b);
    if ~rec.ok
        G0a.(b) = struct('run', false, 'pass', false, ...
            'reason', sprintf('%s: %s', rec.error_id, rec.error_msg));
        fprintf('G0a [%s]: NOT RUN -- %s\n', b, G0a.(b).reason);
        continue;
    end
    Jnu2 = invz_bz_couplings(ion, contract.coupling_opts.(b));
    d2 = invz_exact_numeric_digest(Jnu2(:));
    same = strcmp(d2, rec.digest);
    G0a.(b) = struct('run', true, 'digest1', rec.digest, 'digest2', d2, 'pass', same);
    fprintf('G0a [%s]: digest1=%s... digest2=%s... -> %s\n', b, rec.digest(1:12), d2(1:12), pf(same));
end

% ---- G0b/G0c: need both backends -----------------------------------------------------------
bf = contract.backends.bruteforce;  ew = contract.backends.ewald;
if ~bf.ok || ~ew.ok
    missing = {};
    if ~bf.ok, missing{end+1} = sprintf('bruteforce (%s: %s)', bf.error_id, bf.error_msg); end
    if ~ew.ok, missing{end+1} = sprintf('ewald (%s: %s)', ew.error_id, ew.error_msg); end
    reason = sprintf('backend(s) unavailable: %s', strjoin(missing, '; '));
    G0b = struct('run', false, 'pass', false, 'reason', reason);
    G0c = struct('run', false, 'pass', false, 'reason', reason);
    fprintf('G0b: NOT RUN -- %s\n', reason);
    fprintf('G0c: NOT RUN -- %s\n', reason);
else
    qvec = ew.detail_summary.qvec;
    isG = invz_is_gamma_equiv(qvec, ion.tau);      % same test invz_jq_modes/invz_jq_path use
    n_gamma = nnz(isG);
    Uew = ew.detail_summary.Jnu_unflat(~isG, :);   qk = qvec(~isG, :);
    if n_gamma > 0
        fprintf('   %d Gamma-equivalent rows excluded per the documented exclusion.\n', n_gamma);
    else
        fprintf('   0 Gamma-equivalent rows in this grid: the exclusion is verified vacuous here.\n');
    end
    fprintf('   documented Gamma formulas -- bruteforce: %s\n', contract.qgrid.gamma_formula_bruteforce);
    fprintf('                                 ewald:      %s\n', contract.qgrid.gamma_formula_ewald);

    % --- G0b: dpRng convergence ladder (bruteforce, fresh calls) vs the ONE ewald reference ---
    rel_floor = 1e-6;                              % meV; avoids blow-up at branch zero-crossings
    nd = numel(dpRng_ladder);
    ladder(nd) = struct('dpRng',[],'time_s',[],'max_abs',[],'q_max_abs',[],'branch_max_abs',[], ...
        'max_rel',[],'q_max_rel',[],'branch_max_rel',[],'Jbf_at_max_rel',[],'Jew_at_max_rel',[], ...
        'rms',[]);
    for k = 1:nd
        dp = dpRng_ladder(k);
        bzOpts = struct('grid', contract.qgrid.grid_used, 'dpRng', dp, 'cache', false, 'dipole', 'bruteforce');
        t0 = tic;
        [~, ~, ~, detail] = invz_bz_couplings(ion, bzOpts);
        tk = toc(t0);
        assert(isequal(detail.qvec, qvec), 'invz:wp0GateGridMismatch', ...
            'dpRng=%d: bruteforce/ewald not evaluated on the same q-grid.', dp);
        Ubf = detail.Jnu_unflat(~isG, :);
        absd = abs(Uew - Ubf);  reld = absd ./ max(abs(Ubf), rel_floor);
        [ma, ia] = max(absd(:));  [ra, ca] = ind2sub(size(absd), ia);
        [mr, ir] = max(reld(:));  [rr, cr] = ind2sub(size(reld), ir);
        ladder(k) = struct('dpRng', dp, 'time_s', tk, 'max_abs', ma, 'q_max_abs', qk(ra,:), ...
            'branch_max_abs', ca, 'max_rel', mr, 'q_max_rel', qk(rr,:), 'branch_max_rel', cr, ...
            'Jbf_at_max_rel', Ubf(rr,cr), 'Jew_at_max_rel', Uew(rr,cr), ...
            'rms', sqrt(mean(absd(:).^2)));
        fprintf(['G0b ladder dpRng=%3d (t=%.0fs): max_abs=%.4e meV @ q=[%.3f %.3f %.3f] br%d; ' ...
            'max_rel=%.4f @ q=[%.3f %.3f %.3f] br%d (|Jbf|=%.3e |Jew|=%.3e meV); rms=%.4e meV\n'], ...
            dp, tk, ma, qk(ra,:), ca, mr, qk(rr,:), cr, Ubf(rr,cr), Uew(rr,cr), sqrt(mean(absd(:).^2)));
    end
    rmsv = [ladder.rms];
    scale = abs(bf.Jcc0);  tol_b = 1e-2;
    prod_idx = find([ladder.dpRng] == 30, 1);
    best_idx = numel(ladder);
    has_prod = ~isempty(prod_idx);
    prod_pass = has_prod && (ladder(prod_idx).max_abs/scale) < tol_b;
    best_pass = (ladder(best_idx).max_abs/scale) < tol_b;
    G0b = struct('run', true, 'n_compared', numel(qk)*size(Uew,2), 'n_gamma_excluded', n_gamma, ...
        'ladder', ladder, 'rms_monotone_nonincreasing', all(diff(rmsv) <= 0), ...
        'scale', scale, 'tol', tol_b, 'production_dpRng', 30, 'production_pass', prod_pass, ...
        'best_tested_dpRng', ladder(best_idx).dpRng, 'best_tested_pass', best_pass, ...
        'pass', best_pass);
    if has_prod
        fprintf(['G0b summary: %d (q,branch) pairs/rung, %d Gamma rows excluded; RMS monotone ' ...
            'non-increasing across the ladder = %d; production dpRng=30: max_abs/scale=%.4f%% ' ...
            '(tol %.0f%%) -> %s; best-tested dpRng=%d: max_abs/scale=%.4f%% -> %s (G0b.pass)\n'], ...
            G0b.n_compared, n_gamma, G0b.rms_monotone_nonincreasing, ...
            100*ladder(prod_idx).max_abs/scale, 100*tol_b, pf(prod_pass), ladder(best_idx).dpRng, ...
            100*ladder(best_idx).max_abs/scale, pf(best_pass));
    else
        fprintf(['G0b summary: %d (q,branch) pairs/rung, %d Gamma rows excluded; RMS monotone ' ...
            'non-increasing across the ladder = %d; dpRng=30 not in this ladder (no ' ...
            'production_pass); best-tested dpRng=%d: max_abs/scale=%.4f%% -> %s (G0b.pass)\n'], ...
            G0b.n_compared, n_gamma, G0b.rms_monotone_nonincreasing, ladder(best_idx).dpRng, ...
            100*ladder(best_idx).max_abs/scale, pf(best_pass));
    end

    % --- G0c ---
    Cc = contract.units.constants;
    lorz = 4*pi/(3*ion.Vc)*Cc.gfac;
    Jg_bf = bf.info_summary.Jgamma_cc;  Jg_ew = ew.info_summary.Jgamma_cc;
    dGamma_mat = max(abs(Jg_ew(:) - Jg_bf(:)));
    dJcc0 = ew.Jcc0 - bf.Jcc0;
    tol_c = 0.05;
    G0c = struct('run', true, 'lorz', lorz, 'Jcc0_bruteforce', bf.Jcc0, 'Jcc0_ewald', ew.Jcc0, ...
        'diff_Jcc0', dJcc0, 'ratio_to_lorz', abs(dJcc0)/lorz, ...
        'max_abs_Jgamma_matrix_diff', dGamma_mat, 'matrix_ratio_to_lorz', dGamma_mat/lorz, ...
        'tol_frac_of_lorz', tol_c, ...
        'pass', abs(dJcc0)/lorz < tol_c && dGamma_mat/lorz < tol_c);
    fprintf(['G0c: lorz=%.6g meV; Jcc0 bruteforce=%.8g ewald=%.8g diff=%.4g meV (%.3g%% of lorz); ' ...
        'Jgamma_cc max abs diff=%.4g meV (%.3g%% of lorz); tol=%.0f%% of lorz -> %s\n'], ...
        lorz, bf.Jcc0, ew.Jcc0, dJcc0, 100*G0c.ratio_to_lorz, dGamma_mat, ...
        100*G0c.matrix_ratio_to_lorz, 100*tol_c, pf(G0c.pass));
end

% ---- G0d: declared limit order (RECORD, do not test) -------------------------------------
G0d = struct('run', true, 'declared_order', [ ...
    '(1) solve the stationary equations at NONZERO conjugate (symmetry-breaking) source; ' ...
    '(2) take the thermodynamic limit (coupling grid/momentum sum -> continuum BZ) at FIXED ' ...
    'nonzero source and a declared dipolar shape/boundary convention; (3) only then take the ' ...
    'source to zero. Momentum (q->0+) and zero-frequency (omega->0) limits sit inside a ' ...
    'declared finite-volume sequence taken before step (3), never substituting for it. ' ...
    '[invzp_convg_fix.md Sec.9 WP0 L393-394; Sec.2.5 L74-78 "controlled limits"; ' ...
    'Sec.4.2 L180-186 "essential order of limits"]']);
fprintf('G0d (RECORDED, not tested): %s\n', G0d.declared_order);

ran = G0a.bruteforce.run && G0a.ewald.run && G0b.run;
if ran
    fprintf('OVERALL (gates that ran): %s\n', pf(G0a.bruteforce.pass && G0a.ewald.pass && G0b.pass && G0c.pass));
else
    fprintf('OVERALL: INCOMPLETE -- one or more gates NOT RUN (see reasons above)\n');
end

G = struct('contract', contract, 'G0a', G0a, 'G0b', G0b, 'G0c', G0c, 'G0d', G0d);
outDir = fullfile(root, 'docs','execution','out');
save(fullfile(outDir, 'wp0_contract.mat'), 'contract');
if ~isempty(save_path), save(save_path, 'G'); end
fprintf('contract saved: %s\n', fullfile(outDir, 'wp0_contract.mat'));
end

function s = pf(t)
if t, s = 'PASS'; else, s = 'FAIL'; end
end
