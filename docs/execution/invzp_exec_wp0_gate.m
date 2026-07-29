function G = invzp_exec_wp0_gate(save_path)
%INVZP_EXEC_WP0_GATE WP0 gate (invzp_convg_fix.md Sec.9, L389-401): "two independent
% constructors reproduce all frozen microscopic anchors and the same serialized coupling
% set." Builds the frozen contract (invzf_wp0_contract) and checks:
%   G0a reproducibility -- a SECOND independent invz_bz_couplings call, per backend,
%                          reproduces the contract's exact-byte digest.
%   G0b independent constructors agree at finite q -- bruteforce vs ewald, mode by mode, at
%                          every q the grid retains (Gamma excluded -- the exclusion is
%                          VERIFIED via invz_is_gamma_equiv, not assumed).
%   G0c Gamma-point difference is the documented one -- quantified against the Lorentz term.
%   G0d limit order -- RECORDED (not tested), transcribed from invzp_convg_fix.md Sec.9/2.5/4.2.
%
% Chosen (not repo-frozen) tolerances, stated plainly so they are auditable, never tuned
% after seeing the numbers below:
%   G0b: max|branch diff| < 1% of the bruteforce Jcc0 scale (physical coupling energy scale;
%        an eigenvalue-wise relative error is not used because Jcc(q) branches cross zero).
%   G0c: |Jcc0 diff| and max|Jgamma_cc diff| both < 5% of the Lorentz term itself (a mismatch
%        of ORDER the Lorentz term would mean the correction was not applied to one backend).
if nargin < 1, save_path = ''; end
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
    qvec = bf.detail_summary.qvec;
    assert(isequal(qvec, ew.detail_summary.qvec), 'invz:wp0GateGridMismatch', ...
        'bruteforce/ewald backends were not evaluated on the same q-grid.');
    isG = invz_is_gamma_equiv(qvec, ion.tau);      % same test invz_jq_modes/invz_jq_path use
    n_gamma = nnz(isG);
    Ubf = bf.detail_summary.Jnu_unflat(~isG, :);   ew_full = ew.detail_summary.Jnu_unflat;
    Uew = ew_full(~isG, :);                        qk = qvec(~isG, :);

    % --- G0b ---
    absd = abs(Uew - Ubf);
    rel_floor = 1e-6;                              % meV; avoids blow-up at branch zero-crossings
    reld = absd ./ max(abs(Ubf), rel_floor);
    [maxabs, iabs] = max(absd(:));  [ra, ca] = ind2sub(size(absd), iabs);
    [maxrel, irel] = max(reld(:));  [rr, cr] = ind2sub(size(reld), irel);
    scale = abs(bf.Jcc0);  tol_b = 1e-2;
    G0b = struct('run', true, 'n_compared', numel(absd), 'n_gamma_excluded', n_gamma, ...
        'max_abs_diff', maxabs, 'q_at_max_abs', qk(ra,:), 'branch_at_max_abs', ca, ...
        'max_rel_diff', maxrel, 'rel_floor', rel_floor, 'q_at_max_rel', qk(rr,:), ...
        'branch_at_max_rel', cr, 'scale', scale, 'tol', tol_b, 'pass', (maxabs/scale) < tol_b);
    fprintf(['G0b: %d (q,branch) pairs compared, %d Gamma-equivalent rows excluded; ' ...
        'max abs diff=%.6g meV at q=[%.4g %.4g %.4g] branch %d; max rel diff=%.4g ' ...
        '(floor %.1g) at q=[%.4g %.4g %.4g] branch %d -> %s\n'], G0b.n_compared, n_gamma, ...
        maxabs, qk(ra,:), ca, maxrel, rel_floor, qk(rr,:), cr, pf(G0b.pass));
    if n_gamma > 0
        fprintf('   %d Gamma-equivalent rows excluded per the documented exclusion.\n', n_gamma);
    else
        fprintf('   0 Gamma-equivalent rows in this grid: the exclusion is verified vacuous here.\n');
    end
    fprintf('   documented Gamma formulas -- bruteforce: %s\n', contract.qgrid.gamma_formula_bruteforce);
    fprintf('                                 ewald:      %s\n', contract.qgrid.gamma_formula_ewald);

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
