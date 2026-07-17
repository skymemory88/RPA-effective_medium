function pt = invz_solve_point(ion, T, Bx, Jnu_flat, opts)
%INVZ_SOLVE_POINT Self-consistent 1/z solution at one paramagnetic (T, Bx) point.
% Bx: scalar (transverse, historical) or [Bx By Bz] vector (T).
% Outer loop: EMT K (Task 7) <-> lambda_p (Task 8) <-> Sigma, at fixed single-ion input.
% Inside the ordered phase the paramagnetic EMT fixed point does not exist; outputs may be non-finite and pt.converged false — invz_critical relies on this as the ordered-phase signal. Always check pt.converged.
%
% ODD extension (T1.4, opt-in): opts.odd = false (default) | true. When true,
% opts.odd_blocks = struct('Vca','Vcb','Vcc','Jcc0') (precomputed ONCE by the
% caller from invz_odd_blocks — P0.4: no disk/cache reads inside solver loops;
% Jcc0 is the UNSHIFTED Gamma value) is REQUIRED and Jnu_flat MUST be [] (both
% guarded by error id 'invz:oddArgs'). The cc modes are then rebuilt here from
% Vcc + deltaJ(T,Bx) (E1/E4 via invz_chiperp + invz_odd_deltaJ, evaluated at the
% SAME converged single-ion state as the rest of the solve), and the uniform
% coupling takes the explicit -d (E5): callers keep passing the UNSHIFTED
% info.Jcc0 as opts.J0eff, the shift is applied HERE exactly once. pt gains
% pt.odd = struct('d', 'Xp'). Flag off: byte-identical pre-ODD path.
%
% Retardation (T2.1, opt-in on top of opts.odd): opts.odd_retarded = true
% replaces the static (adiabatic) chi_perp(0) in the mediated coupling by
% chi_perp(i*wn) on the solve's own Matsubara grid, via the scalar surrogate
%   chi_perp(i*wn) ~ r_n * chi_perp(0),  r_n = Xpn(1,1,n)/Xpn(1,1,1) in (0,1],
% which is EXACT in the deltaJ MATRIX (E1 and E4 are linear in Xp, so
% dJ_n = r_n*dJ_0 and d_n = r_n*d follow identically); only the eigenvalue map
% is then approximated to first order around the static-ODD modes:
%   Jnu_n(q,nu) = Jnu_static(q,nu) + (r_n - 1) * u_nu'*dJ(q)*u_nu.
% The per-frequency spectra go to the generalized invz_emt_scalar as a
% [nq*4, nwn] matrix (column n <-> wn(n)). The residual error of the scalar
% surrogate itself (chi_ab and the aa/bb split of Xpn not scaling with r_n) is
% quantified by pt.odd.ab_ratio_max / pt.odd.bb_aa_dev_max. opts.
% odd_retarded_exact = true is the cross-check path: full per-(q,n) eig of
% Vcc + r_n*dJ (exact under the scalar surrogate; used by tests to measure the
% first-order-perturbation error, not in production). Both flags require
% opts.odd ('invz:oddArgs'); opts.odd_rn_override (scalar or [nwn,1], test
% hook) replaces the computed r_n. J0eff ALWAYS keeps the STATIC -d: the
% criticality condition pt.crit lives in the elastic n = 0 sector where
% r_1 = 1 exactly and chi_perp(0) is exact -- retardation never touches the
% uniform q = 0 bookkeeping (this is also why the invz_critical_T0field closed
% form is intrinsically static; ODD plan section 4 mitigation). pt.odd gains
% r_n (plus the smallness diagnostics when computed, and retarded_exact).
if nargin < 5, opts = struct(); end
Ecut  = getf(opts, 'Ecut', 40);
hyp   = getf(opts, 'hyp', true);
J0eff = getf(opts, 'J0eff', ion.J0eff);
Jxx0  = getf(opts, 'Jxx0', ion.Jxx0);
tmf   = getf(opts, 'transverse_mf', 'legacy_x');
mixo  = getf(opts, 'mix_outer', 0.7);
tolo  = getf(opts, 'tol_outer', 1e-8);
maxo  = getf(opts, 'max_outer', 60);
eopts = getf(opts, 'emt', struct());
Bx = invz_field_vec(Bx);                       % scalar -> [Bx 0 0]; 3-vector passes through

% --- ODD diversion (T1.4): strictly additive and opt-in; everything else in
% this function is the pre-ODD code path, byte-untouched when the flag is off
% (regression test test_solve_point_flag_off_bitwise).
oddOn = isfield(opts, 'odd') && ~isempty(opts.odd) && ~isequal(opts.odd, false);
oddRet   = isfield(opts, 'odd_retarded') && ~isempty(opts.odd_retarded) ...
    && ~isequal(opts.odd_retarded, false);
oddRetEx = isfield(opts, 'odd_retarded_exact') && ~isempty(opts.odd_retarded_exact) ...
    && ~isequal(opts.odd_retarded_exact, false);          % wins over odd_retarded
if (oddRet || oddRetEx) && ~oddOn
    error('invz:oddArgs', ['opts.odd_retarded / opts.odd_retarded_exact require ' ...
        'opts.odd = true: retardation modifies the ODD-mediated coupling deltaJ, ' ...
        'which only exists on the ODD path.']);
end
if oddOn
    ob = getf(opts, 'odd_blocks', []);
    if ~(isstruct(ob) && isscalar(ob) && all(isfield(ob, {'Vca', 'Vcb', 'Vcc', 'Jcc0'})))   % Jcc0: contract/audit field, required for caller-side consistency; unused here (J0eff comes via opts)
        error('invz:oddArgs', ['opts.odd = true requires opts.odd_blocks = struct(' ...
            '''Vca'',''Vcb'',''Vcc'',''Jcc0'') precomputed once by the caller from ' ...
            'invz_odd_blocks (P0.4: no disk/cache reads inside solver loops; Jcc0 UNSHIFTED).']);
    end
    if ~isempty(Jnu_flat)
        error('invz:oddArgs', ['opts.odd = true requires Jnu_flat = []: the cc modes are ' ...
            'rebuilt here from odd_blocks + deltaJ, and a caller-supplied baseline Jnu ' ...
            'would silently override the rebuild.']);
    end
    % chi_perp at the SAME single-ion options the solve below resolves (T1.2's
    % same-converged-state requirement); oddXi.si is that state, reused below.
    [Xp, oddXi] = invz_chiperp(ion, T, Bx, struct('hyp', hyp, 'Jxx0', Jxx0, 'transverse_mf', tmf));
    [dJ, d] = invz_odd_deltaJ(ob.Vca, ob.Vcb, Xp);
    nqo = size(ob.Vcc, 3);
    Jnu_odd = zeros(nqo, 4);
    w_odd = zeros(nqo, 4);                     % T2.1 first-order weights u_nu'*dJ*u_nu
    for iq = 1:nqo
        M = ob.Vcc(:,:,iq) + dJ(:,:,iq);
        M = (M + M')/2;                        % both terms Hermitian; cleans rounding only
        Jnu_odd(iq,:) = sort(real(eig(M))).';
        if oddRet && ~oddRetEx
            % Eigvecs from a SEPARATE eig call: the values-only eig above stays
            % bitwise-identical to the static path (LAPACK jobz='N' vs 'V' may
            % differ in last bits), so the elastic column of the retarded solve
            % is exactly the static spectrum. Degenerate pairs make individual
            % u_nu ambiguous; the first-order error that induces is part of
            % what the surrogate-vs-exact gate (AbsTol 1e-3 on Sigma0,
            % test_invz_odd_retarded) measures.
            [U, ev] = eig(M, 'vector');
            [~, is] = sort(real(ev));  U = U(:, is);
            for nu = 1:4
                w_odd(iq, nu) = real(U(:,nu)' * dJ(:,:,iq) * U(:,nu));
            end
        end
    end
    Jnu_flat = Jnu_odd(:);                     % column-major: all q at nu=1, then nu=2, ...
    % E5 uniform shift, applied HERE exactly once (T1.3 bookkeeping rule: the grid
    % matrices' diagonal already carries -d via E4; J0eff carries the explicit -d;
    % NO other q = 0 handling). Callers pass the UNSHIFTED info.Jcc0 as opts.J0eff.
    % That unshifted value comes from invz_odd_blocks' infoB.Jcc0 (or the flag-off
    % invz_jq_modes info.Jcc0) -- NOT from invz_jq_modes' opts.odd path, whose
    % exported info.Jcc0 is ALREADY shifted by -d.
    J0eff = J0eff - d;
    % --- T2.1 retardation: per-frequency mode spectra Jnu_flat [nq*4, nwn] ---
    if oddRet || oddRetEx
        % Same deterministic (T, Ecut) grid as the [wn, wts, beta] call below,
        % so column n of Jnu_flat pairs with the solve's wn(n) (wn(1) = 0).
        wn_r = invz_matsubara(T, Ecut);
        nwn  = numel(wn_r);
        odd_ret_diag = struct();
        rn_ov = getf(opts, 'odd_rn_override', []);
        if ~isempty(rn_ov)
            % Test hook: forced r_n (scalar or [nwn,1]); skips the chi_perp
            % Matsubara solve, so the smallness diagnostics are unavailable.
            if isscalar(rn_ov), r_n = repmat(rn_ov, nwn, 1); else, r_n = rn_ov(:); end
            if numel(r_n) ~= nwn || ~all(isfinite(r_n))
                error('invz:oddArgs', ...
                    'odd_rn_override must be a finite scalar or [%d,1] vector (the Matsubara grid).', nwn);
            end
        else
            Xpn = invz_chiperp(ion, T, Bx, struct('hyp', hyp, 'Jxx0', Jxx0, ...
                'transverse_mf', tmf, 'z', 1i*wn_r));    % [2,2,nwn]
            xaa = squeeze(Xpn(1,1,:));
            r_n = xaa / xaa(1);                          % r_n(1) = 1 exactly (x/x)
            % Smallness diagnostics justifying the scalar surrogate
            % Xp_n ~ r_n*Xp_0: the gyrotropic-symmetrized chi_ab is tiny and
            % chi_bb tracks chi_aa (stored on pt.odd; reported by the tests).
            odd_ret_diag.ab_ratio_max  = max(abs(squeeze(Xpn(1,2,:))) ./ xaa);
            odd_ret_diag.bb_aa_dev_max = max(abs(squeeze(Xpn(2,2,:)) - xaa) ./ xaa);
            assert(abs(r_n(1) - 1) <= 1e-12, 'invz:oddRn', ...
                'r_n(1) = %.17g differs from 1 beyond 1e-12.', r_n(1));
            assert(all(isfinite(r_n)) && all(r_n > 0 & r_n <= 1 + 1e-12), 'invz:oddRn', ...
                'r_n outside (0, 1]: min %.6g, max %.17g.', min(r_n), max(r_n));
            assert(all(diff(r_n) <= 1e-12), 'invz:oddRn', ...
                'r_n not monotone non-increasing along the Matsubara grid.');
        end
        if oddRetEx
            % EXACT cross-check path: E1 and E4 are LINEAR in Xp, so under the
            % scalar surrogate Xp_n = r_n*Xp_0 the retarded coupling matrix is
            % dJ_n = r_n*dJ EXACTLY (and d_n = r_n*d) -- the surrogate is exact
            % in the MATRIX; only the eigenvalue map is approximated by the
            % first-order route below. Full per-(q,n) eig; r_n(1) = 1 makes the
            % elastic column bitwise the static spectrum.
            Jnu_flat = zeros(nqo*4, nwn);
            for n = 1:nwn
                Jn = zeros(nqo, 4);
                for iq = 1:nqo
                    M = ob.Vcc(:,:,iq) + r_n(n)*dJ(:,:,iq);
                    M = (M + M')/2;
                    Jn(iq,:) = sort(real(eig(M))).';
                end
                Jnu_flat(:, n) = Jn(:);                  % same column-major order
            end
        else
            % First-order eigenvalue perturbation around the static-ODD modes.
            % Jnu_odd(:) and w_odd(:) share the [nq,4] column-major flattening
            % (all q at nu=1, then nu=2, ...), so every column keeps the static
            % path's row order; with r_n identically 1 the columns are bitwise
            % the static vector (x + 0*w == x) and the generalized EMT's
            % constant-column path is bitwise-equal to the vector path
            % (test_emt_matrix_column_constant_bitwise).
            Jnu_flat = Jnu_odd(:) + (r_n(:).' - 1) .* w_odd(:);   % [nq*4, nwn]
        end
        % J0eff keeps the STATIC -d applied above: criticality (pt.crit) sits in
        % the elastic n = 0 sector where r_n = 1 and chi_perp(0) is exact, so the
        % n-dependent d_n = r_n*d never enters the uniform q = 0 bookkeeping.
    end
end

[wn, wts, beta] = invz_matsubara(T, Ecut);
if oddOn
    % Bit-identical reuse: invz_chiperp built its si with the SAME (ion, T, Bx)
    % and the SAME (hyp, Jxx0, transverse_mf) options as the call below
    % (invz_field_vec is idempotent), so this skips one 136-dim diagonalization
    % without changing any result.
    si = oddXi.si;
else
si  = invz_single_ion(ion, T, Bx, struct('hyp', hyp, 'Jxx0', Jxx0, 'transverse_mf', tmf));
end
c0  = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
G0  = -real(squeeze(c0(3,3,:)));                 % full (electro)nuclear cc Green function
tl  = invz_twolevel(ion, T, Bx, struct('Jxx0', Jxx0, 'transverse_mf', tmf));   % electronic two-level params for Sigma
g   = real(invz_g(tl, 1i*wn));

% Warm-start seed (finding #6): a converged neighbouring point's Sigma speeds up
% the outer fixed point. It must be the SAME length as wn (same T), and only
% changes the iteration path, not the converged fixed point; ignored otherwise.
% (Only Sigma is seeded: the closed-form EMT solves K directly and ignores any K
% seed, so seeding K had no effect.)
Sigma = zeros(size(wn));
if isfield(opts,'Sigma_seed') && numel(opts.Sigma_seed) == numel(wn), Sigma = opts.Sigma_seed(:); end
converged = false;
for outer = 1:maxo
    med = invz_emt_scalar(G0, Sigma, Jnu_flat, eopts);
    K   = med.K;
    lam = invz_lambdas(K, g, wts, beta, [1 2]);
    sg  = invz_sigma(tl, lam, K, g, beta);
    dS  = max(abs(sg.Sigma - Sigma));
    Sigma = Sigma + mixo*(sg.Sigma - Sigma);
    if dS < tolo, converged = true; break; end
end
pt.Sigma0 = Sigma(1);  pt.alpha = sg.alpha;  pt.lambda = lam;
pt.K = K;  pt.G = med.G;  pt.Sigma = Sigma;  pt.tl = tl;  pt.si = si;
pt.chi0cc0 = -G0(1);
pt.crit = 1 + pt.Sigma0 - J0eff*pt.chi0cc0;
pt.sumrule_rel = abs(sum(wts.*med.G)/beta + si.JzJz_fluct) / si.JzJz_fluct;
pt.converged = converged && med.converged;
pt.outer_iters = outer;
if oddOn
    pt.odd = struct('d', d, 'Xp', Xp);         % T1.4 diagnostics (absent when flag off)
    if oddRet || oddRetEx                      % T2.1 diagnostics
        pt.odd.r_n = r_n;
        pt.odd.retarded_exact = oddRetEx;
        fn = fieldnames(odd_ret_diag);         % ab_ratio_max / bb_aa_dev_max
        for k = 1:numel(fn), pt.odd.(fn{k}) = odd_ret_diag.(fn{k}); end
    end
end
end
