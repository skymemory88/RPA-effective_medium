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
    for iq = 1:nqo
        M = ob.Vcc(:,:,iq) + dJ(:,:,iq);
        M = (M + M')/2;                        % both terms Hermitian; cleans rounding only
        Jnu_odd(iq,:) = sort(real(eig(M))).';
    end
    Jnu_flat = Jnu_odd(:);
    % E5 uniform shift, applied HERE exactly once (T1.3 bookkeeping rule: the grid
    % matrices' diagonal already carries -d via E4; J0eff carries the explicit -d;
    % NO other q = 0 handling). Callers pass the UNSHIFTED info.Jcc0 as opts.J0eff.
    % That unshifted value comes from invz_odd_blocks' infoB.Jcc0 (or the flag-off
    % invz_jq_modes info.Jcc0) -- NOT from invz_jq_modes' opts.odd path, whose
    % exported info.Jcc0 is ALREADY shifted by -d.
    J0eff = J0eff - d;
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
end
end
