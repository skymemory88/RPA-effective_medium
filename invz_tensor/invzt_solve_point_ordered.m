function pt = invzt_solve_point_ordered(ion, T, B, lat, opts)
%INVZT_SOLVE_POINT_ORDERED Tensor a1 solve at one FERROMAGNETIC (T, B) point.
%   pt = INVZT_SOLVE_POINT_ORDERED(ion, T, B, lat, opts) is the ordered-phase
%   counterpart of INVZT_SOLVE_POINT mode 'a1': the single-ion problem is solved
%   with the longitudinal ORDERING mean field (spontaneous m0 = <Jz>, hz =
%   J0z*<Jz> with J0z = lat.info.Jcc0 -- tensor-native threading, NOT ion.J0eff),
%   the self-energy uses the moment form (INVZ_SIGMA_ORDERED, HTML eqs 37-38,
%   lambda_{1,2,3}), and the medium step is the SAME tensor lattice engine as
%   the PM solver (INVZT_GCC_LATTICE 12x12 RPA + LOCKED K bookkeeping).
%
%   MEDIUM (2026-07-20 amendment, WHOLE-CC, no dominant/rest split): the local
%   chi0 fed to the lattice engine is the FULL electronuclear c0 = INVZ_CHI0Z(si,
%   T, i*wn, 'elastic', true), renormalized in one piece (ctil = c0./(1+Sigma)) --
%   mirroring the projected ordered reference (INVZ_SOLVE_POINT_ORDERED), which
%   also renormalizes the full G0, never a dominant-sector subset. The PM
%   solver's E < Esplit content cut (INVZT_CHI0_SPLIT) is NOT used here: measured
%   at T = 0.1 K, Bx = 3.0 T it selects only ndom = 13 states and carries a
%   DRIFTING ~48% of the local weight in "rest" (fdom_cc0 = 0.4762, review P0-3),
%   which made the no-ODD interop parity fail (dSigma0 = 8.6e-3, dalpha_m =
%   5.8e-2 at the T3 dispatch, both > 5e-3) when that rest weight was dropped by
%   chi_rest = false. A fixed content-defined dominant-vertex basis is a3d's
%   concern (Task 7D), not this a1 solver's -- so 'Esplit'/'chi_rest' are REJECTED
%   options here (invzt:orderedSplitKnobs), not silently accepted knobs.
%
%   SCOPE: TRANSVERSE fields only (spontaneous route; invzt:orderedLongitudinal
%   otherwise -- no forced_moment port, 2026-07-19 decision). Modes: 'a1' (this
%   task) and, once Task 7D lands, 'a3d' (full-response dominant-vertex,
%   Matsubara-only); full-dress 'a3' is PERMANENTLY rejected (invzt:orderedMode,
%   136-state vertex budget-refused). nlevels 'std' only
%   (invzt:orderedNlevels). Option-A parity with the projected ordered solver:
%   m0 is the BARE mean-field order parameter -- it onsets at the MF boundary
%   (~5.0 T at 0.1 K, MEASURED 2026-07-19), well above the tensor crit = 0 QPT
%   (4.65-4.70 T). A converged ordered point is therefore NOT evidence of being
%   in the FM phase; phase assignment belongs to INVZT_SOLVE_AUTO's
%   stability-based rule, never to this solver.
%
%   Mixing: plain damped mixing, CHECK-BEFORE-MIX (same loop ordering as the PM
%   solver, review P2-2): on every exit -- converged or not -- the returned
%   (Sigma, K, G, lambda, alpha, alpha_m) describe the SAME evaluated state.
%
%   OPTIONS (getf defaults): Ecut 40, hyp true, transverse_mf 'legacy_x',
%   mix_outer 0.7, tol_outer 1e-8, max_outer 80, m_tol 1e-2, bz_tol 1e-9,
%   odd true, rank_tol 1e-12, mz_seed/mf_maxit/mf_mix forwarded to
%   invz_single_ion. 'Esplit'/'chi_rest' (the PM dominant/rest split knobs) are
%   REJECTED here (invzt:orderedSplitKnobs) -- see MEDIUM above.
%
%   Returns the INVZT_SOLVE_POINT pt surface plus m0, alpha_m, is_ordered,
%   moment_branch ('spontaneous' | 'none'), J0z_used. Early returns (PM
%   relaxation |m0| <= m_tol; MF non-convergence) carry the projected-parity
%   fixed set: m0, is_ordered=false, converged=false, Sigma0=NaN, crit=NaN,
%   si, tl=[], moment_branch. Always check pt.is_ordered AND pt.converged.
%
%   See also INVZT_SOLVE_POINT, INVZ_SOLVE_POINT_ORDERED (projected reference),
%   INVZ_SIGMA_ORDERED, INVZ_TWOLEVEL_ORDERED, INVZT_GCC_LATTICE,
%   INVZT_CRIT_STATIC, INVZT_SOLVE_AUTO.
if nargin < 5, opts = struct(); end
mode   = getf(opts, 'mode', 'a1');
% At Task-3 time only 'a1' exists; Task 7D extends this gate to accept 'a3d'
% (full-response dominant-vertex, Matsubara-only). 'a3' (full dress) stays
% rejected PERMANENTLY: the 136-state vertex is budget-refused (A4 ladder gate).
if ~ismember(char(mode), {'a1'})            % Task 7D: {'a1','a3d'}
    error('invzt:orderedMode', ['invzt_solve_point_ordered implements mode ''a1'' ' ...
        '(and ''a3d'' once Task 7D lands) -- got ''%s''. Full-dress ordered ''a3'' ' ...
        'is permanently out of scope: the 136-state vertex is budget-refused.'], char(mode));
end
Ecut   = getf(opts, 'Ecut', 40);
hyp    = getf(opts, 'hyp', true);
tmf    = getf(opts, 'transverse_mf', 'legacy_x');
mixo   = getf(opts, 'mix_outer', 0.7);
tolo   = getf(opts, 'tol_outer', 1e-8);
maxo   = getf(opts, 'max_outer', 80);
mtol   = getf(opts, 'm_tol', 1e-2);
bztol  = getf(opts, 'bz_tol', 1e-9);
rank_tol = getf(opts, 'rank_tol', 1e-12);
odd = ~isfield(opts, 'odd') || isempty(opts.odd) || ~isequal(opts.odd, false);
% 2026-07-20 amendment (P0-3 fallout): the ordered medium is WHOLE-CC (no
% dominant/rest split) -- 'Esplit'/'chi_rest' are the PM solver's split knobs
% and have no meaning here. Fail loud rather than silently ignore them.
if isfield(opts, 'Esplit') || isfield(opts, 'chi_rest')
    error('invzt:orderedSplitKnobs', ['invzt_solve_point_ordered renormalizes the ' ...
        'WHOLE local cc susceptibility (no dominant/rest split, 2026-07-20 amendment): ' ...
        '''Esplit''/''chi_rest'' are PM-solver-only knobs and are not accepted here. The ' ...
        'E < Esplit content cut drifts to a ~48%% ''rest'' weight at ordered 0.1 K points ' ...
        '(review P0-3), which is ill-defined this deep in the ordered phase; a fixed ' ...
        'content-defined dominant-vertex basis is a3d''s concern (Task 7D), not a1''s.']);
end
nlevels = getf(opts, 'nlevels', 'std');
if ~strcmp(char(nlevels), 'std')
    error('invzt:orderedNlevels', ['invzt_solve_point_ordered supports nlevels = ' ...
        '''std'' only (full electronuclear single ion); got ''%s''.'], char(nlevels));
end

B = invz_field_vec(B);
if abs(B(3)) > bztol
    error('invzt:orderedLongitudinal', ['invzt_solve_point_ordered is the ' ...
        'TRANSVERSE (spontaneous-moment) route only: got Bz = %.3g T > bz_tol = ' ...
        '%.3g T. The forced-moment longitudinal route is not ported to the ' ...
        'tensor branch (2026-07-19 scope decision).'], B(3), bztol);
end
B(3) = 0;                                            % dead band: exactly transverse

J0z  = lat.info.Jcc0;                                % tensor-native uniform cc coupling
Jxx0 = lat.info.Jaa0;                                % v3 threading (parity with PM solver)

[wn, wts, beta] = invz_matsubara(T, Ecut);
nwn = numel(wn);

% --- ordered single-ion mean-field fixed point (full electronuclear space) -------
siopts = struct('hyp', hyp, 'order', true, 'J0z', J0z, 'Jxx0', Jxx0, 'transverse_mf', tmf);
for f = {'mz_seed', 'mf_maxit', 'mf_mix'}
    if isfield(opts, f{1}), siopts.(f{1}) = opts.(f{1}); end
end
si = invz_single_ion(ion, T, B, siopts);
m0 = si.Jexp(3);
if ~si.mf_converged
    pt = early_return(m0, si, 'spontaneous');
    return;
end
if abs(m0) <= mtol
    pt = early_return(m0, si, 'none');               % paramagnetic point: use invzt_solve_point
    return;
end

% --- two-level driver at the converged ordering field hz; WHOLE-CC chi0 (2026-07-20:
%     no dominant/rest split -- mirrors the projected ordered reference's full G0) ---
tl = invz_twolevel_ordered(ion, T, B, si.hz, struct('Jxx0', Jxx0, 'transverse_mf', tmf));
g  = real(invz_g(tl, 1i*wn));
c0 = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
c0_cc = real(squeeze(c0(3,3,:)));
chi0cc0 = c0_cc(1);

% --- odd = false mask on a COPY of lat.Jt (same rule as the PM solver) -----------
lat_eff = lat;
if ~odd
    lat_eff.Jt = invzt_odd_mask(lat.Jt);
end

% --- C2 assertion at Gamma (same as PM solver) -----------------------------------
JG = lat.JtGamma;
odd_ca = JG(3:3:12, 1:3:12);  odd_cb = JG(3:3:12, 2:3:12);
oddmag = max(abs([odd_ca(:); odd_cb(:)]));
if ~(oddmag < 1e-10*abs(lat.info.Jcc0))
    error('invzt:a1OddGamma', ['lat.JtGamma ODD (c<->a,b) blocks do not vanish ' ...
        '(max %.3e >= 1e-10*|Jcc0| = %.3e): C2 symmetry violated at Gamma.'], ...
        oddmag, 1e-10*abs(lat.info.Jcc0));
end

% --- outer moment-form loop: damped mixing, CHECK-BEFORE-MIX (P2-2) --------------
Sigma = zeros(nwn, 1);
denom = @(s) reshape(1 + s, 1, 1, nwn);
converged = false;
Gloc = nan(nwn, 1);  K = nan(nwn, 1);  diag4 = nan(4, nwn);
lam = nan(3, 1);  sg = struct('alpha', NaN, 'alpha_m', NaN, 'gamma', nan(nwn,1), 'Sigma', nan(nwn,1));
for outer = 1:maxo
    ctil = c0 ./ denom(Sigma);                       % [3,3,nwn] WHOLE-CC renormalized chi0
    [Gcc, diag4] = invzt_gcc_lattice(ctil, lat_eff);
    Gloc  = -Gcc(:);
    G0til = -(c0_cc ./ (1 + Sigma));
    K = 1 ./ Gloc - 1 ./ G0til;                      % LOCKED K bookkeeping (same as PM)
    lam = invz_lambdas(K, g, wts, beta, [1 2 3]);    % moment form needs lambda_3
    sg  = invz_sigma_ordered(tl, lam, K, g, beta);   % HTML eqs 37-38
    dS  = max(abs(sg.Sigma - Sigma));
    if dS < tolo, converged = true; break; end       % exit BEFORE mixing: state consistent
    Sigma = Sigma + mixo*(sg.Sigma - Sigma);
end

ctil0 = c0(:,:,1) / (1 + Sigma(1));
[crit, crit_clipped_mass, crit_active_rank] = invzt_crit_static(ctil0, lat.JtGamma, rank_tol);

% --- assemble pt (invzt_solve_point surface + ordered extras) --------------------
pt.m0     = m0;
pt.is_ordered = true;
pt.moment_branch = 'spontaneous';
pt.Sigma0 = Sigma(1);
pt.Sigma  = Sigma;
pt.alpha  = sg.alpha;
pt.alpha_m = sg.alpha_m;
pt.lambda = lam;
pt.K = K;
pt.G = Gloc;
pt.tl = tl;
pt.si = si;
pt.chi0cc0 = chi0cc0;
pt.crit = crit;
pt.crit_clipped_mass = crit_clipped_mass;
pt.crit_active_rank  = crit_active_rank;
pt.sumrule_rel = abs(sum(wts.*Gloc)/beta + si.JzJz_fluct) / max(abs(si.JzJz_fluct), 1e-12);
phys_finite = all(isfinite(Sigma)) && all(isfinite(K)) && all(isfinite(Gloc));
pt.converged = converged && phys_finite;
pt.outer_iters = outer;
pt.diag4_spread = max(max(diag4, [], 1) - min(diag4, [], 1));
pt.mode = 'a1';
pt.odd  = odd;
pt.chi_rest = true;   % interface parity with invzt_solve_point: WHOLE response always
                       % renormalized here (2026-07-20 amendment, no dominant/rest split)
pt.mspec = [];         % no chi0 split performed (mspec is a PM-solver-only diagnostic)
pt.Jxx0_used = Jxx0;
pt.J0z_used  = J0z;
pt.nlevels = 'std';
pt.lat = struct('qvec_hash', hash_vec(lat.qvec(:)), 'conv', lat.conv, ...
    'JtGamma', lat.JtGamma, 'info', lat.info, 'Jxx0_used', Jxx0);
end

% ------------------------------- local helpers ----------------------------------
function pt = early_return(m0, si, branch)
%EARLY_RETURN Complete field set for every non-accepted exit (projected parity).
pt = struct('m0', m0, 'is_ordered', false, 'converged', false, 'Sigma0', NaN, ...
            'crit', NaN, 'si', si, 'tl', [], 'moment_branch', branch);
end

function h = hash_vec(v)
% Weak grid-identity hash, same formula as invz_cache_key / invzt_solve_point.
v = v(:);
h = sprintf('%dv_%08x', numel(v), typecast(single(sum(v.*(1:numel(v))')), 'uint32'));
end
