function [Tc, out] = invz_odd_zero_field(ion, opts)
%INVZ_ODD_ZERO_FIELD Zero-field ordering temperature Tc(0) with/without ODD.
%   [Tc, out] = INVZ_ODD_ZERO_FIELD(ion, opts) is the T1.5 measurement engine and
%   the V4.1 zero-field endpoint source. It builds the ODD geometric blocks on one
%   or more Gamma-excluded uniform BZ meshes, forms the closed-form critical
%   self-energy Sigma_c(0;T) and the uniform coupling J(0;T), and root-finds the
%   1/z ordering temperature Tc via invz_critical_T0field. With ODD off it
%   reproduces the published Rönnow-2007 route bit-for-bit; with ODD on it applies
%   the E1/E4/E5 mediated coupling deltaJ^{cc}(T) and the E5 uniform reduction -d(T).
%
%   opts fields (all optional):
%     grids  ({12, 24})   cell array (or numeric vector) of per-axis grid sizes n;
%                         two grids => Richardson extrapolation (see below). A
%                         single grid returns that grid's value as *_rich.
%     dpRng  (30)         real-space dipole cutoff for invz_odd_blocks.
%     cache  (true)       invz_odd_blocks (odd1_) file cache.
%     mode   ('full')     'off' | 'full' | 'uniform_only' | 'qstruct_only'.
%
%   MODES (the headline Tc of each; mode 'full' additionally emits out.split):
%     'off'          : published route. Sc(T) = invz_sigma_crit(Jcc0, modes(Vcc)),
%                      T-independent; J0 = Jcc0. Reproduces Sigma_c(0) = 0.2980 and
%                      Tc(0) = 1.743 K on the (12^3, 24^3) Richardson pair.
%     'full'   (a)   : (Jcc0 - d(T))*chi0cc = 1 + Sc_odd(T), with
%                      Sc_odd(T) = invz_sigma_crit(Jcc0 - d(T), modes(Vcc + dJ(T))).
%     'uniform_only' : DIAGNOSTIC (source-plan b, kept as reported). Modes of Vcc,
%                      J0 = Jcc0 - d(T). This drives J0 below the peak finite-q Vcc
%                      mode (excluded modes: nex) -- DS2023's naive-MF inconsistency.
%     'qstruct_only' : DIAGNOSTIC (source-plan c). Modes of Vcc + dJ + d*eye(4),
%                      J0 = Jcc0. This is a simultaneous uniform +d shift of (a)'s
%                      couplings AND J0, which leaves the R2007 criticality
%                      invariant => Tc == Tc_full to numerical precision (the
%                      closed-form THEOREM validating the E4/E5 bookkeeping).
%
%   DECOMPOSITION (controller adjudication 2026-07-17, GOVERNING split; the source
%   plan's dTc-space split was ill-posed -- see the module note and ODD-LOG T1.5).
%   out.split (mode 'full') is the sequential condition/Sigma-space factorial in
%   (J0-shift) x (Sigma-source); neither governing leg enters the invalid regime:
%     off    : J0 = Jcc0,      Sc = Sc_off                  (Tc_off)
%     (a)    : J0 = Jcc0 - d,  Sc = Sc_odd                  (Tc = headline)
%     (b)    : J0 = Jcc0 - d,  Sc = Sc_off (FROZEN no-ODD)  (Tc_b, DS2023 MF analog)
%     (c)    : J0 = Jcc0,      Sc = Sc_odd                  (Tc_c, fluctuation piece)
%   where Sc_off = invz_sigma_crit(Jcc0, modes(Vcc)) [no exclusions possible] and
%   Sc_odd(T) = invz_sigma_crit(Jcc0 - d(T), modes(Vcc + dJ(T))) [ODD, at its own
%   consistent config -- dJ lowers the peak so the uniform mode stays critical].
%   out.split also carries: closure_defect = (a) - [(b) + (c) - Tc_off] (the genuine
%   interaction term -- REPORTED, never gated), Tc_b_naive (+ nex.b_naive counts,
%   source-plan b), Tc_c1_literal (== Tc_a, the bookkeeping theorem), and
%   Tc_c_factorial (modes Vcc + dJ, J0 = Jcc0).
%
%   GRID / GAMMA-FILTER / RICHARDSON CONVENTIONS -- copied verbatim from the
%   published Sigma_c benchmark so ODD-off reproduces the number bitwise:
%     invz_projected/tests/test_invz_sigma_crit.m, test_lihof4_sigma_crit
%       line 41: [T,qvec] = evalc("qVec_generator(ion.a,'mode','grid','grid',...
%                            [n n n],'range',[-0.5 0.5])");   % inclusive linspace
%       line 42: qvec = qvec(any(abs(qvec) > 1e-12, 2), :);  % Gamma-exclusion
%       line 46: Sc = 2*S(2) - S(1);                         % Richardson (12,24)
%     The range [-0.5 0.5] inclusive grid has no exact Gamma node for n = 12, 24
%     (linspace step 1/(n-1)), so the filter is a no-op there and 1728 / 13824
%     points survive -- matching the jq4_ / odd1_ production caches. Richardson
%     line 46 is the LINEAR O(1/n) form 2*S_fine - S_coarse (integrable Gamma
%     singularity => sublinear raw convergence); generalized here to an arbitrary
%     pair as X_rich = (n2*X2 - n1*X1)/(n2 - n1), which reduces to line 46 for
%     {12, 24}. NOT the O(1/n^2) squared form.
%
%   OUTPUTS
%     Tc   the Richardson-extrapolated Tc of `mode` (== out.Tc_rich).
%     out  struct: .mode, .grids, .Tc [1,ngrid], .Tc_rich, .Sc_at_Tc [1,ngrid],
%       .Sc_rich, .d_at_Tc (meV, finest grid at Tc_rich; 0 for 'off'), .nex (struct
%       of per-grid excluded-mode counts per computed variant), .timings, and
%       (mode 'full') .split as documented above.
%
%   DESIGN: chi_perp(T) is Van Vleck-dominated and nearly T-flat (invz_chiperp,
%   T1.2); it is rebuilt per root-find iteration but the geometric blocks are
%   loaded ONCE per grid, OUTSIDE the root find (P0.4). Xp uses Jxx0 = info.Jaa0
%   (the live demag-aware transverse coupling), matching the driver seam this
%   engine replaces; at B = 0 <Jx> = 0 so the choice is immaterial to ~1e-4.
%
%   See also INVZ_ODD_BLOCKS, INVZ_CHIPERP, INVZ_ODD_DELTAJ, INVZ_SIGMA_CRIT,
%   INVZ_CRITICAL_T0FIELD.

if nargin < 2, opts = struct(); end
grids = getf(opts, 'grids', {12, 24});
if iscell(grids), grids = cell2mat(grids); end
grids = grids(:).';
dpRng = getf(opts, 'dpRng', 30);
cache = getf(opts, 'cache', true);
mode  = getf(opts, 'mode', 'full');

% Variant recipes: J0 source (const = Jcc0, shift = Jcc0 - d) x Sigma source.
R = recipes();
if ~isfield(R, mode)
    error('invz:oddZeroFieldMode', 'Unknown mode ''%s''.', mode);
end
primary = mode;
% mode 'full' also needs the governing legs + diagnostics for out.split.
if strcmp(mode, 'full')
    variants = {'full', 'off', 'b_cond', 'c_sig', 'b_naive', 'c1_lit', 'c_fact'};
else
    variants = {mode};
end

ng = numel(grids);
nv = numel(variants);
Tcg  = nan(nv, ng);   Scg = nan(nv, ng);   nexg = nan(nv, ng);
tbuild = nan(1, ng);
tall = tic;

finestB = [];
for ig = 1:ng
    n = grids(ig);
    % --- Gamma-excluded uniform mesh: SAME generator + range + Gamma filter as the
    % Sigma_c benchmark (test_invz_sigma_crit.m lines 41-42); 'verbose',false replaces the
    % benchmark's evalc noise-capture -- byte-identical qvec (verbose only gates fprintf),
    % as invz_bz_couplings (the driver's shared grid builder) also does.
    [qvec, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', [n n n], 'range', [-0.5 0.5], 'verbose', false);
    qvec = qvec(any(abs(qvec) > 1e-12, 2), :);
    % --- geometric blocks ONCE per grid, outside every root find (P0.4) -------
    tb = tic;
    [Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', dpRng, 'cache', cache));
    tbuild(ig) = toc(tb);
    Jnu0 = local_modes(Vcc, zeros(4,4,size(Vcc,3)));   % T-independent Vcc modes
    B = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0, ...
        'Jaa0', infoB.Jaa0, 'Jnu0', Jnu0, 'Sc_off', sigma_crit_quiet(infoB.Jcc0, Jnu0(:)));
    if ig == ng, finestB = B; end
    for iv = 1:nv
        [Tcg(iv,ig), Scg(iv,ig), ~, nexg(iv,ig)] = ...
            solve_variant(ion, B, R.(variants{iv}));
    end
end

% --- Richardson per reported quantity (benchmark line 46 form) ---------------
ip = strcmp(variants, primary);
Tc_rich = richardson(grids, Tcg(ip, :));
Sc_rich = richardson(grids, Scg(ip, :));

if strcmp(primary, 'off')
    d_at_Tc = 0;
else
    [~, ~, d_at_Tc] = variant_eval(ion, finestB, R.(primary), Tc_rich);
end

out = struct();
out.mode     = mode;
out.grids    = grids;
out.Tc       = Tcg(ip, :);
out.Tc_rich  = Tc_rich;
out.Sc_at_Tc = Scg(ip, :);
out.Sc_rich  = Sc_rich;
out.d_at_Tc  = d_at_Tc;
out.nex      = struct();
for iv = 1:nv, out.nex.(variants{iv}) = nexg(iv, :); end
out.timings  = struct('build', tbuild, 'total', toc(tall));

if strcmp(mode, 'full')
    Tc_a    = Tc_rich;                              % (a) full == headline
    Tc_off  = richardson(grids, Tcg(strcmp(variants,'off'), :));
    Tc_b    = richardson(grids, Tcg(strcmp(variants,'b_cond'), :));
    Tc_c    = richardson(grids, Tcg(strcmp(variants,'c_sig'), :));
    out.split = struct( ...
        'Tc_a',           Tc_a, ...
        'Tc_b',           Tc_b, ...                 % governing condition-level leg
        'Tc_c',           Tc_c, ...                 % governing Sigma-level leg
        'Tc_off',         Tc_off, ...
        'closure_defect', Tc_a - (Tc_b + Tc_c - Tc_off), ...   % interaction (reported)
        'Tc_b_naive',     richardson(grids, Tcg(strcmp(variants,'b_naive'), :)), ...
        'Tc_c1_literal',  richardson(grids, Tcg(strcmp(variants,'c1_lit'), :)), ...
        'Tc_c_factorial', richardson(grids, Tcg(strcmp(variants,'c_fact'), :)), ...
        'nex', struct('b_naive', nexg(strcmp(variants,'b_naive'), :)));
end

Tc = Tc_rich;
end

% =========================================================================
function R = recipes()
%RECIPES Variant table: j0 in {'const','shift'}; sc in {'off','odd','naive_b','c1','cfact'}.
R.off          = struct('j0', 'const', 'sc', 'off');
R.full         = struct('j0', 'shift', 'sc', 'odd');      % (a)
R.b_cond       = struct('j0', 'shift', 'sc', 'off');      % governing (b): shift J0, FROZEN Sc
R.c_sig        = struct('j0', 'const', 'sc', 'odd');      % governing (c): unshifted J0, ODD Sc
R.b_naive      = struct('j0', 'shift', 'sc', 'naive_b');  % source-plan b (diagnostic)
R.c1_lit       = struct('j0', 'const', 'sc', 'c1');       % source-plan c literal == (a)
R.c_fact       = struct('j0', 'const', 'sc', 'cfact');    % factorial c (diagnostic)
R.uniform_only = R.b_naive;                               % mode alias (source-plan name)
R.qstruct_only = R.c1_lit;                                % mode alias (source-plan name)
end

% =========================================================================
function [Tc, Sc_at, d_at, nex_at] = solve_variant(ion, B, recipe)
%SOLVE_VARIANT Tc for one grid + one recipe via invz_critical_T0field, then
% evaluate Sigma_c / d / excluded-mode-count at the converged Tc.
if strcmp(recipe.j0, 'const') && strcmp(recipe.sc, 'off')
    [Sc0, J00] = variant_eval(ion, B, recipe, NaN);    % pure off: T-independent
    Tc = invz_critical_T0field(ion, Sc0, J00);
else
    % Memoize (Sc, J0) across the paired J0T/ScT calls invz_critical_T0field
    % makes at each bisection T (halves the chi_perp + deltaJ + eig work / step).
    mT = NaN;  mSc = NaN;  mJ0 = NaN;
    Tc = invz_critical_T0field(ion, @scHandle, @j0Handle);
end
[Sc_at, ~, d_at, nex_at] = variant_eval(ion, B, recipe, Tc);

    function s = scHandle(T), ensure(T);  s = mSc;  end
    function j = j0Handle(T), ensure(T);  j = mJ0;  end
    function ensure(T)
        if ~(T == mT)                      % NaN == x is false => first call recomputes
            [mSc, mJ0] = variant_eval(ion, B, recipe, T);
            mT = T;
        end
    end
end

% =========================================================================
function [Sc, J0, d, nex] = variant_eval(ion, B, recipe, T)
%VARIANT_EVAL Closed-form Sigma_c, the criticality-LHS J0, d, and the excluded-mode
% count at (recipe, T). J0 (returned, used on the criticality LHS) is distinct from
% J0sc (the argument inside invz_sigma_crit): the GOVERNING (c) leg has J0 = Jcc0 on
% the LHS but evaluates Sc at the consistent J0sc = Jcc0 - d.
if strcmp(recipe.j0, 'const') && strcmp(recipe.sc, 'off')     % pure off, no ODD
    d = 0;  J0 = B.Jcc0;  Sc = B.Sc_off;
    nex = sum((B.Jcc0 - B.Jnu0(:)) <= 1e-12);
    return;
end
Xp = invz_chiperp(ion, T, [0 0 0], struct('Jxx0', B.Jaa0));
[dJ, d] = invz_odd_deltaJ(B.Vca, B.Vcb, Xp);
if strcmp(recipe.j0, 'shift'), J0 = B.Jcc0 - d; else, J0 = B.Jcc0; end
switch recipe.sc
    case 'off'                                   % FROZEN no-ODD Sc (governing b)
        Sc = B.Sc_off;  J0sc = B.Jcc0;  Jf = B.Jnu0(:);
    case 'odd'                                   % full ODD Sc (a, governing c)
        Jf = local_modes(B.Vcc, dJ);            J0sc = B.Jcc0 - d;
        Sc = sigma_crit_quiet(J0sc, Jf);
    case 'naive_b'                               % source-plan b: Vcc modes vs Jcc0-d
        Jf = B.Jnu0(:);                          J0sc = B.Jcc0 - d;
        Sc = sigma_crit_quiet(J0sc, Jf);
    case 'c1'                                    % source-plan c literal: Vcc+dJ+d*I vs Jcc0
        Jf = local_modes(B.Vcc, dJ + d*eye(4));  J0sc = B.Jcc0;
        Sc = sigma_crit_quiet(J0sc, Jf);
    case 'cfact'                                 % factorial c: Vcc+dJ modes vs Jcc0
        Jf = local_modes(B.Vcc, dJ);             J0sc = B.Jcc0;
        Sc = sigma_crit_quiet(J0sc, Jf);
    otherwise
        error('invz:oddZeroFieldSc', 'Unknown Sc source ''%s''.', recipe.sc);
end
nex = sum((J0sc - Jf) <= 1e-12);                 % honest excluded-mode count (ODD-LOG T1.3)
end

% =========================================================================
function Jf = local_modes(Vcc, dJ)
%LOCAL_MODES Sorted real cc eigenvalues of Vcc + dJ per q, flattened [nq*4, 1].
% Mirrors invz_jq_modes' eig assembly (dJ = 0 reproduces its Jnu bitwise).
nq  = size(Vcc, 3);
Jnu = zeros(nq, 4);
for iq = 1:nq
    M = Vcc(:,:,iq) + dJ(:,:,iq);
    M = (M + M')/2;                              % both Hermitian; cleans rounding only
    Jnu(iq,:) = sort(real(eig(M))).';
end
Jf = Jnu(:);
end

% =========================================================================
function Sc = sigma_crit_quiet(J0, Jf)
%SIGMA_CRIT_QUIET invz_sigma_crit with the known invz:sigmaCritExcluded warning
% suppressed (counted explicitly by the caller). Error-safe restore (house
% pattern, invz_run_phase_diagram T1.4 seam): capture, try, restore-on-error,
% restore-on-success -- never leaves the warning OFF for the session.
ws = warning('off', 'invz:sigmaCritExcluded');
try
    Sc = invz_sigma_crit(J0, Jf);
catch err
    warning(ws);
    rethrow(err);
end
warning(ws);
end

% =========================================================================
function xr = richardson(grids, xg)
%RICHARDSON Linear O(1/n) extrapolation over the coarsest/finest grid pair,
% X_rich = (nf*Xf - nc*Xc)/(nf - nc). Reduces to test_invz_sigma_crit.m line 46
% (2*S_fine - S_coarse) for a 2:1 pair such as {12, 24}. One grid => passthrough.
if isscalar(grids)
    xr = xg(1);  return;
end
[gs, idx] = sort(grids);           % ascending: coarse -> fine
nc = gs(1);      nf = gs(end);
xc = xg(idx(1)); xf = xg(idx(end));
xr = (nf*xf - nc*xc) / (nf - nc);
end
