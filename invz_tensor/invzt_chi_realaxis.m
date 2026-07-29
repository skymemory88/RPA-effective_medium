function out = invzt_chi_realaxis(ion, T, B, pt, w, opts)
%INVZT_CHI_REALAXIS Tensor real-axis spectra: A1 scalar-Sigma continuation.
%   out = INVZT_CHI_REALAXIS(ion, T, B, pt, w, opts) continues the converged
%   A1 fixed point pt (INVZT_SOLVE_POINT, PM, or INVZT_SOLVE_POINT_ORDERED,
%   ordered) to the real frequency axis z = w + 1i*eta (elastic OFF) via the
%   SAME scalar self-energy continuation pattern as the projected
%   INVZ_CHI_REALAXIS: PM continues its non-ordered branch bit-for-bit;
%   ordered continues the projected HTML-eq-37 moment form (both frozen
%   static K, Kw seeded at Gamma):
%
%     PM (pt.is_ordered false/absent):
%       Sigma_w = pt.alpha + pref*(pt.lambda(1) - (1-n01^2)*Kw) .* g(z)
%     ordered (pt.is_ordered true, INVZT_SOLVE_POINT_ORDERED's pt):
%       Sigma_w = (pt.alpha - pt.alpha_m)
%                 + (gamma_w - (2*tl.m^2/tl.M2)*gamma0) .* g(z),
%       gamma_w = pref*(pt.lambda(1) - (1-n01^2)*Kw),
%       gamma0  = pref*(pt.lambda(1) - (1-n01^2)*pt.K(1))
%       (with the frozen static Kw = pt.K(1) seeding used here, gamma_w ==
%       gamma0 identically; both symbols are kept for parity with the
%       projected realaxis_sigma, whose npass/Jfull re-solve makes them
%       differ -- a possible future port, not done here)
%
%   pref = pt.tl.M2/pt.tl.n01^2, Kw = pt.K(1) (the converged static/n=0 K,
%   FROZEN across the whole w grid -- not re-solved on the real axis; see
%   SCOPE below), g(z) = INVZ_G(pt.tl, z). This is bit-for-bit the projected
%   Kw-seeding pattern (INVZ_CHI_REALAXIS's "no opts.Jfull" single-shot path).
%
%   The renormalized LOCAL tensor response ctil(w) [3,3,nw] is REBUILT (not
%   reused from pt, which only carries it on the Matsubara grid), at
%   z = w+1i*eta with elastic OFF, by a branch-dependent rule:
%
%     PM:       [cdom, crest] = INVZT_CHI0_SPLIT(pt.si, T, z, ...);
%               ctil(:,:,k) = cdom(:,:,k)/(1+Sigma_w(k)) + crest(:,:,k)
%               (the SAME dominant/rest rule used inside INVZT_SOLVE_POINT;
%               crest zeroed when chi_rest is false, matching the converged
%               point's own pt.chi_rest by default).
%     ordered:  c0w = INVZ_CHI0Z(pt.si, T, z, struct('elastic', false));
%               ctil(:,:,k) = c0w(:,:,k) / (1+Sigma_w(k))
%               (WHOLE-CC rebuild, no dominant/rest split -- 2026-07-20
%               amendment, matching INVZT_SOLVE_POINT_ORDERED's own WHOLE-CC
%               Matsubara medium; the exact tensor analog of the projected
%               ordered chit = -G0/(1+Sigma_w)).
%
%   ctil is then paged through INVZT_CHI_RPA against a q-selection
%   (opts.qsel) -- this stage is IDENTICAL for both branches:
%
%     'gamma_uniform' (DEFAULT) -- the Task-4 uniform page
%         Jd = kron(ones(4)/4, diag([Jaa0, Jaa0, Jcc0]))
%       built from pt.lat.info (NOT ion.Jxx0/ion.J0eff -- the live lattice-sum
%       values threaded through the solve). The S4-uniform-mode projector
%       u = kron(ones(4,1)/2, eye(3)) recovers the uniform-mode response
%       out.chi_uniform = u'*X*u for ANY chi0, by the SAME exact linear-
%       algebra identity (resolvent
%       identity (I-AB)^-1*A = A*(I-BA)^-1 plus the rank-1 sublattice
%       projector ones(4)/4 = v*v', v = ones(4,1)/2) -- this is what makes the
%       force_sigma0 exact-identity gate below hold to numerical precision,
%       not physics tuning.
%     'gamma' -- uses the ACTUAL pt.lat.JtGamma page (still projected through
%       the same u onto out.chi_uniform). This differs from 'gamma_uniform'
%       whenever JtGamma's other q=0 sublattice eigenmodes (nonzero-coupling
%       "optical" branches, if any) RPA-mix into the S4-uniform channel --
%       'gamma_uniform' isolates that channel exactly (Jd's only nonzero
%       sublattice eigenvalue is the S4-uniform one), 'gamma' does not.
%     explicit [nq,3] qvec -- rebuilds lat via INVZT_JQ_TENSOR at opts.dpRng
%       (own cache namespace; see INVZT_JQ_TENSOR) and ALSO returns the
%       per-q, sublattice-averaged site-diagonal cc response out.chi_cc_q
%       [nq,nw], COMPLEX -- imag() is chi'' (positive up to the frozen-Kw
%       caveat below) (contract parity with the projected INVZ_CHI_REALAXIS's
%       chi_cc_q). The sublattice diag4 mean mirrors INVZT_GCC_LATTICE's
%       pattern but WITHOUT its real(): that projection is Matsubara-axis-only
%       and would delete the whole dissipative part here (fixed 2026-07-18,
%       Codex review F1).
%       out.chi_uniform is still populated via the SAME Jd uniform-mode
%       construction (a q-grid-independent quantity, always available).
%
%   Returns:
%     out.chi_uniform [3,3,nw] : S4-uniform-mode tensor response, ALL qsel.
%     out.chi_cc_q    [nq,nw]  : COMPLEX per-q sublattice-averaged cc response
%                                (imag = chi''; NB the full-Sigma frozen-Kw
%                                continuation carries a known near-resonance
%                                negative-chi'' artifact -- shared with
%                                chi_uniform and the projected reference
%                                scheme, e.g. -313 beside a +652 peak at T=1.6
%                                K/B=2 T -- exact non-negativity holds in the
%                                force_sigma0 bare-RPA limit; see
%                                test_qsel_explicit_q_complex_response);
%                                populated ONLY for an explicit qvec (empty []
%                                for 'gamma_uniform'/'gamma').
%     out.Sigma_w     [nw,1]   : the A1 scalar self-energy continuation.
%
%   B : scalar (transverse-along-a, historical) or [Bx By Bz] (T); validated
%       (INVZ_FIELD_VEC) for interface parity with INVZT_SOLVE_POINT/
%       INVZT_SOLVE_POINT_ORDERED and the projected INVZ_CHI_REALAXIS, but NOT
%       otherwise consumed here (PM or ordered alike): pt.si/pt.tl already
%       encode the field-dependent physics from the original solve, so no
%       fresh field-dependent single-ion rebuild is needed on the real axis.
%
%   OPTIONS (getf defaults):
%     eta          5e-3            : Im(z) = w + 1i*eta.
%     qsel         'gamma_uniform' : 'gamma_uniform' | 'gamma' | [nq,3] qvec.
%     dpRng        30              : explicit-qvec branch only (INVZT_JQ_TENSOR).
%     cache        true            : explicit-qvec branch only (INVZT_JQ_TENSOR).
%     chi_rest     pt.chi_rest     : PM branch only -- false zeroes crest
%                                    (matches the solve by default; override
%                                    only for diagnostics). Ignored on the
%                                    ordered branch (WHOLE-CC rebuild, no
%                                    dominant/rest split; pt.chi_rest is
%                                    always true for an ordered pt, so the
%                                    default read is safe but unused there).
%     dominant_count pt.mspec.ndom : PM branch fixed-rank dominant selector.
%     Esplit       absent          : PM branch LEGACY energy-cut override.
%                                    With neither option, the point's recorded
%                                    selection is reused exactly. Both are ignored
%                                    on the ordered branch (no split).
%     odd          pt.odd          : must equal pt.odd -- mismatch errors
%                                    'invzt:oddMismatch' (Sigma_w/alpha/lambda/K
%                                    already bake in the odd flag used at solve
%                                    time; omitting opts.odd skips the check).
%                                    odd=false ALSO applies INVZT_ODD_MASK's
%                                    Cartesian-off-diagonal mask to the coupling
%                                    tensor consumed at finite q ('gamma' and
%                                    explicit-qvec Jpage) -- the SAME geometric
%                                    rule INVZT_SOLVE_POINT applies to lat.Jt
%                                    before solving, not merely a bookkeeping
%                                    check (R1, 2026-07-18 second review; a
%                                    pre-fix odd=false point got ODD-on
%                                    couplings at finite q). 'gamma_uniform' is
%                                    unaffected: its Jd page is Cartesian-
%                                    diagonal by construction.
%     force_sigma0 false           : TEST HOOK. Forces Sigma_w = 0 identically
%                                    (bypasses pt.alpha/lambda/K entirely),
%                                    isolating the bare-chi0 RPA limit used by
%                                    the exact-identity gate.
%
%   SCOPE (v2, LOCKED): A1 scalar-Sigma continuation, PM AND ordered
%   branches (Task 4, 2026-07-20); A2/A3 continuation remains the open item:
%   continuing the full Vmat(i*omega_n) needs either direct real-axis kernel
%   evaluation (the tensor kernels accept complex frequency arguments) or a
%   fitted continuation -- a separate future work item
%   (historical tensor Tasks 8 and 12; see invzp_convg_diagnosis.md Section 2.1).
%
%   See also INVZT_SOLVE_POINT, INVZT_SOLVE_POINT_ORDERED, INVZT_CHI_RPA,
%   INVZT_CHI0_SPLIT, INVZ_CHI0Z, INVZT_JQ_TENSOR, INVZT_GCC_LATTICE,
%   INVZ_CHI_REALAXIS (projected Kw-seeding reference / interop peak-parity
%   target).
if nargin < 6, opts = struct(); end
eta      = getf(opts, 'eta', 5e-3);
qsel     = getf(opts, 'qsel', 'gamma_uniform');
dpRng    = getf(opts, 'dpRng', 30);
cacheq   = getf(opts, 'cache', true);
chi_rest = getf(opts, 'chi_rest', pt.chi_rest);
odd_req  = getf(opts, 'odd', pt.odd);
force_sigma0 = isfield(opts, 'force_sigma0') && ~isempty(opts.force_sigma0) ...
    && ~isequal(opts.force_sigma0, false);

ptmode = getf(pt, 'mode', 'a1');
if ~strcmp(ptmode, 'a1')
    error('invzt:realaxisMode', ['invzt_chi_realaxis is the A1 scalar-Sigma ' ...
        'continuation ONLY (LOCKED scope; see the SCOPE box in this header): ' ...
        'got pt.mode = ''%s''. Re-solve the point with mode ''a1''. Continuing ' ...
        'the A2/A3 matrix objects is an open item (README section 10).'], ptmode);
end

if ~isequal(odd_req, pt.odd)
    error('invzt:oddMismatch', ['opts.odd (%s) must equal pt.odd (%s): ' ...
        'invzt_chi_realaxis continues the SAME converged A1 point, whose ' ...
        'Sigma/alpha/lambda/K already bake in the odd flag used at solve ' ...
        'time.'], invzt_str(odd_req), invzt_str(pt.odd));
end
B = invz_field_vec(B);   %#ok<NASGU> % validated for interface parity only; see header

w  = w(:);  nw = numel(w);
z  = w + 1i*eta;

% --- A1 scalar-Sigma continuation (Kw-seeding pattern; PM = projected non-ordered
%     branch, ordered = projected HTML-eq-37 moment form, both frozen static K) ---
ordered = isfield(pt, 'is_ordered') && pt.is_ordered;
tl   = pt.tl;
g    = invz_g(tl, z);
pref = tl.M2 / tl.n01^2;
Kw   = pt.K(1) * ones(nw, 1);
if force_sigma0
    Sw = zeros(nw, 1);
elseif ordered
    gamma_w = pref*(pt.lambda(1) - (1 - tl.n01^2)*Kw);
    gamma0  = pref*(pt.lambda(1) - (1 - tl.n01^2)*pt.K(1));
    Sw = (pt.alpha - pt.alpha_m) + (gamma_w - (2*tl.m^2/tl.M2)*gamma0) .* g;
else
    Sw = pt.alpha + pref*(pt.lambda(1) - (1 - tl.n01^2)*Kw) .* g;
end

% --- rebuild ctil(w) on the real axis (elastic OFF) ---------------------------
if ordered
    % WHOLE-CC rebuild (2026-07-20 amendment, matching Task 3's whole-cc medium):
    % no dominant/rest split -- the exact tensor analog of the projected ordered
    % chit = -G0/(1+Sigma_w).
    c0w  = invz_chi0z(pt.si, T, z, struct('elastic', false));   % full local chi0, elastic OFF
    ctil = c0w ./ reshape(1 + Sw, 1, 1, nw);                    % Sw [nw,1] -> [1,1,nw] broadcast
else
    split_opts = struct('elastic', false);
    if isfield(opts, 'dominant_count') && ~isempty(opts.dominant_count)
        split_opts.dominant_count = opts.dominant_count;
    elseif isfield(opts, 'Esplit') && ~isempty(opts.Esplit)
        split_opts.Esplit = opts.Esplit;
    elseif isfield(pt, 'mspec') && strcmp(getf(pt.mspec, 'selection', ''), 'fixed_rank')
        split_opts.dominant_count = pt.mspec.ndom;
    elseif isfield(pt, 'mspec') && strcmp(getf(pt.mspec, 'selection', ''), 'energy')
        split_opts.Esplit = pt.mspec.Esplit;
    end
    [cdom, crest] = invzt_chi0_split(pt.si, T, z, split_opts);
    if ~chi_rest
        crest = zeros(size(crest));
    end
    ctil = cdom ./ reshape(1 + Sw, 1, 1, nw) + crest;       % Sw [nw,1] -> [1,1,nw] broadcast
end

% --- q-selection: page invzt_chi_rpa, project the S4-uniform mode ------------
Jaa0 = pt.lat.info.Jaa0;
Jcc0 = pt.lat.info.Jcc0;
u  = kron(ones(4,1)/2, eye(3));                          % isometric uniform-mode embedding
Jd = kron(ones(4)/4, diag([Jaa0, Jaa0, Jcc0]));          % Task-4 uniform page

out.chi_uniform = zeros(3, 3, nw);
out.chi_cc_q = [];
% Resolve the S4-uniform-mode page ONCE: Jd for 'gamma_uniform' AND the explicit-qvec
% branch (the uniform mode is a q-grid-independent quantity), JtGamma only for 'gamma'.
% The chi_uniform loop then runs ONCE; the explicit-qvec branch adds the extra chi_cc_q.
Jpage = Jd;
explicitq = false;
if ischar(qsel) || isstring(qsel)
    qsel = char(qsel);
    switch qsel
        case 'gamma_uniform'
            Jpage = Jd;                              % Cartesian-diagonal by construction; no mask
        case 'gamma'
            Jpage = pt.lat.JtGamma;
            if ~pt.odd
                % Same odd=false rule INVZT_SOLVE_POINT applies to lat.Jt before solving
                % (R1, 2026-07-18 second review), applied here for uniformity -- a no-op
                % at Gamma to numerical precision: invzt_solve_point asserts JtGamma's ODD
                % (c<->a,b) blocks vanish there by C2 symmetry; the a<->b block has no
                % runtime assert but is measured ~1e-19, so masking changes nothing.
                Jpage = invzt_odd_mask(Jpage);
            end
        otherwise
            error('invzt:qsel', ['opts.qsel must be ''gamma_uniform'', ''gamma'', ' ...
                'or an [nq,3] numeric qvec; got string ''%s''.'], qsel);
    end
else
    if ~(isnumeric(qsel) && ismatrix(qsel) && size(qsel, 2) == 3)
        error('invzt:qsel', ['opts.qsel must be ''gamma_uniform'', ''gamma'', ' ...
            'or an [nq,3] numeric qvec; got a %s.'], invzt_str(qsel));
    end
    explicitq = true;
    qvec = qsel;
    nq   = size(qvec, 1);
    latq = invzt_jq_tensor(ion, qvec, struct('dpRng', dpRng, 'cache', cacheq));
    if ~pt.odd
        % R1 fix (2026-07-18 second Codex review): the pre-fix code consumed
        % latq.Jt UNMASKED here, so an odd=false point got ODD-on couplings at
        % finite q (17.2% response error measured by the review). Apply the
        % SAME Cartesian-off-diagonal mask INVZT_SOLVE_POINT applies to lat.Jt
        % before solving.
        latq.Jt = invzt_odd_mask(latq.Jt);
    end
    out.chi_cc_q = complex(zeros(nq, nw));
end
for k = 1:nw
    X = invzt_chi_rpa(ctil(:,:,k), Jpage);
    out.chi_uniform(:,:,k) = u' * X * u;
    if explicitq
        Xq = invzt_chi_rpa(ctil(:,:,k), latq.Jt);        % [12,12,nq]
        for iq = 1:nq
            acc = 0;
            for s = 1:4
                acc = acc + Xq(3*(s-1)+3, 3*(s-1)+3, iq);
            end
            out.chi_cc_q(iq, k) = acc / 4;
        end
    end
end
out.Sigma_w = Sw;
end
