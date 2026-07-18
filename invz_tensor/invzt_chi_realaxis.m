function out = invzt_chi_realaxis(ion, T, B, pt, w, opts)
%INVZT_CHI_REALAXIS Tensor real-axis spectra: A1 scalar-Sigma continuation ONLY.
%   out = INVZT_CHI_REALAXIS(ion, T, B, pt, w, opts) continues the converged
%   A1 fixed point pt (INVZT_SOLVE_POINT) to the real frequency axis
%   z = w + 1i*eta (elastic OFF) via the SAME scalar self-energy continuation
%   as the projected INVZ_CHI_REALAXIS's paramagnetic branch (A1 tensor solves
%   are PM-only -- INVZT_SOLVE_POINT never returns an is_ordered point -- so
%   there is no ordered-phase branch here, unlike the projected reference):
%
%       Sigma_w = pt.alpha + pref*(pt.lambda(1) - (1-n01^2)*Kw) .* g(z)
%
%   pref = pt.tl.M2/pt.tl.n01^2, Kw = pt.K(1) (the converged static/n=0 K,
%   FROZEN across the whole w grid -- not re-solved on the real axis; see
%   SCOPE below), g(z) = INVZ_G(pt.tl, z). This is bit-for-bit the projected
%   Kw-seeding pattern (INVZ_CHI_REALAXIS's "no opts.Jfull" single-shot path).
%
%   The renormalized LOCAL tensor response ctil(w) [3,3,nw] is REBUILT (not
%   reused from pt, which only carries it on the Matsubara grid) from a fresh
%   INVZT_CHI0_SPLIT(pt.si, T, z, ...) at z = w+1i*eta with elastic OFF, and
%   combined with the SAME dominant/rest rule used inside INVZT_SOLVE_POINT:
%
%       ctil(:,:,k) = cdom(:,:,k)/(1+Sigma_w(k)) + crest(:,:,k)
%
%   (crest zeroed when chi_rest is false, matching the converged point's own
%   pt.chi_rest by default). ctil is then paged through INVZT_CHI_RPA against
%   a q-selection (opts.qsel):
%
%     'gamma_uniform' (DEFAULT) -- the Task-4 uniform page
%         Jd = kron(ones(4)/4, diag([Jaa0, Jaa0, Jcc0]))
%       built from pt.lat.info (NOT ion.Jxx0/ion.J0eff -- the live lattice-sum
%       values threaded through the solve). The S4-uniform-mode projector
%       u = kron(ones(4,1)/2, eye(3)) recovers the uniform-mode response
%       out.chi_uniform = u'*X*u for ANY chi0, by the SAME exact linear-
%       algebra identity proven in interop/test_invzt_rpa_parity.m (resolvent
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
%       [nq,nw] (mirrors INVZT_GCC_LATTICE's per-q diag4 mean, evaluated one
%       q at a time rather than BZ-averaged together). out.chi_uniform is
%       still populated via the SAME Jd uniform-mode construction (a
%       q-grid-independent quantity, always available).
%
%   Returns:
%     out.chi_uniform [3,3,nw] : S4-uniform-mode tensor response, ALL qsel.
%     out.chi_cc_q    [nq,nw]  : per-q sublattice-averaged cc response;
%                                populated ONLY for an explicit qvec (empty
%                                [] for 'gamma_uniform'/'gamma').
%     out.Sigma_w     [nw,1]   : the A1 scalar self-energy continuation.
%
%   B : scalar (transverse-along-a, historical) or [Bx By Bz] (T); validated
%       (INVZ_FIELD_VEC) for interface parity with INVZT_SOLVE_POINT and the
%       projected INVZ_CHI_REALAXIS, but NOT otherwise consumed here: pt.si /
%       pt.tl already encode the field-dependent physics from the original
%       solve (mode 'a1' is PM-only, so there is no ordered-phase branch that
%       would need a fresh field-dependent single-ion rebuild).
%
%   OPTIONS (getf defaults):
%     eta          5e-3            : Im(z) = w + 1i*eta.
%     qsel         'gamma_uniform' : 'gamma_uniform' | 'gamma' | [nq,3] qvec.
%     dpRng        30              : explicit-qvec branch only (INVZT_JQ_TENSOR).
%     cache        true            : explicit-qvec branch only (INVZT_JQ_TENSOR).
%     chi_rest     pt.chi_rest     : false zeroes crest (matches the solve by
%                                    default; override only for diagnostics).
%     Esplit       0.4653          : dominant/rest split energy (meV), passed
%                                    to INVZT_CHI0_SPLIT (matches
%                                    INVZT_SOLVE_POINT's own default).
%     odd          pt.odd          : must equal pt.odd -- mismatch errors
%                                    'invzt:oddMismatch' (Sigma_w/alpha/lambda/K
%                                    already bake in the odd flag used at solve
%                                    time; omitting opts.odd skips the check).
%     force_sigma0 false           : TEST HOOK. Forces Sigma_w = 0 identically
%                                    (bypasses pt.alpha/lambda/K entirely),
%                                    isolating the bare-chi0 RPA limit used by
%                                    the exact-identity gate.
%
%   SCOPE (v2, LOCKED): this is the A1 scalar-Sigma analytic continuation
%   ONLY. It does NOT extend to A3: continuing the full Vmat(i*omega_n) needs
%   either direct real-axis kernel evaluation (the tensor kernels accept
%   complex frequency arguments) or a fitted continuation -- a separate future
%   work item (docs/superpowers/plans/2026-07-17-invz-tensor-full.md Task 8,
%   Task 12).
%
%   See also INVZT_SOLVE_POINT, INVZT_CHI_RPA, INVZT_CHI0_SPLIT,
%   INVZT_JQ_TENSOR, INVZT_GCC_LATTICE, INVZ_CHI_REALAXIS (projected
%   Kw-seeding reference / interop peak-parity target).
if nargin < 6, opts = struct(); end
eta      = getf(opts, 'eta', 5e-3);
qsel     = getf(opts, 'qsel', 'gamma_uniform');
dpRng    = getf(opts, 'dpRng', 30);
cacheq   = getf(opts, 'cache', true);
chi_rest = getf(opts, 'chi_rest', pt.chi_rest);
Esplit   = getf(opts, 'Esplit', 0.4653);
odd_req  = getf(opts, 'odd', pt.odd);
force_sigma0 = isfield(opts, 'force_sigma0') && ~isempty(opts.force_sigma0) ...
    && ~isequal(opts.force_sigma0, false);

if ~isequal(odd_req, pt.odd)
    error('invzt:oddMismatch', ['opts.odd (%s) must equal pt.odd (%s): ' ...
        'invzt_chi_realaxis continues the SAME converged A1 point, whose ' ...
        'Sigma/alpha/lambda/K already bake in the odd flag used at solve ' ...
        'time.'], invzt_str(odd_req), invzt_str(pt.odd));
end
B = invz_field_vec(B);   %#ok<NASGU> % validated for interface parity only; see header

w  = w(:);  nw = numel(w);
z  = w + 1i*eta;

% --- A1 scalar-Sigma continuation (Kw-seeding pattern, mirrors the projected
%     INVZ_CHI_REALAXIS's non-ordered branch exactly) --------------------------
tl   = pt.tl;
g    = invz_g(tl, z);
pref = tl.M2 / tl.n01^2;
Kw   = pt.K(1) * ones(nw, 1);
if force_sigma0
    Sw = zeros(nw, 1);
else
    Sw = pt.alpha + pref*(pt.lambda(1) - (1 - tl.n01^2)*Kw) .* g;
end

% --- rebuild ctil(w) on the real axis (elastic OFF), dominant/rest split ------
[cdom, crest] = invzt_chi0_split(pt.si, T, z, struct('Esplit', Esplit, 'elastic', false));
if ~chi_rest
    crest = zeros(size(crest));
end
ctil = cdom ./ reshape(1 + Sw, 1, 1, nw) + crest;       % Sw [nw,1] -> [1,1,nw] broadcast

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
            Jpage = Jd;
        case 'gamma'
            Jpage = pt.lat.JtGamma;
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
    out.chi_cc_q = zeros(nq, nw);
end
for k = 1:nw
    X = invzt_chi_rpa(ctil(:,:,k), Jpage);
    out.chi_uniform(:,:,k) = u' * X * u;
    if explicitq
        Xq = invzt_chi_rpa(ctil(:,:,k), latq.Jt);        % [12,12,nq]
        for iq = 1:nq
            acc = 0;
            for s = 1:4
                acc = acc + real(Xq(3*(s-1)+3, 3*(s-1)+3, iq));
            end
            out.chi_cc_q(iq, k) = acc / 4;
        end
    end
end
out.Sigma_w = Sw;
end
