function out = invz_chi_realaxis(ion, T, Bx, pt, w, opts)
%INVZ_CHI_REALAXIS 1/z-renormalized cc susceptibility on the real axis.
% Bx: scalar transverse field along the crystallographic a axis (T).
% Paramagnet (HTML eqs 22-23, 29-30):
%   Sigma(w) = alpha + (M2/n01^2)[lambda1 - (1-n01^2)K(w)] g(w)
% Ordered phase (pt.is_ordered = true; HTML eq 37), with gamma(w) as above:
%   Sigma(w) = (alpha - alpha_m) + [gamma(w) - (2 m^2/M^2) gamma(0)] g(w), gamma(0) at static Matsubara K(0)
% chi~0(w) = chi0_cc(w)/(1+Sigma(w));  chi(q,nu,w) = chi~0/(1 - J_nu chi~0)  (HTML eq 29-30).
%
% For an ordered point the single-ion response must come from the ORDERED eigenstates
% (pass opts.si or rely on pt.si) so chi0_cc matches the symmetry-broken state; alpha,
% alpha_m, lambda1, K(0) come from the converged Matsubara solve (pt).
%
% opts.Jxx0   (ion.Jxx0)  transverse MF coupling for the internally built single-ion state.
% opts.hyp    (true)      hyperfine manifold; must match the Matsubara medium's Hilbert space.
% opts.Jshape (0)         strict-uniform demag correction (info.Jshape_cc), applied to chi_cc_q;
%                         leave 0 for finite-q (intrinsic) paths.
% opts.chi0cc_w           precomputed chi0_cc(w) on this exact grid; when supplied, skips the
%                         single-ion diagonalization so a field point can share it across evaluations.
% opts.bare_rpa (false)   build an INDEPENDENT bare-MF state for a Sigma=0 RPA response:
%                         the PM mass 1-J0eff*chi0cc(0) selects the PM/ordered branch, and
%                         the ordered branch solves hz=J0eff*<Jz>. Requires scalar opts.J0eff.
%                         The selector is local to this response owner so projected spectra
%                         cannot accidentally reuse a 1/z-selected state for their RPA leg.
if nargin < 6, opts = struct(); end
if ~(isnumeric(Bx) && isreal(Bx) && isscalar(Bx) && isfinite(Bx))
    error('invz:field', 'Bx must be a finite real scalar transverse field.');
end
Bvec = [Bx 0 0];
eta    = getf(opts, 'eta', 5e-3);
npass  = getf(opts, 'npass', 3);
Jsel   = getf(opts, 'Jsel', ion.J0eff);
Jxx0   = getf(opts, 'Jxx0', ion.Jxx0);
Jshape = getf(opts, 'Jshape', 0);
hyp    = getf(opts, 'hyp', true);
bareRpa = getf(opts, 'bare_rpa', false);
if bareRpa
    if isfield(opts, 'chi0cc_w') && ~isempty(opts.chi0cc_w)
        error('invz:rpaStateReuse', ...
            'opts.bare_rpa cannot reuse chi0cc_w from another (e.g. 1/z-selected) state.');
    end
    [pt, rpaPhase, rpaMassPm] = bare_rpa_point(ion, T, Bx, opts);
    opts.si = pt.si;                                  % override any 1/z-state carrier
end
ordered = isfield(pt,'is_ordered') && pt.is_ordered;

z  = w(:) + 1i*eta;
if isfield(opts,'chi0cc_w') && ~isempty(opts.chi0cc_w)
    G0 = -opts.chi0cc_w(:);                        % reuse a shared bare cc; no si/chi0z needed
else
    if isfield(opts,'si') && ~isempty(opts.si)
        si = opts.si;                              % caller-supplied single-ion state
    elseif ordered && isfield(pt,'si') && ~isempty(pt.si)
        si = pt.si;                                % ordered eigenstates from the solve
    else
        si = invz_single_ion(ion, T, Bvec, struct('hyp', hyp, 'Jxx0', Jxx0));
    end
    c0 = invz_chi0z(si, T, z, struct('elastic', false));
    G0 = -squeeze(c0(3,3,:));
end
tl = pt.tl;
g  = invz_g(tl, z);
pref = tl.M2/tl.n01^2;
if isfield(pt,'K') && ~isempty(pt.K)
    Kw = pt.K(1)*ones(size(z));                    % seed with static Matsubara K
    K0 = pt.K(1);
else
    Kw = zeros(size(z));                           % RPA-limit callers may omit K
    K0 = 0;
end
for pass = 1:npass
    Sw = realaxis_sigma(pt, tl, pref, Kw, K0, g, ordered);
    if isfield(opts,'Jfull') && ~isempty(opts.Jfull)
        med = invz_emt_scalar(G0, Sw, opts.Jfull, struct('max_iter', 100, 'tol', 1e-8));
        Kw = med.K;
    else
        break;                                     % no lattice pass requested: single shot
    end
end
Sw = realaxis_sigma(pt, tl, pref, Kw, K0, g, ordered);
chit = (-G0) ./ (1 + Sw);
out.Sigma_w  = Sw;
out.chi0cc_w = -G0;
out.chi_cc_q = zeros(numel(Jsel), numel(z));
for k = 1:numel(Jsel)
    out.chi_cc_q(k,:) = (chit ./ (1 - Jsel(k)*chit)).';
end
if Jshape ~= 0
    % Sample-shape correction for the STRICT-UNIFORM measured observable only:
    % chi_meas = chi_int/(1 + Jshape*chi_int)  (demag-limited: the soft mode
    % saturates at 1/Jshape instead of diverging). Callers evaluating a finite-q
    % path (intrinsic longitudinal probe) must NOT pass Jshape.
    out.chi_cc_q = out.chi_cc_q ./ (1 + Jshape*out.chi_cc_q);
end
out.Jshape = Jshape;
if bareRpa
    out.phase_rpa = rpaPhase;                         % 1 = bare-MF ordered, 2 = stable PM
    out.rpa_mass_pm = rpaMassPm;                      % PM-side diagnostic at this field
end
end

function Sw = realaxis_sigma(pt, tl, pref, Kw, K0, g, ordered)
%REALAXIS_SIGMA Sigma(w) for the paramagnet (HTML eq 22) or ordered phase (HTML eq 37).
gamma_w = pref*(pt.lambda(1) - (1 - tl.n01^2)*Kw);
if ordered
    gamma0 = pref*(pt.lambda(1) - (1 - tl.n01^2)*K0);
    if tl.M2 == 0
        % Exact M2 cancellation in (2*m^2/M2)*gamma(0). Restrict the
        % reassociation to the exact endpoint so positive-M2 spectra keep
        % their historical arithmetic bit-for-bit.
        Q0 = (2*tl.m^2/tl.n01^2) * ...
            (pt.lambda(1) - (1-tl.n01^2)*K0);
    else
        Q0 = (2*tl.m^2/tl.M2)*gamma0;
    end
    Sw = (pt.alpha - pt.alpha_m) + (gamma_w - Q0) .* g;
else
    Sw = pt.alpha + gamma_w .* g;
end
end

function [pt, phase, massPm] = bare_rpa_point(ion, T, Bx, opts)
%BARE_RPA_POINT Independent bare-MF state selected by the RPA theory's own static mass.
% This deliberately does not inspect a 1/z point or phase label. The PM state decides which
% side of the bare transition the field occupies; only the ordered side pays for the
% longitudinal MF solve. The same uniform scalar J0eff is used for that solve and for the
% Gamma RPA denominator; q-path callers may additionally supply dispersive finite-q Jsel.
J0eff = getf(opts, 'J0eff', []);
if ~(isnumeric(J0eff) && isreal(J0eff) && isscalar(J0eff) && isfinite(J0eff))
    error('invz:rpaJ0eff', 'opts.bare_rpa requires a finite real scalar opts.J0eff.');
end
Jxx0 = getf(opts, 'Jxx0', ion.Jxx0);
hyp   = getf(opts, 'hyp', true);
[si, phase, massPm] = invz_bare_mf_state( ...
    ion, T, Bx, J0eff, Jxx0, hyp);
if phase == 1
    tl = invz_twolevel_ordered(ion, T, Bx, si.hz, struct('Jxx0', Jxx0));
else
    tl = invz_twolevel(ion, T, Bx, struct('Jxx0', Jxx0));
end

pt = struct('alpha', 0, 'alpha_m', 0, 'lambda', [0; 0; 0], ...
    'K', [], 'tl', tl, 'si', si, 'is_ordered', phase == 1);
end
