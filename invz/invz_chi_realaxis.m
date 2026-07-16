function out = invz_chi_realaxis(ion, T, Bx, pt, w, opts)
%INVZ_CHI_REALAXIS 1/z-renormalized cc susceptibility on the real axis.
% Bx: scalar (transverse, historical) or [Bx By Bz] vector (T).
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
if nargin < 6, opts = struct(); end
eta   = 5e-3; if isfield(opts,'eta'),   eta   = opts.eta;   end
npass = 3;    if isfield(opts,'npass'), npass = opts.npass; end
Jsel  = ion.J0eff; if isfield(opts,'Jsel'), Jsel = opts.Jsel; end
Jxx0   = ion.Jxx0; if isfield(opts,'Jxx0'),   Jxx0   = opts.Jxx0;   end
Jshape = 0;        if isfield(opts,'Jshape'), Jshape = opts.Jshape; end
hyp    = true;     if isfield(opts,'hyp'),    hyp    = opts.hyp;    end
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
        si = invz_single_ion(ion, T, invz_field_vec(Bx), struct('hyp', hyp, 'Jxx0', Jxx0));   % paramagnet
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
end

function Sw = realaxis_sigma(pt, tl, pref, Kw, K0, g, ordered)
%REALAXIS_SIGMA Sigma(w) for the paramagnet (HTML eq 22) or ordered phase (HTML eq 37).
gamma_w = pref*(pt.lambda(1) - (1 - tl.n01^2)*Kw);
if ordered
    gamma0 = pref*(pt.lambda(1) - (1 - tl.n01^2)*K0);
    Sw = (pt.alpha - pt.alpha_m) + (gamma_w - (2*tl.m^2/tl.M2)*gamma0) .* g;
else
    Sw = pt.alpha + gamma_w .* g;
end
end
