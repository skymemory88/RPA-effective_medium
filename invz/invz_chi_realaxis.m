function out = invz_chi_realaxis(ion, T, Bx, pt, w, opts)
%INVZ_CHI_REALAXIS 1/z-renormalized cc susceptibility on the real axis.
% Paramagnet (HTML eqs 22-23, 29-30):
%   Sigma(w) = alpha + (M2/n01^2)[lambda1 - (1-n01^2)K(w)] g(w)
% Ordered phase (pt.is_ordered = true; HTML eq 37):
%   Sigma(w) = (alpha - alpha_m) + [gamma(w) - (2 m^2/M^2) gamma(0)] g(w)
% with gamma(w) = (M2/n01^2)[lambda1 - (1-n01^2)K(w)], gamma(0) at the static Matsubara K(0).
% Then chi~0(w) = chi0_cc(w)/(1+Sigma(w));  chi(q,nu,w) = chi~0/(1 - J_nu chi~0)  (HTML eq 29-30).
%
% For an ordered point the single-ion response must come from the ORDERED eigenstates: pass
% opts.si (or rely on pt.si), so chi0_cc is evaluated in the symmetry-broken state. alpha, alpha_m,
% lambda1 and K(0) are fixed by the converged Matsubara solve (pt).
if nargin < 6, opts = struct(); end
eta   = 5e-3; if isfield(opts,'eta'),   eta   = opts.eta;   end
npass = 3;    if isfield(opts,'npass'), npass = opts.npass; end
Jsel  = ion.J0eff; if isfield(opts,'Jsel'), Jsel = opts.Jsel; end
ordered = isfield(pt,'is_ordered') && pt.is_ordered;

if isfield(opts,'si') && ~isempty(opts.si)
    si = opts.si;                                  % caller-supplied single-ion state
elseif ordered && isfield(pt,'si') && ~isempty(pt.si)
    si = pt.si;                                    % ordered eigenstates from the solve
else
    si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', true));   % paramagnet
end
z  = w(:) + 1i*eta;
c0 = invz_chi0z(si, T, z, struct('elastic', false));
G0 = -squeeze(c0(3,3,:));
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
