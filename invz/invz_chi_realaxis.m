function out = invz_chi_realaxis(ion, T, Bx, pt, w, opts)
%INVZ_CHI_REALAXIS 1/z-renormalized cc susceptibility on the real axis.
% Sigma(w) = alpha + (M2/n01^2)*[lambda1 - (1-n01^2)*K(w)] * g(w)   (HTML eqs 22-23),
% alpha and lambda1 fixed by the converged Matsubara solve (pt).
% chi~0(w) = chi0_cc(w)/(1+Sigma(w));  chi(q,nu,w) = chi~0/(1 - J_nu*chi~0)   (HTML eq 29-30).
if nargin < 6, opts = struct(); end
eta   = 5e-3; if isfield(opts,'eta'),   eta   = opts.eta;   end
npass = 3;    if isfield(opts,'npass'), npass = opts.npass; end
Jsel  = ion.J0eff; if isfield(opts,'Jsel'), Jsel = opts.Jsel; end
si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', true));
z  = w(:) + 1i*eta;
c0 = invz_chi0z(si, T, z, struct('elastic', false));
G0 = -squeeze(c0(3,3,:));
tl = pt.tl;
g  = invz_g(tl, z);
pref = tl.M2/tl.n01^2;
if isfield(pt,'K') && ~isempty(pt.K)
    Kw = pt.K(1)*ones(size(z));                    % seed with static Matsubara K
else
    Kw = zeros(size(z));                           % RPA-limit callers may omit K
end
for pass = 1:npass
    Sw  = pt.alpha + pref*(pt.lambda(1) - (1 - tl.n01^2)*Kw).*g;
    if isfield(opts,'Jfull') && ~isempty(opts.Jfull)
        med = invz_emt_scalar(G0, Sw, opts.Jfull, struct('max_iter', 100, 'tol', 1e-8));
        Kw = med.K;
    else
        break;                                     % no lattice pass requested: single shot
    end
end
Sw = pt.alpha + pref*(pt.lambda(1) - (1 - tl.n01^2)*Kw).*g;
chit = (-G0) ./ (1 + Sw);
out.Sigma_w  = Sw;
out.chi0cc_w = -G0;
out.chi_cc_q = zeros(numel(Jsel), numel(z));
for k = 1:numel(Jsel)
    out.chi_cc_q(k,:) = (chit ./ (1 - Jsel(k)*chit)).';
end
end
