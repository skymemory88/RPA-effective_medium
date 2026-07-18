function [cdom, crest, mspec] = invzt_chi0_split(si, T, z, opts)
%INVZT_CHI0_SPLIT Transition-mask split of chi0(z) into a dominant
%(ground-manifold) block and the rest, with a stated elastic convention.
%
%   [cdom, crest, mspec] = INVZT_CHI0_SPLIT(si, T, z, opts) partitions every
%   transition (a,b) of INVZ_CHI0Z's single-ion susceptibility sum into
%   DOMINANT (both si.E(a) < Esplit and si.E(b) < Esplit) or REST, per
%   Jensen's "dominant transition renormalized, weak transitions kept at
%   RPA" rule (PRB 49 Sec. III; odd_implementation_plan.html). Esplit
%   defaults to 0.4653 meV, half the ~0.9306 meV (10.8 K) gap from the
%   LiHoF4 Ho3+ ground Gamma_{3,4} doublet to the first excited
%   crystal-field level, so the dominant group is exactly the ground
%   electronic doublet together with its (near-degenerate, hyperfine-scale
%   split) nuclear sublevels: ndom = 16 with opts.hyp = true (2 electronic
%   x 8 nuclear I=7/2 sublevels), ndom = 2 without hyperfine.
%
%   cdom is a masked re-run of INVZ_CHI0Z's per-transition sum, restricted
%   to dominant-dominant pairs (inelastic transitions AND elastic/degenerate
%   pairs alike -- see MASKED_CHI0 below). crest = full - cdom BY
%   CONSTRUCTION, where
%       full = INVZ_CHI0Z(si, T, z, opts_pass)
%   and opts_pass forwards degtol/ztol/elastic. This guarantees
%   cdom + crest == full to floating-point round-off, independent of any
%   physics judgment call folded into cdom's mask or elastic convention.
%
%   ELASTIC CONVENTION (LOCKED, docs/superpowers/plans/
%   2026-07-17-invz-tensor-full.md Task 5): INVZ_CHI0Z's z~0 elastic term
%   is, per Cartesian pair (mu,nu):
%       el(mu,nu) = beta*( sum_{(a,b) degenerate} M_mu(a,b) M_nu(b,a) p(a) )
%                   - beta*Jexp(mu)*Jexp(nu),      Jexp = si.Jexp
%   i.e. a connected degenerate-pair sum minus a counterterm bilinear in the
%   GLOBAL mean moment. cdom's elastic term restricts the connected sum to
%   degenerate pairs with BOTH endpoints dominant (the same (a,b) mask used
%   for the inelastic sector) and replaces the counterterm with the
%   DOMINANT-GROUP mean moment:
%       Jdom(mu) = sum_{p in dom} si.P(p) * M_mu(p,p)     (real part; NOT
%                   renormalized by sum_{p in dom} si.P(p) -- literal
%                   transcription of the locked box's Sigma_{p in dom} P_p M_pp)
%       cdom_el(mu,nu) = beta*( sum_{(a,b) degenerate, a,b in dom} M_mu(a,b)M_nu(b,a)p(a) )
%                        - beta*Jdom(mu)*Jdom(nu)
%   crest then absorbs, by construction, (i) any degenerate elastic weight
%   with an endpoint outside the dominant group, and (ii) the counterterm
%   residual beta*(Jdom(mu)Jdom(nu) - Jexp(mu)Jexp(nu)) coming purely from
%   using the dominant-group mean instead of the global mean. Both are
%   convention artifacts of WHERE the (exact, rank-1) counterterm bilinear
%   is assigned -- not physics. Their combined size relative to the total
%   elastic term is mspec.elastic_conv_share (below); it is bounded by the
%   excited-state population share (~1e-3 at 1.6 K), since
%   Jexp - Jdom = sum_{p not in dom} P(p)*M_mu(p,p) is itself O(that share).
%
%   mspec fields (fdom_* and elastic_conv_share are evaluated at the STATIC
%   point z=0 regardless of the z argument, since they are properties of
%   si/T alone, not of whichever frequencies were requested for cdom/crest):
%     ndom               number of states with si.E < Esplit
%     fdom_cc0           real(cdom_cc(0)) / real(full_cc(0)), cc = (3,3)
%                        (report; expected close to 1 -- the Ising/c-axis
%                        response is carried almost entirely by the ground
%                        manifold)
%     fdom_perp0         real(cdom_aa(0)) / real(full_aa(0)), aa = (1,1)
%                        (report; expected SMALL -- the transverse response
%                        is Van Vleck, i.e. lives mostly in crest)
%     elastic_conv_share norm(beta*(Jdom*Jdom.' - Jexp*Jexp.'), 'fro')
%                        / norm(full elastic tensor, 'fro')   (report)
%
%   si   : output of INVZ_SINGLE_ION.
%   T    : temperature, K.
%   z    : frequency argument(s) (meV), forwarded to INVZ_CHI0Z.
%   opts : Esplit (meV, default 0.4653), degtol/ztol/elastic (passthrough
%          to INVZ_CHI0Z's opts_pass, same defaults as INVZ_CHI0Z).
%
%   See also INVZ_CHI0Z, INVZ_SINGLE_ION.
if nargin < 4, opts = struct(); end
Esplit = getf(opts, 'Esplit', 0.4653);
degtol = getf(opts, 'degtol', 1e-8);
ztol   = getf(opts, 'ztol', 1e-12);
elast  = getf(opts, 'elastic', true);
opts_pass = struct('degtol', degtol, 'ztol', ztol, 'elastic', elast);

E = si.E;
domVec  = E < Esplit;                  % dominant-group membership, n x 1
domMask = domVec & domVec.';           % (a,b) dominant iff BOTH endpoints are
ndom    = sum(domVec);
Jdom    = dominant_mean(si, domVec);   % [3,1] dominant-group counterterm mean

full  = invz_chi0z(si, T, z, opts_pass);
cdom  = masked_chi0(si, T, z, domMask, Jdom, degtol, ztol, elast);
crest = full - cdom;                   % exact by construction

% ---- mspec: static-point (z=0) diagnostics, independent of the requested z ----
full0 = invz_chi0z(si, T, 0, opts_pass);
cdom0 = masked_chi0(si, T, 0, domMask, Jdom, degtol, ztol, elast);
mspec.ndom       = ndom;
mspec.fdom_cc0   = real(cdom0(3,3)) / real(full0(3,3));
mspec.fdom_perp0 = real(cdom0(1,1)) / real(full0(1,1));

% elastic_conv_share: isolate the elastic-only tensor (elastic minus
% inelastic-at-z=0, both via the reference INVZ_CHI0Z, independent of the
% opts.elastic flag used above) and compare the counterterm residual to it.
C = invz_const();  beta = 1/(C.kB*T);
Jexp = si.Jexp(:);
full_true  = invz_chi0z(si, T, 0, struct('degtol', degtol, 'ztol', ztol, 'elastic', true));
full_false = invz_chi0z(si, T, 0, struct('degtol', degtol, 'ztol', ztol, 'elastic', false));
full_el  = full_true(:,:,1) - full_false(:,:,1);
conv_dep = beta*(Jdom*Jdom.' - Jexp*Jexp.');
mspec.elastic_conv_share = norm(conv_dep, 'fro') / norm(full_el, 'fro');
end

% ---------------------------------------------------------------------
function Jbar = dominant_mean(si, domVec)
%DOMINANT_MEAN <J>_dom(mu) = sum_{p in dom} si.P(p)*M_mu(p,p) (real part),
% the dominant-group mean moment used as the elastic counterterm by the
% locked convention (see INVZT_CHI0_SPLIT header) -- NOT renormalized by
% sum_{p in dom} si.P(p).
M = {si.Mx, si.My, si.Mz};
w = si.P(:) .* domVec(:);
Jbar = zeros(3,1);
for mu = 1:3
    Jbar(mu) = real(diag(M{mu})).' * w;
end
end

% ---------------------------------------------------------------------
function c = masked_chi0(si, T, z, domMask, Jbar, degtol, ztol, elast)
%MASKED_CHI0 Dominant-masked re-run of INVZ_CHI0Z's per-transition sum --
% identical loop / population / matrix-element conventions as INVZ_CHI0Z
% (same dE, dp, inel, Nmn, P2, el, idx0 structure) -- restricted to
% transitions (a,b) with domMask(a,b) true, with the elastic counterterm
% built from Jbar in place of INVZ_CHI0Z's si.Jexp (see the ELASTIC
% CONVENTION note in INVZT_CHI0_SPLIT's header).
C = invz_const();  beta = 1/(C.kB*T);
E = si.E;  p = si.P;  n = numel(E);
dE = E.' - E;                    % dE(a,b) = E(b)-E(a)
dp = p - p.';                    % dp(a,b) = p(a)-p(b)
inel = abs(dE) > degtol;
M = {si.Mx, si.My, si.Mz};
z = z(:); nz = numel(z);
c = zeros(3,3,nz);
% loop-invariants (independent of the mu,nu Cartesian indices): hoist ONCE
keep = inel & domMask;
dEi  = dE(keep);
if elast
    P2    = repmat(p, 1, n);              % P2(a,b) = p(a)
    emask = (~inel) & domMask;
    idx0  = abs(z) < ztol;
end
for mu = 1:3
    for nu = 1:3
        Nmn  = M{mu} .* (M{nu}.');            % M_mu(a,b)*M_nu(b,a)
        wi   = Nmn(keep) .* dp(keep);
        for iz = 1:nz
            c(mu,nu,iz) = sum(wi ./ (dEi - z(iz)));
        end
        if elast
            el = beta*(sum(Nmn(emask).*P2(emask)) - Jbar(mu)*Jbar(nu));
            c(mu,nu,idx0) = c(mu,nu,idx0) + el;
        end
    end
end
end
