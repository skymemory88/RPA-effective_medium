function out = invz_scalar_ewald_ref(q, a, tau, alpha, r_cut, g_cut)
%INVZ_SCALAR_EWALD_REF Gate-B independent scalar-Coulomb Ewald oracle (Step-4 Task 8).
%
%   out = INVZ_SCALAR_EWALD_REF(q, a, tau, alpha, r_cut, g_cut)
%
% A GENUINELY SEPARATE code path from invz_dipole_ewald.m / MF_dipole.m: it
% builds the scalar-Coulomb Ewald LATTICE POTENTIAL phi_lat(x), with the
% phase on the Bravais vector R ONLY (not on the total displacement x=R+d_nm,
% the primitive's own convention), and differentiates it NUMERICALLY
% (finite-difference Hessian + Richardson extrapolation) to obtain the
% dipolar tensor. It never calls invz_dipole_ewald, its private local
% functions (unreachable from another file by MATLAB's own file-scoping, in
% any case), or its geometry builder -- this file depends on nothing but
% (q, a, tau, alpha, r_cut, g_cut).
%
%   qcart = q*B,  B = 2*pi*inv(a)'
%   phi_lat(x) = sum_{R, |R+x|<=r_cut}  exp(-i*qcart.R) erfc(alpha*|R+x|)/|R+x|
%              + (4*pi/Vc) sum_{k=qcart+G, k~=0, |k|<=g_cut}
%                  exp(-|k|^2/(4*alpha^2))/|k|^2 * exp(+i*k.x)
%   dip_nm(q)  = -exp(-i*qcart.d_nm) * Hessian_x[phi_lat](x=d_nm)
%
% with k = qcart+Gcart everywhere (never a bare G): the reciprocal
% denominator, Gaussian, cutoff, AND phase all use k. The scalar k=0 term is
% a constant in x and drops under the Hessian, so it is simply never formed
% (k~=0 is already required by the candidate filter). For n=m, the integer
% real cell R_hkl=[0,0,0] is dropped BY CELL INDEX from every displaced FD
% evaluation (the "regularized non-self" potential -- see prereg sec 4); for
% n~=m no cell is dropped. The "-exp(-i*qcart.d_nm)" prefactor converts the
% phase-on-R convention above into the primitive's total-displacement
% convention (see docs/invzp_ewald_prereg.md sec 4, FROZEN); it is applied to
% each raw per-h tensor before Richardson extrapolation (a linear operation,
% so extrapolating before or after multiplying by this q/pair-fixed complex
% scalar is identical).
%
% Finite-difference stencils (identical family to Task 3's fd_hessian_gab,
% applied to phi_lat instead of the bare screened function), Cartesian steps
% h in {4e-3, 2e-3, 1e-3} Angstrom:
%   diagonal (a=b):  [f(x+h e_a) - 2f(x) + f(x-h e_a)] / h^2
%   mixed (a~=b):    [f(x+h e_a+h e_b) - f(x+h e_a-h e_b)
%                      - f(x-h e_a+h e_b) + f(x-h e_a-h e_b)] / (4h^2)
% (never the mixed formula for a=b -- that silently turns the labelled step h
% into a 2h diagonal stencil). Richardson (both O(h^2)-removed, from adjacent
% step pairs): R12=(4*H_h2-H_h1)/3 (coarse pair), R23=(4*H_h3-H_h2)/3 (fine
% pair).
%
% INPUTS
%   q      [1,3] reduced reciprocal (Miller) coordinates, finite real. Used
%          DIRECTLY as qcart=q*B -- there is no extended-zone K-reduction
%          here (unlike the primitive): the oracle is only ever exercised at
%          the frozen canonical fx.q_int, already inside [-0.5,0.5) on every
%          component, where the primitive's own K=floor(q+0.5) reduction is
%          itself a no-op (K=[0,0,0]), so the two conventions coincide.
%   a      [3,3]    direct lattice, rows = lattice vectors, finite real
%                    nonsingular.
%   tau    [ntau,3] basis positions, fractional coordinates, finite real.
%   alpha, r_cut, g_cut   finite positive scalars (matched generous cutoffs).
%
% OUTPUT out -- machine-readable struct, every field finite:
%   h                   [1,3] the frozen step ladder [4e-3 2e-3 1e-3].
%   raw                 [3,3,ntau,ntau,3] the three raw FD dip_nm(q) estimates,
%                       one per h (already in the primitive's dip convention:
%                       the -exp(-i*qcart.d_nm) prefactor is applied).
%   richardson_coarse   [3,3,ntau,ntau] R12=(4*raw(:,:,:,:,2)-raw(:,:,:,:,1))/3.
%   richardson_fine     [3,3,ntau,ntau] R23=(4*raw(:,:,:,:,3)-raw(:,:,:,:,2))/3.
%   self_analytic       [3,3,ntau,ntau] -delta_nm delta_ab*4*alpha^3/(3*sqrt(pi)).
%   adjacent_residual   [3,3,ntau,ntau] richardson_fine - richardson_coarse.
%   q, alpha, r_cut, g_cut   echoed inputs (provenance only).
%
% See docs/invzp_ewald_prereg.md sec 4 (Gate B, FROZEN 2026-07-24) and
% docs/superpowers/plans/2026-07-24-ewald-dipolar-primitive.md, Task 8.

% ---- light shape/finiteness sanity (a test oracle, not a validated
% production entry point -- deliberately NOT the primitive's invz:* contract)
q = reshape(q, 1, 3);
assert(all(isfinite(q)), 'invz_scalar_ewald_ref: q must be finite.');
assert(isequal(size(a), [3 3]) && all(isfinite(a(:))), ...
    'invz_scalar_ewald_ref: a must be a finite [3,3] matrix.');
assert(size(tau,2) == 3 && size(tau,1) >= 1 && all(isfinite(tau(:))), ...
    'invz_scalar_ewald_ref: tau must be a finite [ntau,3] matrix.');
assert(isscalar(alpha) && isfinite(alpha) && alpha > 0, ...
    'invz_scalar_ewald_ref: alpha must be a finite positive scalar.');
assert(isscalar(r_cut) && isfinite(r_cut) && r_cut > 0, ...
    'invz_scalar_ewald_ref: r_cut must be a finite positive scalar.');
assert(isscalar(g_cut) && isfinite(g_cut) && g_cut > 0, ...
    'invz_scalar_ewald_ref: g_cut must be a finite positive scalar.');

ntau    = size(tau, 1);
taucart = tau*a;
B       = 2*pi*inv(a).';                             %#ok<MINV>
Vc      = abs(det(a));
qcart   = q*B;

h_list = [4e-3 2e-3 1e-3];

% ---- reciprocal candidate list: shared by EVERY pair and EVERY stencil
% point (k=qcart+Gcart depends only on q and alpha, never on x or the pair) --
sb = min(svd(B));
[c1, c2, c3] = ndgrid([-0.5 0.5], [-0.5 0.5], [-0.5 0.5]);
qcorner = [c1(:) c2(:) c3(:)]*B;
qmax    = max(vecnorm(qcorner, 2, 2));
nmax_G  = ceil((g_cut + qmax)/sb) + 2;                % +2 integer candidate slack
rngG    = -nmax_G:nmax_G;
[GH, GK, GL] = ndgrid(rngG, rngG, rngG);
Ghkl  = [GH(:) GK(:) GL(:)];
Gcart = Ghkl*B;
k  = Gcart + qcart;
kk = sum(k.^2, 2);
keepk = (kk <= g_cut^2) & (kk > 0);                   % k~=0 exactly
ksel   = k(keepk, :);
kk2    = kk(keepk);
kernel = (4*pi/Vc) * exp(-kk2/(4*alpha^2)) ./ kk2;    % [NGk,1], independent of x

% ---- real-space geometry + FD assembly, every ordered pair independent ---
sa = min(svd(a));
FD_SLACK = 0.05;    % Angstrom, >> max FD offset magnitude (sqrt(2)*4e-3~5.7e-3)

selfval = 4*alpha^3/(3*sqrt(pi));

raw           = complex(zeros(3, 3, ntau, ntau, 3));
self_analytic = complex(zeros(3, 3, ntau, ntau));
for n = 1:ntau
    self_analytic(:,:,n,n) = -selfval*eye(3);
end

for n = 1:ntau
    for m = 1:ntau
        d = taucart(m,:) - taucart(n,:);
        exclude_R0 = (n == m);

        nmax_r = ceil((r_cut + norm(d) + FD_SLACK)/sa);
        rngR = -nmax_r:nmax_r;
        [RH, RK, RL] = ndgrid(rngR, rngR, rngR);
        Rhkl = [RH(:) RK(:) RL(:)];
        if exclude_R0
            Rhkl = Rhkl(~all(Rhkl == 0, 2), :);        % drop R=[0,0,0] BY CELL INDEX
        end
        Rcart  = Rhkl*a;                               % [NR,3]
        phaseR = exp(-1i*(Rcart*qcart.'));              % [NR,1]; phase on R ONLY

        pref = -exp(-1i*(qcart*d.'));                   % scalar: phase-on-R -> primitive convention

        for hi = 1:3
            Hh = local_fd_hessian(d, h_list(hi), Rcart, phaseR, r_cut, alpha, ksel, kernel);
            raw(:,:,n,m,hi) = pref*Hh;
        end
    end
end

richardson_coarse = (4*raw(:,:,:,:,2) - raw(:,:,:,:,1))/3;
richardson_fine   = (4*raw(:,:,:,:,3) - raw(:,:,:,:,2))/3;
adjacent_residual = richardson_fine - richardson_coarse;

out = struct();
out.h                 = h_list;
out.raw               = raw;
out.richardson_coarse = richardson_coarse;
out.richardson_fine   = richardson_fine;
out.self_analytic     = self_analytic;
out.adjacent_residual  = adjacent_residual;
out.q = q; out.alpha = alpha; out.r_cut = r_cut; out.g_cut = g_cut;

assert(all(isfinite(out.h(:))),                 'invz_scalar_ewald_ref: h is not finite.');
assert(all(isfinite(out.raw(:))),               'invz_scalar_ewald_ref: raw is not finite.');
assert(all(isfinite(out.richardson_coarse(:))), 'invz_scalar_ewald_ref: richardson_coarse is not finite.');
assert(all(isfinite(out.richardson_fine(:))),   'invz_scalar_ewald_ref: richardson_fine is not finite.');
assert(all(isfinite(out.self_analytic(:))),     'invz_scalar_ewald_ref: self_analytic is not finite.');
assert(all(isfinite(out.adjacent_residual(:))), 'invz_scalar_ewald_ref: adjacent_residual is not finite.');
end

% =====================================================================
% phi_lat(x) evaluation (real part uses the precomputed per-pair candidate
% set Rcart/phaseR; reciprocal part uses the shared per-call ksel/kernel).
% =====================================================================
function val = local_phi_eval(x, Rcart, phaseR, r_cut, alpha, ksel, kernel)
x  = reshape(x, 1, 3);
Xr = Rcart + x;
r  = vecnorm(Xr, 2, 2);
keep = (r <= r_cut);
val_real  = sum(phaseR(keep) .* (erfc(alpha*r(keep)) ./ r(keep)));
val_recip = sum(kernel .* exp(1i*(ksel*x.')));
val = val_real + val_recip;
end

% =====================================================================
% Central-difference Hessian of phi_lat at x0, step h -- SAME stencil family
% as Task 3's fd_hessian_gab (test_invz_dipole_ewald.m), applied to the full
% lattice potential rather than the bare screened function. Every (a,b) cell
% of the 3x3 loop is evaluated independently (both (a,b) and (b,a) for
% a~=b), never assembled by symmetry.
% =====================================================================
function H = local_fd_hessian(x0, h, Rcart, phaseR, r_cut, alpha, ksel, kernel)
x0 = reshape(x0, 1, 3);
H  = zeros(3, 3);
f0 = local_phi_eval(x0, Rcart, phaseR, r_cut, alpha, ksel, kernel);
for aa = 1:3
    ea = zeros(1,3); ea(aa) = 1;
    for bb = 1:3
        if aa == bb
            fp = local_phi_eval(x0+h*ea, Rcart, phaseR, r_cut, alpha, ksel, kernel);
            fm = local_phi_eval(x0-h*ea, Rcart, phaseR, r_cut, alpha, ksel, kernel);
            H(aa,bb) = (fp - 2*f0 + fm)/h^2;
        else
            eb = zeros(1,3); eb(bb) = 1;
            fpp = local_phi_eval(x0+h*ea+h*eb, Rcart, phaseR, r_cut, alpha, ksel, kernel);
            fpm = local_phi_eval(x0+h*ea-h*eb, Rcart, phaseR, r_cut, alpha, ksel, kernel);
            fmp = local_phi_eval(x0-h*ea+h*eb, Rcart, phaseR, r_cut, alpha, ksel, kernel);
            fmm = local_phi_eval(x0-h*ea-h*eb, Rcart, phaseR, r_cut, alpha, ksel, kernel);
            H(aa,bb) = (fpp - fpm - fmp + fmm)/(4*h^2);
        end
    end
end
end
