function ts = invzt_threestate(ion, T, B, opts)
%INVZT_THREESTATE  Explicit three-state toy single-ion model for A3 (Task 12).
%
%   ts = INVZT_THREESTATE(ion, T, B, opts) builds a COMPLETE single-ion struct
%   (the same field surface as INVZ_SINGLE_ION, so it drops straight into
%   INVZ_CHI0Z / INVZT_CHI0_SPLIT / the sum-rule reporting) for a minimal
%   3-level model {|1>,|2>,|3>}: a low Ising doublet (|1>,|2>) plus one excited
%   spectator |3>.  It is the single-ion input the A3 tensor 1/z self-energy
%   (INVZT_SIGMA_TENSOR) dresses.
%
%   THREE-STATE MODEL CONTRACT (v3, review Blocking 2, option (a) -- an INDEPENDENT
%   low-doublet tunnelling knob Delta1; rho reserved for the SPECTATOR channel).
%   In the fixed basis {|1>,|2>,|3>}, with gmuB = ion.gL*muB and Bx = B(1):
%       H3 = diag(0, 0, Delta2) + (Delta1/2)*Sx12 - (gmuB*Bx + hx)*Ja0 - hz*Jc0
%       Sx12 = |1><2| + |2><1|                 % DIRECT low-doublet tunnelling (the splitting knob)
%       Jc0  = diag(m0, -m0, 0)                 % Ising doublet moment (cc channel)
%       Ja0(1,3)=Ja0(3,1)=Ja0(2,3)=Ja0(3,2)=rho/sqrt(2)          % SPECTATOR (transverse a) coupling
%       Jb0  = (rho/sqrt(2))*( i|1><3| - i|3><1| - i|2><3| + i|3><2| )   % 90-deg partner (b), Hermitian
%   WHY the direct Delta1 term (v3): with it the hz=0 doublet block is
%   [[0,Delta1/2],[Delta1/2,0]] -- splitting Delta1, transition moment
%   |<+|Jc0|->| = m0 (M2 = m0^2), INDEPENDENT of rho.  So rho carries ONLY the
%   transverse/spectator channel, and rho -> 0 leaves a GENUINE two-level system
%   (the splitting survives as Delta1), making the exact scalar-compatibility
%   identity (INVZT_SIGMA_TENSOR / mode 'a3') consistent.  The field term couples
%   the doublet to |3> ONLY through Ja0, so at rho = 0 the field decouples entirely.
%
%   THE WELL-POSED MATCH.  The constructor NUMERICALLY solves the three knobs
%   (Delta1, m0, rho) against three targets (each knob mainly drives one target --
%   a diagonally-dominant 3x3 Newton solve on the FULL eig(H3), NOT a fragile
%   closed-form splitting expression):
%       (1) eig(H3) doublet splitting  ts.E(2)-ts.E(1) = tl.Delta
%       (2) |Mz(1,2)|^2                                = tl.M2      (Mz = V'*Jc0*V)
%       (3) transverse susceptibility  chi_perp = real chi0_aa(0) = chiperp_target
%   with tl = INVZ_TWOLEVEL(ion, T, Bx) (default Jxx0 = ion.Jxx0 -- the toy is a
%   standalone electronic model targeting the two-level object the scalar chain
%   uses).  Seeds: Delta1 = tl.Delta, m0 = sqrt(tl.M2),
%   rho = sqrt(chiperp_target*Delta2/2) (exact rho=0 doublet + leading spectator
%   estimate chi_perp ~ 2 rho^2 / Delta2).  Errors invzt:threeStateMatch if the
%   residual exceeds 1e-10 after Newton (a bad toy must not silently proceed).
%
%   opts.chiperp_scale (default 1) RESCALES the matched rho -> chiperp_scale*rho and
%   RE-REFINES (Delta1, m0) at that rescaled rho so the doublet ALWAYS reproduces
%   (tl.Delta, tl.M2) and ONLY chi_perp varies (proportional to chiperp_scale^2) --
%   NOT a post-hoc multiply that would leave (Delta1, m0) carrying the old rho^2
%   correction.  At chiperp_scale = 0 the re-refinement gives EXACTLY the two-level
%   (Delta1, m0) = (tl.Delta, sqrt(tl.M2)) with |3> DECOUPLED (Ja0 = Jb0 = 0): the
%   basis of the exact rho->0 scalar gate (exact at ANY Delta2 -- disconnection, not
%   depopulation, makes the limit exact).
%
%   OPTIONS (getf defaults):
%     far_excited     false   : Delta2 = 0.9306 meV (LiHoF4 CF gap); true -> 40 meV.
%     Delta2          [per above] : explicit override of the spectator gap (meV).
%     chiperp_target  11.05   : transverse susceptibility target (meV^-1), in-band.
%     chiperp_scale   1       : rescale rho -> chiperp_scale*rho + re-refine (Delta1,m0).
%     hz              0        : longitudinal field (meV). Paramagnetic rung -> 0.
%     transverse_mf   'none'   : 'legacy_x' iterates the toy transverse MF hx =
%                               Jxx0*<Ja0> (Jxx0 default ion.Jxx0); 'none' (default)
%                               keeps hx = 0 (the matching absorbs the real MF into
%                               Delta1, so the toy needs no MF of its own).
%     tol / maxit     1e-10 / 60 : Newton tolerance (max abs residual) and cap.
%   Hyperfine: EXCLUDED at this rung (a pure 3-state electronic toy).  By/Bz are
%   ignored (a minimal transverse-along-a model; the tests use Bx-only fields).
%
%   RETURNS a COMPLETE single-ion struct ts:
%     E [3x1] ground-shifted (from eig(H3)), V [3x3] eigenvectors, P [3x1] Boltzmann,
%     Mx/My/Mz = V'*(Ja0/Jb0/Jc0)*V (eigenbasis), Jexp [3x1], JzJz_fluct,
%     hx, hz, mf_converged, E0, transverse_mf, and match (achieved Delta/M2/chiperp,
%     targets, residuals, and the solved Delta1/m0/rho -- REPORTED, since chi_perp is
%     matched not imposed-by-luck).
%
%   See also INVZ_SINGLE_ION, INVZ_TWOLEVEL, INVZ_CHI0Z, INVZT_CHI0_SPLIT,
%   INVZT_SIGMA_TENSOR, INVZT_SOLVE_POINT.
if nargin < 4, opts = struct(); end
C = invz_const();
B = invz_field_vec(B);
gmuBx = ion.gL * C.muB * B(1);                 % transverse Zeeman energy scale (meV)

far  = getf(opts, 'far_excited', false);
D2   = getf(opts, 'Delta2', tern(far, 40, 0.9306));
xtar = getf(opts, 'chiperp_target', 11.05);
xscl = getf(opts, 'chiperp_scale', 1);
hz   = getf(opts, 'hz', 0);
tmf  = getf(opts, 'transverse_mf', 'none');
tol  = getf(opts, 'tol', 1e-10);
maxit = getf(opts, 'maxit', 60);

% --- targets from the two-level object the scalar chain uses (default Jxx0) -------
tl = invz_twolevel(ion, T, B(1), struct());
Dtar  = tl.Delta;
M2tar = tl.M2;

% --- toy transverse mean field (optional). Default 'none' -> hx = 0 --------------
if strcmp(tmf, 'legacy_x')
    Jxx0 = getf(opts, 'Jxx0', ion.Jxx0);
else
    Jxx0 = 0;                                  % 'none': hx stays 0
end

% ================= 3-param match: (Delta1, m0, rho) -> (Delta, M2, chi_perp) ======
% Each knob mainly drives one target (diagonally dominant); numeric-Jacobian Newton
% on the FULL eig(H3). hx (if any) is solved as an inner MF loop per H3 build.
buildsi = @(D1, m0, rho) toy_si(D1, m0, rho, D2, gmuBx, hz, T, Jxx0, tmf);
tgt3 = @(x) targets3(buildsi(x(1), x(2), x(3)), T) - [Dtar; M2tar; xtar];
x0 = [Dtar; sqrt(max(M2tar, 0)); sqrt(max(xtar*D2/2, 0))];
[x3, res3] = newton(tgt3, x0, tol, maxit);
D1 = x3(1); m0 = x3(2); rho = x3(3);

% ================= chiperp_scale: rescale rho, re-refine (Delta1, m0) =============
% Rescaling rho changes chi_perp (proportional to scale^2); re-solve (Delta1,m0) so
% the doublet STILL reproduces (Delta, M2). At scale=0: rho=0, and the re-refinement
% collapses to the exact two-level (Delta1,m0)=(Dtar, sqrt(M2tar)) with |3> decoupled.
rho = xscl * rho;
tgt2 = @(x) subsref_targets12(buildsi(x(1), x(2), rho), T) - [Dtar; M2tar];
[x2, res2] = newton(tgt2, [D1; m0], tol, maxit);
D1 = x2(1); m0 = x2(2);

% pick the binding residual: full 3-target for scale=1, doublet 2-target otherwise
if abs(xscl - 1) < 1e-14
    res_bind = res3;
else
    res_bind = res2;
end
if ~(res_bind < 1e-10)
    error('invzt:threeStateMatch', ...
        ['three-state match residual %.3e exceeds 1e-10 (targets: Delta=%.6g, ', ...
         'M2=%.6g, chiperp=%.6g, scale=%.6g) -- a bad toy must not silently proceed.'], ...
        res_bind, Dtar, M2tar, xtar, xscl);
end

% ================= assemble the final toy si + match report ======================
ts = buildsi(D1, m0, rho);
tg = targets3(ts, T);
ts.match = struct('Delta', tg(1), 'M2', tg(2), 'chiperp', tg(3), ...
    'Delta_target', Dtar, 'M2_target', M2tar, 'chiperp_target', xtar, ...
    'chiperp_scale', xscl, 'res_Delta', abs(tg(1) - Dtar), 'res_M2', abs(tg(2) - M2tar), ...
    'res_chiperp', abs(tg(3) - xscl^2*xtar), 'residual', res_bind, ...
    'Delta1', D1, 'm0', m0, 'rho', rho);
ts.mf_residual = res_bind;
end

% ======================================================================= %
%  Build the full toy single-ion struct from the three knobs.  hx is an inner
%  self-consistent transverse MF loop (a no-op when Jxx0 = 0 / tmf = 'none').
% ======================================================================= %
function si = toy_si(D1, m0, rho, D2, gmuBx, hz, T, Jxx0, tmf)
beta = 1/(0.08617333*T);                        % C.kB (meV/K), matches invz_const
r2 = rho/sqrt(2);
Jc0  = diag([m0, -m0, 0]);
Sx12 = [0 1 0; 1 0 0; 0 0 0];
Ja0  = [0 0 r2; 0 0 r2; r2 r2 0];
Jb0  = [0 0 1i*r2; 0 0 -1i*r2; -1i*r2 1i*r2 0];   % Hermitian 90-deg partner
H0   = diag([0 0 D2]) + (D1/2)*Sx12 - gmuBx*Ja0 - hz*Jc0;
hx = 0;
if Jxx0 ~= 0                                    % inner transverse MF fixed point (hx = Jxx0*<Ja0>)
    for it = 1:200
        H = (H0 - hx*Ja0);  H = (H + H')/2;
        [V, D] = eig(H, 'vector');  [E, ix] = sort(real(D));  V = V(:, ix);
        p = exp(-beta*(E - E(1)));  p = p/sum(p);
        jx = real(diag(V'*Ja0*V)).' * p;
        hx_new = Jxx0*jx;
        if abs(hx_new - hx) < 1e-12, hx = hx_new; break; end
        hx = hx + (hx_new - hx);                % undamped (near-linear, converges fast)
    end
end
H = (H0 - hx*Ja0);  H = (H + H')/2;
[V, D] = eig(H, 'vector');  [E, ix] = sort(real(D));  V = V(:, ix);
p = exp(-beta*(E - E(1)));  p = p/sum(p);
si.E  = E - E(1);
si.V  = V;
si.P  = p;
si.Mx = V'*Ja0*V;  si.My = V'*Jb0*V;  si.Mz = V'*Jc0*V;
si.Jexp = [real(diag(si.Mx)).'*p; real(diag(si.My)).'*p; real(diag(si.Mz)).'*p];
si.hx = hx;  si.hy = 0;  si.hz = hz;
jc2 = real(diag(V'*(Jc0*Jc0)*V)).' * p;
si.JzJz_fluct = jc2 - si.Jexp(3)^2;
si.mf_converged = true;
si.mf_iters = 1;
si.E0 = E(1);
si.transverse_mf = tmf;
end

% ----- three toy targets [splitting; M2; chi_perp] -----
function t = targets3(si, T)
c0 = invz_chi0z(si, T, 0, struct('elastic', true));
t = [si.E(2) - si.E(1); abs(si.Mz(1,2))^2; real(c0(1,1))];
end

function t = subsref_targets12(si, T)
t3 = targets3(si, T);
t = t3(1:2);
end

% ----- damped Newton with backtracking (central-difference Jacobian) -----
function [x, res] = newton(fun, x0, tol, maxit)
x = x0(:);  n = numel(x);
f = fun(x);  res = max(abs(f));
for it = 1:maxit
    if res < tol, break; end
    J = zeros(numel(f), n);
    for j = 1:n
        h = 1e-7 * max(1, abs(x(j)));
        xp = x; xp(j) = xp(j) + h;
        xm = x; xm(j) = xm(j) - h;
        J(:, j) = (fun(xp) - fun(xm)) / (2*h);
    end
    step = J \ f;
    lam = 1;  accepted = false;
    for bt = 1:20                               % backtracking line search
        xt = x - lam*step;
        ft = fun(xt);
        if max(abs(ft)) < res
            x = xt;  f = ft;  res = max(abs(ft));  accepted = true;  break;
        end
        lam = lam/2;
    end
    if ~accepted, break; end                    % stalled -> return best so far
end
end

function v = tern(c, a, b)
if c, v = a; else, v = b; end
end
