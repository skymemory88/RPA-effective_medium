function si = invz_single_ion(ion, T, B, opts)
%INVZ_SINGLE_ION Exact single-ion diagonalization with mean field(s).
% B = [Bx By Bz] in T. opts.hyp: include nuclear I=7/2 (dim 136).
%
% Three mean-field modes:
%   Paramagnetic (default): solves transverse hx = Jxx0*<Jx> only; <Jz> stays 0.
%   Ordered (opts.order = true): also solves longitudinal hz = J0z*<Jz> (J0z default
%     ion.J0eff), seeded by opts.mz_seed. si.Jexp(3) = <Jz> is the FM order parameter,
%     relaxing to ~0 if (T,Bx) is not actually ordered.
%   Fixed longitudinal field (opts.hz_fixed = hz, meV): imposes a FIXED -hz*Jz term (only
%     hx is iterated), used to give invz_twolevel_ordered the same hz as the full
%     electronuclear order parameter. Mutually exclusive with opts.order.
%
% Returned diagnostics: si.mf_converged/mf_iters/mf_residual (mean-field loop state),
% si.E0 (unshifted ground energy), si.F_mf (variational MF free energy; NaN with hz_fixed).
% In order mode the mz_seed default is sign-aware: -1 when B(3) < 0 (aligned branch).
%
% opts.transverse_mf selects how the transverse (ab-plane) mean field is handled,
% independently of the longitudinal mode above:
%   'legacy_x' (default): only hx = Jxx0*<Jx> is solved; hy stays exactly 0 (bit-for-bit
%     with the pre-existing single-axis code path).
%   'none': both transverse channels are forced to zero (hx = hy = 0) regardless of Jxx0,
%     i.e. the bare CF + Zeeman problem (equivalent to opts.Jxx0 = 0).
%   'vector_ab': solves both hx = Jxx0*<Jx> and hy = Jxx0*<Jy> self-consistently. Needed
%     for any in-plane field, including along the a/b axes: the B64s crystal-field term
%     breaks the mirror symmetry that would otherwise force <Jy> = 0, so <Jy> ~= 0 even
%     for a field applied purely along x.
% si.hy (meV) is the converged transverse-b mean field (0 unless 'vector_ab' finds
% <Jy> ~= 0); si.transverse_mf echoes the resolved mode string. si.F_mf includes the
% matching 0.5*hy*si.Jexp(2) term.
if nargin < 4, opts = struct(); end
hyp  = isfield(opts,'hyp') && opts.hyp;
Jxx0 = ion.Jxx0;    if isfield(opts,'Jxx0'), Jxx0 = opts.Jxx0; end
tmf = 'legacy_x'; if isfield(opts,'transverse_mf'), tmf = opts.transverse_mf; end
if ~(ischar(tmf) || isstring(tmf)) || ~any(strcmp(tmf, {'legacy_x','none','vector_ab'}))
    error('invz:transverseMF', ...
        'transverse_mf must be ''legacy_x'', ''none'' or ''vector_ab''.');
end
tmf = char(tmf);
vecmf  = strcmp(tmf, 'vector_ab');
nonemf = strcmp(tmf, 'none');
order = isfield(opts,'order') && opts.order;
hzfix = isfield(opts,'hz_fixed');
J0z  = ion.J0eff;   if isfield(opts,'J0z'),   J0z   = opts.J0z;   end
mzsd = 1.0;         if order && B(3) < 0, mzsd = -1.0; end
                    % sign-aware default: an explicit longitudinal field selects the
                    % ALIGNED branch (the +1 seed can trap the metastable mirror state
                    % below the transition; spec 2026-07-16, review finding 2)
if isfield(opts,'mz_seed'), mzsd = opts.mz_seed; end
mix  = 1.0;         if isfield(opts,'mf_mix'), mix = opts.mf_mix; elseif order, mix = 0.6; end
maxit = 200;        if isfield(opts,'mf_maxit'), maxit = opts.mf_maxit; elseif order, maxit = 800; end
C   = invz_const();
oJ  = stevens_ops(ion.J);
B44s = 0;  if isfield(ion,'B44s'), B44s = ion.B44s; end   % sine partner of B44: only set by
                                                          % the rotated-CF route (invz_cfrot)
Hcf = ion.B20*oJ.O20 + ion.B40*oJ.O40 + ion.B44*oJ.O44 + B44s*oJ.O44s ...
    + ion.B60*oJ.O60 + ion.B64c*oJ.O64c + ion.B64s*oJ.O64s;
if hyp
    oI = stevens_ops(ion.I);
    nI = size(oI.Jz,1);
    kJ = @(M) kron(M, eye(nI));
    Hhf = ion.A*(kron(oJ.Jx,oI.Jx) + kron(oJ.Jy,oI.Jy) + kron(oJ.Jz,oI.Jz));
else
    kJ = @(M) M;  Hhf = 0;
end
Jx = kJ(oJ.Jx);  Jy = kJ(oJ.Jy);  Jz = kJ(oJ.Jz);
H0 = kJ(Hcf) + Hhf - ion.gL*C.muB*(B(1)*Jx + B(2)*Jy + B(3)*Jz);
beta = 1/(C.kB*T);
hx = 0;
hy = 0;
hz = 0;
if order, hz = J0z*mzsd; end                   % symmetry-breaking seed for the ordered branch
if hzfix, hz = opts.hz_fixed; end              % imposed longitudinal mean field (held fixed)
converged = false;
for it = 1:maxit                               % mean-field fixed point (hx[,hy], and hz if ordered)
    H = H0 - hx*Jx - hy*Jy - hz*Jz;
    H = (H + H')/2;
    [V, D] = eig(H, 'vector');
    [E, ix] = sort(real(D));  V = V(:, ix);
    p = exp(-beta*(E - E(1)));  p = p/sum(p);
    jx = real(diag(V'*Jx*V)).'*p;
    if nonemf
        hx_new = 0;  hy_new = 0;
    else
        hx_new = Jxx0*jx;
        hy_new = 0;
        if vecmf
            jy = real(diag(V'*Jy*V)).'*p;
            hy_new = Jxx0*jy;
        end
    end
    if order
        jz = real(diag(V'*Jz*V)).'*p;
        hz_new = J0z*jz;
    else
        hz_new = hz;                           % para (0) or held-fixed hz_fixed
    end
    dmf = max([abs(hx_new - hx), abs(hy_new - hy), abs(hz_new - hz)]);
    if dmf < 1e-12
        hx = hx_new;  hy = hy_new;  hz = hz_new;
        converged = true;
        break;
    end
    hx = hx + mix*(hx_new - hx);
    hy = hy + mix*(hy_new - hy);
    hz = hz + mix*(hz_new - hz);
end
if ~converged
    warning('invz:mfNotConverged', 'Mean field not converged after %d iterations: |dmf| = %.3g meV', it, dmf);
end
% Recompute H, eig, populations, and all output fields ONCE from the final converged
% (hx, hy, hz), so the returned struct is exactly self-consistent with si.hx / si.hy / si.hz.
H = H0 - hx*Jx - hy*Jy - hz*Jz;
H = (H + H')/2;
[V, D] = eig(H, 'vector');
[E, ix] = sort(real(D));  V = V(:, ix);
p = exp(-beta*(E - E(1)));  p = p/sum(p);
si.E  = E - E(1);
si.V  = V;
si.P  = p;
si.Mx = V'*Jx*V;  si.My = V'*Jy*V;  si.Mz = V'*Jz*V;
si.Jexp = [real(diag(si.Mx)).'*p; real(diag(si.My)).'*p; real(diag(si.Mz)).'*p];
si.hx = hx;
si.hy = hy;
si.hz = hz;
jz2 = real(diag(V'*(Jz*Jz)*V)).'*p;
si.JzJz_fluct = jz2 - si.Jexp(3)^2;
si.mf_converged = converged;
si.mf_iters     = it;
si.mf_residual  = dmf;
si.E0 = E(1);                                  % unshifted ground energy (si.E stays shifted)
si.transverse_mf = tmf;
% Variational MF free energy (branch diagnostic, spec SR1): the 0.5*h*<J> terms
% restore the -1/2 J <J>^2 double counting; equals Fsite + (hx^2+hy^2)/(2Jxx0) + hz^2/(2J0z)
% at a self-consistent point. Undefined (NaN) under hz_fixed: hz is imposed, not a MF.
Fsite = E(1) - log(sum(exp(-beta*(E - E(1)))))/beta;
if hzfix
    si.F_mf = NaN;
else
    si.F_mf = Fsite + 0.5*(hx*si.Jexp(1) + hy*si.Jexp(2) + hz*si.Jexp(3));
end
end
