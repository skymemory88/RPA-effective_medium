function si = invz_single_ion(ion, T, B, opts)
%INVZ_SINGLE_ION Exact single-ion diagonalization with mean field(s).
% B = [Bx By Bz] in T. opts.hyp: include nuclear I=7/2 (dim 136).
%
% Two mean-field modes:
%   Paramagnetic (default): only the transverse mean field hx = Jxx0*<Jx> is solved; the
%     longitudinal moment stays <Jz> = 0 (caller may check si.Jexp(3) ~ 0).
%   Ordered (opts.order = true): additionally solve the longitudinal ORDERING mean field
%     hz = J0z*<Jz> (J0z defaults to ion.J0eff = Jcc(0)) for its symmetry-broken branch,
%     seeded by opts.mz_seed. This is the ferromagnetic phase: si.Jexp(3) = <Jz> = m0 is the
%     order parameter, and it relaxes back to ~0 if the (T,Bx) point is not actually ordered.
%   Fixed longitudinal field (opts.hz_fixed = hz, meV): apply a FIXED -hz*Jz term without
%     self-solving it (only hx is iterated). Used to impose the full-space ordering mean field
%     on the electronic-doublet solve (invz_twolevel) so the two-level Sigma params see the
%     same hz as the full electronuclear order parameter. Mutually exclusive with opts.order.
if nargin < 4, opts = struct(); end
hyp  = isfield(opts,'hyp') && opts.hyp;
Jxx0 = ion.Jxx0;    if isfield(opts,'Jxx0'), Jxx0 = opts.Jxx0; end
order = isfield(opts,'order') && opts.order;
hzfix = isfield(opts,'hz_fixed');
J0z  = ion.J0eff;   if isfield(opts,'J0z'),   J0z   = opts.J0z;   end
mzsd = 1.0;         if isfield(opts,'mz_seed'), mzsd = opts.mz_seed; end
mix  = 1.0;         if isfield(opts,'mf_mix'), mix = opts.mf_mix; elseif order, mix = 0.6; end
maxit = 200;        if isfield(opts,'mf_maxit'), maxit = opts.mf_maxit; elseif order, maxit = 800; end
C   = invz_const();
oJ  = stevens_ops(ion.J);
Hcf = ion.B20*oJ.O20 + ion.B40*oJ.O40 + ion.B44*oJ.O44 ...
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
hz = 0;
if order, hz = J0z*mzsd; end                   % symmetry-breaking seed for the ordered branch
if hzfix, hz = opts.hz_fixed; end              % imposed longitudinal mean field (held fixed)
converged = false;
for it = 1:maxit                               % mean-field fixed point (hx, and hz if ordered)
    H = H0 - hx*Jx - hz*Jz;
    H = (H + H')/2;
    [V, D] = eig(H, 'vector');
    [E, ix] = sort(real(D));  V = V(:, ix);
    p = exp(-beta*(E - E(1)));  p = p/sum(p);
    jx = real(diag(V'*Jx*V)).'*p;
    hx_new = Jxx0*jx;
    if order
        jz = real(diag(V'*Jz*V)).'*p;
        hz_new = J0z*jz;
    else
        hz_new = hz;                           % para (0) or held-fixed hz_fixed
    end
    dmf = max(abs(hx_new - hx), abs(hz_new - hz));
    if dmf < 1e-12
        hx = hx_new;  hz = hz_new;
        converged = true;
        break;
    end
    hx = hx + mix*(hx_new - hx);
    hz = hz + mix*(hz_new - hz);
end
if ~converged
    warning('invz:mfNotConverged', 'Mean field not converged after %d iterations: |dmf| = %.3g meV', it, dmf);
end
% Recompute H, eig, populations, and all output fields ONCE from the final converged
% (hx, hz), so the returned struct is exactly self-consistent with si.hx / si.hz.
H = H0 - hx*Jx - hz*Jz;
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
si.hz = hz;
jz2 = real(diag(V'*(Jz*Jz)*V)).'*p;
si.JzJz_fluct = jz2 - si.Jexp(3)^2;
end
