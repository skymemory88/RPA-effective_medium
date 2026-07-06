function si = invz_single_ion(ion, T, B, opts)
%INVZ_SINGLE_ION Exact single-ion diagonalization with transverse mean field.
% B = [Bx By Bz] in T. opts.hyp: include nuclear I=7/2 (dim 136). Paramagnetic phase:
% no <Jz> mean field is applied; caller must check si.Jexp(3) ≈ 0.
if nargin < 4, opts = struct(); end
hyp = isfield(opts,'hyp') && opts.hyp;
Jxx0 = ion.Jxx0; if isfield(opts,'Jxx0'), Jxx0 = opts.Jxx0; end
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
converged = false;
for it = 1:200                                   % transverse mean-field fixed point
    H = H0 - hx*Jx;
    H = (H + H')/2;
    [V, D] = eig(H, 'vector');
    [E, ix] = sort(real(D));  V = V(:, ix);
    p = exp(-beta*(E - E(1)));  p = p/sum(p);
    jx = real(diag(V'*Jx*V)).'*p;
    hx_new = Jxx0*jx;
    if abs(hx_new - hx) < 1e-12
        hx = hx_new;
        converged = true;
        break;
    end
    hx = hx_new;
end
if ~converged
    warning('invz:mfNotConverged', 'Transverse mean field not converged after %d iterations: |dhx| = %.3g meV', it, abs(hx_new-hx));
end
% Recompute H, eig, populations, and all output fields ONCE from the final
% converged hx, so the returned struct is exactly self-consistent with si.hx.
H = H0 - hx*Jx;
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
jz2 = real(diag(V'*(Jz*Jz)*V)).'*p;
si.JzJz_fluct = jz2 - si.Jexp(3)^2;
end
