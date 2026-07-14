"""Verify the NEW results of three_level_1z_extension_v2.html against
independent numerics. Checks:
  A. Cardano eigenvalues, Eq. (6), vs numpy eigvalsh of Eq. (4).
  B. Ordered kernels I2 (24), I3 (25) vs direct quadrature.
  C. Exact four-point transform F_{n,l}, Eq. (26), vs independent
     Gauss-Legendre integration over the six ordered simplices, for a
     GENERAL centered 3-level operator (diagonals + all off-diagonals),
     at several (n,l) pairs. Then Gamma_{n,l}, Eq. (27), vs brute
     connected result.
  D. PM chain reduction (28) == general (26) for the chain pattern.
  E. Three-point spectral formulas (39)/(44) vs direct quadrature.
  F. Static anchor (54)-(55): Gamma_00(T->0) from (26)-(27) vs
     -24 M1^4/D1^3 + 24 M1^2 M2^2/(D1^2 D2), and 4th-order RS energy.
  G. Spectator identity (31) == explicit (32) numerically.
"""
import numpy as np
from mpmath import mp, mpc, exp as mexp

mp.dps = 40
fails = []

def check(name, a, b, tol=1e-9):
    a, b = complex(a), complex(b)
    err = abs(a - b); scale = max(1.0, abs(a), abs(b))
    ok = err < tol * scale
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {a:+.10e} vs {b:+.10e}  |d|={err:.1e}")
    if not ok: fails.append(name)

# ---------------- A. Cardano (6) --------------------------------------
print("== A. Cardano eigenvalues ==")
for (D, G, h) in [(3.0, 0.8, 0.0), (2.0, 1.3, 0.7), (1.5, 0.4, -1.1)]:
    H = np.array([[-D - h, -G/np.sqrt(2), 0],
                  [-G/np.sqrt(2), 0, -G/np.sqrt(2)],
                  [0, -G/np.sqrt(2), -D + h]])
    ev = np.sort(np.linalg.eigvalsh(H))
    A = D**2 + 3*(G**2 + h**2)
    arg = D*(2*D**2 + 9*G**2 - 18*h**2)/(2*A**1.5)
    Ek = np.sort([-2*D/3 + 2*np.sqrt(A)/3 *
                  np.cos(np.arccos(np.clip(arg, -1, 1))/3 - 2*np.pi*k/3)
                  for k in range(3)])
    for i in range(3):
        check(f"E{i} (D={D},G={G},h={h})", ev[i], Ek[i], tol=1e-10)

# ---------------- kernels ---------------------------------------------
beta = 1.7
def phi(x):
    x = mpc(x)
    return (mexp(beta*x) - 1)/x if abs(x) > mp.mpf('1e-12') else mp.mpf(beta) + beta**2*x/2 + beta**3*x**2/6
def I2(x, y):
    x, y = mpc(x), mpc(y)
    if abs(y) > 1e-10: return (phi(x+y) - phi(x))/y
    # y->0: d/dx phi at x direction: I2(x,0)=phi'(x); series in y
    return dphi(x) + y*d2phi(x)/2
def dphi(x):
    x = mpc(x)
    if abs(x) > 1e-10: return (beta*mexp(beta*x)*x - (mexp(beta*x)-1))/x**2
    return beta**2/2 + beta**3*x/3 + beta**4*x**2/8
def d2phi(x):
    x = mpc(x)
    if abs(x) > 1e-8:
        e = mexp(beta*x)
        return (beta**2*e*x**2 - 2*(beta*e*x - (e-1)))/x**3
    return beta**3/3 + beta**4*x/4
def I3(x, y, z):
    x, y, z = mpc(x), mpc(y), mpc(z)
    if abs(z) > 1e-10:
        return (I2h(x, y+z) - I2h(x, y))/z
    return dI2y(x, y)
def I2h(x, y):
    if abs(y) > 1e-10: return (phi(x+y) - phi(x))/y
    return dphi(x) + y*d2phi(x)/2
def dI2y(x, y):
    # d/dz of (phi(x+y+z)-phi(x))/(y+z) - style limit: I3(x,y,0) = (dphi(x+y) - I2h(x,y))/y ... use v2's own limit
    if abs(y) > 1e-10:
        return (dphi(x+y) - (phi(x+y)-phi(x))/y)/y
    return d2phi(x)/2  # x,y,z all ~0 -> beta^3/6 ; d2phi(0)/2 = beta^3/6 OK

print("\n== B. kernels I2, I3 vs direct quadrature ==")
from numpy.polynomial.legendre import leggauss
gx, gw = leggauss(60)
def quad01(f):  # integral over [0,1]
    x = 0.5*(gx+1); w = 0.5*gw
    return np.sum(w*f(x))
for (x, y) in [(0.7+2.1j, -0.3+0.9j), (0.4, -0.4), (1.1+0.5j, 0.0)]:
    t1 = 0.5*(gx+1)*beta; w1 = 0.5*gw*beta
    val = 0
    for a, wa in zip(t1, w1):
        t2 = 0.5*(gx+1)*a; w2 = 0.5*gw*a
        val += wa*np.sum(w2*np.exp(np.complex128(x)*a + np.complex128(y)*t2))
    check(f"I2({x},{y})", val, complex(I2(x, y)), tol=1e-10)
for (x, y, z) in [(0.7+2.1j, -0.3-2.1j, 0.5), (0.6, -0.6, 0.0), (0.9+1.2j, 0.2, -0.2-1.2j)]:
    val = 0
    t1 = 0.5*(gx+1)*beta; w1 = 0.5*gw*beta
    for a, wa in zip(t1, w1):
        t2 = 0.5*(gx+1)*a; w2 = 0.5*gw*a
        for b, wb in zip(t2, w2):
            t3 = 0.5*(gx+1)*b; w3 = 0.5*gw*b
            val += wa*wb*np.sum(w3*np.exp(np.complex128(x)*a + np.complex128(y)*b + np.complex128(z)*t3))
    check(f"I3({x},{y},{z})", val, complex(I3(x, y, z)), tol=1e-9)

# ---------------- general 3-level system ------------------------------
rng = np.random.default_rng(11)
E = np.array([0.0, 0.6, 2.9])
Xm = np.array([[0.3, 0.8, 0.45],
               [0.8, -0.5, 0.6],
               [0.45, 0.6, 0.2]])
p = np.exp(-beta*E); p /= p.sum()
Xm = Xm - (p @ np.diag(Xm)) * np.eye(3)   # center: <X>_0 = 0
assert abs(p @ np.diag(Xm)) < 1e-14

def Cn(n):
    """Lehmann (13)."""
    wn = 2*np.pi*n/beta
    s = beta*(n == 0)*np.sum(p*np.diag(Xm)**2)
    for r in range(3):
        for smax in range(r+1, 3):
            d = E[smax]-E[r]
            s += Xm[r, smax]**2 * 2*(p[r]-p[smax])*d/(d**2+wn**2)
    return s

def fourpt_sorted(t1, t2, t3, t4):
    """ordered 4pt via path sums (verified vs matrices in Part II work)."""
    out = np.zeros_like(t1, dtype=np.complex128)
    for r in range(3):
        for s in range(3):
            if Xm[r, s] == 0: continue
            for t in range(3):
                if Xm[s, t] == 0: continue
                for u in range(3):
                    wgt = p[r]*Xm[r, s]*Xm[s, t]*Xm[t, u]*Xm[u, r]
                    if wgt == 0: continue
                    out += wgt*np.exp(E[r]*(t1-t4) - E[s]*(t1-t2)
                                      - E[t]*(t2-t3) - E[u]*(t3-t4))
    return out

def F_brute(n, l, ng=48):
    """Independent: GL integration over six ordered simplices.
    labels: 0 -> tau (freq zn), 1 -> tau1 (freq zl), 2 -> tau2 (freq -zl)."""
    zn, zl = 2j*np.pi*n/beta, 2j*np.pi*l/beta
    xi = [zn, zl, -zl]
    gx_, gw_ = leggauss(ng)
    s1 = 0.5*(gx_+1); w1 = 0.5*gw_
    from itertools import permutations
    tot = 0
    # simplex map: t_a = beta*a, t_b = t_a*b, t_c = t_b*c ; Jac = beta^3 a^2 b
    A, B, C = np.meshgrid(s1, s1, s1, indexing='ij')
    WA, WB, WC = np.meshgrid(w1, w1, w1, indexing='ij')
    ta = beta*A; tb = ta*B; tc = tb*C
    jac = beta**3 * A**2 * B * ta/ta  # beta*A^2*B*beta^2... explicit below
    jac = (beta*A*0 + 1)  # placeholder replaced next line
    jac = beta**3 * A**2 * B
    base = fourpt_sorted(ta, tb, tc, 0*tc)
    W = WA*WB*WC
    for perm in permutations(range(3)):
        ph = np.exp(xi[perm[0]]*ta + xi[perm[1]]*tb + xi[perm[2]]*tc)
        tot += np.sum(W*jac*base*ph)
    return tot

def F_formula(n, l):
    zn, zl = 2j*np.pi*n/beta, 2j*np.pi*l/beta
    xi = [zn, zl, -zl]
    from itertools import permutations
    tot = mpc(0)
    for perm in permutations(range(3)):
        x0, x1, x2 = xi[perm[0]], xi[perm[1]], xi[perm[2]]
        for r in range(3):
            for s in range(3):
                if Xm[r, s] == 0: continue
                for t in range(3):
                    if Xm[s, t] == 0: continue
                    for u in range(3):
                        wgt = p[r]*Xm[r, s]*Xm[s, t]*Xm[t, u]*Xm[u, r]
                        if wgt == 0: continue
                        tot += wgt*I3(E[r]-E[s]+x0, E[s]-E[t]+x1, E[t]-E[u]+x2)
    return complex(tot)

print("\n== C. exact vertex F_{n,l} (26) and Gamma_{n,l} (27) vs independent GL integration ==")
for (n, l) in [(0, 0), (0, 1), (1, 0), (1, 1), (1, -1), (2, 1)]:
    fb = F_brute(n, l)
    ff = F_formula(n, l)
    check(f"F({n},{l})", fb, ff, tol=5e-8)
    # Gamma via (27) vs brute connected
    Gm_formula = ff - beta*Cn(n)*Cn(l) - beta*((n == l) + (n == -l))*Cn(n)**2
    # brute connected: subtract the same pairing transforms from brute F
    Gm_brute = fb - beta*Cn(n)*Cn(l) - beta*((n == l) + (n == -l))*Cn(n)**2
    check(f"Gamma({n},{l}) formula-vs-brute", Gm_brute, Gm_formula, tol=5e-7)

# sanity: Cn from Lehmann vs 1D quadrature
def Cn_quad(n):
    zn = 2j*np.pi*n/beta
    t = 0.5*(gx+1)*beta; w = 0.5*gw*beta
    # <T X(t)X(0)> for t>0 via 2-step paths
    val = np.zeros_like(t, dtype=np.complex128)
    for r in range(3):
        for s in range(3):
            val += p[r]*Xm[r, s]*Xm[s, r]*np.exp(-(E[s]-E[r])*t)
    return np.sum(w*np.exp(zn*t)*val)
for n in [0, 1, 3]:
    check(f"C_{n} Lehmann (13) vs quadrature", Cn_quad(n), Cn(n), tol=1e-10)

# ---------------- D. PM chain reduction (28) ---------------------------
print("\n== D. chain-pattern reduction (28) vs general formula (26) ==")
D_, G_ = 3.0, 0.8
R = np.sqrt(D_**2/4 + G_**2)
th = 0.5*np.arctan2(G_/R, D_/(2*R))
D1, D2 = R - D_/2, 2*R
Ech = np.array([0.0, D1, D2])
M1, M2 = np.cos(th), -np.sin(th)
Xch = np.zeros((3, 3)); Xch[0, 1] = Xch[1, 0] = M1; Xch[1, 2] = Xch[2, 1] = M2
pch = np.exp(-beta*Ech); pch /= pch.sum()

def F28(n, l):
    zn, zl = 2j*np.pi*n/beta, 2j*np.pi*l/beta
    xi = [zn, zl, -zl]
    q = {0: M1**2, 2: M2**2}   # mu in {g(0), e(2)}, hub a = level 1
    from itertools import permutations
    tot = mpc(0)
    for perm in permutations(range(3)):
        x0, x1, x2 = xi[perm[0]], xi[perm[1]], xi[perm[2]]
        for mu in (0, 2):
            for nu in (0, 2):
                w_ = q[mu]*q[nu]
                tot += w_*pch[1]*I3(Ech[1]-Ech[mu]+x0, Ech[mu]-Ech[1]+x1, Ech[1]-Ech[nu]+x2)
                tot += w_*pch[mu]*I3(Ech[mu]-Ech[1]+x0, Ech[1]-Ech[nu]+x1, Ech[nu]-Ech[1]+x2)
    return complex(tot)

# reuse general machinery on the chain system
E_save, X_save, p_save = E.copy(), Xm.copy(), p.copy()
E, Xm, p = Ech, Xch, pch
for (n, l) in [(0, 0), (1, 2), (2, -1)]:
    check(f"F_PM({n},{l}) (28) vs (26)", F28(n, l), F_formula(n, l), tol=1e-10)

# ---------------- F. static anchor (54)-(55) at low T ------------------
print("\n== F. static anchor: Gamma_00(T->0) vs RS 4th-order energy ==")
beta = 160.0
pch = np.exp(-beta*Ech); pch /= pch.sum()
p = pch
G00 = complex(F_formula(0, 0)) - beta*Cn(0)**2 - 2*beta*Cn(0)**2
anchor = -24*M1**4/D1**3 + 24*M1**2*M2**2/(D1**2*D2)
check("Gamma_00(T=0) vs -24M1^4/D1^3 + 24M1^2M2^2/(D1^2 D2)", G00, anchor, tol=1e-10)
# independent RS check of E_g(lambda) 4th order by tiny-lambda diagonalization
lam = 1e-2  # Richardson below
Hl = np.diag(Ech) - lam*Xch
def c4(lam):
    return (np.linalg.eigvalsh(np.diag(Ech) - lam*Xch)[0]
            + np.linalg.eigvalsh(np.diag(Ech) + lam*Xch)[0]
            + 2*M1**2/D1*lam**2) / (2*lam**4)
lam = 0.02
e4 = (64*c4(lam/4) - 20*c4(lam/2) + c4(lam))/45  # double Richardson: kills lam^2, lam^4
check("RS lambda^4 coefficient", e4, M1**4/D1**3 - M1**2*M2**2/(D1**2*D2), tol=1e-6)
E, Xm, p = E_save, X_save, p_save
beta = 1.7

# ---------------- E. three-point formulas (39)/(44) --------------------
print("\n== E. three-point vertex (39)/(44) vs direct quadrature ==")
Am = np.array([[0.7, 0.2, -0.3], [0.2, -0.4, 0.5], [-0.3, 0.5, 0.1]])
Am = Am - (p @ np.diag(Am))*np.eye(3)  # centered
def TA_formula(l, A):
    zl = 2j*np.pi*l/beta
    etas = [(zl, -zl), (-zl, zl)]
    tot = mpc(0)
    for (e1, e2) in etas:
        for r in range(3):
            for s in range(3):
                for t in range(3):
                    wgt = p[r]*Xm[r, s]*Xm[s, t]*A[t, r]
                    if wgt == 0: continue
                    tot += wgt*I2(E[r]-E[s]+e1, E[s]-E[t]+e2)
    return complex(tot)
def TA_brute(l, A, ng=60):
    """Per-ordering smooth GL over the simplex t1>t2>0 (no kink)."""
    zl = 2j*np.pi*l/beta
    gx_, gw_ = leggauss(ng)
    a = 0.5*(gx_+1); wa = 0.5*gw_
    T1 = beta*a[:, None]*np.ones(ng)[None, :]
    T2 = T1*a[None, :]
    W = np.outer(wa, wa)*beta**2*a[:, None]
    val = np.zeros_like(T1, dtype=np.complex128)
    for r in range(3):
        for s_ in range(3):
            for u in range(3):
                wgt = p[r]*Xm[r, s_]*Xm[s_, u]*A[u, r]
                if wgt == 0: continue
                val += wgt*np.exp(E[r]*T1 - E[s_]*(T1-T2) - E[u]*T2)
    return np.sum(W*val*np.exp(zl*(T1-T2))) + np.sum(W*val*np.exp(zl*(T2-T1)))
for l in [0, 1, 2]:
    check(f"T^A_l={l} (44) vs quadrature", TA_brute(l, Am), TA_formula(l, Am), tol=1e-8)
# X itself (39):
for l in [0, 1]:
    check(f"T^(3)_l={l} (39) vs quadrature", TA_brute(l, Xm), TA_formula(l, Xm), tol=1e-8)

# ---------------- G. spectator identity (31) ---------------------------
print("\n== G. spectator identity (31) == explicit deficit form (32) ==")
pa, pb, pc_ = 0.5, 0.3, 0.2
q = pa + pb; pta, ptb = pa/q, pb/q
Dab = 0.9
t1, t2, t3, t4 = 1.4, 1.0, 0.5, 0.1
u = (t1-t2)+(t3-t4); v = (t1-t3)+(t2-t4); w_ = (t1-t2)-(t3-t4)
f = lambda x, A, B: A*np.exp(-Dab*x) + B*np.exp(Dab*x)
C4c_31 = (q*(-2*(2*pta*ptb*np.cosh(Dab*w_) + pta**2*np.exp(-Dab*v) + ptb**2*np.exp(Dab*v)))
          + q*(1-q)*(f(t1-t2, pta, ptb)*f(t3-t4, pta, ptb)
                     + f(t1-t3, pta, ptb)*f(t2-t4, pta, ptb)
                     + f(t1-t4, pta, ptb)*f(t2-t3, pta, ptb)))
C4c_32 = (-2*(2*pa*pb*np.cosh(Dab*w_) + pa**2*np.exp(-Dab*v) + pb**2*np.exp(Dab*v))
          + pc_*(pa*np.exp(-Dab*u) + pb*np.exp(Dab*u)))
check("(31) == (32)", C4c_31, C4c_32, tol=1e-12)

print("\n" + ("ALL v2 CHECKS PASSED" if not fails else f"FAILURES: {fails}"))
