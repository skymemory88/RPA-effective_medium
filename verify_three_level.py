"""
Numerical verification of the three-level (anisotropy-Ising) 1/z formulas.

Checks, against direct 3x3 matrix computation in the MF eigenbasis:
  1. The exact path-sum representation of the ordered four-point function.
  2. The bare-propagator channel functions and the generalized sum rule.
  3. The pure-channel semi-invariants S4^(01), S4^(12) (two-level form +
     spectator-population deficit term).
  4. The cross-channel semi-invariant S4^(cross) for the chain pattern
     (exact form, matched-cancellation collapsed form, and the n2->0 form).
  5. The classical (diagonal/Blume-Capel) fourth cumulant kappa4 = Q(1-3Q).
  6. The S=1 anisotropy + transverse-field example: chain matrix elements,
     M01^2 + M12^2 = 1, M02 = 0.
  7. Delta -> 0 classical limit of Jensen's two-level gamma, alpha:
     Sigma(0) -> beta M^2 K(0).
All tolerances 1e-10 unless noted.
"""
import numpy as np

rng = np.random.default_rng(7)
TOL = 1e-10
fails = []


def check(name, a, b, tol=TOL):
    err = abs(a - b)
    ok = err < tol * max(1.0, abs(a), abs(b))
    print(f"{'PASS' if ok else 'FAIL'}  {name}:  {a:+.12e}  vs  {b:+.12e}   |diff|={err:.2e}")
    if not ok:
        fails.append(name)


# ---------------------------------------------------------------- direct machinery
def make_system(E, J, beta):
    E = np.asarray(E, float)
    w = np.exp(-beta * E)
    Z = w.sum()
    n = w / Z
    return E, np.asarray(J, float), n


def Jtau(E, J, t):
    """Heisenberg operator in imaginary time, matrix elements J_{nm} e^{t(E_n - E_m)}."""
    return J * np.exp(t * (E[:, None] - E[None, :]))


def ordered_4pt(E, J, n, ta, tb, tc, td):
    """<J(ta)J(tb)J(tc)J(td)>_0 for ta>=tb>=tc>=td, direct matrix product."""
    M = Jtau(E, J, ta) @ Jtau(E, J, tb) @ Jtau(E, J, tc) @ Jtau(E, J, td)
    return float(np.sum(n * np.diag(M)))


def path_4pt(E, J, n, ta, tb, tc, td):
    """Exact path-sum formula, eq. (path) of the document."""
    N = len(E)
    s = 0.0
    for nu in range(N):
        for m1 in range(N):
            if J[nu, m1] == 0: continue
            for m2 in range(N):
                if J[m1, m2] == 0: continue
                for m3 in range(N):
                    w = n[nu] * J[nu, m1] * J[m1, m2] * J[m2, m3] * J[m3, nu]
                    if w == 0: continue
                    ex = (E[nu] * (ta - td) - E[m1] * (ta - tb)
                          - E[m2] * (tb - tc) - E[m3] * (tc - td))
                    s += w * np.exp(ex)
    return s


def pair_F(E, J, n, t):
    """<J(t)J(0)>_0 ordered (t>=0), via two-step paths."""
    N = len(E)
    s = 0.0
    for nu in range(N):
        for mu in range(N):
            s += n[nu] * J[nu, mu] * J[mu, nu] * np.exp(-(E[mu] - E[nu]) * t)
    return s


def S4_direct(E, J, n, ta, tb, tc, td):
    """Semi-invariant: ordered 4pt of Jtilde minus the three pairings of Ftilde."""
    Jav = float(np.sum(n * np.diag(J)))
    Jt = J - Jav * np.eye(len(E))
    P = ordered_4pt(E, Jt, n, ta, tb, tc, td)
    Ft = lambda t: pair_F(E, Jt, n, t)
    pair = (Ft(ta - tb) * Ft(tc - td) + Ft(ta - tc) * Ft(tb - td)
            + Ft(ta - td) * Ft(tb - tc))
    return P - pair


# ---------------------------------------------------------------- test 1 & 2: path formula, sum rule
print("== 1. path-sum formula vs direct (random 3-level J incl. diagonals) ==")
beta = 1.3
E, J, n = make_system([0.0, 0.7, 3.1],
                      [[0.15, 0.80, 0.35],
                       [0.80, -0.40, 0.55],
                       [0.35, 0.55, 0.25]], beta)
for _ in range(4):
    ts = np.sort(rng.uniform(0, beta, 4))[::-1]
    check(f"path formula t={np.round(ts,3)}",
          ordered_4pt(E, J, n, *ts), path_4pt(E, J, n, *ts))

print("\n== 2. generalized sum rule: (1/beta) sum_n G0 = -<Jtilde^2>_0 ==")
# channel identity route: sum = -[ sum_{nu<mu} M^2 (n_nu+n_mu) + c_d ]
m_diag = np.diag(J)
Jav = float(np.sum(n * m_diag))
lhs = 0.0
for nu in range(3):
    for mu in range(nu + 1, 3):
        lhs -= J[nu, mu] ** 2 * (n[nu] + n[mu])
lhs -= float(np.sum(n * m_diag ** 2) - Jav ** 2)
J2avg = float(np.sum(n * np.diag(J @ J)))
check("channel decomposition of -<Jt^2>", lhs, -(J2avg - Jav ** 2))
# also check coth identity numerically
for nu in range(3):
    for mu in range(nu + 1, 3):
        d = E[mu] - E[nu]
        check(f"coth identity ({nu}{mu})",
              (n[nu] - n[mu]) / np.tanh(beta * d / 2), n[nu] + n[mu])

# ---------------------------------------------------------------- chain pattern setup
print("\n== 3-4. chain pattern: S4 sector formulas ==")
M1, M2 = 0.83, 0.61
d1, d2 = 0.45, 3.4
db = d2 - d1
beta = 1.1
Jc = np.zeros((3, 3)); Jc[0, 1] = Jc[1, 0] = M1; Jc[1, 2] = Jc[2, 1] = M2
E, Jc, n = make_system([0.0, d1, d2], Jc, beta)
n0, n1, n2 = n
f1 = lambda t: n0 * np.exp(-d1 * t) + n1 * np.exp(d1 * t)
f2 = lambda t: n1 * np.exp(-db * t) + n2 * np.exp(db * t)


def S4_11(u, v, w):
    return -2 * (2 * n0 * n1 * np.cosh(d1 * w) + n0 ** 2 * np.exp(-d1 * v)
                 + n1 ** 2 * np.exp(d1 * v)) + n2 * f1(u)


def S4_22(u, v, w):
    return -2 * (2 * n1 * n2 * np.cosh(db * w) + n1 ** 2 * np.exp(-db * v)
                 + n2 ** 2 * np.exp(db * v)) + n0 * f2(u)


def S4_cross_exact(ta, tb, tc, td):
    tab, tbc, tcd = ta - tb, tb - tc, tc - td
    tac, tbd, tad = ta - tc, tb - td, ta - td
    A0 = n0 * np.exp(-d1 * (tab + tcd) - d2 * tbc)
    A1a = n1 * np.exp(d1 * (tad - tbc) - d2 * tcd)
    A1b = n1 * np.exp(d1 * (tad - tbc) - d2 * tab)
    A2 = n2 * np.exp(d2 * tad - d1 * (tab + tcd))
    pair = sum(f1(x) * f2(y) + f2(x) * f1(y)
               for x, y in [(tab, tcd), (tac, tbd), (tad, tbc)])
    return A0 + A1a + A1b + A2 - pair


def S4_cross_collapsed(ta, tb, tc, td):
    """Matched-cancellation form: (1-n1) on the four path exponentials,
    minus the 20 unmatched pairing terms."""
    tab, tbc, tcd = ta - tb, tb - tc, tc - td
    tac, tbd, tad = ta - tc, tb - td, ta - td
    E0 = np.exp(-d1 * tad - db * tbc)
    E1a = np.exp(d1 * tab - db * tcd)
    E1b = np.exp(d1 * tcd - db * tab)
    E2 = np.exp(db * tad + d1 * tbc)
    s = (1 - n1) * (n0 * E0 + n1 * E1a + n1 * E1b + n2 * E2)
    # all 24 pairing cross terms, then remove the 4 matched ones
    terms = []
    for x, y in [(tab, tcd), (tac, tbd), (tad, tbc)]:
        for (ca, ea) in [(n0, -d1 * x), (n1, d1 * x)]:      # f1(x)
            for (cb, eb) in [(n1, -db * y), (n2, db * y)]:  # f2(y)
                terms.append((ca * cb, ea + eb))
        for (ca, ea) in [(n1, -db * x), (n2, db * x)]:      # f2(x)
            for (cb, eb) in [(n0, -d1 * y), (n1, d1 * y)]:  # f1(y)
                terms.append((ca * cb, ea + eb))
    matched = [(n1 * n1, d1 * tab - db * tcd),
               (n1 * n1, d1 * tcd - db * tab),
               (n0 * n1, -d1 * tad - db * tbc),
               (n2 * n1, db * tad + d1 * tbc)]
    for cm, em in matched:
        for i, (c, e) in enumerate(terms):
            if abs(c - cm) < 1e-15 and abs(e - em) < 1e-12:
                terms.pop(i); break
        else:
            raise RuntimeError("matched term not found")
    assert len(terms) == 20
    s -= sum(c * np.exp(e) for c, e in terms)
    return s


for _ in range(4):
    ts = np.sort(rng.uniform(0, beta, 4))[::-1]
    ta, tb, tc, td = ts
    u = ta - tb + tc - td
    v = ta + tb - tc - td
    w = ta - tb - tc + td
    direct = S4_direct(E, Jc, n, *ts)
    formula = (M1 ** 4 * S4_11(u, v, w) + M2 ** 4 * S4_22(u, v, w)
               + M1 ** 2 * M2 ** 2 * S4_cross_exact(*ts))
    check(f"S4 total (sector sum) t={np.round(ts,3)}", direct, formula)
    check("  cross exact vs collapsed", S4_cross_exact(*ts), S4_cross_collapsed(*ts))

# pure-channel checks in isolation (kill the other channel)
print("\n-- pure channel (01) alone, three levels populated --")
Jp = np.zeros((3, 3)); Jp[0, 1] = Jp[1, 0] = 1.0
for _ in range(3):
    ts = np.sort(rng.uniform(0, beta, 4))[::-1]
    ta, tb, tc, td = ts
    u, v, w = ta - tb + tc - td, ta + tb - tc - td, ta - tb - tc + td
    check(f"S4^(01) deficit formula t={np.round(ts,3)}",
          S4_direct(E, Jp, n, *ts), S4_11(u, v, w))
print("-- pure channel (12) alone --")
Jp = np.zeros((3, 3)); Jp[1, 2] = Jp[2, 1] = 1.0
for _ in range(3):
    ts = np.sort(rng.uniform(0, beta, 4))[::-1]
    ta, tb, tc, td = ts
    u, v, w = ta - tb + tc - td, ta + tb - tc - td, ta - tb - tc + td
    check(f"S4^(12) deficit formula t={np.round(ts,3)}",
          S4_direct(E, Jp, n, *ts), S4_22(u, v, w))

# KMS rewriting: every n2-weighted term is an end-of-circle image of an
# n1- or n0-weighted term. Verify the identities and the pointwise low-T limit.
print("\n-- KMS identities and pointwise T->0 limit of the cross sector --")
for _ in range(3):
    t = rng.uniform(0, beta)
    check(f"n2 e^{{db t}} = n1 e^{{-db(beta-t)}}  (t={t:.3f})",
          n2 * np.exp(db * t), n1 * np.exp(-db * (beta - t)))
ts = np.sort(rng.uniform(0, beta, 4))[::-1]
ta, tb, tc, td = ts
tad, tab, tcd = ta - td, ta - tb, tc - td
check("A2 KMS image: n2 e^{d2 tad - d1(tab+tcd)} = n0 e^{-d2(beta-tad) - d1(tab+tcd)}",
      n2 * np.exp(d2 * tad - d1 * (tab + tcd)),
      n0 * np.exp(-d2 * (beta - tad) - d1 * (tab + tcd)))
check("integrated equality: int n1 e^{-db t} = int n2 e^{+db t} over [0,beta]",
      n1 * (1 - np.exp(-db * beta)) / db, n2 * (np.exp(db * beta) - 1) / db)
# pointwise low-T limit: S4_cross -> + e^{-d1 tad - db tbc} at fixed times
bigbeta = 60.0
Ee, Je, ne = make_system([0.0, d1, d2], Jc, bigbeta)
nn0, nn1, nn2 = ne
ff1 = lambda t: nn0 * np.exp(-d1 * t) + nn1 * np.exp(d1 * t)
ff2 = lambda t: nn1 * np.exp(-db * t) + nn2 * np.exp(db * t)
ta, tb, tc, td = 1.9, 1.2, 0.8, 0.3
tbc = tb - tc
cross_lowT = (S4_direct(Ee, Je, ne, ta, tb, tc, td)
              - M1 ** 4 * (-2 * (2 * nn0 * nn1 * np.cosh(d1 * (ta - tb - tc + td))
                                 + nn0 ** 2 * np.exp(-d1 * (ta + tb - tc - td))
                                 + nn1 ** 2 * np.exp(d1 * (ta + tb - tc - td)))
                           + nn2 * ff1(ta - tb + tc - td))
              - M2 ** 4 * (-2 * (2 * nn1 * nn2 * np.cosh(db * (ta - tb - tc + td))
                                 + nn1 ** 2 * np.exp(-db * (ta + tb - tc - td))
                                 + nn2 ** 2 * np.exp(db * (ta + tb - tc - td)))
                           + nn0 * ff2(ta - tb + tc - td))) / (M1 ** 2 * M2 ** 2)
check("T->0 pointwise: S4_cross -> +e^{-d1 tad - db tbc}  (positive)",
      cross_lowT, np.exp(-d1 * (ta - td) - db * tbc), tol=1e-8)

# ---------------------------------------------------------------- classical sector
print("\n== 5. classical Blume-Capel cumulants ==")
bD = 1.7
Z = 2 + np.exp(-bD)
Q = 2 / Z
sig = np.array([1.0, -1.0, 0.0])
p = np.array([1 / Z, 1 / Z, np.exp(-bD) / Z])
c = float(p @ sig ** 2)
k4 = float(p @ sig ** 4) - 3 * c ** 2
check("Q = 2/(2+e^{-beta D})", c, Q)
check("kappa4 = Q(1-3Q)", k4, Q * (1 - 3 * Q))

# ---------------------------------------------------------------- S=1 worked example
print("\n== 6. S=1, H = -D Sz^2 - Gam Sx: chain pattern ==")
D, Gam = 3.0, 0.8
Sz = np.diag([1.0, 0.0, -1.0])
Sx = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]) / np.sqrt(2)
H = -D * Sz @ Sz - Gam * Sx
ev, U = np.linalg.eigh(H)
Sze = U.T @ Sz @ U
idx = np.argsort(ev)
ev, Sze = ev[idx], Sze[np.ix_(idx, idx)]
check("E0 = -D/2 - sqrt(D^2/4+Gam^2)", ev[0], -D / 2 - np.sqrt(D ** 2 / 4 + Gam ** 2))
check("E1 = -D (state |a>)", ev[1], -D)
check("E2 = -D/2 + sqrt(D^2/4+Gam^2)", ev[2], -D / 2 + np.sqrt(D ** 2 / 4 + Gam ** 2))
check("M02 = 0", Sze[0, 2], 0.0, tol=1e-9)
check("diagonals vanish", float(np.abs(np.diag(Sze)).max()), 0.0, tol=1e-9)
check("M01^2 + M12^2 = 1", Sze[0, 1] ** 2 + Sze[1, 2] ** 2, 1.0)
check("Delta1 -> Gam^2/D check (ratio ~1 small Gam)",
      (np.sqrt(D ** 2 / 4 + 0.01 ** 2) - D / 2) / (0.01 ** 2 / D), 1.0, tol=1e-4)

# ---------------------------------------------------------------- Jensen classical limit
print("\n== 7. Delta->0 limit of Jensen's two-level Sigma(0) ==")
# Sigma(0) = alpha + gamma(0) g(0) with a static-only kernel K(iw') = K0 * beta*delta_{n'0}
# lambda_p = (1/beta) K(0) g(0)^p in the classical limit (K(iw'!=0) -> 0 like G0).
beta = 1.9; M = 1.0; K0 = 0.37
for Delta in [1e-2, 1e-3, 1e-4]:
    n01 = np.tanh(beta * Delta / 2)
    g0 = 2 * n01 / Delta
    lam1 = (1 / beta) * K0 * g0
    lam2 = (1 / beta) * K0 * g0 ** 2
    gam0 = (M ** 2 / n01 ** 2) * (lam1 - (1 - n01 ** 2) * K0)
    alp = (M ** 2 / n01 ** 2) * (lam2 - 0.5 * (g0 + beta * (1 - n01 ** 2)) * lam1)
    Sig0 = alp + gam0 * g0
    check(f"Sigma(0)->beta M^2 K(0)  (Delta={Delta})", Sig0, beta * M ** 2 * K0,
          tol=5 * (beta * Delta) ** 2)

print("\n" + ("ALL CHECKS PASSED" if not fails else f"FAILURES: {fails}"))
