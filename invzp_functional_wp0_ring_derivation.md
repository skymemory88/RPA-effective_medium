# WP0 ring-functional derivation and exact gates

**Recorded:** 2026-07-28

**Scope:** strict scalar, centred local reference, complete unicyclic `O(1/z)` vacuum family

**Status:** derivation record for the isolated `invz_functional/` prototype; it does not modify or
authorize replacement of the Jensen production path

## 1. Scope correction

The provisional WP0 specification described a local `C4` vacuum vertex as if it were the first
non-Gaussian contribution at the same order as the leading ring. That mixes two different
representations and two different orders.

For

\[
\mathcal H_{\rm int}
=-\frac12\sum_{i\ne j}J_{ij}\,\delta X_i\delta X_j,\qquad
\delta X_i=X_i-m_{\rm loc}(h),\qquad
\langle\delta X_i\rangle_0=0,
\]

the Euclidean interaction appearing in the partition function is

\[
A_J=
\frac12\sum_{i\ne j}J_{ij}
\int_0^\beta d\tau\,\delta X_i(\tau)\delta X_j(\tau).
\]

The linked generator is

\[
W-W_0=\log\langle e^{A_J}\rangle_0
=\sum_{p\ge1}\frac1{p!}\langle A_J^p\rangle_{0,c},
\qquad
\Delta F=-\frac1\beta(W-W_0).
\]

Here `m_loc(h)=<X>_0,h` centres the local graph vertices. It is not the independently varied physical
moment `m` in the common `(m,h)` functional of section 4. At pure mean-field stationarity the two are
equal; after the ring correction they differ by `-partial_h f_ring`. The exact Hamiltonian
rearrangement in terms of `X-m` fixes the mean-field source and double counting, but it is not an
off-shell claim that `X-m` remains centred. The ring expansion and every source derivative below use
the same locally centred generator.

With `J_ii=0` and a centred local reference, the `O(J)` term vanishes. At `O(J^2)` the only
surviving connected vacuum graph is the double bond,

\[
W_2=
\frac14\sum_{i\ne j}J_{ij}^2
\sum_{\omega_n}C_{2,i}(i\omega_n)C_{2,j}(-i\omega_n),
\]

and therefore

\[
\Delta F_2=
-\frac1{4\beta}\sum_{i\ne j,n}J_{ij}^2
C_{2,i}(i\omega_n)C_{2,j}(-i\omega_n).
\]

This fixes the sign, factor, and signed-frequency measure.

For a high-density coordination sequence, a connected graph with `V` occupied sites and `E` bonds
scales per site as

\[
z^{V-1}z^{-E}=z^{-(E-V+1)}.
\]

After tadpoles have been removed, a connected unicyclic graph has `E=V` and every occupied local
vertex has degree two. Consequently the complete strict `O(1/z)` vacuum family contains only `C2`
vertices: it is the ring family. Non-Gaussian local cumulants enter at this order through source
derivatives of the ring, not as separate vacuum vertices.

## 2. Complete strict ring functional

Let `C_n` be the diagonal matrix of positive connected local susceptibilities at signed Matsubara
frequency `n`, and let `J` be real symmetric with zero diagonal. The tadpole-subtracted ring
functional is

\[
\Delta F_{\rm ring}
=\frac1{2\beta}\sum_n
\operatorname{Tr}\left[
\log(I-JC_n)+JC_n
\right].
\]

It expands as

\[
\Delta F_{\rm ring}
=-\frac1{4\beta}\sum_n\operatorname{Tr}(JC_n)^2+O(J^3),
\]

and hence reproduces `Delta F_2`. A positive-frequency implementation may replace the signed sum by
weights `w_0=1`, `w_n=2` only after the scalar reality/even-frequency identity has been checked.

Define

\[
L_n=I-JC_n,\qquad
P_n=L_n^{-1}J,\qquad
Q_n=P_n-J.
\]

For an arbitrary diagonal variation `delta C_n`,

\[
\delta\Delta F_{\rm ring}
=-\frac1{2\beta}\sum_n
\operatorname{Tr}(Q_n\,\delta C_n),
\]

because

\[
\delta P_n=P_n\,\delta C_n\,P_n.
\]

The symmetric second variation is therefore

\[
\delta_b\delta_a\Delta F_{\rm ring}
=-\frac1{2\beta}\sum_n\left[
\operatorname{Tr}(P_n\delta_bC_nP_n\delta_aC_n)
+\operatorname{Tr}(Q_n\delta_b\delta_aC_n)
\right].
\]

Cyclicity of the trace makes the first term symmetric under `a <-> b`; the second is symmetric
because all local cumulants are derivatives of the same `W0`.

## 3. One-point and two-point inventory

For a static local source,

\[
\partial_h C_2(i\omega_n)
=C_3(i\omega_n,-i\omega_n,0),
\]

\[
\partial_h^2 C_2(i\omega_n)
=C_4(i\omega_n,-i\omega_n,0,0).
\]

Thus the complete ring one-point correction contains one `C3` insertion:

\[
-\partial_h\Delta F_{\rm ring}
=\frac1{2\beta}\sum_n
\operatorname{Tr}(Q_n C_{3,n}).
\]

Its fixed-reference source curvature contains both allowed classes,

\[
-\partial_h^2\Delta F_{\rm ring}
=\frac1{2\beta}\sum_n\left[
\operatorname{Tr}(P_nC_{3,n}P_nC_{3,n})
+\operatorname{Tr}(Q_nC_{4,n})
\right].
\]

The first term is the two-`C3` insertion class and the second is the one-`C4` insertion class. Neither
may be dropped in the ordered/source-biased problem. At a `Z2`-symmetric paramagnetic point `C3=0`,
but the `C4` term remains.

The physical uniform susceptibility is obtained from the Hessian of the full independently varied
`(m,h)` functional, not by treating the fixed-reference ring curvature above as the complete
response.

## 4. Common strict scalar functional

At uniform source and moment, the strict functional per site is

\[
f(m,h;H)
=f_0(h)+(h-H)m-\frac12J_0m^2+f_{\rm ring}(h).
\]

Its stationarity equations are

\[
0=\partial_m f=h-H-J_0m,
\]

\[
0=\partial_h f=-m_{\rm loc}(h)+m+\partial_hf_{\rm ring}.
\]

Equivalently,

\[
h=H+J_0m,\qquad
m=m_{\rm loc}(h)-\partial_hf_{\rm ring}.
\]

The same scalar function supplies the stationary-state comparison and the response Hessian. No
separate ordered static replacement is used.

For a translational scalar problem with coupling eigenvalues `J_q`,

\[
f_{\rm ring}
=\frac1{2\beta}\sum_nw_n
\left\langle
\log(1-J_qC_n)+J_qC_n
\right\rangle_q .
\]

Writing

\[
q_n=\left\langle\frac{J_q^2C_n}{1-J_qC_n}\right\rangle_q,\qquad
p_n=\left\langle\frac{J_q^2}{(1-J_qC_n)^2}\right\rangle_q ,
\]

gives the directly implementable derivatives

\[
\partial_h f_{\rm ring}
=-\frac1{2\beta}\sum_nw_n q_n\,\partial_hC_n,
\]

\[
\partial_h^2 f_{\rm ring}
=-\frac1{2\beta}\sum_nw_n
\left[
p_n(\partial_hC_n)^2+q_n\partial_h^2C_n
\right].
\]

### 4.1 Temperature derivative and stationary internal energy

For beta differentiation, hold `m`, `h`, physical `H`, `J`, and the integer Matsubara label `n`
fixed. The frequency `omega_n=2*pi*n/beta` therefore changes. For the scalar two-level local object,

\[
C_n=M^2\left[
\frac{\delta^2}{E^2}
\frac{4E\tanh(\beta E)}{4E^2+\omega_n^2}
+\delta_{n0}\,\beta\frac{h^2M^2}{E^2}\operatorname{sech}^2(\beta E)
\right],
\qquad \delta=\Delta/2,
\]

and `partial_beta C_n` is the analytic derivative of this complete expression, including the
frequency denominator and elastic term. With `q_n` defined above,

\[
\partial_\beta[\beta f_{\rm ring}]
=-\frac12\sum_nw_n q_n\,\partial_\beta C_n.
\]

At a stationary solution, the envelope theorem gives the internal energy per site without
differentiating the root path:

\[
u
=u_0(h)+(h-H)m-\frac12J_0m^2
+\partial_\beta[\beta f_{\rm ring}],
\]

\[
u_0(h)=\frac\Delta2-E\tanh(\beta E)
\]

for the prototype convention `Hcf=Delta*|1><1|`. Independent reoptimization at shifted beta verifies
`u=d(beta f_star)/d beta`.

### 4.2 Matsubara-tail control

For `n>0`, the scalar local susceptibility is exactly `C_n=A/(n^2+a^2)`. If
`rho=Jmax*A/(N+1)^2<1`, then

\[
|\log(1-x)+x|\le\frac{x^2}{2(1-|x|)}
\]

gives the rigorous positive-frequency free-energy remainder

\[
|\delta f_{n>N}|
\le
\frac{\mu_2A^2}{6\beta(1-\rho)N^3},
\qquad
\mu_2=\langle J_q^2\rangle_q .
\]

The prototype also performs successive-cutoff Richardson checks of `f`, `u`, the gradient, and the
Hessian, with the expected `N^-3` tail.

## 5. Two-site nonzero-source oracle

Take

\[
\mathcal H_2(h,j)
=\mathcal H_l(h)\otimes I+I\otimes\mathcal H_l(h)-jX\otimes X,
\qquad
\mathcal H_l=\mathcal H_{\rm cf}-hX.
\]

Let `m=<X>_h`, `C_0=partial_h m`, and

\[
G(\tau)=\langle T_\tau X(\tau)X(0)\rangle_h
=m^2+C(\tau).
\]

The exact linked expansion, without assuming that `X` commutes with the local Hamiltonian, is

\[
F_2(h,j)
=2f_0(h)-jm^2-\frac{j^2}{2}B_2+O(j^3),
\]

\[
B_2
=\int_0^\beta d\tau\,[G(\tau)^2-m^4]
=2m^2C_0+\frac1\beta\sum_n C_n^2 .
\]

The stationary mean-field part generates the `-j^2 m^2 C_0` coefficient, while the strict ring
generates `-(j^2/(2 beta)) sum_n C_n^2`. Hence their sum reproduces the complete exact coefficient
through `j^2`.

The corresponding source derivatives are

\[
\mathcal M
=2m+2jmC_0+\frac{j^2}{2}\partial_hB_2+O(j^3),
\]

\[
\chi_{\rm pair}
=2C_0+2j(C_0^2+m\partial_hC_0)
+\frac{j^2}{2}\partial_h^2B_2+O(j^3).
\]

For `Hcf=(Delta/2) sigma_z`, `X=sigma_x`, define

\[
E=\sqrt{(\Delta/2)^2+h^2},\quad
t=\tanh(\beta E),\quad
s^2=\operatorname{sech}^2(\beta E),\quad
a=\frac{h^2}{E^2},\quad
b=\frac{(\Delta/2)^2}{E^2}.
\]

Then

\[
m=\frac hE t,\qquad
C_0=\frac{(\Delta/2)^2}{E^3}t+\beta\frac{h^2}{E^2}s^2,
\]

\[
B_2=
\beta a^2(1-t^4)
+2ab\frac tE
+b^2\left(\frac{\beta s^2}{2}+\frac t{2E}\right).
\]

At `beta=1.7`, `Delta=1.3`, `h=0.37`, the exact pair-series coefficients are:

| quantity | `j^0` | `j^1` | `j^2` |
|---|---:|---:|---:|
| `F` | `-1.58491332420055` | `-0.178566036841724` | `-0.463315074168926` |
| `M` | `0.845141495471190` | `0.824061879107636` | `-0.0982124822695` |
| `chi_pair` | `1.95011577001836` | `0.900465007030457` | `-1.14185578461` |

At `j=0.005`, the exact diagonalization values

\[
F=-1.58581772859600,\quad
\mathcal M=0.849259310666641,\quad
\chi_{\rm pair}=1.95458953899283
\]

differ from the displayed second-order series by `8.7e-9`, `-3.9e-8`, and `-9.7e-9`,
respectively, consistent with an `O(j^3)` remainder.

These displayed absolute energies use the symmetric local convention
`Hcf=(Delta/2)*sigma_z`. The repository Hamiltonian and the prototype use
`Hcf=Delta*|1><1|`, which adds `Delta/2` per site. Consequently the prototype pair values are shifted
by `+Delta` (`F=-0.28581772859600` in this fixture), while every `j` coefficient, moment,
susceptibility, and comparison above is unchanged. Internal energies carry the same declared
constant shift.

## 6. Jensen paramagnetic dynamic gate

At `h=0`, put

\[
n_{01}=\tanh(\beta\Delta/2),\qquad
g_n=\frac{2n_{01}\Delta}{\omega_n^2+\Delta^2}.
\]

For the two-site bond, the leading cavity return is

\[
\frac{K_n}{j^2}=g_n.
\]

The exact local susceptibility coefficient separates as

\[
[j^2]C_{11,n}=g_n^3-g_n\sigma_{2,n}.
\]

The first term is the Gaussian return. The second is the connected-`C4`/self-energy contribution.
With signed Matsubara sums,

\[
L_2=\frac1\beta\sum_ng_n^2
=\frac{n_{01}}{\Delta}
+\frac\beta2(1-n_{01}^2),
\]

\[
L_3=\frac1\beta\sum_ng_n^3
=\frac{3n_{01}^2}{2\Delta^2}
+\frac{3\beta n_{01}(1-n_{01}^2)}{4\Delta}
+\frac{\beta^2(1-n_{01}^2)}4,
\]

and the leading Jensen coefficient is

\[
\sigma_{2,n}
=\frac1{n_{01}^2}\left[
L_3-\frac12\{g_0+\beta(1-n_{01}^2)\}L_2
+g_n\{L_2-(1-n_{01}^2)g_n\}
\right].
\]

It matches the exact two-site Lehmann coefficient. Representative values for `Delta=1` are:

| `beta*Delta` | `n` | exact `[j^2] C11` | Gaussian `g^3` | connected `-g sigma2` |
|---:|---:|---:|---:|---:|
| `0.5` | `0` | `+9.66452772e-4` | `+1.1753186396e-1` | `-1.1656541078e-1` |
| `1.7` | `0` | `+2.193406059e-1` | `+2.6403111406` | `-2.4209705343` |
| `1.7` | `1` | `-7.180523450e-2` | `+8.379573299e-4` | `-7.264319184e-2` |
| `5.0` | `0` | `+2.432447806` | `+7.683024217` | `-5.250576411` |
| `5.0` | `2` | `-2.291213499e-1` | `+1.961614822e-2` | `-2.487374979e-1` |

The exact-minus-Jensen discrepancy of the numerical five-point `j` derivative was below `4.5e-10`
over `beta*Delta={0.5,1.7,5}` and `n=-2,...,2`. The near-cancellations at `n=0` are binding
sign/normalization gates. A one-sided unweighted frequency sum fails this comparison.

The single-frequency term has also been checked at the uncontracted vertex level. The exact two-level
connected `C4` contains

\[
\beta M^2pc\,g_n^2(\delta_{l,n}+\delta_{l,-n}),
\]

in addition to its smooth double-frequency kernel. This anomalous KMS/Hermite contribution collapses
one signed frequency sum and produces the primitive's `-pc*sum K_n^2 g_n^2/2` structure. The complete
formula and normalization are recorded in `invzp_functional_wp0_spec.md`; the dense exact vertex
engine agrees below `2e-13` on the `0,+/-1,+/-2` fixture.

This establishes the strict leading coefficient of the Jensen PM self-energy. It does not derive the
fully resummed production EMT/Jensen solver from a common functional.

## 7. Where an explicit `C4` vacuum really occurs

The first explicit non-Gaussian vacuum graph is `C3-C3` at `J^3`; it vanishes only at a
`Z2`-symmetric reference. A one-`C4` graph first occurs at `J^4`:

\[
W_{422}
=\frac1{8\beta}
\sum_i\sum_{\substack{j\ne i,\;k\ne i\\j\ne k}}\sum_{n,m}
J_{ij}^2J_{ik}^2
C_{4,i}(-n,n,-m,m)C_{2,j}(n)C_{2,k}(m).
\]

Here `W_422` is a dimensionless linked-log contribution and
`Delta F_422=-W_422/beta`. The exact symmetry factor is `1/8` for the ordered distinct-neighbour
sum. Coincident returns `j=k` are excluded: that neighbour carries four legs and requires its own
connected `C4_j` plus Wick pieces. Those terms have fewer embeddings and are `O(1/z^3)`.

In an effective local return line
`K_i(n)=sum_j J_ij^2 C_{2,j}(n)+...`, the same topology is

\[
W_{K^2,C4}
=\frac1{8\beta}\sum_{n,m}
K(n)K(m)C_4(-n,n,-m,m).
\]

It is `O(1/z^2)`, not part of the retained strict ring scope. At that order the complete vacuum also
contains all two-`C3` bicyclic classes: the bare triple bond and its `C2`-decorated
theta/dumbbell/handcuff continuations. Their one-point and two-point source derivatives, and those of
the one-`C4` class, require local `C5` and `C6`. The current repository supplies neither a general
`C3` service nor `C5/C6`, so that extension is explicitly deferred.

## 8. Consequences for implementation

The smallest justified prototype is therefore:

1. scalar two-level local `f0`, `m`, `C2`, `partial_h C2`, and `partial_h^2 C2`;
2. the finite-mode strict ring and its analytic first and second source derivatives;
3. the independently varied `(m,h)` common functional and Hessian;
4. exact two-site coefficient and thermodynamic-derivative tests;
5. no EMT, no dressed self-energy, no ordered `tanh/xi` replacement, and no production dispatch.

The quadratic `Phi_Sigma_PM[K]` remains a useful integrability identity of the resummed PM map. It is
not used as the vacuum functional of this strict prototype. A future EMT/2PI work package must derive
the impurity, lattice-log, Legendre, and double-counting terms together before varying `K` or `G`.
