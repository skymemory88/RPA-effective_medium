## The central point: your low-$h$ block is a spinodal, not a solver failure

Differentiate your outer residual. With $g\mu_B$ absorbed, (43) gives $dm/dh=-G_0(0;h)$, so

$$F'(h)=r(h)+J_{0,\rm eff}G_0(0;h)=\frac{G_0(0;h)}{\widetilde G_0(0;h)}\Big[1+J_{0,\rm eff}\widetilde G_0(0;h)\Big]=\frac{r(h)}{\;\chi_{\rm unif}(h)/\tilde\chi_0(h)\;}.$$

$F'$ carries the sign of the **uniform mass** $1+J_{0,\rm eff}\widetilde G_0$. So the zero of the uniform mass and the stationary point of $F$ are the *same point*. Four things then line up:

1. Your own class-A separator, $J_{0,\rm eff}|G_0^{\rm inel}(0;h)|$ vs $1+\Sigma_0$, is algebraically $1+J_{0,\rm eff}\widetilde G_0<0$ once $\widetilde G_0\simeq G_0^{\rm inel}/(1+\Sigma_0)$ and the measured $\xi\sim10^{-3}$ suppresses the elastic part. It separated class A with zero exceptions over 24 nodes — that is a stability criterion, not a solver diagnostic.
2. Failed nodes form one **contiguous low-$h$ block** at every field and every resolution, with exactly nested grids. That is a domain boundary.
3. The certified 4 T component terminates with the *outer map still contractive* ($\lambda_{\rm dom}\approx-0.20$) while the barrier $1+J_{\rm sup}x\to0$. A fold has $\sigma_{\min}\to0$ with the barrier finite; you observe the opposite signature.
4. $F(0)=0$, $F'<0$ below $h_c$, $F'>0$ above: the classic van der Waals loop. $[0,h_c)$ is the **unstable middle branch**, where the Gaussian fluctuation operator is not positive definite and no admissible effective-medium closure can exist.

Corollary worth checking against your data: this structure permits **exactly one** nonzero root. Your earlier "two candidate sign changes at 1 T" is then a marker of contaminated or pseudo-root nodes, not a real fold.

So: deflation, homotopy, interval methods, better damping, finer $n_H$ — none of these will produce a root below $h_c$. Stop spending effort there.

## Consequence for eq. (45): re-anchor, don't re-derive

The derivation is fine; the *boundary condition* $H_0(0)=0$ is the casualty. In the ordered phase at $H=0$, $h=0$ is the unstable paramagnet, so the path is inadmissible on $[0,h_c)$ by construction.

Equation (45) is really the ODE $dH_0/dh=r(h)$, and any exact anchor will do. In order of what I'd try:

**(a) Saturation anchor.** $\delta H(h)\equiv H_0-h=\int_0^h\Sigma_0\,dh'$, and $\delta H\to0$ as fluctuations quench — this is precisely the boundary condition Jensen already uses for (46). Then

$$F(h)=h-\int_h^{\infty}\Sigma_0(0;h')\,dh'-J_{0,\rm eff}\,m(h),$$

which never touches $[0,h_c)$. Two things to verify first: (i) measure the large-$h$ decay exponent of $\Sigma_0$ on your certified 4 T nodes to confirm the tail converges; (ii) **validate the whole construction in the paramagnetic phase**, where $[0,\infty)$ is admissible — check $\int_0^\infty\Sigma_0\,dh\approx0$ there. That single test decides whether (a) is usable, costs one PM sweep, and is the highest-value thing on this list.

**(b) Boundary-exact anchor + continuation in $B_x$.** Near $B_c$, $h_c\to0$ and the linearised closure $\mathcal J_{\rm MF}=\mathcal J(q_{\rm ord})\widetilde G_0(0;0)/G_0(0;0)$ is exact. Given that your actual physics question is the adjacent-mode softening at the QCP, this may be all you need — deep-ordered $B_x$ may not have to be certified at all.

**(c) Before either: settle whether $h_c$ is even where you think it is.** $J_{\rm sup}=J_{0,\rm eff}$ enters as a *discrete* channel while the $16^3$ mesh max is $0.08\%$ lower. A single $q$-point carries zero measure in the thermodynamic limit; a 3D quadratic band edge gives $\Phi(x)$ a finite value with a $\sqrt{\varepsilon}$ cusp at $\varepsilon\equiv1+J_{\max}x\to0$, not a pole. Your hard gate converts an integrable cusp into a hard stop. Recompute $\Phi(x)=\langle x/(1+J_qx)\rangle$ as a DOS quadrature with the $\Gamma$ neighbourhood resolved analytically, and look at $\Phi$ alone — no solver, 1D quadrature. If $\Phi$ is finite and smooth through the edge, $h_c$ moves and part of the block is artifact; if it genuinely diverges, the spinodal reading stands unmodified. This is your old program item (3) and it is now load-bearing.

## Algorithm changes I'd make regardless

- **Collapse the outer unknown.** From (39)–(42), $\Sigma(i\omega_n)$ is a *function* of $\{\lambda_1,\lambda_2,\lambda_3,K_0,g(0),n_{01},m,M\}$; the free unknowns are $(\lambda_1,\lambda_2,\lambda_3,x)$ — four scalars, not a Matsubara vector. Solve that $4\times4$ system with a globalised Newton (trust region on $\|R\|^2$), not damped Picard. This answers your numerical Q1 and Q3 at once: at dimension 4 you can do exhaustive multistart *and* deflation, and the noncontractive 1 T node ($\lambda_{\rm dom}=1.36$) becomes a non-issue since Newton doesn't care about the Picard spectral radius.
- **Remove the gates by reparametrisation.** Set $x=-\sigma(u)/J_{\rm sup}$ with $\sigma$ a sigmoid, $u\in\mathbb R$. The barrier is then satisfied identically and can never cause a spurious `no_admissible_static_root`; genuine branch termination shows up cleanly as $u\to+\infty$. This separates "solver left the domain" from "branch ended" — currently you cannot distinguish them.
- **Pseudo-arclength in $(h,\lambda,x)$, bordered Newton.** Yes to your Q2, with the caveat above: monitor $\sigma_{\min}$ of the unbordered Jacobian *and* the barrier distance separately. Fold ⇒ $\sigma_{\min}\to0$, barrier finite. Endpoint ⇒ barrier $\to0$, $\sigma_{\min}$ finite.
- **Make the node evaluator a pure function** of $(B_x,h,\text{seed})$. Not transactional rollback — purity. Rollback fixes the symptom; a mutable `Sigma` returned from a rejected node is a defect that will recur. This eliminates failure class 3 by construction and makes your census reproducible.
- **Quadrature (your Q6).** Once the branch is certified, integrate the ODE with an adaptive embedded RK (Dormand–Prince) treating $r(h)$ as an expensive black box, and near the endpoint substitute $h=h_c+t^2$ — the $\sqrt{\varepsilon}$ edge cusp becomes analytic in $t$ and the error control stops collapsing.

## On the representation mismatch

Your two-level $(\Delta,M^2,m,n_{01})$ sector carries only 42–60% of static $|G_0^{\rm bare}|$ at the low-$h$ nodes, yet the whole cc $G_0$ is renormalised by $1/(1+\Sigma)$. This will shift $h_c$ by a non-negligible amount and could move the block boundary appreciably. It won't remove the spinodal, but I'd fix it before quoting any $h_c$ as physical.

---

One question that would sharpen all of the above: **is your $\Phi(x)$ lattice average currently a mesh sum with the exact-$\Gamma$ channel added as a discrete weight, and if so what weight does it carry?** If it carries $1/N_q$ or larger, test (c) first — everything downstream depends on whether the endpoint is physics or discretisation.